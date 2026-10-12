"""STEP / IGES import: entities from a CAD model read through OpenCASCADE.

The optional module ``opengrid._occt`` (backends/occt, built with OG_WITH_OCCT) reads and heals
the file and returns its topology with B-spline geometry (see cad_import.hpp). This module turns
that into OpenGrid entities, created as one undo step:

  * every CAD edge -> a ``nurbs`` curve and a connector on it. The connector's dimension follows
    the target ``spacing`` (default: 1/40 of the model's bounding-box diagonal) and resolves
    curvature: at most ``max_turn_angle`` degrees of tangent turning per segment;
  * every CAD face -> a domain. Planar faces get none (the meshers use the plane); other faces
    get a ``nurbs_surface`` (the face's surface trimmed to its (u, v) bounds). With
    ``domains="unstructured"`` all faces become unstructured domains (outer wire, then holes);
    with ``"auto"``, four-sided faces with one wire become structured domains (TFI, on the
    surface when curved) and the connectors of their opposite sides share one dimension (the
    largest needed, propagated through chains of four-sided faces); ``"none"`` makes no domains.
  * faces that cannot be meshed automatically are reported, not guessed at: faces with a
    degenerate edge (sphere poles, cone apexes), faces without a surface, and faces whose wire
    lost an edge. Edges used by one face (open shells) or by more than two faces (non-manifold)
    are reported with their connectors.

Loop orientation: domains in a plane z = const use the 2D convention (outer loop
counter-clockwise, holes clockwise in (x, y)), so their loops are oriented here; other planar and
surface domains take the wires' orientation (OCCT keeps outer and inner wires opposite), and
their frames orient the outer loop.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from opengrid.validation import Issue, ValidationError

DOMAIN_MODES = ("unstructured", "auto", "none")


@dataclass
class CadImport:
    """What an import created and what it could not handle."""

    source: str
    curves: list[str] = field(default_factory=list)
    surfaces: list[str] = field(default_factory=list)
    connectors: list[str] = field(default_factory=list)
    domains: list[str] = field(default_factory=list)
    structured: list[str] = field(default_factory=list)  # the structured ones among domains
    skipped_faces: list[tuple[int, str]] = field(default_factory=list)  # (face number, reason)
    free_edges: list[str] = field(default_factory=list)  # connectors on edges of one face
    non_manifold_edges: list[str] = field(default_factory=list)  # edges of > 2 faces
    log: list[str] = field(default_factory=list)
    solids: int = 0
    shells: int = 0
    max_tolerance: float = 0.0
    spacing: float = 0.0

    def report(self) -> str:
        lines = [
            f"{self.source}: {self.solids} solid(s), {self.shells} shell(s); created "
            f"{len(self.curves)} curves, {len(self.surfaces)} surfaces, "
            f"{len(self.connectors)} connectors, {len(self.domains)} domains "
            f"({len(self.structured)} structured); spacing {self.spacing:.4g}, "
            f"largest CAD tolerance {self.max_tolerance:.3g}",
        ]
        for line in self.log:
            lines.append(f"  {line}")
        for face, reason in self.skipped_faces:
            lines.append(f"  face {face}: no domain ({reason})")
        if self.free_edges:
            lines.append(
                f"  {len(self.free_edges)} free edge(s) (one face; the shell is open): "
                + ", ".join(self.free_edges)
            )
        if self.non_manifold_edges:
            lines.append(
                f"  {len(self.non_manifold_edges)} non-manifold edge(s) (more than two faces): "
                + ", ".join(self.non_manifold_edges)
            )
        return "\n".join(lines)


def read_cad_model(
    path, tolerance: float, scale: float, heal: bool, sew: bool, split_closed: bool
) -> dict:
    """The healed CAD model from opengrid._occt.read_cad."""
    try:
        from opengrid import _occt
    except ImportError as error:
        raise RuntimeError(
            "CAD import needs OpenGrid built with OpenCASCADE (CMake option OG_WITH_OCCT, "
            "OpenCASCADE found at configure time)"
        ) from error
    path = Path(path)
    if not path.is_file():
        raise ValidationError([Issue(str(path), "file not found")])
    try:
        return _occt.read_cad(
            str(path),
            tolerance=float(tolerance),
            scale=float(scale),
            heal=bool(heal),
            sew=bool(sew),
            split_closed=bool(split_closed),
        )
    except (RuntimeError, ValueError) as error:
        raise ValidationError([Issue(str(path), str(error))]) from error


def check_options(spacing, max_turn_angle, domains, tolerance, scale) -> None:
    issues = []
    if spacing is not None and not (isinstance(spacing, (int, float)) and spacing > 0):
        issues.append(Issue("spacing", f"must be a positive number or None, got {spacing!r}"))
    if not (isinstance(max_turn_angle, (int, float)) and 0 < max_turn_angle <= 180):
        issues.append(Issue("max_turn_angle", f"must be in (0, 180], got {max_turn_angle!r}"))
    if domains not in DOMAIN_MODES:
        issues.append(Issue("domains", f"must be one of {', '.join(DOMAIN_MODES)}"))
    for key, value in (("tolerance", tolerance), ("scale", scale)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
            issues.append(Issue(key, f"must be a positive number, got {value!r}"))
    if issues:
        raise ValidationError(issues)


def _curve(edge: dict):
    from opengrid import _core

    c = edge["curve"]
    return _core.Nurbs(
        int(c["degree"]),
        [float(k) for k in c["knots"]],
        [tuple(float(x) for x in p) for p in c["poles"]],
        [float(w) for w in c["weights"]],
    )


def _turning(curve) -> float:
    """Total tangent turning of a curve in degrees (64 samples)."""
    tangents = np.array([curve.d1(t) for t in np.linspace(0.0, 1.0, 65)])
    lengths = np.linalg.norm(tangents, axis=1)
    tangents = tangents[lengths > 0] / lengths[lengths > 0, None]
    if len(tangents) < 2:
        return 0.0
    cosines = np.clip(np.sum(tangents[:-1] * tangents[1:], axis=1), -1.0, 1.0)
    return float(np.degrees(np.sum(np.arccos(cosines))))


def _dimension(length: float, turning: float, spacing: float, max_turn: float) -> int:
    by_spacing = int(np.ceil(length / spacing - 1e-9)) + 1
    by_turning = int(np.ceil(turning / max_turn - 1e-9)) + 1
    return max(2, by_spacing, by_turning)


def _find(parent: list[int], i: int) -> int:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def _samples(loop, curves) -> np.ndarray:
    """Points along a loop of (edge, reversed) items, closing point not repeated."""
    parts = []
    for edge, reversed_ in loop:
        points = curves[edge].eval_many(np.linspace(0.0, 1.0, 9))
        if reversed_:
            points = points[::-1]
        parts.append(points[:-1])
    return np.concatenate(parts)


def _signed_area_xy(points: np.ndarray) -> float:
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _reverse(loop):
    return [(edge, not reversed_) for edge, reversed_ in reversed(loop)]


def _in_z_plane(loops, curves) -> bool:
    from opengrid.model.frames import PLANARITY

    points = np.concatenate([_samples(loop, curves) for loop in loops])
    size = float(np.max(np.ptp(points, axis=0))) or 1.0
    return float(np.ptp(points[:, 2])) <= PLANARITY * max(size, 1.0)


def _orient_z_loops(loops, curves):
    """Outer loop counter-clockwise and holes clockwise in (x, y)."""
    result = []
    for k, loop in enumerate(loops):
        area = _signed_area_xy(_samples(loop, curves))
        if (k == 0) != (area > 0):
            loop = _reverse(loop)
        result.append(loop)
    return result


def cad_entities(
    model: dict,
    base: str,
    spacing: float | None,
    max_turn_angle: float,
    domains: str,
    generator: dict | None,
    source: str,
) -> tuple[list[tuple[str, dict, str]], CadImport]:
    """Entity items (kind, params, name) for CreateEntities and the import summary."""
    result = CadImport(
        source=source,
        log=list(model["log"]),
        solids=int(model["solids"]),
        shells=int(model["shells"]),
        max_tolerance=float(model["max_tolerance"]),
    )
    edges = model["edges"]
    faces = model["faces"]
    if not edges:
        raise ValidationError([Issue(source, "the file holds no edges")])
    curves = [_curve(edge) for edge in edges]
    if spacing is None:
        poles = np.concatenate([np.asarray(e["curve"]["poles"]) for e in edges])
        diagonal = float(np.linalg.norm(np.ptp(poles, axis=0)))
        spacing = diagonal / 40.0
    result.spacing = float(spacing)

    dims = []
    for edge, curve in zip(edges, curves, strict=True):
        dims.append(_dimension(float(edge["length"]), _turning(curve), spacing, max_turn_angle))

    # Faces as loops of (edge, reversed); decide which can be meshed and how.
    face_loops: list[list[list[tuple[int, bool]]] | None] = []
    structured_faces = set()
    for j, face in enumerate(faces):
        number = j + 1
        loops = []
        for wire in face["wires"]:
            loops.append(list(zip(wire["edges"], wire["reversed"], strict=True)))
        reason = None
        if face["degenerate"]:
            reason = "it has a degenerate edge (pole or apex); mesh it manually"
        elif not loops or not loops[0]:
            reason = "it has no wire"
        elif face["type"] != "plane" and face["surface"] is None:
            reason = "its surface could not be converted to a B-spline"
        if reason is not None:
            result.skipped_faces.append((number, reason))
            face_loops.append(None)
            continue
        face_loops.append(loops)
        four_sided = len(loops) == 1 and len(loops[0]) == 4
        if domains == "auto" and four_sided and len({e for e, _r in loops[0]}) == 4:
            structured_faces.add(j)

    # Opposite sides of structured faces share a dimension (union-find over edges).
    parent = list(range(len(edges)))
    for j in sorted(structured_faces):
        loop = face_loops[j][0]
        for a, b in ((0, 2), (1, 3)):
            ra = _find(parent, loop[a][0])
            rb = _find(parent, loop[b][0])
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    group_dim: dict[int, int] = {}
    for i in range(len(edges)):
        root = _find(parent, i)
        group_dim[root] = max(group_dim.get(root, 2), dims[i])
    dims = [group_dim[_find(parent, i)] for i in range(len(edges))]

    items: list[tuple[str, dict, str]] = []
    curve_names = []
    connector_names = []
    for i, edge in enumerate(edges):
        c = edge["curve"]
        curve_name = f"{base}-curve-{i + 1}"
        params = {
            "degree": int(c["degree"]),
            "knots": [float(k) for k in c["knots"]],
            "control_points": [tuple(float(x) for x in p) for p in c["poles"]],
            "weights": [float(w) for w in c["weights"]],
        }
        items.append(("nurbs", params, curve_name))
        connector_name = f"{base}-con-{i + 1}"
        items.append(("connector", {"segments": curve_name, "dimension": dims[i]}, connector_name))
        curve_names.append(curve_name)
        connector_names.append(connector_name)
        if edge["faces"] == 1:
            result.free_edges.append(connector_name)
        if edge["faces"] > 2:
            result.non_manifold_edges.append(connector_name)
    result.curves = curve_names
    result.connectors = connector_names

    if domains == "none":
        return items, result
    for j, face in enumerate(faces):
        loops = face_loops[j]
        if loops is None:
            continue
        surface_name = None
        if face["type"] != "plane":
            s = face["surface"]
            nu = int(s["nu"])
            nv = int(s["nv"])
            poles = np.asarray(s["poles"], dtype=float).reshape(nu, nv, 3)
            weights = np.asarray(s["weights"], dtype=float).reshape(nu, nv)
            surface_name = f"{base}-surface-{j + 1}"
            params = {
                "degree_u": int(s["degree_u"]),
                "degree_v": int(s["degree_v"]),
                "knots_u": [float(k) for k in s["knots_u"]],
                "knots_v": [float(k) for k in s["knots_v"]],
                "control_points": [[tuple(p) for p in row.tolist()] for row in poles],
                "weights": [list(row) for row in weights.tolist()],
            }
            items.append(("nurbs_surface", params, surface_name))
            result.surfaces.append(surface_name)
        elif _in_z_plane(loops, curves):
            loops = _orient_z_loops(loops, curves)
        domain_name = f"{base}-face-{j + 1}"
        if j in structured_faces:
            edges_param = [
                [{"connector": connector_names[e], "reversed": bool(r)}] for e, r in loops[0]
            ]
            params = {"edges": edges_param, "generator": None, "surface": surface_name}
            items.append(("domain", params, domain_name))
            result.structured.append(domain_name)
        else:
            loops_param = [
                [{"connector": connector_names[e], "reversed": bool(r)} for e, r in loop]
                for loop in loops
            ]
            params = {"loops": loops_param, "generator": generator, "surface": surface_name}
            items.append(("unstructured_domain", params, domain_name))
        result.domains.append(domain_name)
    return items, result
