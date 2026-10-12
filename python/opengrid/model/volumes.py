"""Unstructured 3D volumes: tetrahedra (and pyramids on quad faces) filling the region enclosed
by a set of domains.

A volume lists its boundary ``faces`` (structured or unstructured domains, planar or on
surfaces; any number of closed shells: the outer boundary and inner bodies). The face meshes are
merged into one surface by exact node coordinates (faces share connector nodes, so they merge
bitwise) and checked:
  * watertight: every surface edge is used by exactly two face cells (else the open or
    non-manifold edges are reported with the domains that use them);
  * orientable: each connected shell can be oriented consistently (cells flipped as needed);
  * the outer shell (largest enclosed volume) is oriented outward and the inner shells inward,
    so every surface normal points out of the meshed region.
Quad faces: with ``pyramids`` "auto" (default) quads of faces shared with structured blocks
become pyramids, so the interface stays conforming, and other quads are split into two triangles;
"always" / "never" force either. A pyramid's apex is placed inside
the region at 0.3 times the quad's mean edge length from its centroid along the inward normal
(side faces at ~31 degrees: pyramids on faces meeting at a right angle stay clear of each other;
much sharper concave edges are not supported; along edges where two pyramid faces meet, the
tetrahedra between the pyramids are thin - minimum dihedral angles of a few degrees - which is why
"auto" splits boundary quads), and the
quad's four side triangles replace it in the surface given to the tetrahedral backend; the
pyramid shares the quad with the face, so neighbouring structured blocks stay conforming.
Backends implement ``TetBackend.mesh(points, triangles, max_size, grading)``: ``points`` are the
surface nodes, ``triangles`` the surface oriented out of the region; they return all points
(the surface nodes first, unchanged) and the tetrahedra. Netgen (LGPL-2.1, conda-forge package
``netgen``) is the first backend, used through its Python module at run time; mmg (optimisation)
and an own 3D Delaunay mesher can be added behind the same interface. The backend is a Python
class rather than the C++ interface sketched in ARCHITECTURE: Netgen's Python module is the
supported, dynamically linked way to drive it, and keeps it an optional run-time dependency.
"""

import hashlib
import math
from collections import OrderedDict

import numpy as np

from opengrid import _core
from opengrid.validation import Issue, ValidationError

BACKENDS = ("netgen",)
DEFAULT_VOLUME_GENERATOR = {
    "backend": "netgen",
    "max_size": None,
    "grading": 0.3,
    "pyramids": "auto",
    "inflation": None,
}
PYRAMID_MODES = ("auto", "always", "never")
# Internal settings of the 3D inflation (opengrid._core.inflate_3d; not in projects).
INFLATION_INTERNALS = {"min_jacobian": 0.2, "max_scale": 2.0, "max_slope": 0.5, "proximity": 1.0}
# Apex height of pyramids on quads of an inflation front, in mean edge lengths.
FRONT_APEX = 0.5
DEFAULT_VOLUME_INFLATION = {
    "faces": (),
    "first_height": None,
    "growth_rate": 1.2,
    "layers": 30,
    "total_thickness": None,
    "transition": True,
    "slide": True,
}


def normalize_volume_inflation(raw, faces, store, label: str, issues: list) -> dict | None:
    """Canonical inflation settings of a volume; wall faces (names or ids) become ids and must be
    faces of the volume."""
    if raw is None:
        return None
    if not isinstance(raw, dict):
        issues.append(Issue(label, f"expected a dict or None, got {raw!r}"))
        return None
    for key in sorted(set(raw) - set(DEFAULT_VOLUME_INFLATION)):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    merged = dict(DEFAULT_VOLUME_INFLATION)
    merged.update(raw)
    walls = []
    refs = merged["faces"]
    if isinstance(refs, (str, int)):
        refs = [refs]
    if not isinstance(refs, (list, tuple)) or not refs:
        issues.append(Issue(f"{label}.faces", "at least one wall face is required"))
        refs = []
    for k, ref in enumerate(refs):
        where = f"{label}.faces[{k}]"
        if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
            issues.append(Issue(where, f"no domain {ref!r}"))
            continue
        entity = store.get(ref)
        if entity.id not in faces:
            issues.append(Issue(where, f"{entity.name!r} is not a face of the volume"))
            continue
        if entity.id not in walls:
            walls.append(entity.id)
    out = {"faces": tuple(walls)}
    for key, optional in (
        ("first_height", False),
        ("growth_rate", False),
        ("total_thickness", True),
    ):
        value = merged[key]
        if value is None and optional:
            out[key] = None
            continue
        if hasattr(value, "item") and not isinstance(value, (int, float)):
            value = value.item()
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
            issues.append(Issue(f"{label}.{key}", f"must be positive, got {value!r}"))
            out[key] = None
            continue
        out[key] = float(value)
    if out["growth_rate"] is not None and out["growth_rate"] < 1.0:
        issues.append(Issue(f"{label}.growth_rate", "must be >= 1"))
    layers = merged["layers"]
    if isinstance(layers, bool) or not isinstance(layers, int) or layers < 1:
        issues.append(Issue(f"{label}.layers", f"must be an integer >= 1, got {layers!r}"))
    out["layers"] = layers
    for key in ("transition", "slide"):
        if not isinstance(merged[key], bool):
            issues.append(Issue(f"{label}.{key}", f"must be true or false, got {merged[key]!r}"))
        out[key] = merged[key]
    return out


def normalize_volume(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"faces", "generator"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a volume"))
    faces = []
    refs = raw.get("faces")
    if isinstance(refs, (str, int)):
        refs = [refs]
    if not isinstance(refs, (list, tuple)) or len(refs) < 2:
        issues.append(Issue(f"{where}: faces", "at least two domains are needed"))
        refs = []
    for k, ref in enumerate(refs):
        label = f"{where}: faces[{k}]"
        if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
            issues.append(Issue(label, f"no domain {ref!r}"))
            continue
        entity = store.get(ref)
        if entity.kind not in ("domain", "unstructured_domain"):
            issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a domain"))
            continue
        if entity.id in faces:
            issues.append(Issue(label, f"{entity.name!r} listed twice"))
            continue
        faces.append(entity.id)
    generator = dict(DEFAULT_VOLUME_GENERATOR)
    graw = raw.get("generator") or {}
    if not isinstance(graw, dict):
        issues.append(Issue(f"{where}: generator", "expected a dict"))
        graw = {}
    for key in sorted(set(graw) - set(DEFAULT_VOLUME_GENERATOR)):
        issues.append(Issue(f"{where}: generator.{key}", "unexpected key"))
    generator.update(graw)
    if generator["backend"] not in BACKENDS:
        issues.append(Issue(f"{where}: generator.backend", f"must be one of {', '.join(BACKENDS)}"))
    for key, optional in (("max_size", True), ("grading", False)):
        value = generator[key]
        if value is None and optional:
            continue
        if hasattr(value, "item") and not isinstance(value, (int, float)):
            value = value.item()
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
            issues.append(Issue(f"{where}: generator.{key}", f"must be positive, got {value!r}"))
            continue
        generator[key] = float(value)
    if generator["pyramids"] not in PYRAMID_MODES:
        issues.append(
            Issue(f"{where}: generator.pyramids", f"must be one of {', '.join(PYRAMID_MODES)}")
        )
    generator["inflation"] = normalize_volume_inflation(
        generator["inflation"], faces, store, f"{where}: generator.inflation", issues
    )
    if issues:
        raise ValidationError(issues)
    return {"faces": tuple(faces), "generator": generator}


def volume_script_params(params: dict, store) -> dict:
    result = dict(params)
    result["faces"] = tuple(store.get(f).name for f in params["faces"])
    inflation = params["generator"]["inflation"]
    if inflation is not None:
        generator = dict(params["generator"])
        generator["inflation"] = dict(inflation)
        generator["inflation"]["faces"] = tuple(store.get(f).name for f in inflation["faces"])
        result["generator"] = generator
    return result


# ---------------------------------------------------------------------------- the surface
def face_cells(session, domain_id: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """Points and cells (triangles and quads, as index arrays) of a domain's mesh."""
    entity = session.store.get(domain_id)
    if entity.kind == "unstructured_domain":
        mesh = session.unstructured_mesh(entity.name)
        cells = [np.asarray(t) for t in mesh["triangles"]] + [np.asarray(q) for q in mesh["quads"]]
        return mesh["points"], cells
    grid = session.domain_grid(entity.name)
    nj, ni, _ = grid.shape
    index = np.arange(ni * nj).reshape(nj, ni)
    quads = np.stack(
        [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
    ).reshape(-1, 4)
    return grid.reshape(-1, 3), list(quads)


def _signed_volume(points, cells) -> float:
    total = 0.0
    for cell in cells:
        p = points[cell]
        for k in range(1, len(cell) - 1):
            total += np.dot(p[0], np.cross(p[k], p[k + 1]))
    return total / 6.0


def assemble_shell(session, face_ids, where: str) -> dict:
    """Merged, oriented boundary surface: {"points", "cells", "cell_face", "shells"}."""
    store = session.store
    pool = []
    raw_cells = []
    cell_face = []
    offset = 0
    for f in face_ids:
        points, cells = face_cells(session, f)
        pool.append(points)
        for cell in cells:
            raw_cells.append(cell + offset)
            cell_face.append(f)
        offset += len(points)
    allpoints = np.concatenate(pool)
    points, inverse = np.unique(allpoints, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    cells = [inverse[c] for c in raw_cells]
    # Edges -> cells.
    edges: dict[tuple, list[int]] = {}
    for c, cell in enumerate(cells):
        n = len(cell)
        for k in range(n):
            a, b = int(cell[k]), int(cell[(k + 1) % n])
            edges.setdefault((min(a, b), max(a, b)), []).append(c)
    bad = {e: u for e, u in edges.items() if len(u) != 2}
    if bad:
        names = sorted({store.get(cell_face[c]).name for u in bad.values() for c in u})
        open_edges = sum(1 for u in bad.values() if len(u) == 1)
        non_manifold = len(bad) - open_edges
        raise ValidationError(
            [
                Issue(
                    where,
                    f"the faces do not form closed shells: {open_edges} open and {non_manifold} "
                    f"non-manifold edges (domains {', '.join(names)})",
                )
            ]
        )
    # Shells (connected components) and consistent orientation by breadth-first propagation:
    # neighbours across an edge must traverse it in opposite directions.
    shell_of = [-1] * len(cells)
    shells = []
    for start in range(len(cells)):
        if shell_of[start] >= 0:
            continue
        shell_id = len(shells)
        members = [start]
        shell_of[start] = shell_id
        queue = [start]
        while queue:
            c = queue.pop()
            cell = cells[c]
            n = len(cell)
            for k in range(n):
                a, b = int(cell[k]), int(cell[(k + 1) % n])
                for o in edges[(min(a, b), max(a, b))]:
                    if o == c:
                        continue
                    other = cells[o]
                    m = len(other)
                    same = any(
                        int(other[j]) == a and int(other[(j + 1) % m]) == b for j in range(m)
                    )
                    if shell_of[o] < 0:
                        if same:
                            cells[o] = other[::-1].copy()
                        shell_of[o] = shell_id
                        members.append(o)
                        queue.append(o)
                    elif same:
                        raise ValidationError(
                            [Issue(where, "a shell is not orientable (Moebius-like faces)")]
                        )
        shells.append(members)
    volumes = [_signed_volume(points, [cells[c] for c in m]) for m in shells]
    outer = int(np.argmax([abs(v) for v in volumes]))
    for s, members in enumerate(shells):
        want_positive = s == outer
        if (volumes[s] > 0) != want_positive:
            for c in members:
                cells[c] = cells[c][::-1].copy()
    if len(shells) > 1:
        # Inner shells must lie inside the outer one (by bounding box).
        lo = points[np.unique(np.concatenate([cells[c] for c in shells[outer]]))].min(axis=0)
        hi = points[np.unique(np.concatenate([cells[c] for c in shells[outer]]))].max(axis=0)
        for s, members in enumerate(shells):
            if s == outer:
                continue
            ids = np.unique(np.concatenate([cells[c] for c in members]))
            if np.any(points[ids].min(axis=0) < lo) or np.any(points[ids].max(axis=0) > hi):
                raise ValidationError([Issue(where, "an inner shell lies outside the outer one")])
    return {
        "points": points,
        "cells": cells,
        "cell_face": cell_face,
        "shells": len(shells),
        "volume": abs(volumes[outer]) - sum(abs(v) for s, v in enumerate(volumes) if s != outer),
    }


# ------------------------------------------------------------------------------ backends
class TetBackend:
    name = "abstract"

    def mesh(self, points, triangles, max_size, grading):  # pragma: no cover - interface
        raise NotImplementedError


class NetgenBackend(TetBackend):
    """Netgen's volume mesher from a fixed surface mesh (nglib through the Python module)."""

    name = "netgen"

    def mesh(self, points, triangles, max_size, grading):
        try:
            import netgen.meshing as ng
        except ImportError as error:
            raise ValidationError(
                [Issue("netgen", f"the Netgen backend is not installed ({error})")]
            ) from None
        m = ng.Mesh(dim=3)
        ids = [m.Add(ng.MeshPoint(ng.Pnt(*map(float, p)))) for p in points]
        fd = m.Add(ng.FaceDescriptor(bc=1, domin=1, domout=0, surfnr=1))
        for t in triangles:
            m.Add(ng.Element2D(fd, [ids[int(t[0])], ids[int(t[1])], ids[int(t[2])]]))
        options = {"grading": grading}
        if max_size is not None:
            options["maxh"] = max_size
        m.GenerateVolumeMesh(ng.MeshingParameters(**options))
        coords = np.array([tuple(p) for p in m.Points()], dtype=float)
        tets = np.array([[v.nr - 1 for v in e.vertices] for e in m.Elements3D()], dtype=np.int64)
        if len(tets) == 0:
            raise ValidationError([Issue("netgen", "the volume mesher produced no tetrahedra")])
        return coords, tets


def backend(name: str) -> TetBackend:
    if name == "netgen":
        return NetgenBackend()
    raise ValidationError([Issue("backend", f"unknown backend {name!r}")])


# ------------------------------------------------------------------------------ cell types
# Cell arrays of a volume mesh, in this order everywhere (quality arrays, exports, picking):
# (key, nodes per cell, VTK cell type). Conventions: a tetrahedron's (0, 1, 2) normal points to
# node 3, a pyramid's base normal to the apex, a prism's first triangle normal to the second, a
# hexahedron's first quad normal to the second (as CGNS, Gmsh and VTK).
CELL_TYPES = (("tets", 4, 10), ("pyramids", 5, 14), ("prisms", 6, 13), ("hexes", 8, 12))
# Faces of each type, oriented outward.
CELL_FACES = {
    "tets": ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)),
    "pyramids": ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
    "prisms": ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
    "hexes": ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}
# Corners (node, three neighbours) whose triple product is positive for a valid cell, and the
# factor scaling the ideal cell's corner value to 1.
CELL_CORNERS = {
    "tets": (((0, 1, 2, 3), (1, 2, 0, 3), (2, 0, 1, 3), (3, 0, 2, 1)), math.sqrt(2.0)),
    "pyramids": (((0, 1, 3, 4), (1, 2, 0, 4), (2, 3, 1, 4), (3, 0, 2, 4)), math.sqrt(2.0)),
    "prisms": (
        ((0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0), (4, 3, 5, 1), (5, 4, 3, 2)),
        2.0 / math.sqrt(3.0),
    ),
    "hexes": (
        (
            (0, 1, 3, 4),
            (1, 2, 0, 5),
            (2, 3, 1, 6),
            (3, 0, 2, 7),
            (4, 7, 5, 0),
            (5, 4, 6, 1),
            (6, 5, 7, 2),
            (7, 6, 4, 3),
        ),
        1.0,
    ),
}
# Ideal dihedral angles (smallest, largest) of each type, for the equiangle skewness.
CELL_IDEAL_ANGLES = {
    "tets": (math.degrees(math.acos(1.0 / 3.0)),) * 2,
    "pyramids": (math.degrees(math.atan(math.sqrt(2.0))), 180.0 - math.degrees(math.acos(1 / 3))),
    "prisms": (60.0, 90.0),
    "hexes": (90.0, 90.0),
}


def _edge_face_pairs(faces):
    """(face a, face b) for every edge: the two faces traversing it in opposite directions."""
    owner = {}
    for f, face in enumerate(faces):
        for k in range(len(face)):
            owner[(face[k], face[(k + 1) % len(face)])] = f
    return tuple(sorted({tuple(sorted((owner[(a, b)], owner[(b, a)]))) for (a, b) in owner}))


CELL_EDGE_FACES = {key: _edge_face_pairs(faces) for key, faces in CELL_FACES.items()}


def volume_cells(mesh):
    """[(key, cells (m, n), VTK type)] of a volume mesh, in the canonical order."""
    out = []
    for key, n, vtk in CELL_TYPES:
        cells = np.asarray(mesh.get(key, np.zeros((0, n), dtype=np.int64)), dtype=np.int64)
        out.append((key, cells.reshape(-1, n), vtk))
    return out


def _area_vector(p):
    """Area vectors of triangles (m, 3, 3) or bilinear quads (m, 4, 3)."""
    if p.shape[1] == 3:
        return 0.5 * np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    return 0.5 * np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])


def cell_volumes(points, key, cells) -> np.ndarray:
    """Signed volumes by the divergence theorem over the outward faces (quads bilinear), so a
    face shared by two cells contributes exactly opposite amounts and the volumes of a conforming
    mesh add up to the volume enclosed by its boundary faces."""
    total = np.zeros(len(cells))
    for face in CELL_FACES[key]:
        p = points[cells[:, list(face)]]
        total += np.einsum("ij,ij->i", p.mean(axis=1), _area_vector(p))
    return total / 3.0


def total_volume(mesh) -> float:
    return float(sum(cell_volumes(mesh["points"], k, c).sum() for k, c, _ in volume_cells(mesh)))


# ------------------------------------------------------------------------------ evaluator
def _tet_volume(p):
    return (
        np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
        / 6.0
    )


def _planar(points) -> bool:
    centred = points - points.mean(axis=0)
    size = float(np.ptp(points, axis=0).max())
    if size == 0.0:
        return False
    return float(np.linalg.svd(centred, compute_uv=False)[-1]) <= 1e-9 * size * math.sqrt(
        len(points)
    )


def _boundary_edges(cells):
    """Directed boundary edges of a set of oriented cells (edges whose reverse is not used)."""
    directed = set()
    for cell in cells:
        n = len(cell)
        for k in range(n):
            directed.add((int(cell[k]), int(cell[(k + 1) % n])))
    return sorted(e for e in directed if (e[1], e[0]) not in directed)


class VolumeEvaluator:
    """Volume meshes, cached by a key over the generator settings and the face meshes."""

    CACHE_SIZE = 8

    def __init__(self, session):
        self.session = session
        self._cache: OrderedDict[str, dict] = OrderedDict()

    def mesh(self, entity_id: int) -> dict:
        """{"points", "tets", "pyramids", "prisms", "hexes", "n_surface", "faces": {domain:
        boundary cells}, "shells", "enclosed_volume", "inflation", "warnings"} (points: the
        merged surface nodes first, unchanged)."""
        session = self.session
        store = session.store
        entity = store.get(entity_id)
        generator = entity.params["generator"]
        shell = assemble_shell(session, entity.params["faces"], entity.name)
        digest = hashlib.sha256(repr(generator).encode())
        digest.update(shell["points"].tobytes())
        for cell in shell["cells"]:
            digest.update(np.asarray(cell, dtype=np.int64).tobytes())
        if generator["inflation"] is not None:
            # Sliding faces are re-meshed with their own settings; who else uses them matters.
            for f in entity.params["faces"]:
                digest.update(repr(store.get(f).params.get("generator")).encode())
            digest.update(repr(sorted(self._shared_faces(entity))).encode())
        key = digest.hexdigest()
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        points = shell["points"]
        n_shell = len(points)
        # Faces that are also faces of structured blocks keep their quads (pyramids); with "auto"
        # other quads are split into triangles.
        block_faces = set()
        for block in store.entities("block"):
            block_faces.update(session.blocks.evaluate(block.id)["faces"].values())
        faces: dict[int, list] = {}
        layer = {key: np.zeros((0, n), dtype=np.int64) for key, n, _ in CELL_TYPES}
        info = None
        warnings: list[str] = []
        if generator["inflation"] is None:
            core = [(c, f, False) for c, f in zip(shell["cells"], shell["cell_face"], strict=True)]
        else:
            points, layer, core, faces, info = self._inflate(entity, shell, generator)
            warnings.extend(info["warnings"])
        # The surface of the region to fill with tetrahedra: boundary cells (face: their domain)
        # and the final inflation front (face None; its quads are faces of layer cells, so they
        # always get pyramids).
        triangles = []
        pyramids = []
        apexes = []
        n_points = len(points)
        for cell, face, front in core:
            cell = np.asarray(cell, dtype=np.int64)
            if len(cell) == 3:
                triangles.append(cell)
                if face is not None:
                    faces.setdefault(face, []).append(cell)
                continue
            mode = generator["pyramids"]
            if not front and (mode == "never" or (mode == "auto" and face not in block_faces)):
                for half in (cell[[0, 1, 2]], cell[[0, 2, 3]]):
                    triangles.append(half)
                    faces.setdefault(face, []).append(half)
                continue
            if face is not None:
                faces.setdefault(face, []).append(cell)
            p = points[cell]
            normal = np.cross(p[2] - p[0], p[3] - p[1])  # outward of the region
            normal = normal / np.linalg.norm(normal)
            size = np.mean(np.linalg.norm(p - np.roll(p, -1, axis=0), axis=1))
            # Boundary quads: apex height 0.3 h, so the side faces rise at ~31 degrees and
            # pyramids on two faces meeting at a right angle do not touch (at 0.5 h their sides
            # would coincide). Inflation fronts are smooth: the apex goes to FRONT_APEX h.
            ratio = 0.3
            if front:
                ratio = FRONT_APEX
            apex = p.mean(axis=0) - ratio * size * normal
            index = n_points + len(apexes)
            apexes.append(apex)
            pyramids.append([cell[0], cell[3], cell[2], cell[1], index])
            for k in range(4):
                # Side triangle of the pyramid, oriented out of the remaining region (that is,
                # into the pyramid): the face edge in its own direction, then the apex.
                triangles.append([cell[k], cell[(k + 1) % 4], index])
        if apexes:
            points = np.concatenate([points, np.asarray(apexes)])
        triangles = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
        # The backend gets only the points of the core surface.
        used = np.unique(triangles)
        local = np.full(len(points), -1, dtype=np.int64)
        local[used] = np.arange(len(used))
        coords, tets = backend(generator["backend"]).mesh(
            points[used], local[triangles], generator["max_size"], generator["grading"]
        )
        if len(coords) < len(used) or not np.array_equal(coords[: len(used)], points[used]):
            raise ValidationError(
                [Issue(entity.name, "the backend moved or reordered boundary points")]
            )
        to_global = np.concatenate([used, len(points) + np.arange(len(coords) - len(used))])
        tets = to_global[tets]
        points = np.concatenate([points, coords[len(used) :]])
        volumes = _tet_volume(points[tets])
        flip = volumes < 0
        tets[flip] = tets[flip][:, [0, 2, 1, 3]]
        if np.any(np.abs(volumes) <= 1e-14 * shell["volume"] / max(len(tets), 1)):
            warnings.append("some tetrahedra are degenerate")
        cells = {
            "tets": np.concatenate([layer["tets"], tets]),
            "pyramids": np.concatenate(
                [layer["pyramids"], np.asarray(pyramids, dtype=np.int64).reshape(-1, 5)]
            ),
            "prisms": layer["prisms"],
            "hexes": layer["hexes"],
        }
        # Drop points no cell uses (interior nodes of re-meshed sliding faces), keeping the order.
        keep = np.zeros(len(points), dtype=bool)
        for array in cells.values():
            keep[array.reshape(-1)] = True
        renumber = np.cumsum(keep) - 1
        points = points[keep]
        for k in cells:
            cells[k] = renumber[cells[k]]
        faces = {
            f: [renumber[np.asarray(c, dtype=np.int64)] for c in cs] for f, cs in faces.items()
        }
        if info is not None:
            # Per-node inflation data in the final numbering (non-wall points: 0 layers).
            n_info = len(info["layers"])
            kept = np.flatnonzero(keep[:n_info])
            layers = np.zeros(len(points), dtype=np.int64)
            height = np.zeros(len(points))
            reason = [""] * len(points)
            top = np.arange(len(points))
            layers[renumber[kept]] = info["layers"][kept]
            height[renumber[kept]] = info["height"][kept]
            top[renumber[kept]] = renumber[info["top"][kept]]
            for i in kept:
                reason[renumber[i]] = info["reason"][i]
            info = {**info, "layers": layers, "height": height, "top": top, "reason": reason}
        result = {
            "points": points,
            **cells,
            "n_surface": int(keep[:n_shell].sum()),
            "faces": faces,
            "shells": shell["shells"],
            "enclosed_volume": shell["volume"],
            "inflation": info,
            "warnings": warnings,
        }
        filled = total_volume(result)
        if abs(filled - shell["volume"]) > 1e-9 * abs(shell["volume"]):
            warnings.append(
                f"the cells fill {filled:.12g} of the enclosed {shell['volume']:.12g} "
                "(overlapping or missing cells)"
            )
        self._cache[key] = result
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)
        return result

    def _shared_faces(self, entity) -> set:
        """Faces of this volume also used by another volume or a block."""
        store = self.session.store
        shared = set()
        for other in store.entities("volume"):
            if other.id != entity.id:
                shared.update(other.params["faces"])
        for block in store.entities("block"):
            shared.update(self.session.blocks.evaluate(block.id)["faces"].values())
        return shared & set(entity.params["faces"])

    def _inflate(self, entity, shell, generator):
        """Inflation layers from the wall faces (opengrid._core.inflate_3d), the sliding faces
        re-meshed around them. Returns (points, layer cells, core surface, boundary faces,
        info)."""
        session = self.session
        store = session.store
        inflation = generator["inflation"]
        points = shell["points"]
        cells = shell["cells"]
        cell_face = shell["cell_face"]
        walls = set(inflation["faces"])
        by_face: dict[int, list[int]] = {}
        for c, f in enumerate(cell_face):
            by_face.setdefault(f, []).append(c)
        # Sliding planes: planar unstructured faces used by this volume only.
        sliding = set()
        if inflation["slide"]:
            shared = self._shared_faces(entity)
            for d in entity.params["faces"]:
                e = store.get(d)
                if d in walls or d in shared or e.kind != "unstructured_domain":
                    continue
                if e.params.get("surface") is not None:
                    continue
                ids = np.unique(np.concatenate([cells[c] for c in by_face[d]]))
                if _planar(points[ids]):
                    sliding.add(d)
        node_group = np.full(len(points), -1, dtype=np.int64)
        for d in sorted(sliding):
            ids = np.unique(np.concatenate([cells[c] for c in by_face[d]]))
            clash = node_group[ids] >= 0
            node_group[ids[clash]] = -2
            node_group[ids[~clash]] = d
        node_group[node_group == -2] = -1
        try:
            res = _core.inflate_3d(
                points,
                [[int(v) for v in c] for c in cells],
                [f in walls for f in cell_face],
                [int(f) for f in cell_face],
                node_group.tolist(),
                inflation["first_height"],
                inflation["growth_rate"],
                inflation["layers"],
                inflation["total_thickness"],
                inflation["transition"],
                generator["max_size"],
                **INFLATION_INTERNALS,
            )
        except ValueError as error:
            raise ValidationError([Issue(entity.name, str(error))]) from None
        new_points = [res["points"]]
        n_points = len(res["points"])
        top = np.asarray(res["top"], dtype=np.int64)
        layer = {
            "tets": res["tets"],
            "pyramids": res["pyramids"],
            "prisms": res["prisms"],
            "hexes": res["hexes"],
        }
        sides: dict[int, list] = {}
        for side, c in zip(res["sides"], res["side_face"], strict=True):
            sides.setdefault(cell_face[c], []).append(np.asarray(side, dtype=np.int64))
        faces: dict[int, list] = {}
        core = []
        for c, f in zip(cells, cell_face, strict=True):
            if f in walls:
                faces.setdefault(f, []).append(c)
            elif f not in sides:
                core.append((c, f, False))
        for poly in res["front"]:
            core.append((np.asarray(poly, dtype=np.int64), None, True))
        # Sliding faces with layers on them: the region outside the layers is re-meshed, its
        # boundary nodes along the walls replaced by the column tops.
        for d in sorted(sides):
            own = [cells[c] for c in by_face[d]]
            edges = [(int(top[a]), int(top[b])) for a, b in _boundary_edges(own)]
            nodes = sorted({v for e in edges for v in e})
            index = {v: k for k, v in enumerate(nodes)}
            loop_xyz = np.concatenate(new_points)[nodes]
            segments = np.asarray([(index[a], index[b]) for a, b in edges], dtype=np.int64)
            xyz, tris = session.unstructured.remesh_planar(d, loop_xyz, segments)
            ids = np.concatenate(
                [np.asarray(nodes, dtype=np.int64), n_points + np.arange(len(xyz) - len(nodes))]
            )
            new_points.append(xyz[len(nodes) :])
            n_points += len(xyz) - len(nodes)
            tris = ids[np.asarray(tris, dtype=np.int64).reshape(-1, 3)]
            allp = np.concatenate(new_points)
            outward = sum(np.cross(allp[c[1]] - allp[c[0]], allp[c[2]] - allp[c[0]]) for c in own)
            if len(tris):
                p = allp[tris[0]]
                if np.dot(np.cross(p[1] - p[0], p[2] - p[0]), outward) < 0:
                    tris = tris[:, ::-1].copy()
            # The layer sides are boundary cells here; the triangles are added with the core.
            faces[d] = list(sides[d])
            for t in tris:
                core.append((t, d, False))
        stops = dict(res["stops"])
        info = {
            "faces": tuple(sorted(walls)),
            "sliding": tuple(sorted(sides)),
            "max_layers": int(res["max_layers"]),
            "top": top,
            "height": np.asarray(res["height"]),
            "layers": np.asarray(res["layers"], dtype=np.int64),
            "reason": list(res["reason"]),
            "cells": {key: len(array) for key, array in layer.items()},
            "stops": stops,
            "warnings": list(res["warnings"]),
        }
        return np.concatenate(new_points), layer, core, faces, info


# ------------------------------------------------------------------------------- quality
def volume_quality(mesh) -> dict:
    """Per-cell metrics of a volume mesh (tetrahedra, pyramids, prisms, hexahedra in that order),
    arrays of shape (1, m):
      * scaled_jacobian: the smallest corner det(e1, e2, e3) / (|e1| |e2| |e3|), scaled to 1 for
        the ideal cell (regular tetrahedron, square pyramid with equilateral sides, right prism
        on an equilateral triangle, cube) and capped at 1; negative for inverted cells;
      * aspect_ratio: tetrahedra circumradius / (3 inradius), other cells longest / shortest edge;
      * min_angle / max_angle: dihedral angles between adjacent faces (degrees);
      * skewness: equiangle skewness of those angles against the type's ideal angles;
      * non_orthogonality: largest angle between a face's area vector and the link of the two
        cell centroids (boundary faces: face centre);
      * volume_ratio: largest volume ratio with a face neighbour; volume."""
    pts = mesh["points"]
    blocks: dict[str, list] = {
        k: []
        for k in ("scaled_jacobian", "aspect_ratio", "min_angle", "max_angle", "skewness", "volume")
    }
    all_faces = {3: [], 4: []}
    all_owner = {3: [], 4: []}
    centroids = []
    offset = 0
    for key, cells, _vtk in volume_cells(mesh):
        m = len(cells)
        if m == 0:
            continue
        p = pts[cells]
        vol = cell_volumes(pts, key, cells)
        corners_table, factor = CELL_CORNERS[key]
        corner = []
        for o, a, b, c in corners_table:
            e1 = p[:, a] - p[:, o]
            e2 = p[:, b] - p[:, o]
            e3 = p[:, c] - p[:, o]
            det = np.einsum("ij,ij->i", e1, np.cross(e2, e3))
            lengths = (
                np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
            )
            # Scaled so the ideal cell's corners are 1; capped at 1 (taller pyramid / prism corners
            # than the ideal one are not better).
            corner.append(np.minimum(det / np.maximum(lengths, 1e-300) * factor, 1.0))
        corner = np.stack(corner, axis=1)
        if key == "tets":
            sj = np.min(np.abs(corner), axis=1) * np.sign(vol)
            aspect = _tet_aspect(p, vol)
        else:
            sj = corner.min(axis=1)
            edges = sorted({tuple(sorted(e)) for e in _cell_edges(key)})
            lengths = np.stack(
                [np.linalg.norm(p[:, b] - p[:, a], axis=1) for a, b in edges], axis=1
            )
            aspect = lengths.max(axis=1) / np.maximum(lengths.min(axis=1), 1e-300)
        normals = []
        for face in CELL_FACES[key]:
            normals.append(_area_vector(p[:, list(face)]))
        dihedral = []
        for fa, fb in CELL_EDGE_FACES[key]:
            na = normals[fa]
            nb = normals[fb]
            cos = np.einsum("ij,ij->i", na, nb) / np.maximum(
                np.linalg.norm(na, axis=1) * np.linalg.norm(nb, axis=1), 1e-300
            )
            dihedral.append(180.0 - np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        dihedral = np.stack(dihedral, axis=1)
        lo = dihedral.min(axis=1)
        hi = dihedral.max(axis=1)
        ideal_lo, ideal_hi = CELL_IDEAL_ANGLES[key]
        skew = np.maximum((hi - ideal_hi) / (180.0 - ideal_hi), (ideal_lo - lo) / ideal_lo)
        blocks["scaled_jacobian"].append(sj)
        blocks["aspect_ratio"].append(aspect)
        blocks["min_angle"].append(lo)
        blocks["max_angle"].append(hi)
        blocks["skewness"].append(np.maximum(skew, 0.0))
        blocks["volume"].append(vol)
        centroids.append(p.mean(axis=1))
        for face in CELL_FACES[key]:
            all_faces[len(face)].append(cells[:, list(face)])
            all_owner[len(face)].append(offset + np.arange(m))
        offset += m
    out = {key: np.concatenate([*parts, np.zeros(0)]) for key, parts in blocks.items()}
    centroid = np.concatenate([*centroids, np.zeros((0, 3))])
    n_cells = offset
    non_orth = np.zeros(n_cells)
    ratio = np.ones(n_cells)
    volume = np.abs(out["volume"])
    for size in (3, 4):
        if not all_faces[size]:
            continue
        faces = np.concatenate(all_faces[size])
        owner = np.concatenate(all_owner[size]).astype(np.int64)
        fp = pts[faces]
        area = _area_vector(fp)
        keys = np.sort(faces, axis=1)
        order = np.lexsort(keys.T[::-1])
        sk = keys[order]
        same = np.all(sk[1:] == sk[:-1], axis=1)
        other = np.full(len(faces), -1, dtype=np.int64)
        first = order[:-1][same]
        second = order[1:][same]
        other[first] = owner[second]
        other[second] = owner[first]
        d = fp.mean(axis=1) - centroid[owner]
        paired = other >= 0
        d[paired] = centroid[other[paired]] - centroid[owner[paired]]
        cosine = np.abs(np.einsum("ij,ij->i", area, d)) / np.maximum(
            np.linalg.norm(area, axis=1) * np.linalg.norm(d, axis=1), 1e-300
        )
        np.maximum.at(non_orth, owner, np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))
        va = volume[owner[paired]]
        vb = volume[other[paired]]
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.maximum(va / vb, vb / va)
        np.maximum.at(ratio, owner[paired], np.nan_to_num(r, nan=1.0, posinf=1.0))
    out["volume_ratio"] = ratio
    out["non_orthogonality"] = non_orth
    return {key: value.reshape(1, -1) for key, value in out.items()}


def _cell_edges(key):
    return [
        (face[k], face[(k + 1) % len(face)]) for face in CELL_FACES[key] for k in range(len(face))
    ]


def _tet_aspect(p, vol):
    """Circumradius / (3 inradius) of tetrahedra (1 for the regular one)."""
    faces = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    normals = [np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a]) for a, b, c in faces]
    area_sum = sum(0.5 * np.linalg.norm(n, axis=1) for n in normals)
    inradius = 3.0 * np.abs(vol) / np.maximum(area_sum, 1e-300)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    num = (
        np.einsum("ij,ij->i", c, c)[:, None] * np.cross(a, b)
        + np.einsum("ij,ij->i", a, a)[:, None] * np.cross(b, c)
        + np.einsum("ij,ij->i", b, b)[:, None] * np.cross(c, a)
    )
    circum = np.linalg.norm(num, axis=1) / np.maximum(12.0 * np.abs(vol), 1e-300)
    return circum / np.maximum(3.0 * inradius, 1e-300)


# ------------------------------------------------------------------------- 3D topology
def analyze_3d(session) -> dict:
    """Boundary / interface faces of a 3D model (blocks and volumes) and their patches.

    Returns {"entries": [{"label", "block": (block id, face) or None, "volume": (volume id,
    domain id) or None, "keys": set, "patch": name or None}], "interfaces": [labels],
    "issues": [Issue]}. Keys are ("domain", id) and ("block", id, face): what patch members name.
    """
    from opengrid.model.blocks import Topology3D, analyze_blocks

    store = session.store
    issues: list[Issue] = []
    if store.entities("block"):
        topo = analyze_blocks(store, session.blocks)
    else:
        topo = Topology3D()
    issues.extend(topo.problems)
    volume_uses: dict[int, list[int]] = {}
    for volume in store.entities("volume"):
        for d in volume.params["faces"]:
            volume_uses.setdefault(d, []).append(volume.id)
    entries = []
    interfaces = []

    def name(entity_id):
        return store.get(entity_id).name

    for a, b in topo.interfaces:
        interfaces.append(f"{name(a.block)}.{a.face} | {name(b.block)}.{b.face}")
    block_domains = set()
    for f in topo.boundary:
        label = f"{name(f.block)}.{f.face}"
        if f.domain is not None:
            block_domains.add(f.domain)
        if f.domain is not None and f.domain in volume_uses:
            vols = ", ".join(name(v) for v in volume_uses[f.domain])
            interfaces.append(f"{label} | {vols} ({name(f.domain)})")
            continue
        keys = {("block", f.block, f.face)}
        if f.domain is not None:
            keys.add(("domain", f.domain))
        entries.append({"label": label, "block": (f.block, f.face), "volume": None, "keys": keys})
    for d, vols in sorted(volume_uses.items()):
        count = len(vols) + int(d in block_domains)
        if count == 1:
            entries.append(
                {
                    "label": f"{name(vols[0])}.{name(d)}",
                    "block": None,
                    "volume": (vols[0], d),
                    "keys": {("domain", d)},
                }
            )
        elif count == 2 and d not in block_domains:
            interfaces.append(f"{name(vols[0])} | {name(vols[1])} ({name(d)})")
        elif count > 2:
            issues.append(Issue(name(d), f"bounds {count} blocks / volumes (at most 2)"))
    # Patches.
    by_key: dict[tuple, list[int]] = {}
    for index, entry in enumerate(entries):
        entry["patch"] = None
        for key in entry["keys"]:
            by_key.setdefault(key, []).append(index)
    assigned: dict[int, list[str]] = {}
    for patch in store.entities("patch"):
        for item in patch.params["members"]:
            if "connector" in item:
                issues.append(
                    Issue(
                        name(item["connector"]),
                        f"is a connector in patch {patch.name}; 3D patches contain domains or "
                        "block faces",
                    )
                )
                continue
            if "block" in item:
                key = ("block", item["block"], item["face"])
                label = f"{name(item['block'])}.{item['face']}"
            else:
                key = ("domain", item["domain"])
                label = name(item["domain"])
            hits = by_key.get(key, [])
            if not hits:
                issues.append(Issue(label, f"is in patch {patch.name} but is not a boundary face"))
            for index in hits:
                assigned.setdefault(index, []).append(patch.name)
        if patch.params["periodic"] is not None:
            issues.append(Issue(patch.name, "periodic pairing of 3D faces is not supported yet"))
    for index, entry in enumerate(entries):
        names = sorted(set(assigned.get(index, [])))
        if not names:
            issues.append(Issue(entry["label"], "is a boundary face in no patch"))
        elif len(names) > 1:
            issues.append(Issue(entry["label"], f"belongs to several patches: {', '.join(names)}"))
        else:
            entry["patch"] = names[0]
    return {"entries": entries, "interfaces": interfaces, "issues": issues, "blocks": topo}
