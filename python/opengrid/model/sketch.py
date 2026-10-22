"""Sketch topology (a floating wireframe of blocks): a coarse block lattice drawn loosely
around the geometry, from which ordinary connectors, domains, blocks and a smoothing group are
generated and kept up to date.

A ``sketch`` holds vertices (points) and blocks: 4 vertex indices (2D, a quadrilateral in a
z = const plane) or 8 (3D, a hexahedron: bottom face 0-1-2-3, top face 4-5-6-7 above it, as
in VTK). Its edges and faces are those of the blocks, shared where blocks share vertices.

* **Dimensions** are set on any edge of an edge family (the edges joined by being opposite in
  some block) and hold for the whole family; a family without one takes ``spacing`` (the
  longest edge of the family divided by it). Two different values in one family are reported
  with both edges.
* **Geometry:** edges are assigned to database curves (``curves``) and 3D faces to surfaces
  (``surfaces``). A 2D edge can instead adopt an existing connector (``connectors``, one made
  by hand, not by a sketch): the edge's vertices go to the connector's ends, the domain uses
  that connector itself (so a domain on its other side shares its nodes), its edge family
  takes the connector's point count and the group keeps the connector fixed. A sketch drawn
  by snapping onto connector ends adopts them (the GUI tool; ``matching_connectors``), and
  edges lying on a connector are offered for adoption like curves. A vertex on edges
  assigned to two or more curves goes to the point where those curves meet; on one curve, to
  its projection onto it; on assigned faces only (3D), to its projection onto the surface.
  Everything else stays where it was drawn. Unassigned edges
  and faces lying near geometry are offered for assignment (:func:`suggestions`, also listed
  as warnings of a generation) and never snapped silently.
* **Unstructured blocks** (``unstructured``: block indices, 2D): those blocks become
  unstructured domains on the same connectors (``generator``: their triangulation settings,
  shared by all of them; a wall on one of their edges becomes inflation layers from that
  connector with the wall's first-cell height, layers and growth, ``transition`` (stop at the
  core cell size) and ``slide`` (run into the neighbouring edges)). An edge leaving such a
  wall is clustered at the wall as for structured blocks, except with ``slide`` when it
  borders only unstructured blocks: then it gets points at the layer heights
  (model/layer_spacing.py; between 45 and 135 degrees to the wall; not a periodic one). They
  are not members of the smoothing group; a connector between a structured and an
  unstructured block is pinned so both sides keep the same nodes. The vertex passes treat
  every block as structured.
* **Generation:** an edge assigned to a curve becomes a connector on the part of the curve
  between its vertices; any other edge a straight line (in 3D, between two faces on the same
  surface: a spline through the line's points projected onto it) and a connector on it. Each
  2D block becomes a domain, each 3D face a domain (on its surface when assigned) and each 3D
  block a six-face block; all members of one smoothing group (10a/10b) with the ``walls``
  (edges in 2D, faces in 3D, with a first-cell height; a wall between two blocks is an internal
  wall, e.g. a wake cut) and every other boundary sliding. With ``relax``, the sliding
  vertices - 2D: joining two edges on one curve; 3D: on boundary faces that all lie on one
  surface or plane - are relaxed in passes: the group is smoothed on a grid with half the
  intervals, with the walls sliding and internal walls pinned, the vertices are moved to
  where it put them (Anderson mixing, projected back onto their geometry) and this repeats
  until they move by less than 1e-3 of the local spacing; interior vertices float in the
  group itself. The final group keeps the walls. So the grid depends on the topology and the
  geometry, not on where the floating vertices were drawn. Vertices on assigned curves in
  3D, at curve ends, on pinned edges or internal walls, and on unassigned boundaries in 2D
  stay where drawn (projected).
* **Extrusion** (2D, ``extrude``: ``{"vector", "points", "distribution"}``): every block's
  domain is also extruded along the vector (a translated block of ``points`` layers), so a
  straight 3D grid is the smoothed 2D grid extruded - the same in every section.
* **High-order meshes** (``interval_multiple`` m): every edge family's interval count (and the
  extrusion's) is rounded up to a multiple of m, so m x m (x m) cells form the elements of an
  order-m mesh (``Session.high_order``).
* **Ownership:** the generated entities are listed in ``outputs`` (role key, entity id) and
  are regenerated in place (same names and ids) when the sketch changes; editing, renaming or
  deleting them directly is refused (:func:`check_ownership`, run by every command), and
  deleting the sketch deletes them.

Method: the topology is a floating wireframe (blocks drawn freely, corners projected onto the
geometry, the grid relaxed with its block corners and interfaces free); the smoothing is that
of opengrid.model.groups. References: P. R. Eiseman, "Multiblock grid generation with automatic zoning", NASA CP-3291 (1995); J. F. Thompson, B. K. Soni, N. P. Weatherill (eds.), Handbook of
Grid Generation, CRC Press (1999), ch. 13 (multiblock topology).
"""

import logging
import math
from dataclasses import dataclass, field

import numpy as np

from opengrid.model.connectors import CHAIN_TOLERANCE
from opengrid.model.groups import _count, _growth
from opengrid.validation import Issue, ValidationError

KIND = "sketch"
log = logging.getLogger("opengrid")
DEFAULT_SKETCH = {
    "vertices": (),
    "blocks": (),
    "dimensions": (),
    "distributions": (),
    "spacing": None,
    "curves": (),
    "connectors": (),
    "surfaces": (),
    "unstructured": (),
    "generator": None,
    "walls": (),
    "fixed": (),
    "pinned": (),
    "symmetry": (),
    "relax": True,
    "slide_limit": 0.5,
    "tolerance": 1e-8,
    "max_iterations": 4000,
    "extrude": None,
    "interval_multiple": 1,
    "outputs": (),
}
DEFAULT_SKETCH_WALL = {
    "spacing": None,
    "decay": 0.1,
    "relaxation": 0.1,
    "max_source": 100.0,
    "layers": 0,
    "growth": 1.2,
    "transition": True,  # unstructured blocks: layers stop at the core cell size
    "slide": False,  # unstructured blocks: layers run into the neighbouring edges
}
WALL_KEYS = ("spacing", "decay", "relaxation", "max_source", "layers", "growth")
# Hexahedron (VTK order): edges and faces by local vertex indices.
HEX_EDGES = (
    (0, 1),
    (3, 2),
    (4, 5),
    (7, 6),
    (0, 3),
    (1, 2),
    (4, 7),
    (5, 6),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
)
HEX_FAMILIES = ((0, 1, 2, 3), (4, 5, 6, 7), (8, 9, 10, 11))  # i, j, k edges of one block
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
MAX_PASSES = 60
PASS_TOLERANCE = 1e-3  # largest vertex move per pass, in local spacings, to stop the passes
PASS_SMOOTHING_TOLERANCE = 1e-8
PASS_MAX_ITERATIONS = 1500
PASS_SLIDE_LIMIT = 3.0  # intervals a sliding node may move per pass
SUGGEST_FRACTION = 0.05  # suggestions: geometry within this fraction of the edge length
SURFACE_SPLINE_POINTS = 17


# ----------------------------------------------------------------------------- parameters
def _int(value, label, issues):
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        value = value.item()
    if isinstance(value, bool) or not isinstance(value, int):
        issues.append(Issue(label, f"expected an integer, got {value!r}"))
        return None
    return value


def _number(value, label, issues, positive=False, optional=False):
    if value is None and optional:
        return None
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        value = value.item()
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        issues.append(Issue(label, f"expected a finite number, got {value!r}"))
        return None
    if positive and value <= 0:
        issues.append(Issue(label, f"must be positive, got {value!r}"))
        return None
    return float(value)


def _vertex_tuple(raw, size, count, label, issues):
    """A tuple of ``size`` vertex indices (or None)."""
    if isinstance(raw, (str, bytes)) or not hasattr(raw, "__len__") or len(raw) != size:
        issues.append(Issue(label, f"expected {size} vertex indices, got {raw!r}"))
        return None
    out = []
    for item in raw:
        index = _int(item, label, issues)
        if index is None:
            return None
        if not 0 <= index < count:
            issues.append(Issue(label, f"vertex {index} does not exist ({count} vertices)"))
            return None
        out.append(index)
    return tuple(out)


def edge_key(a: int, b: int) -> tuple[int, int]:
    if a < b:
        return (a, b)
    return (b, a)


def face_key(vertices) -> tuple[int, ...]:
    return tuple(sorted(vertices))


def _items(raw, label, issues):
    """Assignment lists may be given as a dict {edge: value} or a list of dicts."""
    if isinstance(raw, dict):
        return [(key, value) for key, value in raw.items()]
    if isinstance(raw, (str, bytes)) or not hasattr(raw, "__iter__"):
        issues.append(Issue(label, f"expected a dict or a list, got {raw!r}"))
        return []
    return [(None, item) for item in raw]


def _entity_ref(store, raw, kinds, label, issues):
    if isinstance(raw, bool) or not isinstance(raw, (str, int)) or raw not in store:
        issues.append(Issue(label, f"no entity {raw!r}"))
        return None
    entity = store.get(raw)
    if entity.kind not in kinds:
        issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}"))
        return None
    return entity.id


def normalize_sketch(raw: dict, store, where: str) -> dict:
    """Canonical sketch parameters; raises ValidationError listing every problem (unknown
    vertices, malformed blocks, edges and faces that are not in the lattice, conflicting
    dimensions)."""
    from opengrid.model.connectors import CURVE_KINDS, normalize_distribution
    from opengrid.model.entities import SURFACE_KINDS, _Checker

    issues: list[Issue] = []
    for key in sorted(set(raw) - set(DEFAULT_SKETCH)):
        issues.append(Issue(f"{where}: {key}", "unexpected key"))
    merged = dict(DEFAULT_SKETCH)
    merged.update(raw)
    checker = _Checker(where)
    vertices = checker.points(merged["vertices"], "vertices", 3)
    if vertices is None:
        raise ValidationError([*issues, *checker.issues])
    count = len(vertices)
    blocks = []
    size = None
    for index, item in enumerate(merged["blocks"]):
        label = f"{where}: blocks[{index}]"
        n = None
        if hasattr(item, "__len__"):
            n = len(item)
        if n not in (4, 8):
            issues.append(Issue(label, f"expected 4 (2D) or 8 (3D) vertex indices, got {item!r}"))
            continue
        if size is None:
            size = n
        if n != size:
            issues.append(Issue(label, "all blocks must be 2D (4 vertices) or all 3D (8)"))
            continue
        block = _vertex_tuple(item, n, count, label, issues)
        if block is None:
            continue
        if len(set(block)) != n:
            issues.append(Issue(label, f"repeats a vertex: {block}"))
            continue
        blocks.append(block)
    if not blocks and not issues:
        issues.append(Issue(f"{where}: blocks", "at least one block is required"))
    if issues or checker.issues:
        raise ValidationError([*issues, *checker.issues])
    lattice = Lattice.build(vertices, tuple(blocks), where)
    issues.extend(lattice.issues)

    def edge_of(value, label):
        edge = _vertex_tuple(value, 2, count, label, issues)
        if edge is None:
            return None
        if edge_key(*edge) not in lattice.edges:
            issues.append(Issue(label, f"{edge} is not an edge of the sketch"))
            return None
        return edge

    def face_of(value, label):
        face = _vertex_tuple(value, 4, count, label, issues)
        if face is None:
            return None
        if face_key(face) not in lattice.faces:
            issues.append(Issue(label, f"{face} is not a face of the sketch"))
            return None
        return face_key(face)

    dimensions = {}
    for index, (key, item) in enumerate(_items(merged["dimensions"], f"{where}: dimensions", [])):
        label = f"{where}: dimensions[{index}]"
        if key is not None:
            item = {"edge": key, "points": item}
        if not isinstance(item, dict) or set(item) != {"edge", "points"}:
            issues.append(Issue(label, "expected {'edge': (a, b), 'points': n}"))
            continue
        edge = edge_of(item["edge"], label)
        points = _int(item["points"], label, issues)
        if edge is None or points is None:
            continue
        if points < 2:
            issues.append(Issue(label, f"points must be at least 2, got {points}"))
            continue
        dimensions[edge_key(*edge)] = points
    distributions = {}
    for index, (key, item) in enumerate(_items(merged["distributions"], where, issues)):
        label = f"{where}: distributions[{index}]"
        if key is not None:
            item = {"edge": key, "distribution": item}
        if not isinstance(item, dict) or set(item) != {"edge", "distribution"}:
            issues.append(Issue(label, "expected {'edge': (a, b), 'distribution': {...}}"))
            continue
        edge = edge_of(item["edge"], label)
        try:
            distribution = normalize_distribution(item["distribution"], label)
        except ValidationError as error:
            issues.extend(error.issues)
            continue
        if edge is None:
            continue
        distributions[edge_key(*edge)] = (edge, distribution)
    curves = {}
    for index, (key, item) in enumerate(_items(merged["curves"], where, issues)):
        label = f"{where}: curves[{index}]"
        if key is not None:
            item = {"edge": key, "curve": item}
        if not isinstance(item, dict) or set(item) != {"edge", "curve"}:
            issues.append(Issue(label, "expected {'edge': (a, b), 'curve': name}"))
            continue
        edge = edge_of(item["edge"], label)
        curve = _entity_ref(store, item["curve"], CURVE_KINDS, label, issues)
        if edge is not None and curve is not None:
            curves[edge_key(*edge)] = curve
    connectors = {}
    owned = owned_entities(store)
    for index, (key, item) in enumerate(_items(merged["connectors"], where, issues)):
        label = f"{where}: connectors[{index}]"
        if key is not None:
            item = {"edge": key, "connector": item}
        if not isinstance(item, dict) or set(item) != {"edge", "connector"}:
            issues.append(Issue(label, "expected {'edge': (a, b), 'connector': name}"))
            continue
        if lattice.dim == 3:
            issues.append(Issue(label, "connectors are adopted by the edges of 2D sketches"))
            continue
        edge = edge_of(item["edge"], label)
        connector = _entity_ref(store, item["connector"], ("connector",), label, issues)
        if edge is None or connector is None:
            continue
        if connector in owned:
            maker = store.get(owned[connector]).name
            issues.append(
                Issue(
                    label,
                    f"{store.get(connector).name} is generated by sketch {maker!r}: a sketch "
                    "adopts connectors made by hand",
                )
            )
            continue
        key = edge_key(*edge)
        if key in curves:
            issues.append(
                Issue(
                    label, f"edge {key} has a curve assigned: an edge takes a curve or a connector"
                )
            )
            continue
        if key in distributions:
            issues.append(
                Issue(label, f"edge {key} adopts a connector: its distribution is the connector's")
            )
            continue
        connectors[key] = connector
    surfaces = {}
    for index, (key, item) in enumerate(_items(merged["surfaces"], where, issues)):
        label = f"{where}: surfaces[{index}]"
        if key is not None:
            item = {"face": key, "surface": item}
        if not isinstance(item, dict) or set(item) != {"face", "surface"}:
            issues.append(Issue(label, "expected {'face': (a, b, c, d), 'surface': name}"))
            continue
        if lattice.dim == 2:
            issues.append(Issue(label, "2D sketches have no faces to assign"))
            continue
        face = face_of(item["face"], label)
        surface = _entity_ref(store, item["surface"], SURFACE_KINDS, label, issues)
        if face is not None and surface is not None:
            if face not in lattice.boundary_faces:
                issues.append(Issue(label, f"{face} is not a boundary face"))
                continue
            surfaces[face] = surface
    walls = {}
    for index, (key, item) in enumerate(_items(merged["walls"], where, issues)):
        label = f"{where}: walls[{index}]"
        ref_key = "edge"
        if lattice.dim == 3:
            ref_key = "face"
        if key is not None and isinstance(item, dict):
            item = {ref_key: key, **item}  # {edge: {"spacing": h, "layers": n, ...}}
        elif key is not None:
            item = {ref_key: key, "spacing": item}
        if not isinstance(item, dict) or ref_key not in item:
            issues.append(Issue(label, f"expected {{'{ref_key}': ..., 'spacing': h}}"))
            continue
        for extra in sorted(set(item) - set(DEFAULT_SKETCH_WALL) - {ref_key}):
            issues.append(Issue(f"{label}.{extra}", "unexpected key"))
        wall = dict(DEFAULT_SKETCH_WALL)
        wall.update(item)
        if ref_key == "edge":
            ref = edge_of(wall["edge"], label)
            if ref is None:
                continue
            ref = edge_key(*ref)
        else:
            ref = face_of(wall["face"], label)
            if ref is None:
                continue
        walls[ref] = {
            ref_key: ref,
            "spacing": _number(wall["spacing"], f"{label}.spacing", issues, True, True),
            "decay": _number(wall["decay"], f"{label}.decay", issues, True),
            "relaxation": _number(wall["relaxation"], f"{label}.relaxation", issues, True),
            "max_source": _number(wall["max_source"], f"{label}.max_source", issues, True),
            "layers": _count(wall["layers"], f"{label}.layers", issues),
            "growth": _growth(wall["growth"], f"{label}.growth", issues),
            "transition": wall["transition"],
            "slide": wall["slide"],
        }
        for flag in ("transition", "slide"):
            if not isinstance(wall[flag], bool):
                issues.append(
                    Issue(f"{label}.{flag}", f"must be true or false, got {wall[flag]!r}")
                )
    unstructured = set()
    for index, item in enumerate(merged["unstructured"]):
        label = f"{where}: unstructured[{index}]"
        block = _int(item, label, issues)
        if block is None:
            continue
        if lattice.dim == 3:
            issues.append(Issue(label, "unstructured blocks are 2D (an unstructured domain each)"))
            break
        if not 0 <= block < len(blocks):
            issues.append(Issue(label, f"block {block} does not exist"))
            continue
        unstructured.add(block)
    generator = None
    if merged["generator"] is not None:
        from opengrid.model.unstructured import normalize_generator

        try:
            generator = normalize_generator(merged["generator"], where, store)
        except ValidationError as error:
            issues.extend(error.issues)
        else:
            if generator["inflation"] is not None:
                issues.append(
                    Issue(
                        f"{where}: generator.inflation",
                        "the layers of unstructured blocks come from the sketch's walls "
                        "(first-cell height, layers, growth)",
                    )
                )
    fixed = set()
    for index, item in enumerate(merged["fixed"]):
        vertex = _int(item, f"{where}: fixed[{index}]", issues)
        if vertex is not None and not 0 <= vertex < count:
            issues.append(Issue(f"{where}: fixed[{index}]", f"vertex {vertex} does not exist"))
            continue
        if vertex is not None:
            fixed.add(vertex)
    pinned = set()
    for index, item in enumerate(merged["pinned"]):
        label = f"{where}: pinned[{index}]"
        if lattice.dim == 2:
            ref = edge_of(item, label)
            if ref is not None:
                pinned.add(edge_key(*ref))
        else:
            ref = face_of(item, label)
            if ref is not None:
                pinned.add(ref)
    symmetry = set()
    for index, item in enumerate(merged["symmetry"]):
        label = f"{where}: symmetry[{index}]"
        if lattice.dim != 3:
            issues.append(Issue(label, "symmetry planes are faces of 3D sketches"))
            break
        ref = face_of(item, label)
        if ref is None:
            continue
        if ref not in lattice.boundary_faces:
            issues.append(Issue(label, f"{ref} is not a boundary face"))
            continue
        symmetry.add(ref)
    spacing = _number(merged["spacing"], f"{where}: spacing", issues, True, True)
    multiple = merged["interval_multiple"]
    if isinstance(multiple, bool) or not isinstance(multiple, int) or multiple < 1:
        issues.append(
            Issue(f"{where}: interval_multiple", f"expected a whole number >= 1, got {multiple!r}")
        )
        multiple = 1
    relax = merged["relax"]
    if not isinstance(relax, bool):
        issues.append(Issue(f"{where}: relax", f"must be true or false, got {relax!r}"))
    slide_limit = _number(merged["slide_limit"], f"{where}: slide_limit", issues, True)
    tolerance = _number(merged["tolerance"], f"{where}: tolerance", issues, True)
    iterations = _int(merged["max_iterations"], f"{where}: max_iterations", issues)
    if iterations is not None and iterations < 1:
        issues.append(Issue(f"{where}: max_iterations", "must be positive"))
    extrude = _normalize_extrude(merged["extrude"], lattice, where, issues)
    outputs = []
    for index, item in enumerate(merged["outputs"]):
        label = f"{where}: outputs[{index}]"
        if not hasattr(item, "__len__") or len(item) != 2 or not isinstance(item[0], str):
            issues.append(Issue(label, "expected (role, entity id)"))
            continue
        entity_id = _int(item[1], label, issues)
        if entity_id is not None:
            outputs.append((item[0], entity_id))
    if not issues:
        issues.extend(lattice.dimension_conflicts(dimensions))
        issues.extend(_adoption_conflicts(store, lattice, dimensions, connectors, multiple))
    if issues:
        raise ValidationError(issues)
    return {
        "vertices": vertices,
        "blocks": tuple(blocks),
        "dimensions": tuple({"edge": e, "points": n} for e, n in sorted(dimensions.items())),
        "distributions": tuple(
            {"edge": edge, "distribution": d} for _k, (edge, d) in sorted(distributions.items())
        ),
        "spacing": spacing,
        "curves": tuple({"edge": e, "curve": c} for e, c in sorted(curves.items())),
        "connectors": tuple({"edge": e, "connector": c} for e, c in sorted(connectors.items())),
        "surfaces": tuple({"face": f, "surface": s} for f, s in sorted(surfaces.items())),
        "walls": tuple(w for _k, w in sorted(walls.items())),
        "fixed": tuple(sorted(fixed)),
        "pinned": tuple(sorted(pinned)),
        "symmetry": tuple(sorted(symmetry)),
        "unstructured": tuple(sorted(unstructured)),
        "generator": generator,
        "relax": relax,
        "interval_multiple": multiple,
        "slide_limit": slide_limit,
        "tolerance": tolerance,
        "max_iterations": iterations,
        "extrude": extrude,
        "outputs": tuple(sorted(outputs)),
    }


def _adoption_conflicts(store, lattice, dimensions, connectors, multiple: int) -> list[Issue]:
    """Adopted connectors fix their edge families' point counts: two different counts in one
    family (two adopted connectors, or an adopted one and a given dimension), or a count whose
    intervals are not a multiple of ``interval_multiple``, are reported by edge."""
    issues = []
    counts: dict = {}  # family root -> (points, connector name)
    for key in sorted(connectors):
        entity = store.get(connectors[key])
        points = entity.params["dimension"]
        if (points - 1) % multiple:
            issues.append(
                Issue(
                    f"sketch edge {key}",
                    f"adopts {entity.name} with {points} points: its intervals are not a "
                    f"multiple of interval_multiple {multiple}",
                )
            )
        root = lattice.family[key]
        if root in counts and counts[root][0] != points:
            other, name = counts[root]
            issues.append(
                Issue(
                    f"sketch edge {key}",
                    f"adopts {entity.name} with {points} points, but {name} with {other} "
                    "points is in the same edge family",
                )
            )
            continue
        counts.setdefault(root, (points, entity.name))
    for key in sorted(dimensions):
        root = lattice.family[key]
        if root in counts and counts[root][0] != dimensions[key]:
            points, name = counts[root]
            issues.append(
                Issue(
                    f"sketch edge {key}",
                    f"has {dimensions[key]} points, but its edge family adopts {name} with "
                    f"{points} points",
                )
            )
    return issues


def adopted_ends(store, params: dict) -> dict:
    """Edge key -> (start, end) points of the adopted connectors."""
    from opengrid.model.connectors import ConnectorEvaluator

    evaluator = ConnectorEvaluator(store)
    out = {}
    for item in params.get("connectors", ()):
        chain = evaluator.chain(item["connector"])
        out[item["edge"]] = (
            np.asarray(chain.start(), dtype=float),
            np.asarray(chain.end(), dtype=float),
        )
    return out


def _polyline_distance(point, vertices) -> float:
    """Distance from a point to a polyline (its closest segment)."""
    a = np.asarray(vertices[:-1], dtype=float)
    b = np.asarray(vertices[1:], dtype=float)
    ab = b - a
    lengths = np.einsum("ij,ij->i", ab, ab)
    t = np.zeros(len(a))
    nonzero = lengths > 0.0
    t[nonzero] = np.clip(
        np.einsum("ij,ij->i", point - a[nonzero], ab[nonzero]) / lengths[nonzero], 0.0, 1.0
    )
    closest = a + t[:, None] * ab
    return float(np.min(np.linalg.norm(closest - point, axis=1)))


def _connector_candidates(store, only=None) -> list:
    """(id, name, start, end, nodes) of the connectors not generated by a sketch."""
    from opengrid.model.connectors import ConnectorEvaluator

    owned = owned_entities(store)
    evaluator = ConnectorEvaluator(store)
    out = []
    for connector in store.entities("connector"):
        if connector.id in owned or (only is not None and connector.id not in only):
            continue
        chain = evaluator.chain(connector.id)
        nodes = np.asarray(evaluator.nodes(connector.id).points, dtype=float)
        out.append(
            (
                connector.id,
                connector.name,
                np.asarray(chain.start(), dtype=float),
                np.asarray(chain.end(), dtype=float),
                nodes,
            )
        )
    return out


def _connector_under(pa, pb, candidates, end_limit: float, middle_limit: float):
    """The candidate whose ends are the edge's vertices (within end_limit, either way round)
    and whose nodes pass within middle_limit of the edge's midpoint (so a straight edge is
    not matched by an arc with the same ends): (id, name, gap) or None."""
    middle = 0.5 * (np.asarray(pa) + np.asarray(pb))
    best = None
    for connector_id, name, start, end, nodes in candidates:
        straight = max(np.linalg.norm(start - pa), np.linalg.norm(end - pb))
        flipped = max(np.linalg.norm(end - pa), np.linalg.norm(start - pb))
        gap = float(min(straight, flipped))
        if gap > end_limit or _polyline_distance(middle, nodes) > middle_limit:
            continue
        if best is None or gap < best[2]:
            best = (connector_id, name, gap)
    return best


def matching_connectors(store, vertices, edges, tolerance: float | None = None) -> dict:
    """{edge key: connector id} for the edges whose two vertices lie at the ends of a
    connector not generated by a sketch (within ``tolerance``, by default CHAIN_TOLERANCE of
    the model size) and whose middle lies on it (within SUGGEST_FRACTION of the edge length):
    what a sketch drawn by snapping onto connector ends adopts."""
    points = np.asarray(vertices, dtype=float)
    if tolerance is None:
        tolerance = CHAIN_TOLERANCE * max(1.0, float(np.max(np.abs(points))))
    candidates = _connector_candidates(store)
    out = {}
    for edge in edges:
        key = edge_key(*edge)
        pa = points[key[0]]
        pb = points[key[1]]
        length = float(np.linalg.norm(pb - pa))
        found = _connector_under(pa, pb, candidates, tolerance, SUGGEST_FRACTION * length)
        if found is not None:
            out[key] = found[0]
    return out


def regenerate_adopters(store, entity_id: int) -> list[str]:
    """Regenerate the sketches adopting a connector among this entity and its dependents
    (inside the open transaction): their edge families follow the connector's point count
    and their edges its ends. Returns the sketch names."""
    affected = {entity_id}
    pending = [entity_id]
    while pending:
        current = pending.pop()
        for dependent in store.dependents(current):
            if dependent.id not in affected:
                affected.add(dependent.id)
                pending.append(dependent.id)
    names = []
    for sketch in store.entities(KIND):
        adopted = {item["connector"] for item in sketch.params.get("connectors", ())}
        if adopted & affected:
            store.update(sketch.id, params=dict(sketch.params))  # touched: outputs may change
            generate(store, sketch.id, passes=False)
            names.append(sketch.name)
    return names


def _normalize_extrude(raw, lattice, where: str, issues: list) -> dict | None:
    """A 2D sketch's extrusion: {"vector", "points", "distribution"} (spacings as lengths)."""
    from opengrid.model.connectors import normalize_distribution
    from opengrid.model.entities import _Checker

    if raw is None:
        return None
    label = f"{where}: extrude"
    if lattice.dim != 2:
        issues.append(Issue(label, "only 2D sketches are extruded (3D sketches have blocks)"))
        return None
    if not isinstance(raw, dict):
        issues.append(Issue(label, "expected {'vector', 'points', 'distribution'}"))
        return None
    for key in sorted(set(raw) - {"vector", "points", "distribution"}):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    checker = _Checker(label)
    vector = checker.vec3(raw.get("vector", (0.0, 0.0, 1.0)), "vector")
    issues.extend(checker.issues)
    if vector is not None and not any(vector):
        issues.append(Issue(f"{label}.vector", "must not be zero"))
    points = _int(raw.get("points", 11), f"{label}.points", issues)
    if points is not None and points < 2:
        issues.append(Issue(f"{label}.points", f"must be at least 2, got {points}"))
    try:
        distribution = normalize_distribution(raw.get("distribution", "uniform"), label)
    except ValidationError as error:
        issues.extend(error.issues)
        return None
    if distribution["type"] == "explicit":
        issues.append(Issue(f"{label}.distribution", "uniform, geometric, tanh or sine"))
    return {"vector": vector, "points": points, "distribution": distribution}


def _extrusion(params: dict, index: int) -> dict:
    """The extruded block of a 2D sketch's block ``index`` (spacings converted to fractions
    of the extrusion length, as extruded blocks take them)."""
    extrude = params["extrude"]
    length = float(np.linalg.norm(extrude["vector"]))
    distribution = {}
    for key, value in extrude["distribution"].items():
        if key in ("begin_spacing", "end_spacing") and value is not None:
            value = value / length
        distribution[key] = value
    return {
        "domain": f"@d{index}",
        "mode": "translate",
        "vector": extrude["vector"],
        "dimension": _round_points(extrude["points"], params),
        "distribution": distribution,
    }


def sketch_references(params: dict) -> set[int]:
    """The assigned curves, adopted connectors and surfaces (not the outputs: the sketch owns
    those)."""
    refs = {item["curve"] for item in params["curves"]}
    refs.update(item["connector"] for item in params.get("connectors", ()))
    refs.update(item["surface"] for item in params["surfaces"])
    return refs


def sketch_script_params(params: dict, store) -> dict:
    result = dict(params)
    result["curves"] = tuple(
        {"edge": item["edge"], "curve": store.get(item["curve"]).name} for item in params["curves"]
    )
    result["connectors"] = tuple(
        {"edge": item["edge"], "connector": store.get(item["connector"]).name}
        for item in params.get("connectors", ())
    )
    result["surfaces"] = tuple(
        {"face": item["face"], "surface": store.get(item["surface"]).name}
        for item in params["surfaces"]
    )
    return result


# ------------------------------------------------------------------------------- lattice
@dataclass
class Lattice:
    """Edges, faces and edge families of a sketch's blocks."""

    dim: int
    blocks: tuple  # oriented: 2D counter-clockwise, 3D with a positive volume
    edges: dict  # edge key -> [(block, local edge)]
    faces: dict  # 3D: face key -> [(block, local face, ordered vertices)]
    family: dict  # edge key -> family root (edge key)
    boundary_edges: set
    boundary_faces: set
    issues: list = field(default_factory=list)

    @staticmethod
    def build(vertices, blocks, where: str) -> "Lattice":
        points = np.asarray(vertices, dtype=float)
        issues = []
        dim = 2
        if len(blocks[0]) == 8:
            dim = 3
        oriented = []
        for index, block in enumerate(blocks):
            if dim == 2:
                p = points[list(block)]
                area = 0.5 * sum(
                    p[k, 0] * p[(k + 1) % 4, 1] - p[(k + 1) % 4, 0] * p[k, 1] for k in range(4)
                )
                if area == 0.0:
                    issues.append(Issue(f"{where}: blocks[{index}]", "has no area"))
                if area < 0.0:
                    block = (block[0], block[3], block[2], block[1])
            else:
                p = points[list(block)]
                volume = np.dot(np.cross(p[1] - p[0], p[3] - p[0]), p[4] - p[0])
                if volume == 0.0:
                    issues.append(Issue(f"{where}: blocks[{index}]", "has no volume"))
                if volume < 0.0:
                    # Mirrored in i: swap the vertex pairs along the i edges.
                    block = tuple(block[c] for c in (1, 0, 3, 2, 5, 4, 7, 6))
            oriented.append(tuple(block))
        edges: dict = {}
        parent: dict = {}

        def find(e):
            while parent[e] != e:
                parent[e] = parent[parent[e]]
                e = parent[e]
            return e

        def union(a, b):
            ra = find(a)
            rb = find(b)
            if ra < rb:
                parent[rb] = ra
            elif rb < ra:
                parent[ra] = rb

        faces: dict = {}
        for index, block in enumerate(oriented):
            if dim == 2:
                local = [(block[k], block[(k + 1) % 4]) for k in range(4)]
                families = ((0, 2), (1, 3))
            else:
                local = [(block[a], block[b]) for a, b in HEX_EDGES]
                families = HEX_FAMILIES
                for f, corners in enumerate(HEX_FACES):
                    ordered = tuple(block[c] for c in corners)
                    faces.setdefault(face_key(ordered), []).append((index, f, ordered))
            for k, (a, b) in enumerate(local):
                key = edge_key(a, b)
                edges.setdefault(key, []).append((index, k))
                parent.setdefault(key, key)
            for members in families:
                keys = [edge_key(*local[m]) for m in members]
                for key in keys[1:]:
                    union(keys[0], key)
        boundary_edges = set()
        boundary_faces = set()
        if dim == 2:
            for key, uses in edges.items():
                if len(uses) == 1:
                    boundary_edges.add(key)
                elif len(uses) > 2:
                    issues.append(
                        Issue(
                            f"{where}: edge {key}", f"is shared by {len(uses)} blocks (at most 2)"
                        )
                    )
        else:
            for key, uses in faces.items():
                if len(uses) == 1:
                    boundary_faces.add(key)
                    ordered = uses[0][2]
                    for k in range(4):
                        boundary_edges.add(edge_key(ordered[k], ordered[(k + 1) % 4]))
                elif len(uses) > 2:
                    issues.append(
                        Issue(
                            f"{where}: face {key}", f"is shared by {len(uses)} blocks (at most 2)"
                        )
                    )
        return Lattice(
            dim=dim,
            blocks=tuple(oriented),
            edges=edges,
            faces=faces,
            family={key: find(key) for key in edges},
            boundary_edges=boundary_edges,
            boundary_faces=boundary_faces,
            issues=issues,
        )

    def families(self) -> dict:
        """Family root -> sorted edge keys."""
        out: dict = {}
        for key in sorted(self.edges):
            out.setdefault(self.family[key], []).append(key)
        return out

    def dimension_conflicts(self, dimensions: dict) -> list[Issue]:
        issues = []
        seen: dict = {}
        for key in sorted(dimensions):
            root = self.family[key]
            if root in seen and dimensions[seen[root]] != dimensions[key]:
                other = seen[root]
                issues.append(
                    Issue(
                        f"sketch edges {other} and {key}",
                        f"are in one edge family but have {dimensions[other]} and "
                        f"{dimensions[key]} points",
                    )
                )
                continue
            seen.setdefault(root, key)
        return issues

    def vertex_edges(self, vertex: int) -> list:
        return [key for key in sorted(self.edges) if vertex in key]

    def face_edges(self, key) -> list:
        ordered = self.faces[key][0][2]
        return [edge_key(ordered[k], ordered[(k + 1) % 4]) for k in range(4)]


# ------------------------------------------------------------------------------ placement
@dataclass
class Placement:
    positions: np.ndarray  # (n, 3) vertex positions on the geometry
    classes: list  # per vertex: "fixed", "curve" (may slide), "surface", "free"
    curve_of: dict  # vertex -> curve id (class "curve")
    issues: list
    warnings: list


def _scale(store, params) -> float:
    points = np.asarray(params["vertices"], dtype=float)
    return max(1.0, float(np.max(np.abs(points))))


def _curve_candidates(curve, point) -> list[float]:
    t, _p, _d = curve.project(tuple(float(c) for c in point))
    out = [float(t)]
    if _closed(curve) and (t <= 1e-9 or t >= 1.0 - 1e-9):
        out = [0.0, 1.0]
    return out


def place_vertices(store, params: dict, lattice: Lattice, drawn=None) -> Placement:
    """Vertex positions on the assigned geometry (``drawn``: positions to start from, by
    default the sketch's vertices)."""
    if drawn is None:
        drawn = np.asarray(params["vertices"], dtype=float)
    scale = _scale(store, params)
    tolerance = 1e-6 * scale
    curves = {item["edge"]: item["curve"] for item in params["curves"]}
    surfaces = {item["face"]: item["surface"] for item in params["surfaces"]}
    positions = np.array(drawn, dtype=float)
    classes = []
    curve_of = {}
    issues = []
    warnings = []
    # Adopted connectors: each edge's vertices take the connector's ends (the pairing that
    # moves them least), and every connector adopted at a vertex must end there.
    targets: dict[int, list] = {}
    adopted = {item["edge"]: item["connector"] for item in params.get("connectors", ())}
    for key, (start, end) in adopted_ends(store, params).items():
        a, b = key
        straight = np.linalg.norm(start - drawn[a]) + np.linalg.norm(end - drawn[b])
        flipped = np.linalg.norm(end - drawn[a]) + np.linalg.norm(start - drawn[b])
        if straight <= flipped:
            targets.setdefault(a, []).append((start, adopted[key]))
            targets.setdefault(b, []).append((end, adopted[key]))
        else:
            targets.setdefault(a, []).append((end, adopted[key]))
            targets.setdefault(b, []).append((start, adopted[key]))
    for v in range(len(drawn)):
        on_curves = sorted({curves[e] for e in lattice.vertex_edges(v) if e in curves})
        on_surfaces = sorted(
            {surfaces[f] for f in surfaces if v in f and f in lattice.boundary_faces}
        )
        point = drawn[v]
        if v in targets:
            best = targets[v][0][0]
            apart = [c for q, c in targets[v] if float(np.linalg.norm(q - best)) > tolerance]
            if apart:
                names = ", ".join(store.get(c).name for _q, c in targets[v])
                issues.append(Issue(f"sketch vertex {v}", f"connectors {names} do not meet at it"))
            for curve_id in on_curves:
                if store.geometry(curve_id).project(tuple(best))[2] > tolerance:
                    issues.append(
                        Issue(
                            f"sketch vertex {v}",
                            f"the end of connector {store.get(targets[v][0][1]).name} is not on "
                            f"curve {store.get(curve_id).name}",
                        )
                    )
            positions[v] = best
            classes.append("fixed")
            continue
        if len(on_curves) >= 2:
            geometry = [store.geometry(c) for c in on_curves]
            candidates = []
            for curve in geometry:
                candidates.append(np.asarray(curve.eval(0.0)))
                candidates.append(np.asarray(curve.eval(1.0)))
                candidates.append(np.asarray(curve.project(tuple(point))[1]))
            best = None
            best_distance = math.inf
            for candidate in candidates:
                gap = max(curve.project(tuple(candidate))[2] for curve in geometry)
                if gap > tolerance:
                    continue
                distance = float(np.linalg.norm(candidate - point))
                if distance < best_distance:
                    best = candidate
                    best_distance = distance
            if best is None:
                names = ", ".join(store.get(c).name for c in on_curves)
                issues.append(Issue(f"sketch vertex {v}", f"curves {names} do not meet"))
                classes.append("fixed")
                continue
            positions[v] = best
            classes.append("fixed")
            continue
        if len(on_curves) == 1:
            curve = store.geometry(on_curves[0])
            t, p, _d = curve.project(tuple(point))
            positions[v] = np.asarray(p)
            if (t <= 1e-9 or t >= 1.0 - 1e-9) and not _closed(curve):
                classes.append("fixed")
            else:
                classes.append("curve")
                curve_of[v] = on_curves[0]
            continue
        if on_surfaces:
            surface = store.geometry(on_surfaces[0])
            positions[v] = np.asarray(surface.project(tuple(point))[2])
            if len(on_surfaces) > 1:
                names = ", ".join(store.get(s).name for s in on_surfaces)
                warnings.append(
                    f"sketch vertex {v} is on surfaces {names} without a curve assigned at it: "
                    f"projected onto {store.get(on_surfaces[0]).name}"
                )
            classes.append("surface")
            continue
        classes.append("free")
    return Placement(positions, classes, curve_of, issues, warnings)


def _closed(curve) -> bool:
    start = np.asarray(curve.eval(0.0))
    end = np.asarray(curve.eval(1.0))
    return float(np.linalg.norm(start - end)) <= 1e-9 * max(1.0, float(np.max(np.abs(start))))


def _edge_pieces(curve, pa, pb) -> list[tuple[float, float]]:
    """Pieces (t0, t1) of the curve along an edge from pa to pb: one piece; on a closed curve
    the arc that contains the projection of the edge's midpoint - across the curve's seam
    (two pieces) when needed."""
    ta = _curve_candidates(curve, pa)
    tb = _curve_candidates(curve, pb)
    tm = float(curve.project(tuple(0.5 * (np.asarray(pa) + np.asarray(pb))))[0])
    best = None
    best_score = math.inf
    for a in ta:
        for b in tb:
            if a == b:
                continue
            if min(a, b) - 1e-12 <= tm <= max(a, b) + 1e-12 and abs(b - a) < best_score:
                best = (a, b)
                best_score = abs(b - a)
    if best is not None:
        return [best]
    a = ta[0]
    b = tb[0]
    if not _closed(curve) or a == b:
        return [(a, b)]
    # Across the seam: from a away from b to the curve's end, then on to b.
    if a < b:
        return [(a, 0.0), (1.0, b)]
    return [(a, 1.0), (0.0, b)]


def _pieces_length(curve, pieces) -> float:
    return sum(curve.length_between(min(t0, t1), max(t0, t1)) for t0, t1 in pieces)


# ------------------------------------------------------------------------------ generation
@dataclass
class Plan:
    """Everything a generation writes: entity items by role key."""

    items: list  # (role key, kind, raw params)
    group_key: str
    domains: dict  # 2D: block index -> role key
    warnings: list


def _family_dimensions(store, params, lattice, placement) -> tuple[dict, list[Issue]]:
    given = {item["edge"]: item["points"] for item in params["dimensions"]}
    curves = {item["edge"]: item["curve"] for item in params["curves"]}
    adopted = {item["edge"]: item["connector"] for item in params.get("connectors", ())}
    issues = []
    out = {}
    for root, keys in sorted(lattice.families().items()):
        taken = sorted({store.get(adopted[k]).params["dimension"] for k in keys if k in adopted})
        if taken:
            for key in keys:
                out[key] = taken[0]  # the adopted connector's count (checked at normalisation)
            continue
        values = sorted({given[k] for k in keys if k in given})
        if values:
            for key in keys:
                out[key] = values[0]
            continue
        if params["spacing"] is None:
            listed = ", ".join(str(k) for k in keys)
            issues.append(
                Issue(
                    f"sketch edges {listed}",
                    "have no dimension (set one on any of them, or give spacing)",
                )
            )
            continue
        longest = 0.0
        for key in keys:
            pa = placement.positions[key[0]]
            pb = placement.positions[key[1]]
            length = float(np.linalg.norm(pb - pa))
            if key in curves:
                curve = store.geometry(curves[key])
                length = _pieces_length(curve, _edge_pieces(curve, pa, pb))
            longest = max(longest, length)
        points = max(2, math.ceil(longest / params["spacing"] - 1e-9) + 1)
        for key in keys:
            out[key] = points
        del root
    return {key: _round_points(points, params) for key, points in out.items()}, issues


def _round_points(points: int, params: dict) -> int:
    """Point count with the intervals rounded up to a multiple of ``interval_multiple`` (so
    that p x p cells form the elements of an order-p mesh, 11e)."""
    multiple = params.get("interval_multiple", 1)
    intervals = points - 1
    return math.ceil(intervals / multiple) * multiple + 1


def _along_extrusions(store, surfaces, direction) -> bool:
    """Whether every surface is extruded along ``direction``: the edge between them is then
    the straight line they share."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    for surface in surfaces:
        entity = store.get(surface)
        if entity.kind != "extruded_surface":
            return False
        e = np.asarray(entity.params["direction"], dtype=float)
        if abs(abs(float(d @ e)) / float(np.linalg.norm(e)) - 1.0) > 1e-9:
            return False
    return True


def _edge_length(store, segments, pa, pb, curve_ref) -> float:
    """Length of a sketch edge: along its curve's pieces, else the straight line."""
    if curve_ref is None:
        return float(np.linalg.norm(np.asarray(pb) - np.asarray(pa)))
    curve = store.geometry(curve_ref)
    return float(sum(abs(curve.length_between(s["t0"], s["t1"])) for s in segments))


def _reverse_distribution(d: dict) -> dict:
    out = dict(d)
    out["begin_spacing"] = d["end_spacing"]
    out["end_spacing"] = d["begin_spacing"]
    if d["type"] == "explicit":
        out["positions"] = tuple(1.0 - p for p in reversed(d["positions"]))
    if d["type"] == "geometric" and d["growth_ratio"] is not None:
        out["growth_ratio"] = 1.0 / d["growth_ratio"]
    return out


def _wall_vertices(lattice, params) -> dict:
    """Vertex -> first-cell height of the walls at it (the smallest)."""
    out = {}
    for wall in params["walls"]:
        if wall["spacing"] is None:
            continue
        if "edge" in wall:
            vertices = wall["edge"]
        else:
            vertices = wall["face"]
        for v in vertices:
            out[v] = min(out.get(v, math.inf), wall["spacing"])
    return out


def _wall_edges(lattice, params) -> set:
    out = set()
    for wall in params["walls"]:
        if "edge" in wall:
            out.add(wall["edge"])
        else:
            out.update(lattice.face_edges(wall["face"]))
    return out


def coarse_dimension(n: int) -> int:
    """Point count of the vertex passes' coarser grid (every other interval)."""
    if n < 9:
        return n
    return (n - 1) // 2 + 1


def plan(
    store, sketch, placement: Placement, lattice: Lattice, walls_on: bool = True, coarse=False
) -> Plan:
    """The entities a sketch generates for the given vertex placement (``coarse``: every
    family with half the intervals, for the vertex passes)."""
    params = sketch.params
    base = sketch.name
    dimensions, issues = _family_dimensions(store, params, lattice, placement)
    if coarse:
        dimensions = {key: coarse_dimension(n) for key, n in dimensions.items()}
    if issues or placement.issues:
        raise ValidationError([*placement.issues, *issues])
    curves = {item["edge"]: item["curve"] for item in params["curves"]}
    surfaces = {item["face"]: item["surface"] for item in params["surfaces"]}
    distributions = {
        edge_key(*item["edge"]): (item["edge"], item["distribution"])
        for item in params["distributions"]
    }
    wall_at = _wall_vertices(lattice, params)
    wall_edges = _wall_edges(lattice, params)
    items = []
    warnings = list(placement.warnings)
    pos = placement.positions
    # Adopted connectors: used as they are (``direction``: running from the edge's lower
    # vertex to its higher one); the vertex passes get a sketch-owned copy at half resolution.
    adopted = {item["edge"]: item["connector"] for item in params.get("connectors", ())}
    ends = adopted_ends(store, params)
    direction = {}
    unstructured = set()
    if not coarse:
        unstructured = set(params.get("unstructured", ()))
    for key in sorted(lattice.edges):
        a, b = key
        pa = pos[a]
        pb = pos[b]
        if float(np.linalg.norm(pb - pa)) == 0.0:
            raise ValidationError([Issue(f"sketch edge {key}", "has coincident vertices")])
        if key in adopted:
            start = ends[key][0]
            direction[key] = bool(np.linalg.norm(start - pa) <= np.linalg.norm(start - pb))
            if coarse:
                source = store.get(adopted[key]).params
                items.append(
                    (
                        f"c{a}-{b}",
                        "connector",
                        {
                            "segments": source["segments"],
                            "dimension": dimensions[key],
                            "distribution": source["distribution"],
                            "propagate_spacing": False,
                        },
                    )
                )
            continue
        direction[key] = True
        if key in curves:
            curve = store.geometry(curves[key])
            pieces = [(t0, t1) for t0, t1 in _edge_pieces(curve, pa, pb) if t0 != t1]
            if not pieces:
                raise ValidationError(
                    [Issue(f"sketch edge {key}", "both vertices are at one point of its curve")]
                )
            segments = tuple({"curve": curves[key], "t0": t0, "t1": t1} for t0, t1 in pieces)
        else:
            on = sorted(
                {surfaces[f] for f in surfaces if key in lattice.face_edges(f)},
            )
            line_key = f"e{a}-{b}"
            if lattice.dim == 3 and len(on) == 1:
                surface = store.geometry(on[0])
                points = []
                for s in np.linspace(0.0, 1.0, SURFACE_SPLINE_POINTS):
                    q = (1.0 - s) * pa + s * pb
                    points.append(tuple(np.asarray(surface.project(tuple(q))[2]).tolist()))
                points[0] = tuple(pa.tolist())
                points[-1] = tuple(pb.tolist())
                items.append((line_key, "spline", {"points": tuple(points)}))
            else:
                if lattice.dim == 3 and len(on) > 1 and not _along_extrusions(store, on, pb - pa):
                    names = ", ".join(store.get(s).name for s in on)
                    warnings.append(
                        f"sketch edge {key} lies between surfaces {names} and has no curve: "
                        "it stays straight"
                    )
                items.append(
                    (line_key, "line", {"p0": tuple(pa.tolist()), "p1": tuple(pb.tolist())})
                )
            segments = ({"curve": f"@{line_key}", "t0": 0.0, "t1": 1.0},)
        distribution = {"type": "uniform"}
        spaced = None
        if key in distributions:
            edge, given = distributions[key]
            distribution = given
            if edge != key:
                distribution = _reverse_distribution(given)
        elif key not in wall_edges and (a in wall_at or b in wall_at):
            # (An edge between two walls, e.g. an outlet: layers slide in at both ends.)
            spaced, edge_warnings = _leaving_unstructured_wall(
                store,
                base,
                lattice,
                params,
                key,
                unstructured,
                pos,
                _edge_length(store, segments, pa, pb, curves.get(key)),
                dimensions[key],
            )
            warnings.extend(edge_warnings)
            if spaced is None and (a in wall_at) != (b in wall_at):
                # An edge leaving a wall: clustered at the wall (a better initial guess).
                distribution = {"type": "tanh", "begin_spacing": None, "end_spacing": None}
                if a in wall_at:
                    distribution["begin_spacing"] = wall_at[a]
                else:
                    distribution["end_spacing"] = wall_at[b]
        connector = {
            "segments": segments,
            "dimension": dimensions[key],
            "distribution": distribution,
            "propagate_spacing": False,
        }
        if spaced is not None:
            connector.update(spaced)
        items.append((f"c{a}-{b}", "connector", connector))

    def loop(ordered):
        edges = []
        for k in range(4):
            a = ordered[k]
            b = ordered[(k + 1) % 4]
            key = edge_key(a, b)
            connector = f"@c{key[0]}-{key[1]}"
            if key in adopted and not coarse:
                connector = adopted[key]
            backwards = a > b
            if not direction[key]:
                backwards = a < b
            edges.append([{"connector": connector, "reversed": backwards}])
        return edges

    domains = {}
    walls = []
    group = None
    if lattice.dim == 2:

        def connector_ref(key):
            if key in adopted and not coarse:
                return adopted[key]
            return f"@c{key[0]}-{key[1]}"

        structured = [k for k in range(len(lattice.blocks)) if k not in unstructured]
        for index, block in enumerate(lattice.blocks):
            domains[index] = f"d{index}"
            if index in unstructured:
                generator = _block_generator(params, block, connector_ref)
                outer = [item for edge in loop(block) for item in edge]  # one loop, flat
                items.append(
                    (
                        f"d{index}",
                        "unstructured_domain",
                        {"loops": [outer], "generator": generator},
                    )
                )
                if params["extrude"] is not None:
                    items.append((f"x{index}", "volume", {"extrusion": _extrusion(params, index)}))
                continue
            items.append((f"d{index}", "domain", {"edges": loop(block), "generator": None}))
            if params["extrude"] is not None:
                items.append((f"x{index}", "block", {"extrusion": _extrusion(params, index)}))
        structured_edges = {
            edge_key(lattice.blocks[k][m], lattice.blocks[k][(m + 1) % 4])
            for k in structured
            for m in range(4)
        }
        pinned = list(params["pinned"])
        for wall in params["walls"]:
            if wall["edge"] not in structured_edges:
                continue  # a wall of an unstructured block only: its inflation layers
            if walls_on:
                entry = {k: wall[k] for k in WALL_KEYS}
                entry["connector"] = connector_ref(wall["edge"])
                walls.append(entry)
            elif wall["edge"] not in lattice.boundary_edges:
                pinned.append(wall["edge"])  # internal walls stay in place in the passes
        pinned_refs = [f"@c{a}-{b}" for a, b in sorted(pinned)]
        if not coarse:
            # Adopted connectors keep their nodes: a domain on their other side shares them;
            # so do the connectors between a structured and an unstructured block.
            pinned_refs.extend(adopted[key] for key in sorted(adopted))
            for key, uses in sorted(lattice.edges.items()):
                blocks_at = {use[0] for use in uses}
                if blocks_at & unstructured and blocks_at - unstructured:
                    ref = connector_ref(key)
                    if ref not in pinned_refs:
                        pinned_refs.append(ref)
        group = {
            "domains": tuple(f"@d{k}" for k in structured),
            "walls": tuple(walls),
            "pinned": tuple(pinned_refs),
            "joints": "fixed",
            "slide_limit": params["slide_limit"],
            "tolerance": params["tolerance"],
            "max_iterations": params["max_iterations"],
        }
        if params["relax"]:
            group["joints"] = "slide"
        if not structured:
            group = None  # every block unstructured: nothing to smooth
    else:
        face_role = {}
        for key in sorted(lattice.faces):
            ordered = lattice.faces[key][0][2]
            corners = pos[list(ordered)]
            if float(np.ptp(corners[:, 2])) <= CHAIN_TOLERANCE * max(
                1.0, float(np.max(np.abs(corners)))
            ):
                # A face in a plane z = const (as domains test it): counter-clockwise in (x, y).
                area = sum(
                    corners[k, 0] * corners[(k + 1) % 4, 1]
                    - corners[(k + 1) % 4, 0] * corners[k, 1]
                    for k in range(4)
                )
                if area < 0.0:
                    ordered = (ordered[0], ordered[3], ordered[2], ordered[1])
            role = "f" + "-".join(str(v) for v in key)
            face_role[key] = role
            items.append(
                (
                    role,
                    "domain",
                    {"edges": loop(ordered), "generator": None, "surface": surfaces.get(key)},
                )
            )
        for index, block in enumerate(lattice.blocks):
            keys = [face_key(tuple(block[c] for c in corners)) for corners in HEX_FACES]
            items.append((f"b{index}", "block", {"faces": [f"@{face_role[k]}" for k in keys]}))
        pinned = list(params["pinned"])
        for wall in params["walls"]:
            if walls_on:
                entry = {k: wall[k] for k in WALL_KEYS}
                entry["domain"] = f"@{face_role[wall['face']]}"
                walls.append(entry)
            elif wall["face"] not in lattice.boundary_faces:
                pinned.append(wall["face"])  # internal walls stay in place in the passes
        group = {
            "blocks": tuple(f"@b{k}" for k in range(len(lattice.blocks))),
            "walls": tuple(walls),
            "pinned": tuple(f"@{face_role[f]}" for f in sorted(pinned)),
            "symmetry": tuple(f"@{face_role[f]}" for f in params["symmetry"]),
            # Face corners inside one surface or plane slide on it (on a plane they are
            # symmetry-plane nodes: interior nodes of the sections).
            "joints": "slide",
            "slide_limit": params["slide_limit"],
            "tolerance": params["tolerance"],
            "max_iterations": params["max_iterations"],
        }
    if group is not None:
        items.append(("group", "smoothing_group", group))
    del base
    return Plan(items=items, group_key="group", domains=domains, warnings=warnings)


def _block_generator(params: dict, block, connector_ref) -> dict:
    """The generator of an unstructured block: the sketch's ``generator`` (or the default)
    with inflation layers from the walls on the block's edges (the smallest first-cell
    height, the most layers, the first wall's growth)."""
    from opengrid.model.unstructured import DEFAULT_INFLATION

    generator = {}
    if params.get("generator") is not None:
        generator = dict(params["generator"])
    edges = {edge_key(block[k], block[(k + 1) % 4]) for k in range(4)}
    connectors = []
    heights = []
    layers = []
    growth = None
    transition = True
    slide = False
    for wall in params["walls"]:
        if "edge" not in wall or wall["edge"] not in edges:
            continue
        transition = transition and wall.get("transition", True)
        slide = slide or wall.get("slide", False)
        if wall["spacing"] is None:
            raise ValidationError(
                [
                    Issue(
                        f"sketch edge {wall['edge']}",
                        "is a wall of an unstructured block: give its first-cell height (the "
                        "first height of the inflation layers)",
                    )
                ]
            )
        connectors.append(connector_ref(wall["edge"]))
        heights.append(wall["spacing"])
        if wall["layers"] > 0:
            layers.append(wall["layers"])
        if growth is None:
            growth = wall["growth"]
    generator["inflation"] = None
    if connectors:
        count = DEFAULT_INFLATION["layers"]
        if layers:
            count = max(layers)
        generator["inflation"] = {
            "connectors": tuple(connectors),
            "first_height": min(heights),
            "layers": count,
            "growth_rate": growth,
            "transition": transition,
            "slide": slide,
        }
    return generator


def _leaving_unstructured_wall(
    store, base, lattice, params, key, unstructured, pos, length, dimension
):
    """The spacing of an edge with one end at a wall whose layers slide into it (``slide`` on
    the walls of the unstructured blocks it borders, and no structured block uses it): points
    at the layer heights of the block's inflation, 1 / sin(angle) farther apart along it
    (layer_spacing.layered_params). None otherwise (clustered at the wall as before), and for
    a periodic edge (its partner keeps its nodes). Returns (connector parameters to update,
    warnings)."""
    from opengrid.model.layer_spacing import (
        SLIDE_MAX_ANGLE,
        SLIDE_MIN_ANGLE,
        layer_heights,
        layered_params,
    )
    from opengrid.model.patches import partner_of

    blocks_at = {use[0] for use in lattice.edges[key]}
    if not blocks_at or not blocks_at <= unstructured:
        return None, []
    warnings = []
    stacks = {}
    growth = None
    a, b = key
    for k in sorted(blocks_at):
        block = lattice.blocks[k]
        inflation = _block_generator(params, block, lambda edge: edge)["inflation"]
        if inflation is None or not inflation["slide"]:
            continue
        for v, other in ((a, b), (b, a)):
            walls = [e for e in inflation["connectors"] if v in e and e != key]
            if not walls or v in stacks:
                continue
            wall = walls[0]
            far = wall[0]
            if far == v:
                far = wall[1]
            u1 = pos[far] - pos[v]
            u2 = pos[other] - pos[v]
            cosine = float(np.dot(u1, u2) / (np.linalg.norm(u1) * np.linalg.norm(u2)))
            angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
            if not SLIDE_MIN_ANGLE <= angle <= SLIDE_MAX_ANGLE:
                warnings.append(
                    f"sketch edge {key} meets the wall {wall} at {angle:.0f} degrees (layers "
                    f"slide between {SLIDE_MIN_ANGLE:.0f} and {SLIDE_MAX_ANGLE:.0f}): the "
                    "layers taper there"
                )
                continue
            scale = 1.0 / math.sin(math.radians(angle))
            heights = layer_heights(
                inflation["first_height"], inflation["growth_rate"], inflation["layers"]
            )
            stacks[v] = [h * scale for h in heights]
            growth = inflation["growth_rate"]
    if not stacks:
        return None, warnings
    name = f"{base}-c{a}-{b}"
    if name in store:
        connector_id = store.get(name).id
        for patch in store.entities("patch"):
            members = {m.get("connector") for m in patch.params["members"]}
            if connector_id in members and (
                patch.params["patch_type"] == "periodic" or partner_of(patch, store) is not None
            ):
                warnings.append(
                    f"{name} is periodic: it keeps its points (the layers taper towards it)"
                )
                return None, warnings
    try:
        spaced = layered_params(
            {"dimension": dimension}, length, growth, begin=stacks.get(a), end=stacks.get(b)
        )
    except ValueError as error:
        warnings.append(f"sketch edge {key}: {error}; it keeps its spacing (the layers taper)")
        return None, warnings
    return spaced, warnings


def _resolve(value, ids: dict):
    """Replace "@role" references by entity ids."""
    if isinstance(value, str) and value.startswith("@"):
        return ids[value[1:]]
    if isinstance(value, dict):
        return {k: _resolve(v, ids) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(_resolve(v, ids) for v in value)
    return value


def write_plan(store, sketch_id: int, the_plan: Plan) -> dict:
    """Create, update or delete the sketch's outputs to match ``the_plan`` (inside the open
    transaction); entities keep their ids and names across regenerations. Returns role -> id."""
    from opengrid.model.entities import normalize_params

    sketch = store.get(sketch_id)
    old = dict(sketch.params["outputs"])
    ids = {}
    wanted = {role for role, _kind, _raw in the_plan.items}
    # Outputs no longer generated: removed first (dependents before what they use).
    order = {"smoothing_group": 0, "block": 1, "domain": 2, "connector": 3}
    stale = [(role, entity_id) for role, entity_id in old.items() if role not in wanted]
    stale = [(r, i) for r, i in stale if i in store]
    stale.sort(key=lambda item: (order.get(store.get(item[1]).kind, 4), item[0]))
    from opengrid.model.membership import MEMBERSHIP_KINDS, drop_members

    for role, entity_id in stale:
        outside = [
            d.name
            for d in store.dependents(entity_id)
            if d.id not in old.values() and d.kind != KIND and d.kind not in MEMBERSHIP_KINDS
        ]
        if outside:
            raise ValidationError(
                [
                    Issue(
                        store.get(entity_id).name,
                        f"is no longer generated by sketch {sketch.name!r} but is used by "
                        f"{', '.join(outside)}",
                    )
                ]
            )
        # Patches, zones and groups listing it lose the member.
        for group_id in drop_members(store, {entity_id}):
            if group_id in store and group_id not in old.values():
                store.remove(group_id)
        store.remove(entity_id)
        del role
    for role, kind, raw in the_plan.items:
        resolved = _resolve(raw, ids)
        name = f"{sketch.name}-{role}"
        if role in old and old[role] in store and store.get(old[role]).kind == kind:
            entity = store.get(old[role])
            params = normalize_params(kind, resolved, entity.name, store)
            if params != entity.params:
                store.update(entity.id, params=params)
            ids[role] = entity.id
            continue
        if role in old and old[role] in store and store.get(old[role]).name == name:
            # The role now makes another kind of entity (e.g. a line became a spline on a
            # surface): replace it; what referred to it is regenerated after it.
            store.remove(old[role])
        if name in store:
            raise ValidationError(
                [Issue(name, f"is already in use (sketch {sketch.name!r} generates it)")]
            )
        params = normalize_params(kind, resolved, name, store)
        ids[role] = store.add(kind, params, name).id
    params = dict(sketch.params)
    params["outputs"] = tuple(sorted(ids.items()))
    store.update(sketch_id, params=params)
    return ids


def _evaluators(store):
    """A fresh group evaluator (with its connector, domain and block evaluators) for the
    store as it is inside the transaction."""
    from opengrid.model.blocks import BlockEvaluator
    from opengrid.model.connectors import ConnectorEvaluator
    from opengrid.model.domains import DomainEvaluator
    from opengrid.model.groups import GroupEvaluator

    connectors = ConnectorEvaluator(store)
    domains = DomainEvaluator(store, connectors)
    groups = GroupEvaluator(store, connectors, domains)
    domains.groups = groups
    blocks = BlockEvaluator(store, domains)
    blocks.groups = groups
    groups.blocks = blocks
    return groups


def _corner(grid: np.ndarray, k: int) -> np.ndarray:
    if k == 0:
        return grid[0, 0]
    if k == 1:
        return grid[0, -1]
    if k == 2:
        return grid[-1, -1]
    return grid[-1, 0]


def _local_spacing(positions, lattice, dimensions) -> np.ndarray:
    out = np.full(len(positions), math.inf)
    for key in lattice.edges:
        a, b = key
        h = float(np.linalg.norm(positions[b] - positions[a])) / (dimensions[key] - 1)
        out[a] = min(out[a], h)
        out[b] = min(out[b], h)
    return out


@dataclass
class Generation:
    names: list
    warnings: list
    passes: int
    moved: float  # last vertex move in local spacings (passes)


def generate(store, sketch_id: int, passes: bool = True) -> Generation:
    """Regenerate a sketch's outputs (inside the open transaction). ``passes``: run the
    vertex passes (when the sketch's ``relax`` is on); previews skip them."""
    sketch = store.get(sketch_id)
    params = sketch.params
    lattice = Lattice.build(params["vertices"], params["blocks"], sketch.name)
    placement = place_vertices(store, params, lattice)
    pinned = set(params["fixed"])
    held = _held_vertices(lattice, params)
    joints = _joints(store, params, lattice, placement, pinned | held)
    run = passes
    passes = 0
    moved = 0.0
    if run and params["relax"] and joints:
        passes, moved = _vertex_passes(store, sketch_id, lattice, placement, joints, pinned | held)
        sketch = store.get(sketch_id)
    the_plan = plan(store, sketch, placement, lattice, walls_on=True)
    ids = write_plan(store, sketch_id, the_plan)
    names = [store.get(i).name for _r, i in sorted(ids.items())]
    warnings = list(the_plan.warnings)
    if run:
        warnings.extend(
            f"{s['what']}: {s['message']}" for s in suggestions(store, store.get(sketch_id))
        )
    if passes >= MAX_PASSES and moved >= PASS_TOLERANCE:
        warnings.append(
            f"sketch {sketch.name}: vertices still moving after {passes} passes "
            f"({moved:.2g} local spacings)"
        )
    if passes:
        log.info(
            "sketch %s: %d vertex passes (last move %.2g local spacings)",
            sketch.name,
            passes,
            moved,
        )
    return Generation(names, warnings, passes, moved)


def _vertex_passes(store, sketch_id, lattice, placement, joints: dict, held: set):
    """Move the sliding vertices (``joints``: vertex -> projection onto its curve, plane or
    surface) to where the group smoothing puts them, pass after pass on the coarser grid,
    until they stop (Anderson mixing of their positions); interior vertices follow too.
    Updates ``placement``; returns (passes, last move in local spacings)."""
    sketch = store.get(sketch_id)
    params = sketch.params
    dimensions, _issues = _family_dimensions(store, params, lattice, placement)
    boundary_vertices = _boundary_vertices(lattice)
    free = [
        v
        for v in range(len(placement.positions))
        if placement.classes[v] == "free" and v not in held and v not in boundary_vertices
    ]
    order = sorted(joints)
    x = np.concatenate([placement.positions[v] for v in order])
    history_x: list[np.ndarray] = []
    history_r: list[np.ndarray] = []
    previous = None
    absolute = None
    passes = 0
    moved = 0.0
    while passes < MAX_PASSES:
        for k, v in enumerate(order):
            placement.positions[v] = joints[v](x[3 * k : 3 * k + 3])
        the_plan = plan(store, sketch, placement, lattice, walls_on=False, coarse=True)
        group_raw = dict(the_plan.items[-1][2])
        group_raw["tolerance"] = PASS_SMOOTHING_TOLERANCE
        group_raw["max_iterations"] = PASS_MAX_ITERATIONS
        group_raw["slide_limit"] = PASS_SLIDE_LIMIT
        the_plan.items[-1] = (the_plan.items[-1][0], "smoothing_group", group_raw)
        ids = write_plan(store, sketch_id, the_plan)
        sketch = store.get(sketch_id)
        groups = _evaluators(store)
        result, previous, nodes = _smooth_pass(
            groups, ids[the_plan.group_key], previous, absolute, _member_ids(lattice, ids)
        )
        if absolute is None:
            absolute = _absolute_tolerance(result)
        passes += 1
        smoothed = _smoothed_vertices(lattice, placement.positions, nodes, previous)
        h = _local_spacing(placement.positions, lattice, dimensions)
        moves = [0.0]
        for v in free:
            moves.append(float(np.linalg.norm(smoothed[v] - placement.positions[v])) / h[v])
            placement.positions[v] = smoothed[v]
        image = np.concatenate([joints[v](smoothed[v]) for v in order])
        for k, v in enumerate(order):
            moves.append(
                float(np.linalg.norm(image[3 * k : 3 * k + 3] - x[3 * k : 3 * k + 3])) / h[v]
            )
        moved = max(moves)
        if moved < PASS_TOLERANCE:
            x = image
            break
        size = _block_size(placement.positions, lattice)
        far = [
            float(np.linalg.norm(image[3 * k : 3 * k + 3] - x[3 * k : 3 * k + 3])) > 2.0 * h[v]
            for k, v in enumerate(order)
        ]
        if any(far):
            # Still far from the fixed point: the passes' slide limit clips the moves, so the
            # map is not smooth there - damped plain steps until close.
            history_x.clear()
            history_r.clear()
            x = x + 0.5 * (image - x)
            continue
        x = _anderson_step(x, image, history_x, history_r, 0.25 * size[order])
    for k, v in enumerate(order):
        placement.positions[v] = joints[v](x[3 * k : 3 * k + 3])
    return passes, moved


def _block_size(positions, lattice) -> np.ndarray:
    """Shortest sketch edge at every vertex."""
    out = np.full(len(positions), math.inf)
    for a, b in lattice.edges:
        length = float(np.linalg.norm(positions[b] - positions[a]))
        out[a] = min(out[a], length)
        out[b] = min(out[b], length)
    return out


def _boundary_vertices(lattice) -> set:
    out = set()
    if lattice.dim == 2:
        for edge in lattice.boundary_edges:
            out.update(edge)
    else:
        for face in lattice.boundary_faces:
            out.update(face)
    return out


def _anderson_step(x, image, history_x, history_r, limit, depth: int = 5) -> np.ndarray:
    """Next sliding vertex positions for the fixed point x = image(x): Anderson mixing over
    the last ``depth`` passes, falling back to the plain step when the mixed step is much
    longer or moves a vertex further than ``limit`` (per vertex; the caller projects the
    result back onto the geometry)."""
    residual = image - x
    history_x.append(x.copy())
    history_r.append(residual.copy())
    del history_x[: -depth - 1]
    del history_r[: -depth - 1]
    if len(history_r) < 2:
        return image
    d_r = np.stack([history_r[k + 1] - history_r[k] for k in range(len(history_r) - 1)], axis=1)
    d_x = np.stack([history_x[k + 1] - history_x[k] for k in range(len(history_x) - 1)], axis=1)
    gamma, *_rest = np.linalg.lstsq(d_r, residual, rcond=None)
    mixed = x + residual - (d_x + d_r) @ gamma
    step = float(np.max(np.abs(mixed - x)))
    per_vertex = np.linalg.norm((mixed - x).reshape(-1, 3), axis=1)
    if (
        not np.all(np.isfinite(mixed))
        or step > 4.0 * float(np.max(np.abs(residual)))
        or np.any(per_vertex > limit)
    ):
        history_x.clear()
        history_r.clear()
        return image
    return mixed


def _smooth_pass(groups, group_id: int, previous, absolute, members):
    """One pass's group smoothing, warm-started from the previous pass where the node
    numbering is the same (its interior nodes start where they ended). ``members``: block
    index -> member entity id. Returns (result, final node positions, block index ->
    (corner node ids, their initial positions))."""
    from opengrid.model.groups import run_group

    _key, layout, params = groups.prepare(group_id)
    params = dict(params)
    if absolute is not None:
        params["absolute_tolerance"] = absolute
    if previous is not None and previous.shape == layout.points.shape:
        # Interior nodes only: sliding nodes start from their connectors' new distribution,
        # in order between the joints.
        free = layout.fixed == 0
        for item in layout.slides:
            free[item[0]] = False
        layout.points[free] = previous[free]
    corners = {}
    for block, member in members.items():
        ids = layout.ids[member]
        if ids.ndim == 2:
            nodes = ids[[0, 0, -1, -1], [0, -1, 0, -1]]
        else:
            nodes = np.array([ids[k, j, i] for k in (0, -1) for j in (0, -1) for i in (0, -1)])
        corners[block] = (nodes, layout.points[nodes].copy())
    raw = run_group(layout, params)
    result = groups.finish(params, layout, raw)
    return result, np.asarray(raw[0], dtype=float), corners


def _member_ids(lattice, ids) -> dict:
    """Block index -> the id of its domain (2D) or block (3D) in a pass."""
    if lattice.dim == 2:
        return {b: ids[f"d{b}"] for b in range(len(lattice.blocks))}
    return {b: ids[f"b{b}"] for b in range(len(lattice.blocks))}


def _smoothed_vertices(lattice, positions, corners, final) -> np.ndarray:
    """Where the smoothing put every vertex: the block corner node that started at it."""
    out = np.array(positions, dtype=float)
    for v in range(len(positions)):
        block, _corner = _first_corner(lattice, v)
        nodes, initial = corners[block]
        k = int(np.argmin(np.linalg.norm(initial - positions[v], axis=1)))
        out[v] = final[nodes[k]]
    return out


def _absolute_tolerance(result) -> float:
    """The first pass's relative tolerance as an absolute one for the warm-started passes."""
    return PASS_SMOOTHING_TOLERANCE * float(result.residuals[0])


def _first_corner(lattice, v) -> tuple[int, int]:
    for index, block in enumerate(lattice.blocks):
        if v in block:
            return index, block.index(v)
    raise ValueError(f"vertex {v} is in no block")


def _held_vertices(lattice, params) -> set:
    """Vertices on pinned edges (3D: faces), internal walls or adopted connectors: they stay
    where drawn (placed)."""
    out = set()
    for ref in params["pinned"]:
        out.update(ref)
    for item in params.get("connectors", ()):
        out.update(item["edge"])
    for wall in params["walls"]:
        if "edge" in wall and wall["edge"] not in lattice.boundary_edges:
            out.update(wall["edge"])
        if "face" in wall and wall["face"] not in lattice.boundary_faces:
            out.update(wall["face"])
    return out


def _joints(store, params, lattice, placement, held) -> dict:
    """Vertices that slide during the vertex passes -> their projection (point -> point).
    2D: a vertex joining exactly two boundary edges assigned to one curve (inside it). 3D: a
    boundary vertex without assigned curves whose boundary faces all lie on one assigned
    surface, or are unassigned and in one plane (vertices on curves stay where projected)."""
    out = {}
    if lattice.dim == 2:
        for v, cls in enumerate(placement.classes):
            if cls != "curve" or v in held:
                continue
            boundary = [e for e in lattice.vertex_edges(v) if e in lattice.boundary_edges]
            if len(boundary) == 2:
                out[v] = _curve_projector(store.geometry(placement.curve_of[v]))
        return out
    surfaces = {item["face"]: item["surface"] for item in params["surfaces"]}
    scale = max(1.0, float(np.max(np.abs(placement.positions))))
    for v, cls in enumerate(placement.classes):
        if cls not in ("free", "surface") or v in held:
            continue
        faces = [f for f in lattice.boundary_faces if v in f]
        if not faces:
            continue
        assigned = {surfaces.get(f) for f in faces}
        if len(assigned) != 1:
            continue
        surface = assigned.pop()
        if surface is not None:
            out[v] = _surface_projector(store.geometry(surface))
            continue
        corners = placement.positions[sorted({c for f in faces for c in f})]
        centre = corners.mean(axis=0)
        _u, _s, vt = np.linalg.svd(corners - centre)
        normal = vt[2]
        if float(np.max(np.abs((corners - centre) @ normal))) <= 1e-9 * scale:
            out[v] = _plane_projector(centre, normal)
    return out


def _curve_projector(curve):
    def project(point):
        return np.asarray(curve.project(tuple(float(c) for c in point))[1])

    return project


def _surface_projector(surface):
    def project(point):
        return np.asarray(surface.project(tuple(float(c) for c in point))[2])

    return project


def _plane_projector(origin, normal):
    def project(point):
        point = np.asarray(point, dtype=float)
        return point - float((point - origin) @ normal) * normal

    return project


# ------------------------------------------------------------------------------ suggestions
def suggestions(
    store,
    sketch,
    tolerance: float | None = None,
    only=None,
    fraction: float = SUGGEST_FRACTION,
    middle_fraction: float | None = None,
) -> list[dict]:
    """Unassigned boundary edges (and 3D faces) lying near geometry: both ends (corners)
    within ``tolerance`` - by default ``fraction`` (5 %) of the edge length (face size) - and
    the middle (centre) within ``middle_fraction`` of it (default: ``fraction``; a chord
    of a curve sags in the middle). The nearest curve (surface) is offered; ``only``: the ids
    to consider (default: all that the sketch does not generate). Never applied silently."""
    from opengrid.model.connectors import CURVE_KINDS
    from opengrid.model.entities import SURFACE_KINDS

    if middle_fraction is None:
        middle_fraction = fraction
    params = sketch.params
    lattice = Lattice.build(params["vertices"], params["blocks"], sketch.name)
    points = np.asarray(params["vertices"], dtype=float)
    assigned = {item["edge"] for item in params["curves"]}
    adopted = {item["edge"] for item in params.get("connectors", ())}
    assigned_faces = {item["face"] for item in params["surfaces"]}
    owned = set(owned_entities(store))
    curves = [c for c in store.entities() if c.kind in CURVE_KINDS and c.id not in owned]
    if only is not None:
        curves = [c for c in curves if c.id in only]
    connectors = []
    if lattice.dim == 2:
        connectors = _connector_candidates(store, only)
    out = []
    for key in sorted(lattice.boundary_edges):
        if key in adopted:
            continue
        pa = points[key[0]]
        pb = points[key[1]]
        length = float(np.linalg.norm(pb - pa))
        end_limit = fraction * length
        middle_limit = middle_fraction * length
        if tolerance is not None:
            end_limit = tolerance
            middle_limit = tolerance
        best = _connector_under(pa, pb, connectors, end_limit, middle_limit)
        if best is not None:
            out.append(
                {
                    "what": f"sketch {sketch.name} edge {key}",
                    "edge": key,
                    "connector": best[1],
                    "message": f"lies on connector {best[1]} (not adopted: adopt it to share "
                    "its nodes with the domain on its other side)",
                }
            )
    for key in sorted(lattice.boundary_edges):
        if key in assigned or key in adopted or any(s.get("edge") == key for s in out):
            continue
        pa = points[key[0]]
        pb = points[key[1]]
        length = float(np.linalg.norm(pb - pa))
        end_limit = fraction * length
        middle_limit = middle_fraction * length
        if tolerance is not None:
            end_limit = tolerance
            middle_limit = tolerance
        best = None
        for curve in curves:
            geometry = store.geometry(curve.id)
            ends = max(geometry.project(tuple(q))[2] for q in (pa, pb))
            middle = geometry.project(tuple(0.5 * (pa + pb)))[2]
            gap = max(ends, middle)
            if ends <= end_limit and middle <= middle_limit and (best is None or gap < best[1]):
                best = (curve.name, gap)
        if best is not None:
            out.append(
                {
                    "what": f"sketch {sketch.name} edge {key}",
                    "edge": key,
                    "curve": best[0],
                    "message": f"lies within {best[1]:.3g} of curve {best[0]} (not assigned)",
                }
            )
    if lattice.dim == 3:
        surfaces = [e for e in store.entities() if e.kind in SURFACE_KINDS]
        if only is not None:
            surfaces = [e for e in surfaces if e.id in only]
        for key in sorted(lattice.boundary_faces):
            if key in assigned_faces:
                continue
            corners = points[list(key)]
            centre = corners.mean(axis=0)
            size = float(np.max(np.linalg.norm(corners - centre, axis=1)))
            end_limit = fraction * size
            middle_limit = middle_fraction * size
            if tolerance is not None:
                end_limit = tolerance
                middle_limit = tolerance
            best = None
            for surface in surfaces:
                geometry = store.geometry(surface.id)
                ends = max(geometry.project(tuple(q))[3] for q in corners)
                middle = geometry.project(tuple(centre))[3]
                gap = max(ends, middle)
                if ends <= end_limit and middle <= middle_limit and (best is None or gap < best[1]):
                    best = (surface.name, gap)
            if best is not None:
                out.append(
                    {
                        "what": f"sketch {sketch.name} face {key}",
                        "face": key,
                        "surface": best[0],
                        "message": f"lies within {best[1]:.3g} of surface {best[0]} (not assigned)",
                    }
                )
    return out


# ------------------------------------------------------------------------------ ownership
def owned_entities(store) -> dict[int, int]:
    """Generated entity id -> id of the sketch that owns it."""
    out = {}
    for sketch in store.entities(KIND):
        for _role, entity_id in sketch.params["outputs"]:
            out[entity_id] = sketch.id
    return out


def check_ownership(store, touched: dict) -> None:
    """Refuse a transaction that changed a sketch's outputs without touching the sketch
    (``touched``: entity id -> state before the transaction)."""
    owners = owned_entities(store)
    for entity_id, before in touched.items():
        if before is not None and before.kind == KIND:
            for _role, output in before.params["outputs"]:
                owners.setdefault(output, entity_id)
    issues = []
    for entity_id, before in sorted(touched.items()):
        owner = owners.get(entity_id)
        if owner is None or owner in touched:
            continue
        after = None
        if entity_id in store:
            after = store.get(entity_id)
        if after == before:
            continue
        shown = after
        if shown is None:
            shown = before
        sketch = store.get(owner).name
        issues.append(
            Issue(
                shown.name,
                f"is generated by sketch {sketch!r}: edit the sketch instead "
                f"(og.edit_sketch({sketch!r}, ...))",
            )
        )
    if issues:
        raise ValidationError(issues)


# ------------------------------------------------------------------------- edit functions
def _canonical_input(value):
    """User input (NumPy arrays, dicts with tuple keys) as plain Python values."""
    if isinstance(value, np.ndarray):
        return _canonical_input(value.tolist())
    if isinstance(value, dict):
        return {_canonical_input(k): _canonical_input(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return tuple(_canonical_input(v) for v in value)
    if hasattr(value, "item") and not isinstance(value, (int, float, str)):
        return value.item()
    return value


def create_sketch(store, name=None, **raw):
    """Create a sketch and generate its outputs. Returns ([sketch name], warnings)."""
    from opengrid.model.entities import normalize_params

    raw = {k: _canonical_input(v) for k, v in raw.items()}
    raw["outputs"] = ()
    where = name
    if where is None:
        where = "new sketch"
    params = normalize_params(KIND, raw, where, store)
    entity = store.add(KIND, params, name)
    result = generate(store, entity.id)
    return [entity.name], result.warnings


def edit_sketch(store, name, passes: bool = True, **changes):
    """Change sketch parameters and regenerate. Returns ([sketch name], warnings)."""
    from opengrid.model.entities import normalize_params

    entity = store.get(name)
    if entity.kind != KIND:
        raise ValidationError([Issue(str(name), f"is a {entity.kind}, not a sketch")])
    if "outputs" in changes:
        raise ValidationError([Issue(entity.name, "outputs are maintained by the sketch")])
    merged = dict(sketch_script_params(entity.params, store))
    merged.update({k: _canonical_input(v) for k, v in changes.items()})
    params = normalize_params(KIND, merged, entity.name, store)
    store.update(entity.id, params=params)
    result = generate(store, entity.id, passes=passes)
    return [entity.name], result.warnings


def regenerate_sketch(store, name, passes: bool = True):
    entity = store.get(name)
    if entity.kind != KIND:
        raise ValidationError([Issue(str(name), f"is a {entity.kind}, not a sketch")])
    store.update(entity.id, params=dict(entity.params))
    result = generate(store, entity.id, passes=passes)
    return [entity.name], result.warnings


def report(store, sketch) -> str:
    """Vertices (class and position), edge families and their dimensions, suggestions."""
    params = sketch.params
    lattice = Lattice.build(params["vertices"], params["blocks"], sketch.name)
    placement = place_vertices(store, params, lattice)
    held = _held_vertices(lattice, params)
    joints = _joints(store, params, lattice, placement, set(params["fixed"]) | held)
    lines = [
        f"sketch {sketch.name}: {lattice.dim}D, {len(params['vertices'])} vertices, "
        f"{len(params['blocks'])} blocks, {len(lattice.edges)} edges, "
        f"{len(lattice.families())} edge families"
    ]
    for v, cls in enumerate(placement.classes):
        shown = cls
        if v in params["fixed"]:
            shown = "pinned"
        elif v in held and cls != "fixed":
            shown = "stays where drawn (on a pinned edge or internal wall)"
        elif v in joints and v in placement.curve_of:
            shown = f"slides on {store.get(placement.curve_of[v]).name}"
        elif v in joints:
            shown = "slides on its boundary plane or surface"
        elif cls == "free":
            boundary = any(v in e for e in lattice.boundary_edges)
            if boundary:
                shown = "fixed where drawn (boundary, no geometry)"
            else:
                shown = "floats (interior)"
        p = placement.positions[v]
        lines.append(f"  vertex {v}: ({p[0]:.4g}, {p[1]:.4g}, {p[2]:.4g}) {shown}")
    dimensions, issues = _family_dimensions(store, params, lattice, placement)
    for _root, keys in sorted(lattice.families().items()):
        value = dimensions.get(keys[0], "?")
        lines.append(f"  family {', '.join(str(k) for k in keys)}: {value} points")
    for item in params.get("connectors", ()):
        connector = store.get(item["connector"])
        lines.append(
            f"  edge {item['edge']}: adopts connector {connector.name} "
            f"({connector.params['dimension']} points, shared with its other side)"
        )
    if params.get("unstructured"):
        listed = ", ".join(str(k) for k in params["unstructured"])
        lines.append(f"  unstructured blocks: {listed} (walls on them: inflation layers)")
    for issue in [*placement.issues, *issues]:
        lines.append(f"  ERROR {issue.where}: {issue.message}")
    for warning in placement.warnings:
        lines.append(f"  WARNING {warning}")
    for item in suggestions(store, sketch):
        lines.append(f"  suggestion: {item['what']} {item['message']}")
    return "\n".join(lines)


__all__ = [
    "KIND",
    "Lattice",
    "check_ownership",
    "create_sketch",
    "edit_sketch",
    "generate",
    "matching_connectors",
    "normalize_sketch",
    "owned_entities",
    "regenerate_adopters",
    "regenerate_sketch",
    "report",
    "suggestions",
]


def output_owner(store) -> dict:
    """{output entity id: (sketch entity, role)} for every output of every sketch."""
    owners = {}
    for sketch in store.entities("sketch"):
        for role, output in sketch.params["outputs"]:
            owners[output] = (sketch, role)
    return owners


def _item_edge(item):
    """The edge key an edge-keyed sketch item refers to (a wall dict, a (edge, value) pair
    or a bare edge), or None."""
    if isinstance(item, dict):
        return item.get("edge")
    if isinstance(item, tuple) and len(item) == 2:
        if isinstance(item[0], tuple):
            return item[0]
        if all(isinstance(v, int) for v in item):
            return item
    return None


def remove_blocks(store, name, blocks):
    """Blocks (indices) out of a 2D sketch, in one regeneration: the remaining blocks are
    renumbered (a domain's name follows its block index, so names can shift), the walls,
    pinned edges, dimensions, distributions, curves and adopted connectors of edges left in
    no block are dropped, the unstructured indices follow; the vertices stay. Refused for the
    last block (delete the sketch instead) and for 3D sketches. The GUI deletes a sketch's
    domain this way."""
    entity = store.get(name)
    if entity.kind != "sketch":
        raise ValidationError([Issue(name, f"is a {entity.kind}, not a sketch")])
    params = entity.params
    if isinstance(blocks, int):
        blocks = [blocks]
    issues: list[Issue] = []
    indices = set()
    for item in blocks:
        if isinstance(item, bool) or not isinstance(item, int):
            issues.append(Issue(name, f"block index must be an integer, got {item!r}"))
        elif not 0 <= item < len(params["blocks"]):
            issues.append(Issue(name, f"block {item} does not exist"))
        else:
            indices.add(item)
    if issues:
        raise ValidationError(issues)
    if any(len(block) != 4 for block in params["blocks"]):
        raise ValidationError([Issue(name, "blocks can only be removed from a 2D sketch")])
    if not indices:
        raise ValidationError([Issue(name, "no block given")])
    if len(indices) >= len(params["blocks"]):
        raise ValidationError([Issue(name, "would remove every block: delete the sketch instead")])
    keep = [k for k in range(len(params["blocks"])) if k not in indices]
    remap = {old: new for new, old in enumerate(keep)}
    new_blocks = [params["blocks"][k] for k in keep]
    edges = {edge_key(b[i], b[(i + 1) % 4]) for b in new_blocks for i in range(4)}
    changes = {
        "blocks": new_blocks,
        "unstructured": tuple(remap[b] for b in params["unstructured"] if b in remap),
    }
    for key in ("walls", "pinned", "dimensions", "distributions", "curves", "connectors"):
        changes[key] = tuple(
            item for item in params[key] if _item_edge(item) is None or _item_edge(item) in edges
        )
    return edit_sketch(store, name, **changes)
