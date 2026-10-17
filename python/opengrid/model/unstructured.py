"""Unstructured 2D domains: an outer loop and optional inner loops (holes) of connectors,
triangulated by the C++ meshers (opengrid._core.triangulate): constrained Delaunay with
refinement (algorithm "delaunay", the default) or advancing front ("advancing_front"), then
optionally recombined into quad-dominant or all-quad meshes, smoothed, edge-swapped and untangled
(opengrid._core.mesh_domain). Cells are numbered triangles first, then quads.

Boundary layers: ``generator["inflation"]`` (None or a dict) inflates quad layers from the given
wall connectors (part of the domain's loops) before the core is meshed: first height, growth
rate, at most ``layers`` layers or ``total_thickness``, and with ``transition`` columns stop when
cells become isotropic (core size).

Loops are closed chains of connectors (each possibly reversed). The outer loop runs
counter-clockwise and inner loops clockwise, i.e. the domain lies on the left of every loop.
Boundary nodes are exactly the connector nodes, so neighbouring domains (structured or
unstructured) sharing a connector are conforming.
"""

import hashlib
import math
from collections import OrderedDict

import numpy as np

from opengrid import _core
from opengrid.model.connectors import ConnectorEvaluator
from opengrid.model.domains import _edge_item, _signed_area, _tolerance, edge_points
from opengrid.validation import Issue, ValidationError

ALGORITHMS = ("delaunay", "advancing_front")
ELEMENT_TYPES = ("triangles", "quad_dominant", "quads")
SMOOTHING_METHODS = ("none", "laplacian", "length_weighted", "angle_based")

DEFAULT_GENERATOR = {
    "algorithm": "delaunay",
    "refine": True,
    "min_angle": 25.0,
    "growth_rate": 1.2,
    "max_size": None,
    "sources": (),
    "max_points": 2_000_000,
    "element_type": "triangles",
    "smoothing": "none",
    "smoothing_iterations": 5,
    "edge_swap": True,
    "inflation": None,
}

DEFAULT_INFLATION = {
    "connectors": (),
    "first_height": None,
    "growth_rate": 1.2,
    "layers": 30,
    "total_thickness": None,
    "transition": True,
}


def normalize_inflation(raw, label: str, store, issues: list) -> dict | None:
    """Canonical inflation settings; connectors (names or ids) become ids."""
    if raw is None:
        return None
    if not isinstance(raw, dict):
        issues.append(Issue(label, f"expected a dict or None, got {raw!r}"))
        return None
    for key in sorted(set(raw) - set(DEFAULT_INFLATION)):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    merged = dict(DEFAULT_INFLATION)
    merged.update(raw)
    connectors = []
    refs = merged["connectors"]
    if isinstance(refs, (str, int)):
        refs = [refs]
    if not isinstance(refs, (list, tuple)) or not refs:
        issues.append(Issue(f"{label}.connectors", "at least one wall connector is required"))
        refs = []
    for k, ref in enumerate(refs):
        where = f"{label}.connectors[{k}]"
        if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
            issues.append(Issue(where, f"no connector {ref!r}"))
            continue
        entity = store.get(ref)
        if entity.kind != "connector":
            issues.append(Issue(where, f"{entity.name!r} is a {entity.kind}, not a connector"))
            continue
        if entity.id not in connectors:
            connectors.append(entity.id)
    first = _number(merged["first_height"], f"{label}.first_height", issues)
    growth = _number(merged["growth_rate"], f"{label}.growth_rate", issues)
    if growth is not None and growth < 1.0:
        issues.append(Issue(f"{label}.growth_rate", "must be >= 1"))
    layers = merged["layers"]
    if isinstance(layers, bool) or not isinstance(layers, int) or layers < 1:
        issues.append(Issue(f"{label}.layers", f"must be an integer >= 1, got {layers!r}"))
    thickness = _number(
        merged["total_thickness"], f"{label}.total_thickness", issues, optional=True
    )
    transition = merged["transition"]
    if not isinstance(transition, bool):
        issues.append(Issue(f"{label}.transition", f"must be true or false, got {transition!r}"))
    return {
        "connectors": tuple(connectors),
        "first_height": first,
        "growth_rate": growth,
        "layers": layers,
        "total_thickness": thickness,
        "transition": transition,
    }


def _number(value, label, issues, positive=True, optional=False):
    if value is None and optional:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        issues.append(Issue(label, f"expected a finite number, got {value!r}"))
        return None
    if positive and value <= 0:
        issues.append(Issue(label, f"must be positive, got {value!r}"))
        return None
    return float(value)


def _point(value, label, issues):
    if isinstance(value, (str, bytes)) or not hasattr(value, "__len__") or len(value) not in (2, 3):
        issues.append(Issue(label, f"expected 2 or 3 coordinates, got {value!r}"))
        return None
    coords = [float(v) for v in value]
    if len(coords) == 2:
        coords.append(0.0)
    return tuple(coords)


def normalize_generator(raw, where: str, store=None) -> dict:
    issues: list[Issue] = []
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValidationError([Issue(f"{where}: generator", f"expected a dict, got {raw!r}")])
    for key in sorted(set(raw) - set(DEFAULT_GENERATOR)):
        issues.append(Issue(f"{where}: generator.{key}", "unexpected key"))
    merged = dict(DEFAULT_GENERATOR)
    merged.update(raw)
    label = f"{where}: generator"
    algorithm = merged["algorithm"]
    if algorithm not in ALGORITHMS:
        issues.append(
            Issue(
                f"{label}.algorithm", f"must be one of {', '.join(ALGORITHMS)}, got {algorithm!r}"
            )
        )
    element_type = merged["element_type"]
    if element_type not in ELEMENT_TYPES:
        issues.append(
            Issue(
                f"{label}.element_type",
                f"must be one of {', '.join(ELEMENT_TYPES)}, got {element_type!r}",
            )
        )
    smoothing = merged["smoothing"]
    if smoothing not in SMOOTHING_METHODS:
        issues.append(
            Issue(
                f"{label}.smoothing",
                f"must be one of {', '.join(SMOOTHING_METHODS)}, got {smoothing!r}",
            )
        )
    iterations = merged["smoothing_iterations"]
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 0:
        issues.append(
            Issue(f"{label}.smoothing_iterations", f"must be an integer >= 0, got {iterations!r}")
        )
    inflation = None
    if merged["inflation"] is not None:
        if store is None:
            issues.append(Issue(f"{label}.inflation", "needs the model to resolve connectors"))
        else:
            inflation = normalize_inflation(
                merged["inflation"], f"{label}.inflation", store, issues
            )
    edge_swap = merged["edge_swap"]
    if not isinstance(edge_swap, bool):
        issues.append(Issue(f"{label}.edge_swap", f"must be true or false, got {edge_swap!r}"))
    refine = merged["refine"]
    if not isinstance(refine, bool):
        issues.append(Issue(f"{label}.refine", f"must be true or false, got {refine!r}"))
    min_angle = _number(merged["min_angle"], f"{label}.min_angle", issues)
    if min_angle is not None and not 0 < min_angle < 34:
        issues.append(Issue(f"{label}.min_angle", "must lie in (0, 34) degrees"))
    growth = _number(merged["growth_rate"], f"{label}.growth_rate", issues)
    if growth is not None and growth <= 1.0:
        issues.append(Issue(f"{label}.growth_rate", "must be greater than 1"))
    max_size = _number(merged["max_size"], f"{label}.max_size", issues, optional=True)
    max_points = merged["max_points"]
    if isinstance(max_points, bool) or not isinstance(max_points, int) or max_points < 3:
        issues.append(Issue(f"{label}.max_points", f"must be an integer >= 3, got {max_points!r}"))
    sources = []
    for index, item in enumerate(merged["sources"]):
        sl = f"{label}.sources[{index}]"
        if not isinstance(item, dict) or "p0" not in item or "size" not in item:
            issues.append(Issue(sl, "needs 'p0', 'size' (and optionally 'p1', 'radius')"))
            continue
        for key in sorted(set(item) - {"p0", "p1", "size", "radius"}):
            issues.append(Issue(f"{sl}.{key}", "unexpected key"))
        p0 = _point(item["p0"], f"{sl}.p0", issues)
        p1 = _point(item.get("p1", item["p0"]), f"{sl}.p1", issues)
        size = _number(item["size"], f"{sl}.size", issues)
        radius = _number(item.get("radius", 0.0), f"{sl}.radius", issues, positive=False)
        if radius is not None and radius < 0:
            issues.append(Issue(f"{sl}.radius", "must be >= 0"))
        sources.append({"p0": p0, "p1": p1, "size": size, "radius": radius})
    if issues:
        raise ValidationError(issues)
    return {
        "algorithm": algorithm,
        "refine": refine,
        "min_angle": min_angle,
        "growth_rate": growth,
        "max_size": max_size,
        "sources": tuple(sources),
        "max_points": max_points,
        "element_type": element_type,
        "smoothing": smoothing,
        "smoothing_iterations": iterations,
        "edge_swap": edge_swap,
        "inflation": inflation,
    }


def normalize_unstructured(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"loops", "generator", "surface"}):
        issues.append(Issue(f"{where}: {key}", "unexpected key"))
    surface = raw.get("surface")
    if surface is not None:
        from opengrid.model.entities import SURFACE_KINDS

        if isinstance(surface, bool) or not isinstance(surface, (str, int)) or surface not in store:
            issues.append(Issue(f"{where}: surface", f"no surface {surface!r}"))
            surface = None
        elif store.get(surface).kind not in SURFACE_KINDS:
            entity = store.get(surface)
            issues.append(Issue(f"{where}: surface", f"{entity.name!r} is a {entity.kind}"))
            surface = None
        else:
            surface = store.get(surface).id
    loops_raw = raw.get("loops")
    if not isinstance(loops_raw, (list, tuple)) or not loops_raw:
        raise ValidationError(
            [*issues, Issue(f"{where}: loops", "at least one (outer) loop is required")]
        )
    loops = []
    for k, loop_raw in enumerate(loops_raw):
        if isinstance(loop_raw, (str, int, dict)):
            loop_raw = [loop_raw]
        items = []
        for m, item in enumerate(loop_raw):
            items.append(_edge_item(item, f"{where}: loops[{k}][{m}]", store, issues))
        loops.append(tuple(items))
    try:
        generator = normalize_generator(raw.get("generator"), where, store)
    except ValidationError as error:
        issues.extend(error.issues)
        generator = None
    if issues:
        raise ValidationError(issues)
    if generator["inflation"] is not None:
        in_loops = {item["connector"] for loop in loops for item in loop}
        for ref in generator["inflation"]["connectors"]:
            if ref not in in_loops:
                issues.append(
                    Issue(
                        f"{where}: generator.inflation.connectors",
                        f"{store.get(ref).name!r} is not part of the domain's loops",
                    )
                )
        if issues:
            raise ValidationError(issues)
    params = {"loops": tuple(loops), "generator": generator, "surface": surface}
    check_loops(params, store, ConnectorEvaluator(store), where)
    return params


def loop_points(params, store, connectors, where) -> list[np.ndarray]:
    """Node arrays of the loops (closing node not repeated)."""
    issues: list[Issue] = []
    arrays = []
    for k, loop in enumerate(params["loops"]):
        points = edge_points(loop, connectors, store, f"{where}: loops[{k}]", issues)
        arrays.append(points)
    if issues:
        raise ValidationError(issues)
    result = []
    for k, points in enumerate(arrays):
        tolerance = _tolerance(points)
        gap = float(np.linalg.norm(points[-1] - points[0]))
        if gap > tolerance:
            issues.append(Issue(f"{where}: loops[{k}]", f"is not closed (gap {gap:.3g})"))
            continue
        result.append(points[:-1])
    if issues:
        raise ValidationError(issues)
    return result


def domain_frame(params, store, loops, where: str = ""):
    """2D frame of an unstructured domain (opengrid.model.frames)."""
    from opengrid.model.frames import make_frame

    surface = None
    entity = None
    if params.get("surface") is not None:
        entity = store.get(params["surface"])
        surface = store.geometry(entity.name)
    return make_frame(loops, surface, where, surface_entity=entity, store=store)


def check_loops(params, store, connectors, where) -> None:
    loops3 = loop_points(params, store, connectors, where)
    frame = domain_frame(params, store, loops3, where)
    loops = [frame.to2d(loop) for loop in loops3]
    issues = []
    for k, points in enumerate(loops):
        if len(points) < 3:
            issues.append(Issue(f"{where}: loops[{k}]", "needs at least 3 nodes"))
            continue
        area = _signed_area(points)
        if k == 0 and area <= 0:
            issues.append(Issue(f"{where}: loops[0]", "the outer loop must run counter-clockwise"))
        if k > 0 and area >= 0:
            issues.append(Issue(f"{where}: loops[{k}]", "inner loops (holes) must run clockwise"))
        if k > 0 and not _inside(points[0], loops[0]):
            issues.append(Issue(f"{where}: loops[{k}]", "lies outside the outer loop"))
    if issues:
        raise ValidationError(issues)


def _inside(point, polygon) -> bool:
    x = polygon[:, 0]
    y = polygon[:, 1]
    xn = np.roll(x, -1)
    yn = np.roll(y, -1)
    crosses = (y > point[1]) != (yn > point[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        at = x + (point[1] - y) * (xn - x) / (yn - y)
    return bool(np.count_nonzero(crosses & (point[0] < at)) % 2)


def assemble_loops(
    store, connectors: ConnectorEvaluator, names: list[str], surface=None
) -> list[list[dict]]:
    """Closed loops from a set of connectors: the loop enclosing the largest area becomes the
    outer loop (counter-clockwise), the others holes (clockwise)."""
    ids = []
    for name in names:
        entity = store.get(name)
        if entity.kind != "connector":
            raise ValidationError([Issue(name, f"is a {entity.kind}, not a connector")])
        ids.append(entity.id)
    ids = sorted(set(ids))
    polylines = {i: connectors.nodes(i).points for i in ids}
    labels = {i: store.get(i).name for i in ids}
    surface_geometry = None
    if surface is not None:
        surface_geometry = store.geometry(surface)
    return assemble_polylines(polylines, labels, surface_geometry)


def assemble_polylines(polylines: dict, labels: dict, surface_geometry=None) -> list[list[dict]]:
    """Closed loops from polylines keyed by connector (an id, or the name of a connector about
    to be created): ``[{"connector": key, "reversed": bool}, ...]`` per loop, the loop
    enclosing the largest area first and counter-clockwise, the others (holes) clockwise.
    Ends coincide within the chain tolerance; every end must meet exactly one other end, else
    the error names the connectors (``labels``) and where their open ends are."""
    from opengrid.model.domains import _walk

    keys = []

    def key_of(point):
        tolerance = _tolerance(np.asarray([point]))
        for index, known in enumerate(keys):
            if np.linalg.norm(known - point) <= tolerance:
                return index
        keys.append(np.asarray(point, dtype=float))
        return len(keys) - 1

    nodes = {i: np.asarray(points, dtype=float) for i, points in polylines.items()}
    ends = {i: (key_of(nodes[i][0]), key_of(nodes[i][-1])) for i in nodes}
    degree = {}
    for a, b in ends.values():
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    issues = []
    for i in nodes:
        a, b = ends[i]
        if degree[a] == 1 and degree[b] == 1 and a != b:
            # Both ends loose and close together: a loop that does not quite close.
            gap = float(np.linalg.norm(keys[a] - keys[b]))
            length = float(np.sum(np.linalg.norm(np.diff(nodes[i], axis=0), axis=1)))
            if gap <= 0.1 * length:
                issues.append(
                    Issue(
                        str(labels[i]),
                        f"its two ends are {gap:.3g} apart: close the loop (drag one end onto "
                        "the other, or release a freehand stroke at its start)",
                    )
                )
                continue
        for which, end in (("start", a), ("end", b)):
            if degree[end] == 2:
                continue
            x, y, z = keys[end]
            issues.append(
                Issue(
                    str(labels[i]),
                    f"its {which} at ({x:.4g}, {y:.4g}, {z:.4g}) meets {degree[end] - 1} "
                    "other end(s), not exactly one",
                )
            )
    if issues:
        issues.insert(
            0,
            Issue(
                "unstructured domain",
                "connectors must form closed loops (every end shared by exactly two "
                "connector ends)",
            ),
        )
        raise ValidationError(issues)
    remaining = list(nodes)
    loops = []
    while remaining:
        order = _walk(remaining, ends, remaining[0])
        loops.append(order)
        used = {i for i, _r in order}
        remaining = [i for i in remaining if i not in used]

    def oriented(order):
        parts = []
        for i, reversed_ in order:
            points = nodes[i]
            if reversed_:
                points = points[::-1]
            parts.append(points[:-1])
        return np.concatenate(parts)

    from opengrid.model.frames import make_frame

    frame = make_frame(
        [oriented(order) for order in loops], surface_geometry, "unstructured domain"
    )
    areas = [_signed_area(frame.to2d(oriented(order))) for order in loops]
    outer = int(np.argmax([abs(a) for a in areas]))
    result = []
    for k in [outer] + [k for k in range(len(loops)) if k != outer]:
        order = loops[k]
        want_positive = k == outer
        if (areas[k] > 0) != want_positive:
            order = [(i, not r) for i, r in reversed(order)]
        result.append([{"connector": i, "reversed": r} for i, r in order])
    return result


def _sources(generator) -> list[tuple]:
    sources = []
    for item in generator["sources"]:
        sources.append(
            (
                item["p0"][0],
                item["p0"][1],
                item["p1"][0],
                item["p1"][1],
                item["size"],
                item["radius"],
            )
        )
    return sources


class UnstructuredEvaluator:
    """Triangulations of unstructured domains, cached by a key over the generator settings and
    the boundary nodes."""

    CACHE_SIZE = 32

    def __init__(self, store, connectors: ConnectorEvaluator):
        self.store = store
        self.connectors = connectors
        self._cache: OrderedDict[str, dict] = OrderedDict()

    def invalidate(self) -> None:
        """Nothing to drop: the cache key covers every input."""

    def wall_flags(self, entity_id: int) -> list[bool]:
        """One flag per boundary segment (loop order): True on inflation wall connectors."""
        entity = self.store.get(entity_id)
        inflation = entity.params["generator"]["inflation"]
        walls = set()
        if inflation is not None:
            walls = set(inflation["connectors"])
        flags = []
        for loop in entity.params["loops"]:
            for item in loop:
                count = len(self.connectors.nodes(item["connector"]).points) - 1
                flags.extend([item["connector"] in walls] * count)
        return flags

    def boundary(self, entity_id: int):
        entity = self.store.get(entity_id)
        loops = loop_points(entity.params, self.store, self.connectors, entity.name)
        points = np.concatenate(loops)
        segments = []
        start = 0
        for loop in loops:
            n = len(loop)
            index = np.arange(start, start + n)
            segments.append(np.stack([index, np.roll(index, -1)], axis=1))
            start += n
        return points, np.concatenate(segments).astype(np.int64)

    def remesh_planar(self, entity_id: int, points, segments) -> tuple[np.ndarray, np.ndarray]:
        """Triangles filling the region of a planar domain bounded by other loops (``points``
        (n, 3) and ``segments`` (m, 2); e.g. around inflation layers that slide on the domain), in
        the domain's own frame and with its generator settings (triangles only, no inflation).
        Returns (points: the given points first, unchanged, then interior points; triangles)."""
        entity = self.store.get(entity_id)
        generator = entity.params["generator"]
        loops3 = loop_points(entity.params, self.store, self.connectors, entity.name)
        frame = domain_frame(entity.params, self.store, loops3, entity.name)
        points = np.asarray(points, dtype=float)
        uv = frame.to2d(points)
        segments = np.asarray(segments, dtype=np.int64).reshape(-1, 2)
        # Outer loops counter-clockwise, holes clockwise: the loops' signed areas add up to the
        # (positive) area of the region.
        a = uv[segments[:, 0]]
        b = uv[segments[:, 1]]
        if np.sum(a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1]) < 0:
            segments = segments[:, ::-1]
        try:
            xy, triangles, _quads, _rejected, _warnings, _layers = _core.mesh_domain(
                uv,
                np.ascontiguousarray(segments),
                algorithm=generator["algorithm"],
                element_type="triangles",
                refine=generator["refine"],
                min_angle=generator["min_angle"],
                growth_rate=generator["growth_rate"],
                max_size=generator["max_size"],
                sources=_sources(generator),
                max_points=generator["max_points"],
                smoothing=generator["smoothing"],
                smoothing_iterations=generator["smoothing_iterations"],
                edge_swap=generator["edge_swap"],
                wall=[],
            )
        except ValueError as error:
            raise ValidationError(
                [Issue(entity.name, f"re-meshing around layers: {error}")]
            ) from None
        xyz = np.zeros((len(xy), 3))
        if len(xy) > len(points):
            xyz[len(points) :] = frame.to3d(xy[len(points) :])
        xyz[: len(points)] = points
        return xyz, np.asarray(triangles, dtype=np.int64)

    def mesh(self, entity_id: int) -> dict:
        entity = self.store.get(entity_id)
        generator = entity.params["generator"]
        points, segments = self.boundary(entity_id)
        digest = hashlib.sha256(repr(generator).encode())
        digest.update(points.tobytes())
        if entity.params.get("surface") is not None:
            from opengrid.model.domains import surface_signature

            digest.update(repr(surface_signature(self.store, entity.params["surface"])).encode())
        key = digest.hexdigest()
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        sources = _sources(generator)
        try:
            inflation = generator["inflation"]
            wall = []
            options = {}
            if inflation is not None:
                wall = self.wall_flags(entity_id)
                options = {
                    "first_height": inflation["first_height"],
                    "inflation_growth": inflation["growth_rate"],
                    "layers": inflation["layers"],
                    "total_thickness": inflation["total_thickness"],
                    "transition": inflation["transition"],
                }
            loops3 = loop_points(entity.params, self.store, self.connectors, entity.name)
            frame = domain_frame(entity.params, self.store, loops3, entity.name)
            uv = frame.to2d(points)
            # A loop node where the surface closes on itself (a pole, an apex) is opened into
            # the degenerate side it stands for, so the triangles there fan around it.
            uv, points_in, segments, wall, origin = _open_poles(
                frame, loops3, points, uv, segments, wall
            )
            xy, triangles, quads, rejected, warnings, layers = _core.mesh_domain(
                uv,
                segments,
                algorithm=generator["algorithm"],
                element_type=generator["element_type"],
                refine=generator["refine"],
                min_angle=generator["min_angle"],
                growth_rate=generator["growth_rate"],
                max_size=generator["max_size"],
                sources=sources,
                max_points=generator["max_points"],
                smoothing=generator["smoothing"],
                smoothing_iterations=generator["smoothing_iterations"],
                edge_swap=generator["edge_swap"],
                wall=wall,
                **options,
            )
        except ValueError as error:
            raise ValidationError([Issue(entity.name, str(error))]) from None
        xyz = np.zeros((len(xy), 3))
        if len(xy) > len(points_in):
            xyz[len(points_in) :] = frame.to3d(xy[len(points_in) :])
        # Boundary points are the connector nodes themselves (bitwise), so they merge exactly.
        xyz[: len(points_in)] = points_in
        if origin is not None:
            xyz, triangles, quads = _close_poles(xyz, triangles, quads, origin, len(points))
            if len(quads) == 0:
                target = generator["max_size"]
                if target is None:
                    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
                    target = float(np.median(lengths[lengths > 0]))
                xyz, triangles = _improve_on_surface(
                    xyz, triangles, len(points), frame.surface, target
                )
        result = {
            "points": xyz,
            "triangles": triangles,
            "quads": quads,
            "layers": int(layers),
            "rejected": int(rejected),
            "warnings": list(warnings),
        }
        self._cache[key] = result
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)
        return result


def _degenerate_direction(surface, u: float, v: float, point, tolerance: float):
    """ "u" when the whole u line through (u, v) collapses onto ``point`` (the surface closes
    there along u, e.g. a pole of a surface of revolution), "v" likewise, else None."""
    target = np.asarray(point, dtype=float)
    samples = (0.0, 0.25, 0.5, 0.75, 1.0)
    along_v = max(float(np.linalg.norm(np.asarray(surface.eval(u, t)) - target)) for t in samples)
    if along_v <= tolerance:
        return "v"
    along_u = max(float(np.linalg.norm(np.asarray(surface.eval(t, v)) - target)) for t in samples)
    if along_u <= tolerance:
        return "u"
    return None


def _open_poles(frame, loops3, points, uv, segments, wall):
    """Loop nodes at a pole of the surface (its parameter line collapses there) become two
    2D nodes joined along the degenerate side, taking the neighbours' other coordinate, so the
    2D polygon is the surface's full parameter region. Returns (uv, points, segments, wall,
    origin): the expanded arrays and, when anything was opened, the original boundary index
    of every expanded boundary node (else None)."""
    if getattr(frame, "kind", "") != "surface":
        return uv, points, segments, wall, None
    tolerance = _tolerance(points)
    poles = {}  # boundary index -> "u" | "v"
    for k, (x, y) in enumerate(uv):
        u = x / frame.su
        v = frame.flip_v * y / frame.sv
        direction = _degenerate_direction(frame.surface, u, v, points[k], tolerance)
        if direction is not None:
            poles[k] = direction
    if not poles:
        return uv, points, segments, wall, None
    flags = list(wall)
    if not flags:
        flags = [False] * len(segments)
    new_uv = []
    new_points = []
    new_segments = []
    new_wall = []
    origin = []
    start = 0
    for loop in loops3:
        n = len(loop)
        first = len(new_uv)
        for j in range(n):
            k = start + j
            previous = start + (j - 1) % n
            following = start + (j + 1) % n
            if k not in poles:
                new_uv.append(tuple(uv[k]))
                new_points.append(points[k])
                origin.append(k)
                new_wall.append(flags[k])
                continue
            if poles[k] == "v":
                before = (uv[k][0], uv[previous][1])
                after = (uv[k][0], uv[following][1])
            else:
                before = (uv[previous][0], uv[k][1])
                after = (uv[following][0], uv[k][1])
            new_uv.append(before)
            new_points.append(points[k])
            origin.append(k)
            new_wall.append(False)  # the degenerate side carries no wall
            new_uv.append(after)
            new_points.append(points[k])
            origin.append(k)
            new_wall.append(flags[k])
        count = len(new_uv) - first
        index = np.arange(first, first + count)
        new_segments.append(np.stack([index, np.roll(index, -1)], axis=1))
        start += n
    uv_out = np.ascontiguousarray(np.asarray(new_uv, dtype=float))
    points_out = np.asarray(new_points, dtype=float)
    segments_out = np.concatenate(new_segments).astype(np.int64)
    if not wall:
        new_wall = []
    return uv_out, points_out, segments_out, new_wall, np.asarray(origin, dtype=np.int64)


def _close_poles(xyz, triangles, quads, origin, n_original: int):
    """Undo the pole opening on the mesh: the duplicated boundary nodes fold back onto their
    original node (the first ``n_original`` nodes are the loop nodes again, in loop order),
    the interior nodes follow, cells left with a repeated node go (quads with one become
    triangles)."""
    n_expanded = len(origin)
    remap = np.empty(len(xyz), dtype=np.int64)
    remap[:n_expanded] = origin
    remap[n_expanded:] = np.arange(n_original, n_original + len(xyz) - n_expanded)
    points = np.concatenate([np.zeros((n_original, 3)), xyz[n_expanded:]])
    for j in range(n_expanded):
        points[origin[j]] = xyz[j]
    tris = []
    for cell in np.asarray(triangles, dtype=np.int64).reshape(-1, 3):
        mapped = remap[cell]
        if len(set(mapped.tolist())) == 3:
            tris.append(mapped)
    quad_cells = []
    for cell in np.asarray(quads, dtype=np.int64).reshape(-1, 4):
        mapped = remap[cell]
        unique = []
        for index in mapped.tolist():
            if index not in unique:
                unique.append(index)
        if len(unique) == 4:
            quad_cells.append(mapped)
        elif len(unique) == 3:
            tris.append(np.asarray(unique, dtype=np.int64))
    out_tris = np.asarray(tris, dtype=np.int64).reshape(-1, 3)
    out_quads = np.asarray(quad_cells, dtype=np.int64).reshape(-1, 4)
    return points, out_tris, out_quads


def _triangle_normals(points, triangles):
    a = points[triangles[:, 0]]
    b = points[triangles[:, 1]]
    c = points[triangles[:, 2]]
    return np.cross(b - a, c - a)


def _min_angles(points, triangles):
    """The smallest angle (degrees) of each triangle."""
    out = np.full(len(triangles), 180.0)
    for k in range(3):
        p = points[triangles[:, k]]
        q = points[triangles[:, (k + 1) % 3]]
        r = points[triangles[:, (k + 2) % 3]]
        u = q - p
        v = r - p
        cosine = np.einsum("ij,ij->i", u, v) / (
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) + 1e-300
        )
        out = np.minimum(out, np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
    return out


def _improve_on_surface(points, triangles, n_boundary: int, surface, target: float):
    """Clean up a triangle mesh on a surface whose parametrisation is distorted (a pole):
    interior edges much shorter than ``target`` are collapsed, interior nodes smoothed and put
    back on the surface, interior edges flipped where that opens the angles; the boundary
    nodes (the first ``n_boundary``) never move. Returns (points, triangles)."""
    points = np.array(points, dtype=float)
    triangles = np.array(triangles, dtype=np.int64).reshape(-1, 3)
    if len(triangles) == 0:
        return points, triangles
    alive = np.ones(len(points), dtype=bool)
    for _pass in range(12):
        normals = _triangle_normals(points, triangles)
        edges = np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
        lengths = np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1)
        order = np.argsort(lengths)
        touched = np.zeros(len(points), dtype=bool)
        collapsed = 0
        for index in order:
            if lengths[index] >= 0.55 * target:
                break
            a, b = int(edges[index, 0]), int(edges[index, 1])
            if touched[a] or touched[b]:
                continue
            if a < n_boundary and b < n_boundary:
                continue
            if a < n_boundary:
                a, b = b, a  # the interior node a folds onto b
            uses_a = np.any(triangles == a, axis=1)
            uses_b = np.any(triangles == b, axis=1)
            # The link condition: the nodes next to both a and b must be exactly the apexes
            # of the edge's triangles, else the collapse folds the mesh (a duplicate triangle).
            link_a = set(np.unique(triangles[uses_a]).tolist()) - {a, b}
            link_b = set(np.unique(triangles[uses_b]).tolist()) - {a, b}
            apexes = set(np.unique(triangles[uses_a & uses_b]).tolist()) - {a, b}
            if link_a & link_b != apexes:
                continue
            trial = triangles[uses_a].copy()
            trial[trial == a] = b
            keep = np.array([len(set(row.tolist())) == 3 for row in trial])
            kept = trial[keep]
            if len(kept) == 0:
                continue
            new_normals = _triangle_normals(points, kept)
            old_normals = normals[uses_a][keep]
            if np.any(np.einsum("ij,ij->i", new_normals, old_normals) <= 0.0):
                continue
            triangles = np.concatenate([triangles[~uses_a], kept])
            normals = _triangle_normals(points, triangles)
            alive[a] = False
            touched[a] = True
            touched[b] = True
            collapsed += 1
        if collapsed == 0:
            break
    for _pass in range(4):
        # Laplacian smoothing of the interior nodes, back onto the surface, kept only when no
        # triangle turns over.
        neighbours = [set() for _ in range(len(points))]
        for tri in triangles:
            for k in range(3):
                neighbours[tri[k]].add(int(tri[(k + 1) % 3]))
                neighbours[tri[k]].add(int(tri[(k + 2) % 3]))
        normals = _triangle_normals(points, triangles)
        for node in range(n_boundary, len(points)):
            if not alive[node] or not neighbours[node]:
                continue
            centroid = points[list(neighbours[node])].mean(axis=0)
            _u, _v, projected, _d = surface.project(tuple(float(c) for c in centroid))
            old = points[node].copy()
            points[node] = projected
            uses = np.any(triangles == node, axis=1)
            new_normals = _triangle_normals(points, triangles[uses])
            if np.any(np.einsum("ij,ij->i", new_normals, normals[uses]) <= 0.0):
                points[node] = old
            else:
                normals[uses] = new_normals
        # Edge flips where the two triangles' smallest angle grows.
        edge_triangles: dict = {}
        for t_index, tri in enumerate(triangles):
            for k in range(3):
                key = (min(tri[k], tri[(k + 1) % 3]), max(tri[k], tri[(k + 1) % 3]))
                edge_triangles.setdefault(key, []).append(t_index)
        flipped = 0
        touched_triangles: set = set()
        existing = set(edge_triangles)
        for (a, b), owners in list(edge_triangles.items()):
            if len(owners) != 2 or (a, b) not in existing:
                continue
            t0, t1 = owners
            if t0 in touched_triangles or t1 in touched_triangles:
                continue
            c = [v for v in triangles[t0].tolist() if v not in (a, b)]
            d = [v for v in triangles[t1].tolist() if v not in (a, b)]
            if len(c) != 1 or len(d) != 1:
                continue
            c, d = c[0], d[0]
            if (min(c, d), max(c, d)) in existing:
                continue  # the flipped edge exists already
            before = np.asarray([triangles[t0], triangles[t1]])
            after = np.asarray([[a, c, d], [b, d, c]], dtype=np.int64)
            reference = _triangle_normals(points, before).mean(axis=0)
            new_normals = _triangle_normals(points, after)
            if np.dot(new_normals[0], reference) < 0:
                after = after[:, ::-1]
                new_normals = _triangle_normals(points, after)
            if np.any(new_normals @ reference <= 0.0):
                continue
            if _min_angles(points, after).min() <= _min_angles(points, before).min() + 1e-9:
                continue
            triangles[t0] = after[0]
            triangles[t1] = after[1]
            existing.discard((a, b))
            existing.add((min(c, d), max(c, d)))
            touched_triangles.update((t0, t1))
            flipped += 1
        if flipped == 0 and _pass > 0:
            break
    used = np.zeros(len(points), dtype=bool)
    used[:n_boundary] = True
    used[np.unique(triangles)] = True
    index = np.cumsum(used) - 1
    return points[used], index[triangles]
