"""Unstructured 2D domains: an outer loop and optional inner loops (holes) of connectors,
triangulated by the C++ meshers (opengrid._core.triangulate): constrained Delaunay with
refinement (algorithm "delaunay", the default) or advancing front ("advancing_front").

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

DEFAULT_GENERATOR = {
    "algorithm": "delaunay",
    "refine": True,
    "min_angle": 25.0,
    "growth_rate": 1.2,
    "max_size": None,
    "sources": (),
    "max_points": 2_000_000,
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


def normalize_generator(raw, where: str) -> dict:
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
    }


def normalize_unstructured(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"loops", "generator"}):
        issues.append(Issue(f"{where}: {key}", "unexpected key"))
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
        generator = normalize_generator(raw.get("generator"), where)
    except ValidationError as error:
        issues.extend(error.issues)
        generator = None
    if issues:
        raise ValidationError(issues)
    params = {"loops": tuple(loops), "generator": generator}
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


def check_loops(params, store, connectors, where) -> None:
    loops = loop_points(params, store, connectors, where)
    issues = []
    everything = np.concatenate(loops)
    if np.ptp(everything[:, 2]) > _tolerance(everything):
        raise ValidationError(
            [Issue(where, "unstructured domains must lie in a plane z = const (2D)")]
        )
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


def assemble_loops(store, connectors: ConnectorEvaluator, names: list[str]) -> list[list[dict]]:
    """Closed loops from a set of connectors: the loop enclosing the largest area becomes the
    outer loop (counter-clockwise), the others holes (clockwise)."""
    from opengrid.model.domains import _walk

    ids = []
    for name in names:
        entity = store.get(name)
        if entity.kind != "connector":
            raise ValidationError([Issue(name, f"is a {entity.kind}, not a connector")])
        ids.append(entity.id)
    ids = sorted(set(ids))
    nodes = {i: connectors.nodes(i).points for i in ids}
    keys = []

    def key_of(point):
        tolerance = _tolerance(np.asarray([point]))
        for index, known in enumerate(keys):
            if np.linalg.norm(known - point) <= tolerance:
                return index
        keys.append(np.asarray(point))
        return len(keys) - 1

    ends = {i: (key_of(nodes[i][0]), key_of(nodes[i][-1])) for i in ids}
    degree = {}
    for a, b in ends.values():
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    if any(d != 2 for d in degree.values()):
        raise ValidationError(
            [
                Issue(
                    "unstructured domain",
                    "connectors must form closed loops "
                    "(every end shared by exactly two connector ends)",
                )
            ]
        )
    remaining = list(ids)
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

    areas = [_signed_area(oriented(order)) for order in loops]
    outer = int(np.argmax([abs(a) for a in areas]))
    result = []
    for k in [outer] + [k for k in range(len(loops)) if k != outer]:
        order = loops[k]
        want_positive = k == outer
        if (areas[k] > 0) != want_positive:
            order = [(i, not r) for i, r in reversed(order)]
        result.append([{"connector": i, "reversed": r} for i, r in order])
    return result


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

    def mesh(self, entity_id: int) -> dict:
        entity = self.store.get(entity_id)
        generator = entity.params["generator"]
        points, segments = self.boundary(entity_id)
        digest = hashlib.sha256(repr(generator).encode())
        digest.update(points.tobytes())
        key = digest.hexdigest()
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
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
        try:
            xy, triangles, rejected, warnings = _core.triangulate(
                np.ascontiguousarray(points[:, :2]),
                segments,
                refine=generator["refine"],
                min_angle=generator["min_angle"],
                growth_rate=generator["growth_rate"],
                max_size=generator["max_size"],
                sources=sources,
                max_points=generator["max_points"],
                algorithm=generator["algorithm"],
            )
        except ValueError as error:
            raise ValidationError([Issue(entity.name, str(error))]) from None
        xyz = np.zeros((len(xy), 3))
        xyz[:, :2] = xy
        xyz[:, 2] = points[0, 2]
        # Boundary points are the connector nodes themselves (bitwise), so they merge exactly.
        xyz[: len(points)] = points
        result = {
            "points": xyz,
            "triangles": triangles,
            "rejected": int(rejected),
            "warnings": list(warnings),
        }
        self._cache[key] = result
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)
        return result
