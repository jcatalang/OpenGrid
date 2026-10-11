"""Structured domains: edge loops of connectors, validation, automatic edge assembly and
grid generation (TFI / Hermite TFI, optional elliptic smoothing).

Planar domains in a plane z = const support every generator. Other domains (faces of 3D
blocks) use linear TFI in space, or, with ``surface`` set, parametric TFI on that surface (the
boundary nodes are projected onto it and the interior lies on it).

Edges are listed as a counter-clockwise loop: edge 0 runs along j = 0 with i increasing,
edge 1 along i = imax with j increasing, edge 2 along j = jmax with i *decreasing* and edge 3
along i = 0 with j decreasing. Each edge is a chain of connectors, each possibly reversed.
Opposite edges must have equal point counts. Grids are arrays of shape (nj, ni, 3).
"""

import hashlib
import itertools
import math
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.model.connectors import CHAIN_TOLERANCE, ConnectorEvaluator
from opengrid.validation import Issue, ValidationError

EDGE_NAMES = ("edge 0 (j = 0)", "edge 1 (i = imax)", "edge 2 (j = jmax)", "edge 3 (i = 0)")
METHODS = ("tfi", "hermite_tfi")
CONTROLS = ("laplace", "thomas_middlecoff", "ttm")
SOLVERS = ("point_sor", "line_sor")
WALL_EDGES = (0, 2)  # walls with spacing/orthogonality control: j = 0 and j = jmax
ATTRACTION_KINDS = ("xi_line", "eta_line", "point")

DEFAULT_ELLIPTIC = {
    "control": "thomas_middlecoff",
    "attractions": (),
    "walls": (),
    "solver": "line_sor",
    "omega": 1.0,
    "max_iterations": 2000,
    "tolerance": 1e-6,
}
DEFAULT_WALL = {"edge": 0, "spacing": None, "decay": 0.1, "relaxation": 0.1, "max_source": 100.0}
DEFAULT_GENERATOR = {
    "method": "tfi",
    "hermite_walls": (0,),
    "hermite_strength": 1.0,
    "elliptic": None,
}


# ----------------------------------------------------------------------------- parameters
def _float(value, where, name, issues, positive=False, optional=False):
    if value is None and optional:
        return None
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        value = value.item()
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        issues.append(Issue(f"{where}: {name}", f"expected a finite number, got {value!r}"))
        return None
    if positive and value <= 0:
        issues.append(Issue(f"{where}: {name}", f"must be positive, got {value!r}"))
        return None
    return float(value)


def _unknown_keys(raw: dict, allowed, where: str, issues: list[Issue]) -> None:
    for key in sorted(set(raw) - set(allowed)):
        issues.append(Issue(f"{where}: {key}", "unexpected key"))


def normalize_generator(raw, where: str) -> dict:
    issues: list[Issue] = []
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValidationError([Issue(f"{where}: generator", f"expected a dict, got {raw!r}")])
    _unknown_keys(raw, DEFAULT_GENERATOR, f"{where}: generator", issues)
    merged = dict(DEFAULT_GENERATOR)
    merged.update(raw)
    if merged["method"] not in METHODS:
        issues.append(Issue(f"{where}: generator.method", f"must be one of {', '.join(METHODS)}"))
    walls = tuple(merged["hermite_walls"])
    for edge in walls:
        if edge not in WALL_EDGES:
            issues.append(
                Issue(f"{where}: generator.hermite_walls", f"edges must be 0 or 2, got {edge!r}")
            )
    strength = _float(
        merged["hermite_strength"], where, "generator.hermite_strength", issues, positive=True
    )
    elliptic = merged["elliptic"]
    if elliptic is not None:
        elliptic = _normalize_elliptic(elliptic, where, issues)
    if issues:
        raise ValidationError(issues)
    return {
        "method": merged["method"],
        "hermite_walls": tuple(sorted(set(walls))),
        "hermite_strength": strength,
        "elliptic": elliptic,
    }


def _normalize_elliptic(raw, where: str, issues: list[Issue]) -> dict | None:
    label = f"{where}: generator.elliptic"
    if raw is True:
        raw = {}
    if not isinstance(raw, dict):
        issues.append(Issue(label, f"expected a dict, None or True, got {raw!r}"))
        return None
    _unknown_keys(raw, DEFAULT_ELLIPTIC, label, issues)
    merged = dict(DEFAULT_ELLIPTIC)
    merged.update(raw)
    if merged["control"] not in CONTROLS:
        issues.append(Issue(f"{label}.control", f"must be one of {', '.join(CONTROLS)}"))
    if merged["solver"] not in SOLVERS:
        issues.append(Issue(f"{label}.solver", f"must be one of {', '.join(SOLVERS)}"))
    omega = _float(merged["omega"], label, "omega", issues, positive=True)
    if omega is not None and omega >= 2.0:
        issues.append(Issue(f"{label}.omega", "must be below 2"))
    iterations = merged["max_iterations"]
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 1:
        issues.append(
            Issue(f"{label}.max_iterations", f"must be a positive integer, got {iterations!r}")
        )
    tolerance = _float(merged["tolerance"], label, "tolerance", issues, positive=True)
    walls = []
    seen = set()
    for index, item in enumerate(merged["walls"]):
        wl = f"{label}.walls[{index}]"
        if isinstance(item, int) and not isinstance(item, bool):
            item = {"edge": item}
        if not isinstance(item, dict):
            issues.append(Issue(wl, f"expected an edge number or a dict, got {item!r}"))
            continue
        _unknown_keys(item, DEFAULT_WALL, wl, issues)
        wall = dict(DEFAULT_WALL)
        wall.update(item)
        if wall["edge"] not in WALL_EDGES:
            issues.append(Issue(wl, f"edge must be 0 or 2, got {wall['edge']!r}"))
            continue
        if wall["edge"] in seen:
            issues.append(Issue(wl, f"edge {wall['edge']} listed twice"))
            continue
        seen.add(wall["edge"])
        walls.append(
            {
                "edge": wall["edge"],
                "spacing": _float(
                    wall["spacing"], wl, "spacing", issues, positive=True, optional=True
                ),
                "decay": _float(wall["decay"], wl, "decay", issues, positive=True),
                "relaxation": _float(wall["relaxation"], wl, "relaxation", issues, positive=True),
                "max_source": _float(wall["max_source"], wl, "max_source", issues, positive=True),
            }
        )
    attractions = []
    for index, item in enumerate(merged["attractions"]):
        al = f"{label}.attractions[{index}]"
        if not isinstance(item, dict) or item.get("kind") not in ATTRACTION_KINDS:
            issues.append(Issue(al, f"needs 'kind' in {', '.join(ATTRACTION_KINDS)}"))
            continue
        _unknown_keys(item, ("kind", "xi", "eta", "amplitude", "decay"), al, issues)
        kind = item["kind"]
        xi = None
        eta = None
        if kind in ("xi_line", "point"):
            xi = _float(item.get("xi"), al, "xi", issues)
        if kind in ("eta_line", "point"):
            eta = _float(item.get("eta"), al, "eta", issues)
        attractions.append(
            {
                "kind": kind,
                "xi": xi,
                "eta": eta,
                "amplitude": _float(item.get("amplitude"), al, "amplitude", issues),
                "decay": _float(item.get("decay", 0.5), al, "decay", issues, positive=True),
            }
        )
    walls.sort(key=lambda w: w["edge"])
    return {
        "control": merged["control"],
        "attractions": tuple(attractions),
        "walls": tuple(walls),
        "solver": merged["solver"],
        "omega": omega,
        "max_iterations": iterations,
        "tolerance": tolerance,
    }


def _edge_item(raw, where: str, store, issues: list[Issue]) -> dict | None:
    reversed_ = False
    if isinstance(raw, dict):
        ref = raw.get("connector")
        reversed_ = raw.get("reversed", False)
    elif isinstance(raw, (list, tuple)) and len(raw) == 2:
        ref, reversed_ = raw
    else:
        ref = raw
    if isinstance(ref, bool) or not isinstance(ref, (str, int)):
        issues.append(Issue(where, f"expected a connector name or id, got {ref!r}"))
        return None
    if ref not in store:
        issues.append(Issue(where, f"no entity {ref!r}"))
        return None
    entity = store.get(ref)
    if entity.kind != "connector":
        issues.append(Issue(where, f"{entity.name!r} is a {entity.kind}, not a connector"))
        return None
    if not isinstance(reversed_, bool):
        issues.append(Issue(where, f"'reversed' must be true or false, got {reversed_!r}"))
        return None
    return {"connector": entity.id, "reversed": reversed_}


def normalize_domain(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    _unknown_keys(raw, ("edges", "generator", "surface"), where, issues)
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
    edges_raw = raw.get("edges")
    if not isinstance(edges_raw, (list, tuple)) or len(edges_raw) != 4:
        raise ValidationError([*issues, Issue(f"{where}: edges", "exactly 4 edges are required")])
    edges = []
    for k, edge_raw in enumerate(edges_raw):
        if isinstance(edge_raw, (str, int, dict)):
            edge_raw = [edge_raw]
        if not edge_raw:
            issues.append(Issue(f"{where}: {EDGE_NAMES[k]}", "is empty"))
            continue
        items = []
        for m, item in enumerate(edge_raw):
            normalized = _edge_item(item, f"{where}: {EDGE_NAMES[k]}[{m}]", store, issues)
            items.append(normalized)
        edges.append(tuple(items))
    try:
        generator = normalize_generator(raw.get("generator"), where)
    except ValidationError as error:
        issues.extend(error.issues)
        generator = None
    if issues:
        raise ValidationError(issues)
    params = {"edges": tuple(edges), "generator": generator, "surface": surface}
    check_domain(params, store, ConnectorEvaluator(store), where)
    return params


# ------------------------------------------------------------------------------ geometry
def _tolerance(points: np.ndarray) -> float:
    return CHAIN_TOLERANCE * max(1.0, float(np.max(np.abs(points))))


def edge_points(
    edge, connectors: ConnectorEvaluator, store, where: str, issues: list[Issue]
) -> np.ndarray | None:
    """Concatenated nodes of one edge (junction nodes once)."""
    parts = []
    for m, item in enumerate(edge):
        nodes = connectors.nodes(item["connector"]).points
        if item["reversed"]:
            nodes = nodes[::-1]
        if parts:
            gap = float(np.linalg.norm(parts[-1][-1] - nodes[0]))
            if gap > _tolerance(nodes):
                name = store.get(item["connector"]).name
                issues.append(
                    Issue(
                        where,
                        f"connector {m} ({name}) does not start where the previous one ends "
                        f"(gap {gap:.3g})",
                    )
                )
                return None
            nodes = nodes[1:]
        parts.append(nodes)
    return np.concatenate(parts)


def _edge_label(edge, store) -> str:
    names = []
    for item in edge:
        name = store.get(item["connector"]).name
        if item["reversed"]:
            name += " (reversed)"
        names.append(name)
    return ", ".join(names)


def loop_edges(params: dict, store, connectors: ConnectorEvaluator, where: str) -> list[np.ndarray]:
    """Edge node arrays of a domain in loop order; raises ValidationError."""
    issues: list[Issue] = []
    arrays = []
    for k, edge in enumerate(params["edges"]):
        arrays.append(edge_points(edge, connectors, store, f"{where}: {EDGE_NAMES[k]}", issues))
    if issues:
        raise ValidationError(issues)
    return arrays


def check_domain(params: dict, store, connectors: ConnectorEvaluator, where: str) -> None:
    arrays = loop_edges(params, store, connectors, where)
    issues: list[Issue] = []
    tolerance = _tolerance(np.concatenate(arrays))
    for k in range(4):
        gap = float(np.linalg.norm(arrays[k][-1] - arrays[(k + 1) % 4][0]))
        if gap > tolerance:
            issues.append(
                Issue(
                    f"{where}: {EDGE_NAMES[k]}",
                    f"does not end where {EDGE_NAMES[(k + 1) % 4]} starts (gap {gap:.3g})",
                )
            )
    for k in (0, 1):
        n_a = len(arrays[k])
        n_b = len(arrays[k + 2])
        if n_a != n_b:
            issues.append(
                Issue(
                    f"{where}",
                    f"{EDGE_NAMES[k]} [{_edge_label(params['edges'][k], store)}] has {n_a} "
                    f"points but the opposite {EDGE_NAMES[k + 2]} "
                    f"[{_edge_label(params['edges'][k + 2], store)}] has {n_b}",
                )
            )
    if issues:
        raise ValidationError(issues)
    loop = np.concatenate([a[:-1] for a in arrays])
    if np.ptp(loop[:, 2]) > tolerance or params.get("surface") is not None:
        # Face of a 3D block: linear TFI in space or on a surface; no orientation requirement.
        # Surface domains may be smoothed (Laplace / Thomas-Middlecoff, nodes kept on it).
        generator = params["generator"]
        elliptic = generator["elliptic"]
        if (
            params.get("surface") is not None
            and generator["method"] == "tfi"
            and elliptic is not None
            and elliptic["control"] in ("laplace", "thomas_middlecoff")
            and not elliptic["walls"]
            and not elliptic["attractions"]
        ):
            return
        if generator["method"] != "tfi" or elliptic is not None:
            raise ValidationError(
                [
                    Issue(
                        where,
                        "Hermite TFI and elliptic smoothing need a planar domain in a plane "
                        "z = const; domains in 3D use method 'tfi' (on a surface, elliptic "
                        "smoothing with Laplace or Thomas-Middlecoff control and no walls)",
                    )
                ]
            )
        return
    area = 0.5 * float(
        np.sum(loop[:, 0] * np.roll(loop[:, 1], -1) - np.roll(loop[:, 0], -1) * loop[:, 1])
    )
    if area <= 0.0:
        raise ValidationError(
            [
                Issue(
                    where,
                    "the edge loop runs clockwise (or encloses no area); "
                    "list the edges counter-clockwise",
                )
            ]
        )


def boundary_arrays(
    arrays: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(bottom, right, top, left), each ordered by increasing grid index."""
    return arrays[0], arrays[1], arrays[2][::-1].copy(), arrays[3][::-1].copy()


# ------------------------------------------------------------------- automatic assembly
def _turning_angle(incoming: np.ndarray, outgoing: np.ndarray) -> float:
    a = incoming[-1] - incoming[-2]
    b = outgoing[1] - outgoing[0]
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(a, b) / (na * nb))))))


def assemble_edges(store, connectors: ConnectorEvaluator, names: list[str]) -> list[list[dict]]:
    """Edges for a domain from a set of connectors.

    Supported: a single closed loop of >= 4 connectors (corners chosen among the junctions so
    that opposite edges have equal point counts, preferring the sharpest turns), and an
    O-type set of two closed loops joined by one connector (used on both sides of the cut).
    The loop is oriented counter-clockwise; edge 0 contains the lowest-id connector (simple
    loop) or is the inner loop (O-type).
    """
    if not names:
        raise ValidationError([Issue("domain", "no connectors given")])
    ids = []
    for name in names:
        entity = store.get(name)
        if entity.kind != "connector":
            raise ValidationError([Issue(name, f"is a {entity.kind}, not a connector")])
        ids.append(entity.id)
    ids = sorted(set(ids))
    nodes = {i: connectors.nodes(i).points for i in ids}
    # Cluster end points into graph nodes.
    keys: list[np.ndarray] = []

    def key_of(point) -> int:
        tol = _tolerance(np.asarray([point]))
        for index, known in enumerate(keys):
            if np.linalg.norm(known - point) <= tol:
                return index
        keys.append(np.asarray(point))
        return len(keys) - 1

    ends = {i: (key_of(nodes[i][0]), key_of(nodes[i][-1])) for i in ids}
    degree: dict[int, int] = {}
    for a, b in ends.values():
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1

    if all(d == 2 for d in degree.values()):
        return _assemble_loop(store, ids, nodes, ends)
    cut = [
        i
        for i in ids
        if degree[ends[i][0]] == 3 and degree[ends[i][1]] == 3 and ends[i][0] != ends[i][1]
    ]
    if (
        len(cut) == 1
        and all(d in (2, 3) for d in degree.values())
        and sum(d == 3 for d in degree.values()) == 2
    ):
        return _assemble_o_type(store, ids, nodes, ends, cut[0])
    raise ValidationError(
        [
            Issue(
                "domain",
                "connectors must form one closed loop, or two closed loops joined by one connector",
            )
        ]
    )


def _walk(ids, ends, start_id, start_reversed=False):
    """Walk a closed chain of connectors starting with start_id; returns [(id, reversed)]."""
    order = [(start_id, start_reversed)]
    used = {start_id}
    if start_reversed:
        current = ends[start_id][0]
        first = ends[start_id][1]
    else:
        current = ends[start_id][1]
        first = ends[start_id][0]
    while current != first:
        for i in ids:
            if i in used:
                continue
            if ends[i][0] == current:
                order.append((i, False))
                current = ends[i][1]
                break
            if ends[i][1] == current:
                order.append((i, True))
                current = ends[i][0]
                break
        else:
            raise ValidationError([Issue("domain", "connectors do not close into a loop")])
        used.add(order[-1][0])
    return order


def _oriented(nodes, item):
    points = nodes[item[0]]
    if item[1]:
        return points[::-1]
    return points


def _signed_area(points: np.ndarray) -> float:
    return 0.5 * float(
        np.sum(points[:, 0] * np.roll(points[:, 1], -1) - np.roll(points[:, 0], -1) * points[:, 1])
    )


def _assemble_loop(store, ids, nodes, ends):
    order = _walk(ids, ends, ids[0])
    if len(order) != len(ids):
        raise ValidationError([Issue("domain", "connectors form more than one loop")])
    loop = np.concatenate([_oriented(nodes, item)[:-1] for item in order])
    if _signed_area(loop) < 0.0:
        order = [(i, not r) for i, r in reversed(order)]
    n = len(order)
    if n < 4:
        raise ValidationError(
            [
                Issue(
                    "domain",
                    f"a loop of {n} connector(s) cannot form 4 edges; "
                    "split connectors at the corners",
                )
            ]
        )
    sizes = [len(_oriented(nodes, item)) for item in order]
    turning = [
        _turning_angle(_oriented(nodes, order[k - 1]), _oriented(nodes, order[k])) for k in range(n)
    ]
    best = None
    for corners in itertools.combinations(range(n), 4):
        counts = []
        for c in range(4):
            start = corners[c]
            stop = corners[(c + 1) % 4]
            members = [(start + m) % n for m in range((stop - start) % n or n)]
            counts.append(sum(sizes[k] for k in members) - (len(members) - 1))
        if counts[0] != counts[2] or counts[1] != counts[3]:
            continue
        score = sum(turning[c] for c in corners)
        if best is None or score > best[0] + 1e-9:
            best = (score, corners)
    if best is None:
        raise ValidationError(
            [
                Issue(
                    "domain",
                    "no choice of 4 corners gives matching opposite point counts "
                    f"(connector sizes {sizes})",
                )
            ]
        )
    corners = best[1]
    edges = []
    for c in range(4):
        start = corners[c]
        stop = corners[(c + 1) % 4]
        members = [(start + m) % n for m in range((stop - start) % n or n)]
        edges.append([{"connector": order[k][0], "reversed": order[k][1]} for k in members])
    # Edge 0 contains the lowest-id connector.
    lowest = min(ids)
    for shift in range(4):
        if any(item["connector"] == lowest for item in edges[shift]):
            return edges[shift:] + edges[:shift]
    return edges


def _assemble_o_type(store, ids, nodes, ends, cut):
    rest = [i for i in ids if i != cut]
    a, b = ends[cut]
    loops = []
    for anchor in (a, b):
        start = [i for i in rest if anchor in ends[i]]
        if not start:
            raise ValidationError(
                [Issue("domain", "the cut connector does not join two closed loops")]
            )
        first = start[0]
        order = _walk(rest, ends, first, ends[first][1] == anchor and ends[first][0] != anchor)
        loops.append(order)
    if sum(len(loop) for loop in loops) != len(rest):
        raise ValidationError([Issue("domain", "connectors outside the two loops")])
    areas = [
        _signed_area(np.concatenate([_oriented(nodes, item)[:-1] for item in loop]))
        for loop in loops
    ]
    # Edge 0: the inner loop (smaller enclosed area), traversed clockwise; edge 2: the outer loop.
    inner = 0
    if abs(areas[1]) < abs(areas[0]):
        inner = 1
    outer = 1 - inner
    inner_loop = loops[inner]
    outer_loop = loops[outer]
    if areas[inner] > 0:
        inner_loop = _reverse_closed(inner_loop)
    # Edge 2 is listed in loop order (i decreasing), i.e. counter-clockwise along the outer loop
    # (with i increasing it runs clockwise, like edge 0).
    if areas[outer] < 0:
        outer_loop = _reverse_closed(outer_loop)
    cut_from_inner = (a, b)[inner] == ends[cut][0]
    edge1 = [{"connector": cut, "reversed": not cut_from_inner}]
    edge3 = [{"connector": cut, "reversed": cut_from_inner}]
    edges = [
        [{"connector": i, "reversed": r} for i, r in inner_loop],
        edge1,
        [{"connector": i, "reversed": r} for i, r in outer_loop],
        edge3,
    ]
    return edges


def _reverse_closed(order):
    """Reverse a closed walk keeping its start node."""
    return [(i, not r) for i, r in reversed(order)]


# ---------------------------------------------------------------------------- evaluation
@dataclass
class SmoothingResult:
    grid: np.ndarray
    residuals: np.ndarray
    iterations: int
    converged: bool
    cancelled: bool
    warnings: list[str] = field(default_factory=list)


def elliptic_arguments(elliptic: dict) -> dict:
    lines = []
    points = []
    for item in elliptic["attractions"]:
        if item["kind"] == "xi_line":
            lines.append((True, item["xi"], item["amplitude"], item["decay"]))
        elif item["kind"] == "eta_line":
            lines.append((False, item["eta"], item["amplitude"], item["decay"]))
        else:
            points.append((item["xi"], item["eta"], item["amplitude"], item["decay"]))
    walls = {0: None, 2: None}
    for wall in elliptic["walls"]:
        walls[wall["edge"]] = (
            wall["spacing"],
            wall["decay"],
            wall["relaxation"],
            wall["max_source"],
        )
    return {
        "control": elliptic["control"],
        "lines": lines,
        "points": points,
        "bottom_wall": walls[0],
        "top_wall": walls[2],
        "solver": elliptic["solver"],
        "omega": elliptic["omega"],
        "max_iterations": elliptic["max_iterations"],
        "tolerance": elliptic["tolerance"],
    }


def surface_signature(store, surface_id: int) -> tuple:
    """Parameters of a surface and of the curves it references (for cache keys)."""
    from opengrid.model.entities import SURFACE_CURVES

    entity = store.get(surface_id)
    parts = [entity.kind, repr(entity.params)]
    for key in SURFACE_CURVES.get(entity.kind, ()):
        curve = store.get(entity.params[key])
        parts.append(repr(curve.params))
    return tuple(parts)


class DomainEvaluator:
    """Grids of structured domains. TFI grids are recomputed on demand; smoothed grids are
    cached by a key over the generator settings and the boundary nodes, so a smoothing run
    stays valid until something that affects it changes."""

    CACHE_SIZE = 32

    def __init__(self, store, connectors: ConnectorEvaluator):
        self.store = store
        self.connectors = connectors
        self._boundaries: dict[int, tuple] = {}
        self._smoothed: OrderedDict[str, SmoothingResult] = OrderedDict()

    def invalidate(self) -> None:
        self._boundaries.clear()

    def boundary(self, entity_id: int):
        if entity_id not in self._boundaries:
            entity = self.store.get(entity_id)
            arrays = loop_edges(entity.params, self.store, self.connectors, entity.name)
            self._boundaries[entity_id] = boundary_arrays(arrays)
        return self._boundaries[entity_id]

    def dimensions(self, entity_id: int) -> tuple[int, int]:
        bottom, right, _top, _left = self.boundary(entity_id)
        return len(bottom), len(right)

    def key(self, entity_id: int) -> str:
        params = self.store.get(entity_id).params
        digest = hashlib.sha256(repr(params["generator"]).encode())
        if params.get("surface") is not None:
            digest.update(repr(surface_signature(self.store, params["surface"])).encode())
        for array in self.boundary(entity_id):
            digest.update(np.ascontiguousarray(array).tobytes())
        return digest.hexdigest()

    def initial_grid(self, entity_id: int) -> np.ndarray:
        params = self.store.get(entity_id).params
        generator = params["generator"]
        bottom, right, top, left = self.boundary(entity_id)
        if params.get("surface") is not None:
            surface = self.store.geometry(params["surface"])
            grid, _distance = _core.tfi_on_surface(
                surface,
                np.ascontiguousarray(bottom),
                np.ascontiguousarray(right),
                np.ascontiguousarray(top),
                np.ascontiguousarray(left),
            )
            return grid
        if generator["method"] == "hermite_tfi":
            return _core.tfi(
                bottom,
                right,
                top,
                left,
                method="hermite",
                hermite_bottom=0 in generator["hermite_walls"],
                hermite_top=2 in generator["hermite_walls"],
                strength=generator["hermite_strength"],
            )
        return _core.tfi(bottom, right, top, left)

    def smoothing_result(self, entity_id: int) -> SmoothingResult | None:
        return self._smoothed.get(self.key(entity_id))

    def store_result(self, key: str, result: SmoothingResult) -> None:
        self._smoothed[key] = result
        self._smoothed.move_to_end(key)
        while len(self._smoothed) > self.CACHE_SIZE:
            self._smoothed.popitem(last=False)

    def smooth(self, entity_id: int, progress=None) -> SmoothingResult:
        """Run (or reuse) the elliptic smoother; synchronous."""
        generator = self.store.get(entity_id).params["generator"]
        if generator["elliptic"] is None:
            raise ValueError(
                f"{self.store.get(entity_id).name!r} has no elliptic smoothing configured"
            )
        key = self.key(entity_id)
        if key in self._smoothed:
            return self._smoothed[key]
        surface = self.store.get(entity_id).params.get("surface")
        if surface is not None:
            e = generator["elliptic"]
            grid, residuals, iterations, converged, cancelled, warnings = _core.smooth_on_surface(
                self.store.geometry(surface),
                np.ascontiguousarray(self.initial_grid(entity_id)),
                e["control"],
                e["omega"],
                e["max_iterations"],
                e["tolerance"],
                progress,
            )
            result = SmoothingResult(
                grid, np.asarray(residuals), iterations, converged, cancelled, list(warnings)
            )
        else:
            result = run_smoother(self.initial_grid(entity_id), generator["elliptic"], progress)
        if not result.cancelled:
            self.store_result(key, result)
        return result

    def grid(self, entity_id: int, compute: bool = True) -> np.ndarray:
        """Final grid: TFI, or the smoothed grid (computed now if ``compute`` and missing;
        otherwise the TFI grid is returned while smoothing is pending)."""
        generator = self.store.get(entity_id).params["generator"]
        if generator["elliptic"] is None:
            return self.initial_grid(entity_id)
        cached = self.smoothing_result(entity_id)
        if cached is not None:
            return cached.grid
        if compute:
            return self.smooth(entity_id).grid
        return self.initial_grid(entity_id)

    def status(self, entity_id: int) -> str:
        generator = self.store.get(entity_id).params["generator"]
        if generator["elliptic"] is None:
            return generator["method"]
        cached = self.smoothing_result(entity_id)
        if cached is None:
            return "smoothing pending (showing TFI)"
        if cached.converged:
            return f"smoothed, converged in {cached.iterations} iterations"
        return f"smoothed, NOT converged after {cached.iterations} iterations"


def run_smoother(initial: np.ndarray, elliptic: dict, progress=None) -> SmoothingResult:
    grid, residuals, iterations, converged, cancelled, warnings = _core.smooth_elliptic(
        initial, progress=progress, **elliptic_arguments(elliptic)
    )
    return SmoothingResult(grid, residuals, iterations, converged, cancelled, list(warnings))
