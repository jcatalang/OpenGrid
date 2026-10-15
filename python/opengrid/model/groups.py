"""Smoothing groups: planar structured domains smoothed together as one grid.

A ``smoothing_group`` lists structured domains. Their grids are no longer smoothed one by one:
the TFI grids of all members are merged into one node set (nodes shared through connectors
are one node) and relaxed together by ``opengrid._core.smooth_multiblock``, so interfaces and
block corners move and grid lines pass smoothly from one domain to the next. Each member's
grid is its part of the result; its own ``elliptic`` setting is ignored while it is in a group.

Connectors, by their uses:

* used by two member domains (or twice by one, a cut): an *interface*, its nodes float;
* used by one member only: a *boundary*. Walls (``walls``) keep their nodes and get
  first-cell height and orthogonality control; other boundary connectors either slide along
  their curves (``slide``: "auto" slides every non-wall boundary connector, or a tuple of
  connectors, or () for none) or stay fixed. The end nodes of boundary connectors stay fixed;
* used by any domain outside the group: fixed (the neighbour needs those nodes).

Members must be planar (z = const) structured domains without a surface. Parameters never
use the keys id/name/type (project record keys).
"""

import hashlib
import math
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.model.domains import MAX_ACCELERATION
from opengrid.model.topology import _position
from opengrid.validation import Issue, ValidationError

KIND = "smoothing_group"
DEFAULT_GROUP = {
    "domains": (),
    "walls": (),
    "slide": "auto",
    "omega": 1.0,
    "max_iterations": 2000,
    "tolerance": 1e-6,
    # Deeper than for single domains: Laplace smoothing of a whole multi-block grid has slow
    # global modes (a 68k-node O-grid: 1367 iterations at depth 5, 456 at depth 10).
    "acceleration": 10,
}
DEFAULT_GROUP_WALL = {
    "connector": None,
    "spacing": None,
    "decay": 0.1,
    "relaxation": 0.1,
    "max_source": 100.0,
}
NODE_CLASSES = ("fixed", "regular", "irregular", "sliding")


# ----------------------------------------------------------------------------- parameters
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


def _ref(store, raw, kind, label, issues):
    if isinstance(raw, bool) or not isinstance(raw, (str, int)) or raw not in store:
        issues.append(Issue(label, f"no entity {raw!r}"))
        return None
    entity = store.get(raw)
    if entity.kind != kind:
        issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a {kind}"))
        return None
    return entity.id


def member_connectors(store, domain_ids) -> set[int]:
    return {
        item["connector"]
        for d in domain_ids
        for edge in store.get(d).params["edges"]
        for item in edge
    }


def normalize_group(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - set(DEFAULT_GROUP)):
        issues.append(Issue(f"{where}: {key}", "unexpected key"))
    merged = dict(DEFAULT_GROUP)
    merged.update(raw)
    domains_raw = merged["domains"]
    if isinstance(domains_raw, (str, int)):
        domains_raw = [domains_raw]
    domains = []
    for index, item in enumerate(domains_raw):
        ref = _ref(store, item, "domain", f"{where}: domains[{index}]", issues)
        if ref is None:
            continue
        if ref in domains:
            issues.append(Issue(f"{where}: domains[{index}]", "listed twice"))
            continue
        if store.get(ref).params.get("surface") is not None:
            issues.append(
                Issue(f"{where}: domains[{index}]", f"{store.get(ref).name!r} lies on a surface")
            )
            continue
        domains.append(ref)
    if not domains_raw:
        issues.append(Issue(f"{where}: domains", "at least one domain is required"))
    connectors = member_connectors(store, domains)
    walls = []
    seen = set()
    for index, item in enumerate(merged["walls"]):
        label = f"{where}: walls[{index}]"
        if isinstance(item, (str, int)) and not isinstance(item, bool):
            item = {"connector": item}
        if not isinstance(item, dict):
            issues.append(Issue(label, f"expected a connector or a dict, got {item!r}"))
            continue
        for key in sorted(set(item) - set(DEFAULT_GROUP_WALL)):
            issues.append(Issue(f"{label}.{key}", "unexpected key"))
        wall = dict(DEFAULT_GROUP_WALL)
        wall.update(item)
        ref = _ref(store, wall["connector"], "connector", f"{label}.connector", issues)
        if ref is None:
            continue
        if ref not in connectors:
            issues.append(Issue(label, f"{store.get(ref).name!r} is not used by a member domain"))
            continue
        if ref in seen:
            issues.append(Issue(label, f"{store.get(ref).name!r} listed twice"))
            continue
        seen.add(ref)
        walls.append(
            {
                "connector": ref,
                "spacing": _number(
                    wall["spacing"], f"{label}.spacing", issues, positive=True, optional=True
                ),
                "decay": _number(wall["decay"], f"{label}.decay", issues, positive=True),
                "relaxation": _number(
                    wall["relaxation"], f"{label}.relaxation", issues, positive=True
                ),
                "max_source": _number(
                    wall["max_source"], f"{label}.max_source", issues, positive=True
                ),
            }
        )
    slide = merged["slide"]
    if slide != "auto":
        if isinstance(slide, (str, int)):
            slide = [slide]
        refs = []
        for index, item in enumerate(slide):
            label = f"{where}: slide[{index}]"
            ref = _ref(store, item, "connector", label, issues)
            if ref is None:
                continue
            if ref not in connectors:
                issues.append(
                    Issue(label, f"{store.get(ref).name!r} is not used by a member domain")
                )
                continue
            if ref in seen:
                issues.append(Issue(label, f"{store.get(ref).name!r} is a wall (walls stay fixed)"))
                continue
            if ref not in refs:
                refs.append(ref)
        slide = tuple(sorted(refs))
    omega = _number(merged["omega"], f"{where}: omega", issues, positive=True)
    if omega is not None and omega >= 2.0:
        issues.append(Issue(f"{where}: omega", "must be below 2"))
    iterations = merged["max_iterations"]
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 1:
        issues.append(
            Issue(f"{where}: max_iterations", f"must be a positive integer, got {iterations!r}")
        )
    tolerance = _number(merged["tolerance"], f"{where}: tolerance", issues, positive=True)
    acceleration = merged["acceleration"]
    if (
        isinstance(acceleration, bool)
        or not isinstance(acceleration, int)
        or not 0 <= acceleration <= MAX_ACCELERATION
    ):
        issues.append(
            Issue(
                f"{where}: acceleration",
                f"must be an integer from 0 (off) to {MAX_ACCELERATION}, got {acceleration!r}",
            )
        )
    if issues:
        raise ValidationError(issues)
    walls.sort(key=lambda w: w["connector"])
    return {
        "domains": tuple(domains),
        "walls": tuple(walls),
        "slide": slide,
        "omega": omega,
        "max_iterations": iterations,
        "tolerance": tolerance,
        "acceleration": acceleration,
    }


def group_references(params: dict) -> set[int]:
    refs = set(params["domains"])
    refs.update(w["connector"] for w in params["walls"])
    if params["slide"] != "auto":
        refs.update(params["slide"])
    return refs


def group_script_params(params: dict, store) -> dict:
    result = dict(params)
    result["domains"] = tuple(store.get(d).name for d in params["domains"])
    walls = []
    for wall in params["walls"]:
        named = dict(wall)
        named["connector"] = store.get(wall["connector"]).name
        walls.append(named)
    result["walls"] = tuple(walls)
    if params["slide"] != "auto":
        result["slide"] = tuple(store.get(c).name for c in params["slide"])
    return result


# ---------------------------------------------------------------------------- assembly
@dataclass
class GroupLayout:
    """Global numbering of a group: node ids per member grid and the solver input."""

    ids: dict[int, np.ndarray]  # domain id -> (nj, ni) global node ids
    points: np.ndarray  # (n, 3) initial positions
    quads: np.ndarray  # (m, 4) counter-clockwise
    fixed: np.ndarray  # (n,) uint8
    chains: list
    slides: list[tuple[int, int]]
    walls: list[tuple[int, float, float, float, float]]
    roles: dict[int, str]  # connector id -> interface / wall / sliding / fixed / outside
    ends: dict[int, tuple[int, int]]  # connector id -> global ids of its first and last node


@dataclass
class GroupResult:
    grids: dict[int, np.ndarray]
    residuals: np.ndarray
    iterations: int
    converged: bool
    cancelled: bool
    node_counts: dict[str, int]
    lines: int
    colours: int
    roles: dict[int, str]
    singular: frozenset  # (connector, 0 | 1): interface ends at irregular points
    max_kink: float  # largest interface kink (degrees) away from irregular points
    warnings: list[str] = field(default_factory=list)


class _UnionFind:
    """Union-find over the few keys that are ever merged (boundary nodes); smaller roots win."""

    def __init__(self):
        self.parent: dict[int, int] = {}

    def find(self, a: int) -> int:
        root = a
        while self.parent.get(root, root) != root:
            root = self.parent[root]
        while self.parent.get(a, a) != root:
            self.parent[a], a = root, self.parent[a]
        return root

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra < rb:
            self.parent[rb] = ra
        elif rb < ra:
            self.parent[ra] = rb


def _signed_area(grid: np.ndarray) -> float:
    loop = np.concatenate([grid[0, :, :2], grid[:, -1, :2], grid[-1, ::-1, :2], grid[::-1, 0, :2]])
    x = loop[:, 0]
    y = loop[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def connector_roles(store, params: dict) -> dict[int, str]:
    members = set(params["domains"])
    uses: dict[int, int] = {}
    for d in params["domains"]:
        for edge in store.get(d).params["edges"]:
            for item in edge:
                uses[item["connector"]] = uses.get(item["connector"], 0) + 1
    walls = {w["connector"] for w in params["walls"]}
    slide = params["slide"]
    roles = {}
    for c in sorted(uses):
        outside = [
            e
            for e in store.dependents(c)
            if e.kind in ("domain", "unstructured_domain") and e.id not in members
        ]
        if outside:
            roles[c] = "outside"
        elif uses[c] >= 2:
            roles[c] = "interface"
        elif c in walls:
            roles[c] = "wall"
        elif slide == "auto" or c in slide:
            roles[c] = "sliding"
        else:
            roles[c] = "fixed"
    return roles


def build_layout(store, connectors, grids: dict[int, np.ndarray], params: dict) -> GroupLayout:
    """Merge the member grids (initial positions in ``grids``) into one node set."""
    domains = params["domains"]
    roles = connector_roles(store, params)
    # Keys: member grid nodes first (domain order, row-major), then connector nodes.
    offsets = {}
    total = 0
    for d in domains:
        offsets[d] = total
        total += grids[d].shape[0] * grids[d].shape[1]
    connector_offset = {}
    for c in sorted(roles):
        connector_offset[c] = total
        total += store.get(c).params["dimension"]
    uf = _UnionFind()
    for d in domains:
        nj, ni = grids[d].shape[:2]
        for k, edge in enumerate(store.get(d).params["edges"]):
            p = 0
            for item in edge:
                c = item["connector"]
                n = store.get(c).params["dimension"]
                for m in range(n):
                    i, j = _position(k, p + m, ni, nj)
                    node = m
                    if item["reversed"]:
                        node = n - 1 - m
                    uf.union(offsets[d] + j * ni + i, connector_offset[c] + node)
                p += n - 1
    # Every connector node is merged into a member grid node (keys of grid nodes are smaller),
    # so the roots of the grid keys are all the distinct nodes; ids follow first appearance.
    member_keys = offsets[domains[-1]] + grids[domains[-1]].shape[0] * grids[domains[-1]].shape[1]
    roots = np.arange(member_keys, dtype=np.int64)
    for key in uf.parent:
        if key < member_keys:
            roots[key] = uf.find(key)
    unique, first, inverse = np.unique(roots, return_index=True, return_inverse=True)
    order = np.argsort(first, kind="stable")
    rank = np.empty(len(unique), dtype=np.int64)
    rank[order] = np.arange(len(unique))
    flat_ids = rank[inverse]
    all_points = np.concatenate([grids[d].reshape(-1, 3) for d in domains])
    points = all_points[first[order]]
    ids = {}
    for d in domains:
        nj, ni = grids[d].shape[:2]
        ids[d] = flat_ids[offsets[d] : offsets[d] + nj * ni].reshape(nj, ni)

    def node_of(c: int, m: int) -> int:
        root = uf.find(connector_offset[c] + m)
        return int(rank[np.searchsorted(unique, root)])

    quads = []
    for d in domains:
        g = ids[d]
        q = np.stack([g[:-1, :-1], g[:-1, 1:], g[1:, 1:], g[1:, :-1]], axis=-1).reshape(-1, 4)
        if _signed_area(grids[d]) < 0.0:
            q = q[:, [0, 3, 2, 1]]
        quads.append(q)
    n = len(points)
    fixed = np.zeros(n, dtype=np.uint8)
    chains = []
    slides = []
    walls = []
    wall_nodes = set()
    wall_settings = {w["connector"]: w for w in params["walls"]}
    for c in sorted(roles):
        role = roles[c]
        dim = store.get(c).params["dimension"]
        if role == "interface":
            continue
        # Boundary and outside connectors: end nodes always fixed.
        fixed[node_of(c, 0)] = 1
        fixed[node_of(c, dim - 1)] = 1
        if role == "sliding":
            chain_index = len(chains)
            chains.append(connectors.chain(c))
            for m in range(1, dim - 1):
                slides.append((node_of(c, m), chain_index))
            continue
        for m in range(dim):
            fixed[node_of(c, m)] = 1
        if role == "wall":
            w = wall_settings[c]
            spacing = w["spacing"]
            if spacing is None:
                spacing = 0.0
            for m in range(dim):
                node = node_of(c, m)
                if node in wall_nodes:
                    continue
                wall_nodes.add(node)
                walls.append((node, spacing, w["decay"], w["relaxation"], w["max_source"]))
    # A node that both slides and is fixed (end of another connector) stays fixed.
    slides = [(node, chain) for node, chain in slides if not fixed[node]]
    return GroupLayout(
        ids=ids,
        points=np.ascontiguousarray(points, dtype=float),
        quads=np.ascontiguousarray(np.concatenate(quads).astype(np.int64)),
        fixed=fixed,
        chains=chains,
        slides=slides,
        walls=walls,
        roles=roles,
        ends={c: (node_of(c, 0), node_of(c, store.get(c).params["dimension"] - 1)) for c in roles},
    )


def _inward(edge: int, i: int, j: int, step: int = 1) -> tuple[int, int]:
    """(i, j) of the node ``step`` lines inside the domain from boundary node (i, j)."""
    if edge == 0:
        return (i, j + step)
    if edge == 1:
        return (i - step, j)
    if edge == 2:
        return (i, j - step)
    return (i + step, j)


def _signed_turn(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Signed turn (degrees) of the polyline a -> b -> c at b."""
    u = b - a
    v = c - b
    return math.degrees(math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]))


def interface_kinks(
    store, params: dict, grids: dict[int, np.ndarray], singular=frozenset(), margin: int = 3
) -> dict[int, float]:
    """Largest kink (degrees) of the grid lines crossing each interface connector of a group:
    the turn of the line at the interface node minus the mean of its turns at the neighbouring
    nodes on either side (zero for a line that crosses smoothly, whatever its curvature). Over
    the connector's interior nodes, except the ``margin`` nodes next to a connector end listed
    in ``singular`` as (connector id, 0 or 1) - irregular points, where lines must turn.
    ``grids``: member domain id -> grid (any state)."""
    uses: dict[int, list] = {}
    for d in params["domains"]:
        nj, ni = grids[d].shape[:2]
        for k, edge in enumerate(store.get(d).params["edges"]):
            p = 0
            for item in edge:
                c = item["connector"]
                n = store.get(c).params["dimension"]
                positions = []
                for m in range(n):
                    node = m
                    if item["reversed"]:
                        node = n - 1 - m
                    positions.append((node, *_position(k, p + m, ni, nj)))
                positions.sort()
                uses.setdefault(c, []).append((d, k, positions))
                p += n - 1
    kinks = {}
    for c, found in sorted(uses.items()):
        if len(found) != 2:
            continue
        (da, ka, pa), (db, kb, pb) = found
        worst = 0.0
        n = len(pa)
        for (m, ia, ja), (_n, ib, jb) in zip(pa[1:-1], pb[1:-1], strict=True):
            if ((c, 0) in singular and m <= margin) or ((c, 1) in singular and m >= n - 1 - margin):
                continue
            ga = grids[da]
            gb = grids[db]
            if min(ga.shape[0], ga.shape[1], gb.shape[0], gb.shape[1]) < 3:
                continue
            g = ga[ja, ia, :2]
            a1 = ga[_inward(ka, ia, ja, 1)[::-1]][:2]
            a2 = ga[_inward(ka, ia, ja, 2)[::-1]][:2]
            b1 = gb[_inward(kb, ib, jb, 1)[::-1]][:2]
            b2 = gb[_inward(kb, ib, jb, 2)[::-1]][:2]
            here = _signed_turn(a1, g, b1)
            around = 0.5 * (_signed_turn(a2, a1, g) + _signed_turn(g, b1, b2))
            worst = max(worst, abs(here - around))
        kinks[c] = worst
    return kinks


# ---------------------------------------------------------------------------- evaluation
def run_group(layout: GroupLayout, params: dict, progress=None):
    """The solver call alone (touches no model state: safe on a worker thread)."""
    return _core.smooth_multiblock(
        layout.points,
        layout.quads,
        layout.fixed,
        layout.chains,
        layout.slides,
        layout.walls,
        omega=params["omega"],
        max_iterations=params["max_iterations"],
        tolerance=params["tolerance"],
        acceleration=params["acceleration"],
        progress=progress,
    )


class GroupEvaluator:
    """Results of smoothing groups, cached by a key over the group settings, the member
    domains' initial grids and the connector roles."""

    CACHE_SIZE = 8

    def __init__(self, store, connectors, domains):
        self.store = store
        self.connectors = connectors
        self.domains = domains
        self._results: OrderedDict[str, GroupResult] = OrderedDict()

    def group_of(self, domain_id: int) -> int | None:
        found = []
        for entity in self.store.entities(KIND):
            if domain_id in entity.params["domains"]:
                found.append(entity)
        if len(found) > 1:
            names = ", ".join(e.name for e in found)
            raise ValidationError(
                [Issue(self.store.get(domain_id).name, f"is in several smoothing groups ({names})")]
            )
        if found:
            return found[0].id
        return None

    def _check(self, group_id: int) -> None:
        entity = self.store.get(group_id)
        issues = []
        for d in entity.params["domains"]:
            others = [
                e.name
                for e in self.store.entities(KIND)
                if e.id != group_id and d in e.params["domains"]
            ]
            if others:
                issues.append(
                    Issue(
                        entity.name,
                        f"domain {self.store.get(d).name!r} is also in {', '.join(others)}",
                    )
                )
            grid = self.domains.initial_grid(d)
            z = grid[..., 2]
            scale = max(1.0, float(np.max(np.abs(grid))))
            if float(np.ptp(z)) > 1e-9 * scale:
                issues.append(
                    Issue(
                        entity.name,
                        f"domain {self.store.get(d).name!r} is not in a z = const plane",
                    )
                )
        zs = {float(self.domains.initial_grid(d)[0, 0, 2]) for d in entity.params["domains"]}
        if len(zs) > 1:
            issues.append(Issue(entity.name, "member domains lie in different planes"))
        if issues:
            raise ValidationError(issues)

    def key(self, group_id: int) -> str:
        entity = self.store.get(group_id)
        digest = hashlib.sha256(repr(entity.params).encode())
        for d in entity.params["domains"]:
            digest.update(self.domains.key(d).encode())
        digest.update(repr(sorted(connector_roles(self.store, entity.params).items())).encode())
        return digest.hexdigest()

    def result(self, group_id: int) -> GroupResult | None:
        return self._results.get(self.key(group_id))

    def prepare(self, group_id: int) -> tuple[str, GroupLayout, dict]:
        """Cache key, solver input and parameters (on the thread that owns the model)."""
        self._check(group_id)
        params = self.store.get(group_id).params
        grids = {d: self.domains.initial_grid(d) for d in params["domains"]}
        return self.key(group_id), build_layout(self.store, self.connectors, grids, params), params

    def finish(self, params: dict, layout: GroupLayout, raw) -> GroupResult:
        """GroupResult from the output of :func:`run_group`."""
        (
            points,
            classes,
            residuals,
            iterations,
            converged,
            cancelled,
            lines,
            colours,
            warnings,
        ) = raw
        counts = np.bincount(np.asarray(classes, dtype=np.int64), minlength=4)
        grids = {d: points[layout.ids[d]] for d in params["domains"]}
        singular = set()
        for c, (first, last) in layout.ends.items():
            if classes[first] == 2:
                singular.add((c, 0))
            if classes[last] == 2:
                singular.add((c, 1))
        singular = frozenset(singular)
        kinks = interface_kinks(self.store, params, grids, singular)
        return GroupResult(
            grids=grids,
            residuals=np.asarray(residuals),
            iterations=iterations,
            converged=converged,
            cancelled=cancelled,
            node_counts={name: int(counts[k]) for k, name in enumerate(NODE_CLASSES)},
            lines=lines,
            colours=colours,
            roles=layout.roles,
            singular=singular,
            max_kink=max(kinks.values(), default=0.0),
            warnings=list(warnings),
        )

    def store_result(self, key: str, result: GroupResult) -> None:
        if result.cancelled:
            return
        self._results[key] = result
        self._results.move_to_end(key)
        while len(self._results) > self.CACHE_SIZE:
            self._results.popitem(last=False)

    def smooth(self, group_id: int, progress=None) -> GroupResult:
        """Run (or reuse) the group's smoothing; synchronous."""
        key = self.key(group_id)
        if key in self._results:
            return self._results[key]
        key, layout, params = self.prepare(group_id)
        result = self.finish(params, layout, run_group(layout, params, progress))
        self.store_result(key, result)
        return result

    def moved_connector_points(self, connector_id: int) -> np.ndarray | None:
        """Nodes of an interface or sliding connector as placed by a cached group result (in
        the connector's own order), or None when no smoothed group moves it."""
        for entity in self.store.entities(KIND):
            if not any(
                item["connector"] == connector_id
                for d in entity.params["domains"]
                for edge in self.store.get(d).params["edges"]
                for item in edge
            ):
                continue
            result = self.result(entity.id)
            if result is None or result.roles.get(connector_id) not in ("interface", "sliding"):
                continue
            n = self.store.get(connector_id).params["dimension"]
            for d in entity.params["domains"]:
                grid = result.grids[d]
                nj, ni = grid.shape[:2]
                for k, edge in enumerate(self.store.get(d).params["edges"]):
                    p = 0
                    for item in edge:
                        m = self.store.get(item["connector"]).params["dimension"]
                        if item["connector"] == connector_id:
                            points = np.array(
                                [grid[_position(k, p + a, ni, nj)[::-1]] for a in range(m)]
                            )
                            if item["reversed"]:
                                points = points[::-1]
                            return np.ascontiguousarray(points[:n])
                        p += m - 1
        return None

    def status(self, group_id: int) -> str:
        entity = self.store.get(group_id)
        cached = self.result(group_id)
        count = len(entity.params["domains"])
        if cached is None:
            return f"{count} domains, smoothing pending (showing TFI)"
        nodes = cached.node_counts
        text = (
            f"{count} domains, {sum(nodes.values())} nodes ({nodes['regular']} regular, "
            f"{nodes['irregular']} irregular, {nodes['sliding']} sliding, {nodes['fixed']} fixed), "
            f"interface kink {cached.max_kink:.2g} deg"
        )
        if cached.converged:
            return f"{text}, converged in {cached.iterations} iterations"
        return f"{text}, NOT converged after {cached.iterations} iterations"
