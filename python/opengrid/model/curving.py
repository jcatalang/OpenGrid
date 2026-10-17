"""Curved boundary elements (Phase 11g): the high-order nodes of the boundary moved onto the
geometry, the displacement extended into the mesh, and every changed element checked.

Method (a posteriori curving of the straight-sided elements of ``model/elements.py``):

1. Targets on the geometry. 2D: the nodes on each boundary connector's edges are projected
   onto the connector's curves (closest point). 3D: the nodes on edges along a face domain's
   connectors go onto those connectors' curves; the other nodes of its boundary faces onto
   the face domain's surface, or, for a structured face domain or block face without a
   surface (linear TFI in space), to the degree-p interpolation of its own grid (which lies
   on the surface it was built on) at the node's fractional grid index; planar unstructured
   faces stay. Periodic patches stay straight (their node correspondence is kept).
2. The displacement is extended into the mesh (``core/src/highorder/curving.cpp``): the
   rest of the boundary is fixed, free nodes within RADIUS_FACTOR boundary-element lengths of
   a moved node follow by a weighted graph Laplacian, so boundary-layer stacks bend with the
   wall; elements far from curved boundaries stay straight-sided.
3. Validity: the Jacobian of every element with a moved node at the order-(p + 2) lattice
   (``model/shape.py``), relative to its straight-sided orientation; the scaled Jacobian is
   min / max. Elements whose Jacobian is not positive everywhere are reported (entity, kind,
   index, location), not fixed.

References: S. J. Sherwin, J. Peiro, "Mesh generation in curvilinear domains using high-order
elements", Int. J. Numer. Meth. Engng 53 (2002) 207-223; P.-O. Persson, J. Peraire, AIAA
2009-949; F. Bassi, S. Rebay, J. Comput. Phys. 138 (1997) 251-285 (why curved walls).
"""

from functools import cache

import numpy as np

from opengrid import _core
from opengrid.io import lagrange
from opengrid.model import shape
from opengrid.validation import Issue, ValidationError

RADIUS_FACTOR = 3.0
MAX_REPORTED = 10
BEZIER_RATIO = 0.3  # rigorous validity for elements whose sampled scaled Jacobian is below


# ------------------------------------------------------------------ mesh topology helpers
@cache
def _neighbour_pairs(kind: str, p: int) -> np.ndarray:
    """Local index pairs of neighbouring lattice nodes (one lattice step; simplex-type
    lattices also along the diagonal directions)."""
    x = lagrange.as_array(lagrange.nodes("gmsh", kind, p))
    limit = np.sqrt(2.0) / p
    if kind in ("line", "quad", "hexa"):
        limit = 1.0 / p
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    a, b = np.nonzero((d > 0) & (d <= limit * (1 + 1e-9)))
    keep = a < b
    return np.column_stack([a[keep], b[keep]])


@cache
def _face_local_nodes(kind: str, p: int) -> list:
    """Per face of a 3D kind (lagrange.GMSH order): local indices of the nodes on it."""
    nodes = lagrange.nodes("gmsh", kind, p)
    support = [frozenset(k for k, w in enumerate(lagrange.weights(kind, n)) if w) for n in nodes]
    out = []
    for face in lagrange.GMSH[kind]["faces"]:
        vertices = frozenset(face)
        out.append(np.array([j for j, s in enumerate(support) if s <= vertices], dtype=np.int64))
    return out


@cache
def _edge_local_nodes(kind: str, p: int) -> list:
    nodes = lagrange.nodes("gmsh", kind, p)
    support = [frozenset(k for k, w in enumerate(lagrange.weights(kind, n)) if w) for n in nodes]
    out = []
    for edge in lagrange.GMSH[kind]["edges"]:
        vertices = frozenset(edge)
        out.append(np.array([j for j, s in enumerate(support) if s <= vertices], dtype=np.int64))
    return out


def lattice_edges(mesh) -> np.ndarray:
    """(m, 2) neighbouring node pairs over every element (duplicates possible)."""
    parts = []
    for kind, conn in mesh.cells.items():
        pairs = _neighbour_pairs(kind, mesh.order)
        parts.append(np.stack([conn[:, pairs[:, 0]], conn[:, pairs[:, 1]]], axis=-1).reshape(-1, 2))
    return np.concatenate(parts)


def boundary_nodes(mesh) -> np.ndarray:
    """Ids of the nodes on the mesh boundary (sides used by one element)."""
    p = mesh.order
    keys, nodes = [], []
    for kind, conn in mesh.cells.items():
        if mesh.dim == 2:
            sides = lagrange.GMSH[kind]["edges"]
            local = _edge_local_nodes(kind, p)
        else:
            sides = lagrange.GMSH[kind]["faces"]
            local = _face_local_nodes(kind, p)
        for index, side in enumerate(sides):
            ids = np.sort(conn[:, list(side)], axis=1)
            if mesh.dim == 3 and ids.shape[1] == 3:
                ids = np.column_stack([ids, np.full(len(ids), -1, dtype=np.int64)])
            keys.append(ids)
            nodes.append(conn[:, local[index]])
    _unique, inverse, counts = np.unique(
        np.vstack(keys), axis=0, return_inverse=True, return_counts=True
    )
    once = counts[inverse.reshape(-1)] == 1
    found = []
    at = 0
    for key, side_nodes in zip(keys, nodes, strict=True):
        found.append(side_nodes[once[at : at + len(key)]].reshape(-1))
        at += len(key)
    return np.unique(np.concatenate(found))


# ------------------------------------------------------------------ targets
class _Targets:
    """Target positions of moved nodes (the first target wins unless overridden)."""

    def __init__(self):
        self.position = {}
        self.radius = {}

    def set(self, node: int, point, radius: float, override: bool = False):
        if node < 0:
            return
        if node in self.position and not override:
            return
        self.position[node] = np.asarray(point, dtype=float)
        self.radius[node] = max(radius, self.radius.get(node, 0.0))


def _closest(geometries, point) -> np.ndarray:
    best = None
    for geometry in geometries:
        result = geometry.project(tuple(point))
        target = result[-2]
        distance = result[-1]
        if best is None or distance < best[1]:
            best = (target, distance)
    return np.asarray(best[0], dtype=float)


def _periodic_connectors(store) -> set:
    out = set()
    for patch in store.entities("patch"):
        if patch.params["patch_type"] != "periodic":
            continue
        for item in patch.params["members"]:
            if "connector" in item:
                out.add(item["connector"])
    return out


def _periodic_sources(store, assembly) -> set[int]:
    """Indices of the 3D boundary face sources in periodic patches (block faces, volume faces,
    the ends and side walls of extruded volumes): they are not curved, so the paired nodes stay
    images of each other."""
    keys = set()
    for patch in store.entities("patch"):
        if patch.params["patch_type"] != "periodic":
            continue
        for item in patch.params["members"]:
            if "block" in item:
                keys.add(("block", item["block"], item["face"]))
            elif "volume" in item:
                keys.add(("volume", item["volume"], item["face"]))
            elif "domain" in item:
                keys.add(("domain", item["domain"]))
            else:
                keys.add(("connector", item["connector"]))
    out = set()
    for index, source in enumerate(assembly.boundary_sources):
        if source["block"] is not None and ("block", *source["block"]) in keys:
            out.add(index)
        if source["domain"] is not None and ("domain", source["domain"]) in keys:
            out.add(index)
        if source["volume"] is not None:
            volume, key = source["volume"]
            if isinstance(key, tuple):
                if key[0] == "top" and ("volume", volume, "top") in keys:
                    out.add(index)
                if key[0] == "side" and ("connector", key[1]) in keys:
                    out.add(index)
            elif ("volume", volume, "bottom") in keys:
                out.add(index)
    return out


def _connector_curves(store, connector_id) -> list:
    return [store.geometry(s["curve"]) for s in store.get(connector_id).params["segments"]]


def _connector_edges(mesh, points, ids, curves, targets: _Targets):
    """Project the edge nodes along consecutive connector nodes onto the curves."""
    a = np.asarray(ids[:-1], dtype=np.int64)
    b = np.asarray(ids[1:], dtype=np.int64)
    nodes = mesh.edge_nodes(a, b)
    length = np.linalg.norm(points[a] - points[b], axis=1)
    for k in range(len(a)):
        for node in nodes[k]:
            if node >= 0:
                targets.set(
                    int(node), _closest(curves, points[node]), RADIUS_FACTOR * length[k], True
                )


def _targets_2d(session, assembly, mesh, targets: _Targets):
    store = session.store
    periodic = _periodic_connectors(store)
    points = mesh.points
    for connector_id, ids in assembly.boundary_connectors:
        if connector_id in periodic:
            continue
        _connector_edges(mesh, points, ids, _connector_curves(store, connector_id), targets)


def _domain_connectors(entity) -> list:
    params = entity.params
    out = []
    for side in params.get("edges", ()):
        out.extend(item["connector"] for item in side)
    for loop in params.get("loops", ()):
        out.extend(item["connector"] for item in loop if "connector" in item)
    return out


def _grid_interpolate(grid, index, p: int) -> np.ndarray:
    """Degree-min(p, n - 1) Lagrange interpolation of a structured grid (n2, n1, 3) at
    fractional indices (m, 2) = (i, j), local stencils of consecutive grid nodes."""
    n2, n1, _ = grid.shape
    out = np.zeros((len(index), 3))
    weights = []
    starts = []
    for axis, n in ((0, n1), (1, n2)):
        q = min(p, n - 1)
        t = index[:, axis]
        start = np.clip(np.round(t - q / 2).astype(np.int64), 0, n - 1 - q)
        local = t - start
        w = np.ones((len(t), q + 1))
        for j in range(q + 1):
            for k in range(q + 1):
                if k != j:
                    w[:, j] *= (local - k) / (j - k)
        weights.append(w)
        starts.append(start)
    wi, wj = weights
    si, sj = starts
    for b in range(wj.shape[1]):
        for a in range(wi.shape[1]):
            out += (wi[:, a] * wj[:, b])[:, None] * grid[sj + b, si + a]
    return out


def _face_node_ids(mesh, kind: str, cells: np.ndarray) -> np.ndarray:
    """(k, n) ids of all nodes of boundary faces given by their vertices (Gmsh order)."""
    p = mesh.order
    ref = lagrange.nodes("gmsh", kind, p)
    nv = len(lagrange.VERTICES[kind])
    out = np.full((len(cells), len(ref)), -1, dtype=np.int64)
    out[:, :nv] = cells
    col = nv
    for a, b in lagrange.GMSH[kind]["edges"]:
        out[:, col : col + p - 1] = mesh.edge_nodes(cells[:, a], cells[:, b])
        col += p - 1
    if col < len(ref):
        from opengrid.model.elements import _lookup

        scale = p * p
        for j in range(col, len(ref)):
            w = [int(v * scale) for v in lagrange.weights(kind, ref[j])]
            ids = cells
            weights = np.broadcast_to(np.array(w, dtype=np.int64), ids.shape)
            order = np.argsort(ids, axis=1, kind="stable")
            ids_sorted = np.take_along_axis(ids, order, axis=1)
            w_sorted = np.take_along_axis(weights, order, axis=1)
            if ids.shape[1] == 3:
                ids_sorted = np.column_stack([ids_sorted, np.full(len(ids), -1, dtype=np.int64)])
                w_sorted = np.column_stack([w_sorted, np.zeros(len(ids), dtype=np.int64)])
            at = _lookup(mesh.sorted_keys("face"), np.column_stack([ids_sorted, w_sorted]))
            found = at >= 0
            out[found, j] = mesh.vertices + len(mesh.edge_keys) + at[found]
    return out


def _targets_3d(session, assembly, mesh, targets: _Targets):
    store = session.store
    p = mesh.order
    points = mesh.points
    periodic = _periodic_connectors(store)
    periodic_sources = _periodic_sources(store, assembly)
    on_connectors = set()
    # Connector edges first (they override face targets): a face domain's border edges (used
    # by one of its faces) lie on its connectors; each goes to the connector whose curves are
    # nearest its midpoint (smoothing groups slide nodes, so coordinates need not match). Edges
    # only between planar faces stay: any position in the plane is on the geometry, and their
    # grid may have slid off the connector on purpose.
    sources_of: dict[tuple[int, int], list[int]] = {}
    oriented: dict[tuple[int, int], tuple[int, int]] = {}
    for index, source in enumerate(assembly.boundary_sources):
        for a, b in _border_edges(source["cells"]):
            key = (min(a, b), max(a, b))
            sources_of.setdefault(key, []).append(index)
            oriented.setdefault(key, (a, b))
    flat = [_planar_source(session, store, source) for source in assembly.boundary_sources]
    candidates = []
    for index, source in enumerate(assembly.boundary_sources):
        domain = source["domain"]
        found = []
        if domain is not None and index not in periodic_sources:
            found = [
                _connector_curves(store, c)
                for c in _domain_connectors(store.get(domain))
                if c not in periodic
            ]
        candidates.append(found)
    for key in sorted(sources_of):
        users = sources_of[key]
        if all(flat[u] for u in users):
            continue
        options = [curves for u in users for curves in candidates[u]]
        if not options:
            continue
        a, b = oriented[key]
        middle = 0.5 * (points[a] + points[b])
        best = None
        for curves in options:
            d = min(curve.project(tuple(middle))[-1] for curve in curves)
            if best is None or d < best[0]:
                best = (d, curves)
        if best[0] > 0.25 * float(np.linalg.norm(points[a] - points[b])):
            continue  # not on any of the faces' connectors
        _connector_edges(mesh, points, [a, b], best[1], targets)
        on_connectors.update(int(n) for n in mesh.edge_nodes([a], [b]).reshape(-1) if n >= 0)
    for index, source in enumerate(assembly.boundary_sources):
        domain = source["domain"]
        if index in periodic_sources:
            continue
        surface = None
        entity = None
        if domain is not None:
            entity = store.get(domain)
            if entity.params.get("surface") is not None:
                surface = store.geometry(entity.params["surface"])
        grid = source["grid"]
        if grid is None and surface is None and entity is not None and entity.kind == "domain":
            grid = session.domain_grid(entity.name)
        if surface is None and grid is None:
            continue  # planar unstructured face
        if surface is None and _planar(grid):
            continue  # a planar grid: its straight-sided faces already lie on the plane
        cells = source["cells"]
        for size, kind in ((3, "triangle"), (4, "quad")):
            rows = [k for k, c in enumerate(cells) if len(c) == size]
            if not rows:
                continue
            verts = np.array([cells[k] for k in rows], dtype=np.int64)
            ids = _face_node_ids(mesh, kind, verts)
            ref = lagrange.nodes("gmsh", kind, p)
            nv = len(lagrange.VERTICES[kind])
            edge_len = np.max(
                [
                    np.linalg.norm(points[verts[:, a]] - points[verts[:, (a + 1) % nv]], axis=1)
                    for a in range(nv)
                ],
                axis=0,
            )
            index = None
            if surface is None:
                index = _grid_index(grid, points, verts, source, rows)
            for j in range(nv, len(ref)):
                w = np.array([float(v) for v in lagrange.weights(kind, ref[j])])
                for r in range(len(verts)):
                    node = int(ids[r, j])
                    if node < 0 or node in on_connectors:
                        continue
                    if surface is not None:
                        target = _closest([surface], points[node])
                    else:
                        frac = w @ index[r]
                        target = _grid_interpolate(grid, frac[None, :], p)[0]
                    targets.set(node, target, RADIUS_FACTOR * edge_len[r])


def _planar(grid) -> bool:
    """True if the grid's points lie in one plane (relative tolerance 1e-10)."""
    x = grid.reshape(-1, 3)
    x = x - x.mean(axis=0)
    extent = float(np.abs(x).max())
    if extent == 0.0:
        return True
    return float(np.linalg.svd(x, compute_uv=False)[-1]) <= 1e-10 * extent * np.sqrt(len(x))


def _planar_source(session, store, source) -> bool:
    """True if a boundary face source is planar: no surface, and a planar grid or an
    unstructured face (planar when it has no surface)."""
    domain = source["domain"]
    entity = None
    if domain is not None:
        entity = store.get(domain)
        if entity.params.get("surface") is not None:
            return False
    grid = source["grid"]
    if grid is None and entity is not None and entity.kind == "domain":
        grid = session.domain_grid(entity.name)
    if grid is None:
        return True
    return _planar(grid)


def _border_edges(cells) -> list[tuple[int, int]]:
    """Edges used by one face of the set (in the faces' own direction)."""
    count: dict[tuple[int, int], int] = {}
    first: dict[tuple[int, int], tuple[int, int]] = {}
    for cell in cells:
        n = len(cell)
        for k in range(n):
            a = int(cell[k])
            b = int(cell[(k + 1) % n])
            key = (min(a, b), max(a, b))
            count[key] = count.get(key, 0) + 1
            first.setdefault(key, (a, b))
    return [first[key] for key in sorted(count) if count[key] == 1]


def _grid_index(grid, points, verts, source, rows) -> np.ndarray:
    """(k, 4, 2) fractional grid indices (i, j) of the vertices of quads on a structured grid."""
    if source["cell_ij"] is not None:
        base = np.asarray(source["cell_ij"])[rows]
        offsets = np.array([(0, 0), (1, 0), (1, 1), (0, 1)])
        return base[:, None, :] + offsets[None, :, :]
    n2, n1, _ = grid.shape
    lookup = {grid[j, i].tobytes(): (i, j) for j in range(n2) for i in range(n1)}
    out = np.zeros((len(verts), verts.shape[1], 2))
    for r, row in enumerate(verts):
        for c, v in enumerate(row):
            out[r, c] = lookup[points[v].tobytes()]
    return out


# ------------------------------------------------------------------ curving
def _assess(mesh, straight, changed) -> dict:
    """Per kind: (rows of the elements with a changed node, their scaled Jacobian curved and
    straight-sided), sampled at the order-(p + 2) lattice."""
    out = {}
    for kind, conn in mesh.cells.items():
        rows = np.flatnonzero(changed[conn].any(axis=1))
        if len(rows) == 0:
            continue
        det = shape.jacobians(kind, mesh.order, mesh.points[conn[rows]], mesh.dim)
        before = shape.jacobians(kind, mesh.order, straight[conn[rows]], mesh.dim)
        sign = np.sign(shape.jacobians_at_centre(kind, mesh.order, straight[conn[rows]], mesh.dim))
        signed = det * sign[:, None]
        signed_before = before * sign[:, None]
        ratio = signed.min(axis=1) / np.maximum(signed.max(axis=1), 1e-300)
        ratio_straight = signed_before.min(axis=1) / np.maximum(signed_before.max(axis=1), 1e-300)
        out[kind] = (rows, ratio, ratio_straight)
    return out


def curve_boundary(session, assembly, mesh, repair: bool = True) -> dict:
    """Move the boundary nodes of ``mesh`` (in place) onto the geometry and extend the
    displacement; ``repair``: optimise the elements that fold or lose most of their quality
    (``model/repair.py``). Returns the report also stored in ``mesh.curving``."""
    """Move the boundary nodes of ``mesh`` (in place) onto the geometry and extend the
    displacement; ``repair``: optimise the elements that fold or lose most of their quality
    (``model/repair.py``). Returns the report also stored in ``mesh.curving``."""
    report = {
        "moved": 0,
        "max_displacement": 0.0,
        "solved": 0,
        "iterations": 0,
        "curved_elements": 0,
        "min_scaled_jacobian": None,  # of the curved elements (None: none curved)
        "min_scaled_jacobian_straight": None,
    }
    if mesh.order < 2:
        mesh.curving = report
        return report
    targets = _Targets()
    if mesh.dim == 2:
        _targets_2d(session, assembly, mesh, targets)
    else:
        _targets_3d(session, assembly, mesh, targets)
    n = len(mesh.points)
    fixed = np.zeros(n, dtype=np.uint8)
    fixed[boundary_nodes(mesh)] = 1
    displacement = np.zeros((n, 3))
    radius = np.zeros(n)
    for node, target in targets.position.items():
        fixed[node] = 1
        displacement[node] = target - mesh.points[node]
        radius[node] = targets.radius[node]
    moved = np.linalg.norm(displacement, axis=1) > 0
    report["moved"] = int(moved.sum())
    if not moved.any():
        mesh.curving = report
        return report
    report["max_displacement"] = float(np.linalg.norm(displacement, axis=1).max())
    edges = lattice_edges(mesh)
    result, solved, _seeds, iterations = _core.extend_displacement(
        np.ascontiguousarray(mesh.points),
        np.ascontiguousarray(edges),
        fixed,
        np.ascontiguousarray(displacement),
        radius,
    )
    report["solved"] = int(solved)
    report["iterations"] = int(iterations)
    straight = mesh.points.copy()
    changed = np.linalg.norm(result, axis=1) > 0
    mesh.points = mesh.points + result
    assessment = _assess(mesh, straight, changed)
    report["curved_elements"] = int(sum(len(a[0]) for a in assessment.values()))
    if repair and assessment:
        from opengrid.model import repair as repair_module

        flagged = {}
        for kind, (rows, ratio, ratio_straight) in assessment.items():
            mask = repair_module.flagged_elements(ratio, ratio_straight)
            if mask.any():
                flagged[kind] = rows[mask]
        if flagged:
            report["repair"] = repair_module.repair(mesh, straight, flagged, fixed)
            changed |= np.linalg.norm(mesh.points - straight, axis=1) > 0
            assessment = _assess(mesh, straight, changed)
    worst, worst_straight, issues = 1.0, 1.0, []
    checked = undetermined = 0
    for kind, (rows, ratio, ratio_straight) in assessment.items():
        worst = min(worst, float(ratio.min()))
        worst_straight = min(worst_straight, float(ratio_straight.min()))
        conn = mesh.cells[kind]
        low = np.flatnonzero(ratio < BEZIER_RATIO)
        if len(low) == 0:
            continue
        from opengrid.model import repair as repair_module

        sign = np.sign(
            shape.jacobians_at_centre(kind, mesh.order, straight[conn[rows[low]]], mesh.dim)
        )
        status, bound = repair_module.validity(
            kind, mesh.order, mesh.points[conn[rows[low]]], mesh.dim, sign
        )
        checked += len(low)
        undetermined += int((status == 0).sum())
        for k in np.flatnonzero(status <= 0):
            if len(issues) >= MAX_REPORTED:
                break
            element = int(rows[low[k]])
            c = mesh.points[conn[element]].mean(axis=0)
            where = f"({c[0]:.4g}, {c[1]:.4g}, {c[2]:.4g})"
            if status[k] < 0:
                issues.append(
                    Issue(
                        f"{kind} {element}",
                        f"folds after curving (scaled Jacobian {ratio[low[k]]:.3f}, Bezier "
                        f"bound {bound[k]:.3f}) near {where}",
                    )
                )
            else:
                issues.append(
                    Issue(
                        f"{kind} {element}",
                        f"validity undetermined after curving (sampled scaled Jacobian "
                        f"{ratio[low[k]]:.3f}, Bezier bound {bound[k]:.3f}) near {where}",
                    )
                )
    report["min_scaled_jacobian"] = worst
    report["min_scaled_jacobian_straight"] = worst_straight  # the same elements, straight-sided
    report["validity"] = {"checked": checked, "undetermined": undetermined}
    mesh.curving = report
    if issues:
        raise ValidationError(issues)
    return report
