"""High-order (curved) elements from the structured grids (Phase 11e).

Method: a structured domain (2D) or block (3D) whose interval counts are multiples of the
order p is the equispaced node set of order-p Lagrange quadrilaterals (hexahedra): p x p
(x p) grid cells form one element, so the smoothing that placed the grid placed every
high-order node, boundary nodes lie on the geometry and neighbours share their edge and face
nodes (core: ``core/src/highorder/elements.cpp``). Other node sets are the element
mappings evaluated at their points: Gauss-Lobatto-Legendre (continuous Galerkin / spectral
elements, DGSEM with boundary points) and Gauss-Legendre (discontinuous Galerkin solution
points). With ``snap`` the GLL nodes on boundary sides are projected onto their curves (2D:
sides not shared with another domain) or the surfaces of their face domains (3D, face
interiors), so they lie on the geometry as the grid nodes do. Validity: the Jacobian of
each element mapping at (p + 2)^d GLL points; scaled Jacobian = min / max (not positive:
invalid).

References: S. Sherwin, J. Peiro, Int. J. Numer. Meth. Engng 53 (2002) 207-223; J. S.
Hesthaven, T. Warburton, Nodal Discontinuous Galerkin Methods (2008); D. A. Kopriva,
Implementing Spectral Methods for PDEs (2009).
"""

from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.validation import Issue, ValidationError

NODE_SETS = ("equispaced", "gll", "gauss")
SIDES_2D = ("bottom", "right", "top", "left")


@dataclass
class HighOrderMesh:
    dim: int
    order: int
    nodes: str  # node set
    xyz: np.ndarray  # (elements, (order + 1)^dim, 3), tensor order, i fastest
    min_jacobian: np.ndarray  # (elements,)
    max_jacobian: np.ndarray
    owners: list = field(
        default_factory=list
    )  # (entity name, first element, count, (ne_i, ne_j, ne_k))
    snapped: int = 0  # boundary nodes projected onto the geometry
    max_snap: float = 0.0  # largest projection distance

    @property
    def elements(self) -> int:
        return len(self.xyz)

    @property
    def scaled_jacobian(self) -> np.ndarray:
        return self.min_jacobian / self.max_jacobian

    def merged(self, tolerance: float = 1e-10):
        """(points (n, 3), connectivity (elements, (order + 1)^dim)): coincident nodes of
        neighbouring elements merged (shared edge and face nodes)."""
        flat = self.xyz.reshape(-1, 3)
        scale = max(float(np.ptp(flat, axis=0).max()), 1.0)
        keys = np.round(flat / (tolerance * scale)).astype(np.int64)
        _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
        order = np.argsort(first, kind="stable")
        rank = np.empty(len(first), dtype=np.int64)
        rank[order] = np.arange(len(first))
        ids = rank[inverse.reshape(-1)].reshape(self.xyz.shape[:2])
        return flat[first[order]], ids


def _fail(where: str, message: str):
    raise ValidationError([Issue(where, message)])


def _side_curves(store, domain) -> list[list]:
    """Curves (C++ geometry) of the connectors on each boundary side (bottom, right, top,
    left); none for a side shared with another domain (an interface: its nodes float in a
    smoothing group, off the connector's original curve)."""
    out = []
    for side in domain.params["edges"]:
        curves = []
        for item in side:
            connector = item["connector"]
            users = [
                e
                for e in store.dependents(connector)
                if e.kind in ("domain", "unstructured_domain") and e.id != domain.id
            ]
            if users:
                curves = []
                break
            for segment in store.get(connector).params["segments"]:
                curves.append(store.geometry(segment["curve"]))
        out.append(curves)
    return out


def _snap_2d(xyz, ne, order, curves) -> tuple[int, float]:
    """Project the nodes on the domain's boundary sides onto their curves."""
    n1 = order + 1
    count = 0
    worst = 0.0
    nodes = xyz.reshape(ne[1], ne[0], n1, n1, 3)  # (ej, ei, b, a, 3)
    for side, side_curves in enumerate(curves):
        if not side_curves:
            continue
        if side == 0:
            view = nodes[0, :, 0, :]  # (ei, a)
        elif side == 1:
            view = nodes[:, -1, :, -1]  # (ej, b)
        elif side == 2:
            view = nodes[-1, :, -1, :]
        else:
            view = nodes[:, 0, :, 0]
        flat = view.reshape(-1, 3)
        for k in range(len(flat)):
            best = None
            for curve in side_curves:
                _t, point, distance = curve.project(tuple(flat[k]))
                if best is None or distance < best[1]:
                    best = (point, distance)
            worst = max(worst, best[1])
            flat[k] = best[0]
            count += 1
        if side == 0:
            nodes[0, :, 0, :] = flat.reshape(view.shape)
        elif side == 1:
            nodes[:, -1, :, -1] = flat.reshape(view.shape)
        elif side == 2:
            nodes[-1, :, -1, :] = flat.reshape(view.shape)
        else:
            nodes[:, 0, :, 0] = flat.reshape(view.shape)
    return count, worst


def _face_view(nodes, face: str):
    """Element nodes (ek, ej, ei, c, b, a, 3) on a block face, without the face's border
    nodes (those on block edges stay as interpolated)."""
    if face == "imin":
        return nodes[:, :, 0, :, :, 0]
    if face == "imax":
        return nodes[:, :, -1, :, :, -1]
    if face == "jmin":
        return nodes[:, 0, :, :, 0, :]
    if face == "jmax":
        return nodes[:, -1, :, :, -1, :]
    if face == "kmin":
        return nodes[0, :, :, 0, :, :]
    return nodes[-1, :, :, -1, :, :]


def _snap_3d(store, xyz, ne, order, faces) -> tuple[int, float]:
    n1 = order + 1
    nodes = xyz.reshape(ne[2], ne[1], ne[0], n1, n1, n1, 3)
    count = 0
    worst = 0.0
    for face, domain_id in faces.items():
        if domain_id is None:
            continue
        surface_id = store.get(domain_id).params.get("surface")
        if surface_id is None:
            continue
        surface = store.geometry(surface_id)
        view = _face_view(nodes, face)  # (e1, e2, q1, q2, 3)
        e1, e2, q1, q2 = view.shape[:4]
        for i in range(e1):
            for j in range(e2):
                for a in range(q1):
                    for b in range(q2):
                        on_border = (i == 0 and a == 0) or (i == e1 - 1 and a == q1 - 1)
                        on_border = (
                            on_border or (j == 0 and b == 0) or (j == e2 - 1 and b == q2 - 1)
                        )
                        if on_border:
                            continue
                        _u, _v, point, distance = surface.project(tuple(view[i, j, a, b]))
                        view[i, j, a, b] = point
                        worst = max(worst, distance)
                        count += 1
    return count, worst


def high_order(session, names=None, order: int = 2, nodes: str = "equispaced", snap: bool = True):
    """Order-``order`` curved elements of structured domains (2D) or blocks (3D); see the
    module docstring. ``names``: the domains or blocks (default: all of them)."""
    store = session.store
    where = "high_order"
    if isinstance(order, bool) or not isinstance(order, int) or order < 1:
        _fail(where, f"order must be a whole number >= 1, got {order!r}")
    if nodes not in NODE_SETS:
        _fail(where, f"nodes must be one of {', '.join(NODE_SETS)}, got {nodes!r}")
    three_d = session.is_3d()
    kind = "domain"
    if three_d:
        kind = "block"
    if names is None:
        names = session.names(kind)
    if isinstance(names, str):
        names = [names]
    if not names:
        _fail(where, f"no {kind}s to convert")
    parts = []
    owners = []
    snapped = 0
    worst = 0.0
    first = 0
    for name in names:
        entity = store.get(name)
        if entity.kind != kind:
            _fail(where, f"{entity.name!r} is a {entity.kind}, not a {kind}")
        if three_d:
            grid = np.ascontiguousarray(session.block_grid(name), dtype=float)
        else:
            grid = np.ascontiguousarray(session.domain_grid(name), dtype=float)
        try:
            xyz, lo, hi, ne = _core.high_order_elements(grid, order, nodes)
        except ValueError as error:
            _fail(f"{where}: {entity.name}", str(error))
        xyz = np.array(xyz)
        if snap and nodes == "gll" and order > 1:
            if three_d:
                faces = session.blocks.evaluate(entity.id, raw=True)["faces"]
                count, far = _snap_3d(store, xyz, ne, order, faces)
            else:
                count, far = _snap_2d(xyz, ne, order, _side_curves(store, entity))
            snapped += count
            worst = max(worst, far)
        parts.append((xyz, np.asarray(lo), np.asarray(hi)))
        owners.append((entity.name, first, len(xyz), tuple(int(v) for v in ne)))
        first += len(xyz)
    dim = 2
    if three_d:
        dim = 3
    return HighOrderMesh(
        dim=dim,
        order=order,
        nodes=nodes,
        xyz=np.concatenate([p[0] for p in parts]),
        min_jacobian=np.concatenate([p[1] for p in parts]),
        max_jacobian=np.concatenate([p[2] for p in parts]),
        owners=owners,
        snapped=snapped,
        max_snap=worst,
    )
