"""Straight-sided high-order elements from the mesh cells (Phase 11f).

Method: every cell of the linear mesh (2D: triangles and quadrilaterals; 3D: tetrahedra,
pyramids, prisms and hexahedra) becomes one order-p Lagrange element with the linear cell's
geometry. Its nodes are the equispaced lattice points of the reference element placed by the
linear map: barycentric for triangles and tetrahedra, bilinear / trilinear for quadrilaterals
and hexahedra, triangle x segment for prisms, and the collapsed map for pyramids
(``io/lagrange.py``). Edges stay straight and triangular faces planar. A quadrilateral face is
the bilinear surface of its four vertices, which is planar when they are coplanar.

Neighbouring elements share their edge and face nodes exactly. A node is identified by the
vertices of the edge or face it lies on and its integer lattice weights over them, so the
merge needs no tolerance and a shared node's coordinates are computed once, from that key.
Boundary lines (2D) and faces (3D) of each patch become order-p elements on the same nodes.
Periodic pairs carry their node correspondences over to the high-order nodes, through the same
keys. Node numbering: the linear mesh's vertices first (same indices), then the edge nodes,
the face nodes and the element interiors, each in a fixed sorted order.

Connectivity is held in Gmsh's node order (the canonical one here); ``connectivity(kind, fmt)``
gives VTK's or CGNS's order.
"""

from dataclasses import dataclass, field

import numpy as np

from opengrid.io import lagrange

# VTK linear cell type -> element kind (the assemblies' vertex orders are the formats' own).
VTK_KIND = {5: "triangle", 9: "quad", 10: "tetra", 14: "pyramid", 13: "prism", 12: "hexa"}
KIND_ORDER = ("line", "triangle", "quad", "tetra", "pyramid", "prism", "hexa")
MAX_ORDER = 6  # verified against Gmsh and VTK; CGNS defines orders up to 4


@dataclass
class ElementMesh:
    dim: int
    order: int
    points: np.ndarray  # (n, 3)
    vertices: int  # points[:vertices] are the linear mesh's points, same indices
    cells: dict  # kind -> (count, nodes) int64 node ids, Gmsh node order
    cell_zone: dict  # kind -> (count,) zone index
    faces: dict  # kind -> (count, nodes): boundary elements (2D: lines; 3D: triangles, quads)
    face_patch: dict  # kind -> (count,) patch index
    zone_names: list
    patch_names: list
    patch_types: list
    # (source patch, target patch, transform, target node ids, matching source node ids)
    periodic: list = field(default_factory=list)

    @property
    def elements(self) -> int:
        return int(sum(len(c) for c in self.cells.values()))

    @property
    def boundary_elements(self) -> int:
        return int(sum(len(c) for c in self.faces.values()))

    def connectivity(self, kind: str, fmt: str = "gmsh", boundary: bool = False) -> np.ndarray:
        """Node ids of the elements of one kind in the node order of ``fmt`` ("gmsh", "vtk",
        "cgns")."""
        table = self.cells
        if boundary:
            table = self.faces
        return table[kind][:, lagrange.permutation(fmt, kind, self.order)]

    def element_nodes(self, kind: str, nodes: str = "equispaced") -> np.ndarray:
        """(count, n, 3) node coordinates of each element of one kind: "equispaced" (Gmsh
        order, the nodes above), or "gll" / "gauss" (Gauss-Lobatto-Legendre / Gauss-Legendre
        points of the reference square / cube mapped by the bilinear / trilinear map; tensor
        order, i fastest; quadrilaterals and hexahedra only)."""
        conn = self.cells[kind]
        if nodes == "equispaced":
            return self.points[conn]
        if kind not in ("quad", "hexa"):
            raise ValueError(
                f"the {nodes!r} node set is defined for quadrilaterals and hexahedra (tensor "
                f"products), not for {kind} elements: export the equispaced nodes"
            )
        from opengrid import _core

        n = self.order + 1
        t = (np.asarray(_core.nodes_1d(nodes, n)) + 1.0) / 2.0  # on [0, 1]
        d = lagrange.DIM[kind]
        grid = np.meshgrid(*([t] * d), indexing="ij")
        ref = np.column_stack([g.transpose().reshape(-1) for g in grid])  # i fastest
        weights = np.array([[float(w) for w in lagrange.weights(kind, x)] for x in ref])
        nv = len(lagrange.VERTICES[kind])
        return np.einsum("jk,ekd->ejd", weights, self.points[conn[:, :nv]])

    def counts(self) -> dict:
        out = {kind: len(c) for kind, c in self.cells.items()}
        for kind, c in self.faces.items():
            out[f"boundary {kind}"] = len(c)
        return out


def _gather(assembly, mask) -> np.ndarray:
    """(count, vertices) connectivity of the selected cells (one cell type)."""
    index = np.flatnonzero(mask)
    sizes = assembly.offsets[index + 1] - assembly.offsets[index]
    size = int(sizes[0])
    return assembly.connectivity[assembly.offsets[index][:, None] + np.arange(size)].astype(
        np.int64
    )


def _integral(value) -> int:
    """An exact lattice weight scaled to an integer (the node keys rely on it)."""
    if value.denominator != 1:
        raise ValueError(f"lattice weight {value} is not integral")
    return int(value)


class _Reference:
    """One kind's order-p reference element in a mesh of dimension ``dim``: lattice weights and
    node classes (a node inside a 3D mesh's boundary triangle or quad is a face node)."""

    def __init__(self, kind: str, p: int, dim: int):
        self.kind = kind
        nodes = lagrange.nodes("gmsh", kind, p)
        weights = [lagrange.weights(kind, x) for x in nodes]
        self.float_weights = np.array([[float(w) for w in row] for row in weights])
        self.vertex = []  # (node, local vertex)
        self.edge = []  # (node, a, b, p * weight on a)
        self.face = []  # (node, local vertices, p^2 * weights)
        self.interior = []  # nodes
        scale = p * p
        for j, row in enumerate(weights):
            support = tuple(k for k, w in enumerate(row) if w != 0)
            if len(support) == 1:
                self.vertex.append((j, support[0]))
            elif len(support) == 2:
                a, b = support
                self.edge.append((j, a, b, _integral(row[a] * p)))
            elif dim == 3 and (len(support) < len(row) or lagrange.DIM[kind] < dim):
                self.face.append((j, support, tuple(_integral(row[k] * scale) for k in support)))
            else:
                self.interior.append(j)
        self.count = len(nodes)


def _edge_keys(ref, conn, p) -> list:
    out = []
    for _, a, b, k in ref.edge:
        va = conn[:, a]
        vb = conn[:, b]
        swap = va > vb
        lo = np.where(swap, vb, va)
        hi = np.where(swap, va, vb)
        weight = np.where(swap, p - k, k)
        out.append(np.column_stack([lo, hi, weight]))
    return out


def _face_keys(ref, conn) -> list:
    out = []
    for _, support, weights in ref.face:
        ids = conn[:, list(support)]
        w = np.broadcast_to(np.asarray(weights, dtype=np.int64), ids.shape)
        order = np.argsort(ids, axis=1, kind="stable")
        ids = np.take_along_axis(ids, order, axis=1)
        w = np.take_along_axis(w, order, axis=1)
        if ids.shape[1] == 3:
            ids = np.column_stack([ids, np.full(len(ids), -1, dtype=np.int64)])
            w = np.column_stack([w, np.zeros(len(w), dtype=np.int64)])
        out.append(np.column_stack([ids, w]))
    return out


def build_elements(assembly, order: int) -> ElementMesh:
    """Order-``order`` straight-sided elements of a linear mesh assembly
    (``io/assembly.py`` / ``io/assembly3d.py``)."""
    if not 1 <= order <= MAX_ORDER:
        raise ValueError(f"order must be 1 .. {MAX_ORDER}, got {order}")
    p = order
    dim = assembly.dimension
    points = np.asarray(assembly.points, dtype=float)
    nv = len(points)
    # Element groups: (role, kind, linear connectivity, tag).
    groups = []
    for kind in KIND_ORDER:
        for vtk_type, k in VTK_KIND.items():
            if k != kind:
                continue
            mask = assembly.cell_type == vtk_type
            if np.any(mask):
                groups.append(("cell", kind, _gather(assembly, mask), assembly.cell_zone[mask]))
    if dim == 2:
        if len(assembly.lines):
            groups.append(
                (
                    "face",
                    "line",
                    np.asarray(assembly.lines, dtype=np.int64),
                    np.asarray(assembly.line_patch),
                )
            )
    else:
        cells = assembly.boundary_cells
        patch = np.asarray(assembly.boundary_patch)
        for size, kind in ((3, "triangle"), (4, "quad")):
            rows = [k for k, c in enumerate(cells) if len(c) == size]
            if rows:
                conn = np.array([cells[k] for k in rows], dtype=np.int64)
                groups.append(("face", kind, conn, patch[rows]))
    refs = {kind: _Reference(kind, p, dim) for _, kind, _, _ in groups}

    # Shared nodes: unique edge and face keys over all groups.
    edge_rows = [k for _, kind, conn, _ in groups for k in _edge_keys(refs[kind], conn, p)]
    face_rows = [k for _, kind, conn, _ in groups for k in _face_keys(refs[kind], conn)]
    edge_keys = np.zeros((0, 3), dtype=np.int64)
    edge_inverse = np.zeros(0, dtype=np.int64)
    if edge_rows:
        edge_keys, edge_inverse = np.unique(np.vstack(edge_rows), axis=0, return_inverse=True)
    face_keys = np.zeros((0, 8), dtype=np.int64)
    face_inverse = np.zeros(0, dtype=np.int64)
    if face_rows:
        face_keys, face_inverse = np.unique(np.vstack(face_rows), axis=0, return_inverse=True)
    edge_inverse = edge_inverse.reshape(-1)
    face_inverse = face_inverse.reshape(-1)
    first_edge = nv
    first_face = first_edge + len(edge_keys)
    first_interior = first_face + len(face_keys)

    # Coordinates of the shared nodes, from their keys.
    lo, hi, k = edge_keys.T
    edge_xyz = (k[:, None] * points[lo] + (p - k)[:, None] * points[hi]) / p
    face_xyz = np.zeros((len(face_keys), 3))
    for m in range(4):
        ids = face_keys[:, m]
        used = ids >= 0
        face_xyz[used] += face_keys[used, 4 + m][:, None] * points[ids[used]]
    face_xyz /= p * p

    cells, cell_zone, faces, face_patch = {}, {}, {}, {}
    interior_xyz = []
    next_interior = first_interior
    edge_at = 0
    face_at = 0
    for role, kind, conn, tag in groups:
        ref = refs[kind]
        ne = len(conn)
        out = np.empty((ne, ref.count), dtype=np.int64)
        for j, v in ref.vertex:
            out[:, j] = conn[:, v]
        for j, _, _, _ in ref.edge:
            out[:, j] = first_edge + edge_inverse[edge_at : edge_at + ne]
            edge_at += ne
        for j, _, _ in ref.face:
            out[:, j] = first_face + face_inverse[face_at : face_at + ne]
            face_at += ne
        if ref.interior:
            n_int = len(ref.interior)
            ids = next_interior + np.arange(ne * n_int).reshape(ne, n_int)
            out[:, ref.interior] = ids
            next_interior += ne * n_int
            w = ref.float_weights[ref.interior]
            interior_xyz.append(np.einsum("jk,ekd->ejd", w, points[conn]).reshape(-1, 3))
        if role == "cell":
            cells[kind] = out
            cell_zone[kind] = np.asarray(tag, dtype=np.int32)
        else:
            faces[kind] = out
            face_patch[kind] = np.asarray(tag, dtype=np.int32)
    all_points = np.vstack([points, edge_xyz, face_xyz, *interior_xyz])
    mesh = ElementMesh(
        dim=dim,
        order=p,
        points=all_points,
        vertices=nv,
        cells=cells,
        cell_zone=cell_zone,
        faces=faces,
        face_patch=face_patch,
        zone_names=list(assembly.zone_names),
        patch_names=list(assembly.patch_names),
        patch_types=list(assembly.patch_types),
    )
    mesh.periodic = _periodic(assembly, mesh, refs, edge_keys, face_keys)
    return mesh


def _periodic(assembly, mesh, refs, edge_keys, face_keys) -> list:
    """Node correspondences of each periodic pair extended to the high-order nodes of the
    target patch's boundary elements, through the vertex correspondence and the node keys."""
    if not assembly.periodic:
        return []
    p = mesh.order
    first_edge = mesh.vertices
    first_face = first_edge + len(edge_keys)
    edge_index = {tuple(row): first_edge + k for k, row in enumerate(edge_keys.tolist())}
    face_index = {tuple(row): first_face + k for k, row in enumerate(face_keys.tolist())}
    out = []
    for source, target, transform, target_ids, source_ids in assembly.periodic:
        vmap = dict(
            zip(np.asarray(target_ids).tolist(), np.asarray(source_ids).tolist(), strict=True)
        )
        pairs = dict(vmap)
        for kind, conn in mesh.faces.items():
            ref = refs[kind]
            rows = conn[mesh.face_patch[kind] == target]
            for row in rows.tolist():
                verts = [row[j] for j, _ in sorted(ref.vertex, key=lambda t: t[1])]
                for j, a, b, k in ref.edge:
                    pairs[row[j]] = _mapped_edge(edge_index, vmap, verts[a], verts[b], k, p)
                for j, support, weights in ref.face:
                    pairs[row[j]] = _mapped_face(
                        face_index, vmap, [verts[s] for s in support], weights
                    )
        targets = np.array(sorted(pairs), dtype=np.int64)
        sources = np.array([pairs[t] for t in targets.tolist()], dtype=np.int64)
        out.append((source, target, transform, targets, sources))
    return out


def _mapped_edge(index, vmap, a, b, k, p) -> int:
    ma = vmap[a]
    mb = vmap[b]
    if ma < mb:
        return index[(ma, mb, k)]
    return index[(mb, ma, p - k)]


def _mapped_face(index, vmap, verts, weights) -> int:
    mapped = sorted(zip((vmap[v] for v in verts), weights, strict=True))
    ids = [m for m, _ in mapped]
    w = [x for _, x in mapped]
    if len(ids) == 3:
        ids.append(-1)
        w.append(0)
    return index[(*ids, *w)]
