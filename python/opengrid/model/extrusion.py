"""Extrusion of an unstructured 2D mesh into layers of prisms and hexahedra: a volume swept
from a planar domain (what "Extrude domain" does with an unstructured domain; a structured
domain becomes a block instead, opengrid.model.blocks).

Method: the domain's mesh points are copied layer by layer along the extrusion rule -
translate: p + s v; rotate: p turned about an axis by s * angle (Rodrigues' formula); normal:
p + s * distance * n with n the plane normal of the domain - with s running through the layer
fractions of a connector-style distribution (opengrid._core.distribute). Every triangle
becomes a prism and every quad a hexahedron between consecutive layers, in the node orders of
opengrid.model.volumes.CELL_TYPES (the first face's normal points into the cell); the base
cells are mirrored when the extrusion runs against their normal, so the cell volumes come out
positive. The boundary is returned as outward faces: the bottom (the domain's cells), the top
(their copy at the last layer) and one strip of quads per connector of the domain's loops.
The mesh's boundary nodes come first, in loop order, with each loop's connectors laid end to
end and junction nodes once (opengrid.model.unstructured.loop_points), which is how the side
strips find their edges.
"""

import numpy as np


def layer_points(points, extrusion: dict, fractions, normal=None) -> np.ndarray:
    """(nk, n, 3): the base points at every layer fraction."""
    base = np.asarray(points, dtype=float).reshape(-1, 3)
    mode = extrusion["mode"]
    out = np.empty((len(fractions), len(base), 3))
    if mode == "translate":
        vector = np.asarray(extrusion["vector"], dtype=float)
        for k, s in enumerate(fractions):
            out[k] = base + s * vector
        return out
    if mode == "normal":
        direction = np.asarray(normal, dtype=float) * extrusion["distance"]
        for k, s in enumerate(fractions):
            out[k] = base + s * direction
        return out
    if mode == "rotate":
        axis = np.asarray(extrusion["vector"], dtype=float)
        length = float(np.linalg.norm(axis))
        if length == 0.0:
            raise ValueError("the rotation axis (vector) must not be zero")
        axis = axis / length
        origin = np.asarray(extrusion["axis_point"], dtype=float)
        p = base - origin
        along = np.outer(p @ axis, axis)
        perpendicular = p - along
        cross = np.cross(axis, p)
        for k, s in enumerate(fractions):
            angle = s * extrusion["angle"]
            out[k] = origin + along + perpendicular * np.cos(angle) + cross * np.sin(angle)
        return out
    raise ValueError(f"unknown extrusion mode {mode!r}")


def stacked_cells(cells: np.ndarray, layers: int, n: int) -> np.ndarray:
    """The base cells (m, c) between every pair of consecutive layers: (layers * m, 2 c), the
    lower nodes first, then the same nodes one layer up (n points per layer)."""
    cells = np.asarray(cells, dtype=np.int64)
    m, c = cells.shape
    out = np.empty((layers * m, 2 * c), dtype=np.int64)
    for k in range(layers):
        out[k * m : (k + 1) * m, :c] = cells + k * n
        out[k * m : (k + 1) * m, c:] = cells + (k + 1) * n
    return out


def connector_edges(loops, counts: dict) -> dict[int, list[tuple[int, int]]]:
    """Boundary edges (node a, node b) of the base mesh per connector id, from the domain's
    loops ([{"connector", "reversed"}, ...] each) and every connector's node count. Nodes are
    numbered as the mesh numbers its boundary: loops end to end, each loop its connectors end
    to end with junction nodes once and the closing node not repeated."""
    out: dict[int, list[tuple[int, int]]] = {}
    offset = 0
    for loop in loops:
        total = sum(counts[item["connector"]] - 1 for item in loop)
        position = 0
        for item in loop:
            edges = out.setdefault(item["connector"], [])
            steps = counts[item["connector"]] - 1
            for p in range(steps):
                a = offset + position + p
                b = offset + (position + p + 1) % total
                edges.append((a, b))
            position += steps
        offset += total
    return out


def extrude_mesh(
    points, triangles, quads, loops, counts: dict, extrusion: dict, fractions, normal=None
) -> dict:
    """{"points" (nk n, 3), "prisms" (., 6), "hexes" (., 8), "bottom", "top" (lists of outward
    cells), "sides" {connector id: [outward quads]}, "mirrored", "boundary_points"}."""
    base = np.asarray(points, dtype=float).reshape(-1, 3)
    n = len(base)
    nk = len(fractions)
    layers = layer_points(base, extrusion, fractions, normal)
    tris = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    qds = np.asarray(quads, dtype=np.int64).reshape(-1, 4)
    # The base cells' normals must point to the next layer (the first face of a prism or
    # hexahedron faces into the cell); otherwise the base is mirrored.
    mirrored = False
    if nk > 1:
        up = layers[1] - layers[0]
        sign = 0.0
        if len(tris):
            p = base[tris]
            normals = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
            sign += float(np.einsum("ij,ij->", normals, up[tris].mean(axis=1)))
        if len(qds):
            p = base[qds]
            normals = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
            sign += float(np.einsum("ij,ij->", normals, up[qds].mean(axis=1)))
        mirrored = sign < 0.0
    if mirrored:
        tris = tris[:, [0, 2, 1]]
        qds = qds[:, [0, 3, 2, 1]]
    top_offset = (nk - 1) * n
    bottom = [cell[::-1] for cell in tris] + [cell[::-1] for cell in qds]
    top = [cell + top_offset for cell in tris] + [cell + top_offset for cell in qds]
    sides: dict[int, list[np.ndarray]] = {}
    edges_of = connector_edges(loops, counts)
    for connector, edges in edges_of.items():
        strip = []
        for k in range(nk - 1):
            lo = k * n
            hi = (k + 1) * n
            for a, b in edges:
                if mirrored:
                    quad = (b + lo, a + lo, a + hi, b + hi)
                else:
                    quad = (a + lo, b + lo, b + hi, a + hi)
                strip.append(np.asarray(quad, dtype=np.int64))
        sides[connector] = strip
    boundary = sum(len(edges) for edges in edges_of.values())
    boundary_points = 2 * n + max(nk - 2, 0) * boundary
    return {
        "points": layers.reshape(-1, 3),
        "prisms": stacked_cells(tris, nk - 1, n),
        "hexes": stacked_cells(qds, nk - 1, n),
        "bottom": bottom,
        "top": top,
        "sides": sides,
        "mirrored": mirrored,
        "boundary_points": int(boundary_points),
    }
