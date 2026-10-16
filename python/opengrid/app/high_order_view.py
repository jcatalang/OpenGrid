"""Viewport overlay of a high-order mesh (``model/elements.py``): the element edges sampled
along their curves and the high-order nodes, so curved elements can be judged before export.

2D: every distinct edge of the cells; 3D: the edges of the boundary elements (the faces of the
patches), or of all cells when there are no patches. Each edge is the degree-p curve through
its p + 1 nodes (Gmsh line order), sampled at ``SAMPLES`` points.
"""

import numpy as np
import pyvista as pv

from opengrid.io import lagrange
from opengrid.model import shape

SAMPLES = 9


def element_edges(mesh) -> np.ndarray:
    """(k, p + 1) node ids of the distinct edges to draw: the two ends, then the interior
    nodes from the first end."""
    p = mesh.order
    tables = mesh.cells
    if mesh.dim == 3 and mesh.faces:
        tables = mesh.faces
    rows = []
    for kind, conn in tables.items():
        if kind == "line":
            continue
        nv = len(lagrange.VERTICES[kind])
        for k, (a, b) in enumerate(lagrange.GMSH[kind]["edges"]):
            interior = conn[:, nv + k * (p - 1) : nv + (k + 1) * (p - 1)]
            rows.append(np.column_stack([conn[:, a], conn[:, b], interior]))
    if not rows:
        return np.zeros((0, p + 1), dtype=np.int64)
    edges = np.vstack(rows)
    _unique, first = np.unique(np.sort(edges[:, :2], axis=1), axis=0, return_index=True)
    return edges[np.sort(first)]


def edge_polydata(mesh, samples: int = SAMPLES) -> pv.PolyData | None:
    """The edges as polylines sampled along their curves (None without edges)."""
    edges = element_edges(mesh)
    if len(edges) == 0:
        return None
    values, _grads = shape.basis("line", mesh.order, np.linspace(0.0, 1.0, samples)[:, None])
    points = np.einsum("qn,enc->eqc", values, mesh.points[edges]).reshape(-1, 3)
    index = np.arange(len(edges) * samples).reshape(len(edges), samples)
    lines = np.column_stack([np.full(len(edges), samples), index]).reshape(-1)
    return pv.PolyData(points, lines=lines)


def node_polydata(mesh) -> pv.PolyData | None:
    """The high-order nodes (those beyond the linear mesh's vertices)."""
    nodes = mesh.points[mesh.vertices :]
    if len(nodes) == 0:
        return None
    return pv.PolyData(np.ascontiguousarray(nodes))
