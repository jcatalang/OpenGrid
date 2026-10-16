"""Writers for high-order meshes (Phase 11e): VTKHDF Lagrange cells and HDF5 node sets.

* VTKHDF: an UnstructuredGrid of VTK_LAGRANGE_QUADRILATERAL (70) or VTK_LAGRANGE_HEXAHEDRON
  (72) cells; VTK's Lagrange cells interpolate equispaced nodes, so only the equispaced node
  set is written (shared nodes merged). Point order within a cell as VTK's
  vtkHigherOrderQuadrilateral / vtkHigherOrderHexahedron::PointIndexFromIJK: corners, edges,
  faces, interior (checked against VTK in the tests). Cell data: scaled_jacobian,
  min_jacobian, block_id; FieldData block_names.
* HDF5 (``.h5``): any node set, for solvers that take element node coordinates directly:
  ``nodes`` (elements, (p + 1)^d, 3) in tensor order (i fastest), ``min_jacobian``,
  ``max_jacobian``, ``block_id``; attributes order, dim, node_set, and the merged
  ``points`` / ``connectivity`` (tensor order) for continuous Galerkin numbering.
"""

from pathlib import Path

import h5py
import numpy as np

from opengrid.io.vtkhdf import VERSION, _group, _write_unstructured

VTK_LAGRANGE_QUADRILATERAL = 70
VTK_LAGRANGE_HEXAHEDRON = 72


def _vertex(i: int, j: int) -> int:
    if i:
        if j:
            return 2
        return 1
    if j:
        return 3
    return 0


def vtk_quad_index(i: int, j: int, p: int) -> int:
    """VTK point index of tensor node (i, j) of an order-p Lagrange quadrilateral."""
    ibdy = i in (0, p)
    jbdy = j in (0, p)
    if ibdy and jbdy:
        return _vertex(i, j)
    offset = 4
    if jbdy:  # an edge along i: bottom, then (after the right edge) top
        top = 0
        if j:
            top = (p - 1) + (p - 1)
        return (i - 1) + top + offset
    if ibdy:  # an edge along j: right, then left
        side = 2 * (p - 1) + (p - 1)
        if i:
            side = p - 1
        return (j - 1) + side + offset
    offset += 4 * (p - 1)
    return offset + (i - 1) + (p - 1) * (j - 1)


def vtk_hex_index(i: int, j: int, k: int, p: int) -> int:
    """VTK point index of tensor node (i, j, k) of an order-p Lagrange hexahedron."""
    ibdy = i in (0, p)
    jbdy = j in (0, p)
    kbdy = k in (0, p)
    count = int(ibdy) + int(jbdy) + int(kbdy)
    m = p - 1
    if count == 3:
        upper = 0
        if k:
            upper = 4
        return _vertex(i, j) + upper
    offset = 8
    if count == 2:
        layer = 0
        if k:
            layer = 2 * (m + m)
        if not ibdy:
            row = 0
            if j:
                row = m + m
            return (i - 1) + row + layer + offset
        if not jbdy:
            column = 2 * m + m
            if i:
                column = m
            return (j - 1) + column + layer + offset
        offset += 4 * m + 4 * m
        return (k - 1) + m * _vertex(i, j) + offset  # edges along k: in vertex order
    offset += 4 * (m + m + m)
    if count == 1:
        if ibdy:
            far = 0
            if i:
                far = m * m
            return (j - 1) + m * (k - 1) + far + offset
        offset += 2 * m * m
        if jbdy:
            far = 0
            if j:
                far = m * m
            return (i - 1) + m * (k - 1) + far + offset
        offset += 2 * m * m
        far = 0
        if k:
            far = m * m
        return (i - 1) + m * (j - 1) + far + offset
    offset += 2 * (m * m + m * m + m * m)
    return offset + (i - 1) + m * ((j - 1) + m * (k - 1))


def vtk_permutation(dim: int, p: int) -> np.ndarray:
    """perm[v] = tensor index (i fastest) of VTK point v of a Lagrange cell."""
    n = p + 1
    if dim == 2:
        perm = np.empty(n * n, dtype=np.int64)
        for j in range(n):
            for i in range(n):
                perm[vtk_quad_index(i, j, p)] = j * n + i
        return perm
    perm = np.empty(n * n * n, dtype=np.int64)
    for k in range(n):
        for j in range(n):
            for i in range(n):
                perm[vtk_hex_index(i, j, k, p)] = (k * n + j) * n + i
    return perm


def _block_ids(mesh) -> np.ndarray:
    ids = np.zeros(mesh.elements, dtype=np.int64)
    for index, (_name, first, count, _ne) in enumerate(mesh.owners):
        ids[first : first + count] = index
    return ids


def write_vtkhdf_high_order(path, mesh) -> dict:
    """VTKHDF of Lagrange cells (equispaced node set only)."""
    if mesh.nodes != "equispaced":
        raise ValueError(
            "VTK Lagrange cells interpolate equispaced nodes: export the 'equispaced' node set "
            f"to VTKHDF (or {mesh.nodes!r} to .h5)"
        )
    points, ids = mesh.merged()
    perm = vtk_permutation(mesh.dim, mesh.order)
    cells = ids[:, perm]
    per = cells.shape[1]
    offsets = np.arange(0, per * (len(cells) + 1), per, dtype=np.int64)
    cell_type = VTK_LAGRANGE_QUADRILATERAL
    if mesh.dim == 3:
        cell_type = VTK_LAGRANGE_HEXAHEDRON
    types = np.full(len(cells), cell_type, dtype=np.uint8)
    cell_data = {
        "scaled_jacobian": mesh.scaled_jacobian,
        "min_jacobian": mesh.min_jacobian,
        "block_id": _block_ids(mesh),
    }
    with h5py.File(Path(path), "w") as f:
        root = _group(f, "VTKHDF")
        _write_unstructured(root, points, offsets, cells.reshape(-1), types, cell_data, VERSION)
        names = [name.encode("ascii") for name, *_rest in mesh.owners]
        root["FieldData"].create_dataset(
            "block_names", data=np.array(names, dtype=object), dtype=h5py.string_dtype("ascii")
        )
    return {"points": len(points), "cells": len(cells), "order": mesh.order}


def write_high_order_h5(path, mesh) -> dict:
    """HDF5 with the element nodes of any node set (see the module docstring)."""
    points, ids = mesh.merged()
    with h5py.File(Path(path), "w") as f:
        f.attrs["order"] = mesh.order
        f.attrs["dim"] = mesh.dim
        f.attrs["node_set"] = mesh.nodes
        f.attrs["layout"] = "tensor, i fastest; elements in grid order per block"
        f.create_dataset("nodes", data=mesh.xyz)
        f.create_dataset("min_jacobian", data=mesh.min_jacobian)
        f.create_dataset("max_jacobian", data=mesh.max_jacobian)
        f.create_dataset("block_id", data=_block_ids(mesh))
        f.create_dataset("points", data=points)
        f.create_dataset("connectivity", data=ids)
        names = [name.encode("ascii") for name, *_rest in mesh.owners]
        f.create_dataset(
            "block_names", data=np.array(names, dtype=object), dtype=h5py.string_dtype("ascii")
        )
        f.create_dataset(
            "elements_per_direction", data=np.array([ne for *_rest, ne in mesh.owners])
        )
    return {"points": len(points), "elements": mesh.elements, "order": mesh.order}
