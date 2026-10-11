"""VTKHDF export (UnstructuredGrid).

Layout checked against the VTKHDF files written by the installed VTK 9.7 (vtkHDFWriter, which
writes format version 2.8): a root group ``VTKHDF`` with attributes ``Version`` and ``Type``
(fixed-length ASCII), per-partition ``NumberOfPoints``/``NumberOfCells``/
``NumberOfConnectivityIds``, ``Points`` (n, 3), ``Offsets``, ``Connectivity``, ``Types``
(uint8) and the groups ``CellData``, ``PointData`` and ``FieldData``. Only features defined
since format version 1.0 are used, so the file is declared as version 1.0 for the widest
reader compatibility.

Structured blocks are written as quads in a single partition; points shared by several blocks
(identical coordinates along shared connectors and cuts) are merged so the mesh is conforming.
"""

from pathlib import Path

import h5py
import numpy as np

VTK_QUAD = 9
VERSION = (1, 0)


def write_vtkhdf(path, blocks: list[dict]) -> dict:
    """Write blocks ``{"name", "grid" (nj, ni, 3), "cell_data": {name: (nj-1, ni-1)}}``.

    Cell data present in every block is written; ``block_id`` (int32), ``cell_i`` and
    ``cell_j`` (int32) are added, and the block names go to FieldData ``block_names``.
    Returns a summary: number of points before/after merging, number of cells and degenerate
    quads (repeated corner points).
    """
    if not blocks:
        raise ValueError("nothing to export")
    all_points = []
    cells = []
    offset = 0
    for block in blocks:
        grid = np.asarray(block["grid"], dtype=float)
        nj, ni, _ = grid.shape
        all_points.append(grid.reshape(-1, 3))
        index = np.arange(ni * nj).reshape(nj, ni) + offset
        quads = np.stack([index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1)
        cells.append(quads.reshape(-1, 4))
        offset += ni * nj
    points = np.concatenate(all_points)
    merged, inverse = np.unique(points, axis=0, return_inverse=True)
    connectivity = inverse.reshape(-1)[np.concatenate(cells)]
    n_cells = len(connectivity)
    degenerate = int(np.sum([len(set(row)) < 4 for row in connectivity.tolist()]))

    cell_data: dict[str, list[np.ndarray]] = {"block_id": [], "cell_i": [], "cell_j": []}
    common = set(blocks[0].get("cell_data", {}))
    for block in blocks[1:]:
        common &= set(block.get("cell_data", {}))
    for name in sorted(common):
        cell_data[name] = []
    for number, block in enumerate(blocks):
        nj, ni, _ = np.asarray(block["grid"]).shape
        jj, ii = np.meshgrid(np.arange(nj - 1), np.arange(ni - 1), indexing="ij")
        cell_data["block_id"].append(np.full((nj - 1) * (ni - 1), number, dtype=np.int32))
        cell_data["cell_i"].append(ii.reshape(-1).astype(np.int32))
        cell_data["cell_j"].append(jj.reshape(-1).astype(np.int32))
        for name in sorted(common):
            cell_data[name].append(np.asarray(block["cell_data"][name], dtype=float).reshape(-1))

    with h5py.File(Path(path), "w") as f:
        root = f.create_group("VTKHDF")
        root.attrs["Version"] = np.array(VERSION, dtype=np.int64)
        root.attrs.create(
            "Type", np.bytes_("UnstructuredGrid"), dtype=h5py.string_dtype("ascii", 16)
        )
        root.create_dataset("NumberOfPoints", data=np.array([len(merged)], dtype=np.int64))
        root.create_dataset("NumberOfCells", data=np.array([n_cells], dtype=np.int64))
        root.create_dataset("NumberOfConnectivityIds", data=np.array([4 * n_cells], dtype=np.int64))
        root.create_dataset("Points", data=merged)
        root.create_dataset("Offsets", data=np.arange(0, 4 * n_cells + 1, 4, dtype=np.int64))
        root.create_dataset("Connectivity", data=connectivity.reshape(-1).astype(np.int64))
        root.create_dataset("Types", data=np.full(n_cells, VTK_QUAD, dtype=np.uint8))
        group = root.create_group("CellData")
        for name, parts in cell_data.items():
            group.create_dataset(name, data=np.concatenate(parts))
        root.create_group("PointData")
        fields = root.create_group("FieldData")
        # Variable-length ASCII strings, as VTK writes them (VTK cannot read UTF-8 here).
        names = [block["name"].encode("ascii") for block in blocks]
        fields.create_dataset(
            "block_names", data=np.array(names, dtype=object), dtype=h5py.string_dtype("ascii")
        )
    return {
        "points_before_merge": len(points),
        "points": len(merged),
        "cells": n_cells,
        "degenerate_cells": degenerate,
    }
