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
        root = f.create_group("VTKHDF", track_order=True)
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
        group = root.create_group("CellData", track_order=True)
        for name, parts in cell_data.items():
            group.create_dataset(name, data=np.concatenate(parts))
        root.create_group("PointData", track_order=True)
        fields = root.create_group("FieldData", track_order=True)
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


VTK_LINE = 3
COLLECTION_VERSION = (2, 1)  # PartitionedDataSetCollection with an Assembly needs format 2.1
RESERVED_NAMES = {"Assembly"}


def _write_unstructured(
    group, points: np.ndarray, cells: np.ndarray, cell_type: int, cell_data: dict
) -> None:
    used, local = np.unique(cells.reshape(-1), return_inverse=True)
    connectivity = local.reshape(cells.shape)
    n_cells, per_cell = cells.shape
    group.attrs["Version"] = np.array(COLLECTION_VERSION, dtype=np.int64)
    group.attrs.create("Type", np.bytes_("UnstructuredGrid"), dtype=h5py.string_dtype("ascii", 16))
    group.create_dataset("NumberOfPoints", data=np.array([len(used)], dtype=np.int64))
    group.create_dataset("NumberOfCells", data=np.array([n_cells], dtype=np.int64))
    group.create_dataset(
        "NumberOfConnectivityIds", data=np.array([per_cell * n_cells], dtype=np.int64)
    )
    group.create_dataset("Points", data=points[used])
    group.create_dataset(
        "Offsets", data=np.arange(0, per_cell * n_cells + 1, per_cell, dtype=np.int64)
    )
    group.create_dataset("Connectivity", data=connectivity.reshape(-1).astype(np.int64))
    group.create_dataset("Types", data=np.full(n_cells, cell_type, dtype=np.uint8))
    cells_group = group.create_group("CellData", track_order=True)
    for name, values in cell_data.items():
        cells_group.create_dataset(name, data=values)
    group.create_group("PointData", track_order=True)
    group.create_group("FieldData", track_order=True)


def _write_periodic_fields(group, assembly, patch_index: int) -> None:
    """Periodic pairing as field data: on the source block the partner name and the 4 x 4
    affine transform (source -> partner); on the target block the source name."""
    from opengrid.model.periodic import affine_matrix

    fields = group["FieldData"]
    ascii = h5py.string_dtype("ascii")
    for source, target, transform, _target_ids, _source_ids in assembly.periodic:
        if source == patch_index:
            partner = assembly.patch_names[target].encode("ascii")
            fields.create_dataset(
                "periodic_partner", data=np.array([partner], dtype=object), dtype=ascii
            )
            fields.create_dataset("periodic_affine", data=affine_matrix(transform).reshape(-1))
        if target == patch_index:
            origin = assembly.patch_names[source].encode("ascii")
            fields.create_dataset(
                "periodic_source", data=np.array([origin], dtype=object), dtype=ascii
            )


def write_vtkhdf_collection(path, assembly) -> dict:
    """PartitionedDataSetCollection: one quad block per zone (Assembly/Volume) and one line
    block per patch (Assembly/Patches), with zone_id / block_id / quality and patch_id cell
    data. Layout as written by VTK 9.7's vtkHDFWriter (see module docstring)."""
    names = list(assembly.zone_names) + list(assembly.patch_names)
    clash = sorted(RESERVED_NAMES.intersection(names))
    if clash:
        raise ValueError(f"reserved VTKHDF block name(s): {', '.join(clash)}")
    index = 0
    summary = {"points": len(assembly.points), "volume_cells": 0, "boundary_cells": 0, "blocks": []}
    # VTK's reader iterates links by creation order: groups must track and index it (VTK's own
    # writer does; without it vtkHDFReader does not return).
    with h5py.File(Path(path), "w") as f:
        root = f.create_group("VTKHDF", track_order=True)
        root.attrs["Version"] = np.array(COLLECTION_VERSION, dtype=np.int64)
        root.attrs.create(
            "Type", np.bytes_("PartitionedDataSetCollection"), dtype=h5py.string_dtype("ascii", 28)
        )
        tree = root.create_group("Assembly", track_order=True)
        volume = tree.create_group("Volume", track_order=True)
        patches = tree.create_group("Patches", track_order=True)
        for zone_index, zone in enumerate(assembly.zone_names):
            mask = assembly.quad_zone == zone_index
            if not np.any(mask):
                continue
            data = {
                "zone_id": assembly.quad_zone[mask],
                "block_id": assembly.quad_domain[mask],
                "cell_i": assembly.quad_ij[mask, 0],
                "cell_j": assembly.quad_ij[mask, 1],
            }
            for key, values in assembly.cell_data.items():
                data[key] = values[mask]
            group = root.create_group(zone, track_order=True)
            group.attrs["Index"] = np.int64(index)
            _write_unstructured(group, assembly.points, assembly.quads[mask], VTK_QUAD, data)
            volume[zone] = h5py.SoftLink(f"/VTKHDF/{zone}")
            summary["volume_cells"] += int(np.count_nonzero(mask))
            summary["blocks"].append(zone)
            index += 1
        for patch_index, patch in enumerate(assembly.patch_names):
            mask = assembly.line_patch == patch_index
            if not np.any(mask):
                continue
            group = root.create_group(patch, track_order=True)
            group.attrs["Index"] = np.int64(index)
            _write_unstructured(
                group,
                assembly.points,
                assembly.lines[mask],
                VTK_LINE,
                {"patch_id": assembly.line_patch[mask]},
            )
            _write_periodic_fields(group, assembly, patch_index)
            patches[patch] = h5py.SoftLink(f"/VTKHDF/{patch}")
            summary["boundary_cells"] += int(np.count_nonzero(mask))
            summary["blocks"].append(patch)
            index += 1
    return summary
