"""VTKHDF export.

Layout checked against the files written by the installed VTK 9.7 (vtkHDFWriter, format 2.8): a
root group ``VTKHDF`` with attributes ``Version`` and ``Type`` (fixed-length ASCII),
per-partition ``NumberOfPoints``/``NumberOfCells``/``NumberOfConnectivityIds``, ``Points``
(n, 3), ``Offsets``, ``Connectivity``, ``Types`` (uint8) and the groups ``CellData``,
``PointData`` and ``FieldData``; string field data as variable-length ASCII. All groups track
link creation order (VTK's reader iterates links in creation order and does not return
otherwise).

Two forms:
  * a single UnstructuredGrid (format 1.0) of all domains (no patches or zones in the model);
  * a PartitionedDataSetCollection (format 2.1) with one cell block per zone and one line block
    per patch, grouped in the Assembly as Volume/<zone> and Patches/<patch>.
Cells are quads (structured domains) and triangles (unstructured domains), or hexahedra with
boundary quads for 3D block models; points shared by several domains / blocks are merged so the
mesh is conforming.
"""

from pathlib import Path

import h5py
import numpy as np

VTK_LINE = 3
VTK_QUAD = 9
VERSION = (1, 0)
COLLECTION_VERSION = (2, 1)  # PartitionedDataSetCollection with an Assembly needs format 2.1
RESERVED_NAMES = {"Assembly"}


def _group(parent, name):
    return parent.create_group(name, track_order=True)


def _write_unstructured(group, points, offsets, connectivity, types, cell_data, version) -> None:
    """One partition; ``connectivity`` refers to ``points`` (unused points are dropped)."""
    used, local = np.unique(connectivity, return_inverse=True)
    n_cells = len(types)
    group.attrs["Version"] = np.array(version, dtype=np.int64)
    group.attrs.create("Type", np.bytes_("UnstructuredGrid"), dtype=h5py.string_dtype("ascii", 16))
    group.create_dataset("NumberOfPoints", data=np.array([len(used)], dtype=np.int64))
    group.create_dataset("NumberOfCells", data=np.array([n_cells], dtype=np.int64))
    group.create_dataset(
        "NumberOfConnectivityIds", data=np.array([len(connectivity)], dtype=np.int64)
    )
    group.create_dataset("Points", data=points[used])
    group.create_dataset("Offsets", data=np.asarray(offsets, dtype=np.int64))
    group.create_dataset("Connectivity", data=local.reshape(-1).astype(np.int64))
    group.create_dataset("Types", data=np.asarray(types, dtype=np.uint8))
    cells_group = _group(group, "CellData")
    for name, values in cell_data.items():
        cells_group.create_dataset(name, data=values)
    _group(group, "PointData")
    _group(group, "FieldData")


def _subset(assembly, mask):
    """Offsets / connectivity / types of the selected cells."""
    cells = assembly.cells(mask)
    lengths = np.array([len(c) for c in cells], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)])
    connectivity = np.concatenate(cells)
    return offsets, connectivity, assembly.cell_type[mask]


def _cell_data(assembly, mask) -> dict:
    data = {
        "zone_id": assembly.cell_zone[mask],
        "block_id": assembly.cell_domain[mask],
        "cell_i": assembly.cell_ij[mask, 0],
        "cell_j": assembly.cell_ij[mask, 1],
    }
    for key, values in assembly.cell_data.items():
        data[key] = values[mask]
    return data


def write_vtkhdf_grid(path, assembly) -> dict:
    """Single UnstructuredGrid of all cells, with block (domain) ids, zone ids, cell indices and
    quality metrics as cell data and the domain names as FieldData ``block_names``."""
    mask = np.ones(len(assembly.cell_type), dtype=bool)
    offsets, connectivity, types = _subset(assembly, mask)
    used = np.unique(connectivity)
    degenerate = 0
    for cell in assembly.cells():
        if len(set(cell.tolist())) < len(cell):
            degenerate += 1
    with h5py.File(Path(path), "w") as f:
        root = _group(f, "VTKHDF")
        _write_unstructured(
            root, assembly.points, offsets, connectivity, types, _cell_data(assembly, mask), VERSION
        )
        names = [name.encode("ascii") for name in assembly.domain_names]
        root["FieldData"].create_dataset(
            "block_names", data=np.array(names, dtype=object), dtype=h5py.string_dtype("ascii")
        )
    return {
        "points_before_merge": assembly.points_before_merge,
        "points": len(used),
        "cells": len(types),
        "degenerate_cells": degenerate,
    }


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
    """PartitionedDataSetCollection: one cell block per zone (Assembly/Volume) and one line
    block per patch (Assembly/Patches)."""
    names = list(assembly.zone_names) + list(assembly.patch_names)
    clash = sorted(RESERVED_NAMES.intersection(names))
    if clash:
        raise ValueError(f"reserved VTKHDF block name(s): {', '.join(clash)}")
    index = 0
    summary = {"points": len(assembly.points), "volume_cells": 0, "boundary_cells": 0, "blocks": []}
    with h5py.File(Path(path), "w") as f:
        root = _group(f, "VTKHDF")
        root.attrs["Version"] = np.array(COLLECTION_VERSION, dtype=np.int64)
        root.attrs.create(
            "Type", np.bytes_("PartitionedDataSetCollection"), dtype=h5py.string_dtype("ascii", 28)
        )
        tree = _group(root, "Assembly")
        volume = _group(tree, "Volume")
        patches = _group(tree, "Patches")
        for zone_index, zone in enumerate(assembly.zone_names):
            mask = assembly.cell_zone == zone_index
            if not np.any(mask):
                continue
            group = _group(root, zone)
            group.attrs["Index"] = np.int64(index)
            offsets, connectivity, types = _subset(assembly, mask)
            _write_unstructured(
                group,
                assembly.points,
                offsets,
                connectivity,
                types,
                _cell_data(assembly, mask),
                COLLECTION_VERSION,
            )
            volume[zone] = h5py.SoftLink(f"/VTKHDF/{zone}")
            summary["volume_cells"] += int(np.count_nonzero(mask))
            summary["blocks"].append(zone)
            index += 1
        for patch_index, patch in enumerate(assembly.patch_names):
            mask = assembly.line_patch == patch_index
            if not np.any(mask):
                continue
            group = _group(root, patch)
            group.attrs["Index"] = np.int64(index)
            count = int(np.count_nonzero(mask))
            width = assembly.lines.shape[1]
            boundary_type = VTK_LINE
            if width == 4:
                boundary_type = VTK_QUAD
            _write_unstructured(
                group,
                assembly.points,
                np.arange(0, width * count + 1, width),
                assembly.lines[mask].reshape(-1),
                np.full(count, boundary_type, dtype=np.uint8),
                {"patch_id": assembly.line_patch[mask]},
                COLLECTION_VERSION,
            )
            _write_periodic_fields(group, assembly, patch_index)
            patches[patch] = h5py.SoftLink(f"/VTKHDF/{patch}")
            summary["boundary_cells"] += count
            summary["blocks"].append(patch)
            index += 1
    return summary
