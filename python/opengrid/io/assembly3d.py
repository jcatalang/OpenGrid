"""Unstructured view of a 3D mesh (structured blocks and unstructured volumes), for the VTKHDF
and Gmsh writers.

Cells: hexahedra of blocks (VTK_HEXAHEDRON, nodes (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k) and
the same at k+1), and the tetrahedra, pyramids, prisms and hexahedra of volumes (VTK types and
node orders as in opengrid.model.volumes.CELL_TYPES);
points shared by several blocks / volumes (interfaces carry identical nodes) are merged.
Boundary faces in patches become quads (block faces) or triangles / quads (volume faces) with
outward normals. Blocks and volumes in no zone are collected in the implicit zone "default".
"""

import numpy as np

from opengrid import _core
from opengrid.io.assembly import DEFAULT_ZONE, MeshAssembly
from opengrid.model.blocks import face_array
from opengrid.model.volumes import volume_cells
from opengrid.validation import Issue, ValidationError

VTK_QUAD = 9
VTK_HEXAHEDRON = 12
VTK_TRIANGLE = 5
# Faces whose (first index x second index) direction points into the block.
REVERSED_FACES = ("imin", "kmin", "jmax")


def block_zones(session) -> tuple[list[str], dict[int, int]]:
    zone_names = []
    block_zone = {}
    for zone in session.store.entities("zone"):
        zone_names.append(zone.name)
        for item in zone.params["members"]:
            if session.store.get(item["domain"]).kind in ("block", "volume"):
                block_zone[item["domain"]] = len(zone_names) - 1
    blocks = session.store.entities("block") + session.store.entities("volume")
    if any(b.id not in block_zone for b in blocks):
        if DEFAULT_ZONE in zone_names:
            raise ValidationError(
                [Issue(DEFAULT_ZONE, "is reserved for the blocks outside every zone")]
            )
        zone_names.append(DEFAULT_ZONE)
        for b in blocks:
            block_zone.setdefault(b.id, len(zone_names) - 1)
    return zone_names, block_zone


def assemble_blocks(session, topology=None, quality: bool = True) -> MeshAssembly:
    """``topology``: the result of opengrid.model.volumes.analyze_3d (patched boundary faces),
    or None for cells only."""
    store = session.store
    blocks = store.entities("block")
    volumes = store.entities("volume")
    if not blocks and not volumes:
        raise ValueError("the model has no blocks or volumes")
    zone_names, block_zone = block_zones(session)
    pool = []
    cells = []
    types = []
    cell_block = []
    cell_zone = []
    cell_ijk = []
    data: dict[str, list] = {}
    offset = 0
    grids = {}
    starts = {}
    number = 0
    for block in blocks:
        grid = session.block_grid(block.name)
        grids[block.id] = grid
        starts[block.id] = offset
        nk, nj, ni, _ = grid.shape
        pool.append(grid.reshape(-1, 3))
        index = np.arange(ni * nj * nk).reshape(nk, nj, ni) + offset
        hexes = np.stack(
            [
                index[:-1, :-1, :-1],
                index[:-1, :-1, 1:],
                index[:-1, 1:, 1:],
                index[:-1, 1:, :-1],
                index[1:, :-1, :-1],
                index[1:, :-1, 1:],
                index[1:, 1:, 1:],
                index[1:, 1:, :-1],
            ],
            axis=-1,
        ).reshape(-1, 8)
        cells.extend(hexes)
        types.append(np.full(len(hexes), VTK_HEXAHEDRON, dtype=np.uint8))
        cell_block.append(np.full(len(hexes), number, dtype=np.int32))
        cell_zone.append(np.full(len(hexes), block_zone[block.id], dtype=np.int32))
        kk, jj, ii = np.meshgrid(
            np.arange(nk - 1), np.arange(nj - 1), np.arange(ni - 1), indexing="ij"
        )
        cell_ijk.append(np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1).astype(np.int32))
        if quality:
            for key, values in _core.quality_3d(grid).items():
                data.setdefault(key, []).append(values.reshape(-1))
        offset += ni * nj * nk
        number += 1
    meshes = {}
    for volume in volumes:
        mesh = session.volume_mesh(volume.name)
        meshes[volume.id] = mesh
        starts[volume.id] = offset
        pool.append(mesh["points"])
        count = 0
        for _key, array, vtk_type in volume_cells(mesh):
            cells.extend(array + offset)
            types.append(np.full(len(array), vtk_type, dtype=np.uint8))
            count += len(array)
        cell_block.append(np.full(count, number, dtype=np.int32))
        cell_zone.append(np.full(count, block_zone[volume.id], dtype=np.int32))
        cell_ijk.append(
            np.stack([np.arange(count), np.full(count, -1), np.full(count, -1)], axis=1).astype(
                np.int32
            )
        )
        if quality:
            for key, values in session.volume_quality(volume.name).items():
                data.setdefault(key, []).append(values.reshape(-1))
        offset += len(mesh["points"])
        number += 1
    points = np.concatenate(pool)
    merged, inverse = np.unique(points, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    connectivity = [inverse[np.asarray(c)] for c in cells]
    ijk = np.concatenate(cell_ijk)
    cell_data = {}
    for key, parts in data.items():
        if len(parts) == number:  # metrics every block / volume provides
            cell_data[key] = np.concatenate(parts)
    cell_data["cell_k"] = ijk[:, 2]
    patch_names = [p.name for p in store.entities("patch")]
    patch_types = [p.params["patch_type"] for p in store.entities("patch")]
    boundary_cells = []
    boundary_patch = []
    boundary_sources = []  # every boundary face source, patched or not
    if topology is not None:
        for entry in topology["entries"]:
            patch = None
            if entry["patch"] is not None:
                patch = patch_names.index(entry["patch"])
            domains = [key[1] for key in entry["keys"] if key[0] == "domain"]
            domain = None
            if domains:
                domain = domains[0]
            if entry["block"] is not None:
                block_id, face = entry["block"]
                grid = grids[block_id]
                nk, nj, ni, _ = grid.shape
                index = np.arange(ni * nj * nk).reshape(nk, nj, ni) + starts[block_id]
                ids = face_array(index[..., None], face)[..., 0]
                q = np.stack([ids[:-1, :-1], ids[:-1, 1:], ids[1:, 1:], ids[1:, :-1]], axis=-1)
                n2, n1 = ids.shape
                jj, ii = np.meshgrid(np.arange(n2 - 1), np.arange(n1 - 1), indexing="ij")
                cell_ij = np.stack([ii.reshape(-1), jj.reshape(-1)], axis=1)
                q = q.reshape(-1, 4)
                boundary_sources.append(
                    {
                        "block": (block_id, face),
                        "volume": None,
                        "domain": domain,
                        "cells": inverse[q],  # face-grid order: (i, j), (i+1, j), ...
                        "grid": face_array(grid, face),
                        "cell_ij": cell_ij,
                    }
                )
                if face in REVERSED_FACES:
                    q = q[:, ::-1]
                if patch is not None:
                    for row in q:
                        boundary_cells.append(inverse[row])
                        boundary_patch.append(patch)
            else:
                volume_id, domain = entry["volume"]
                cells = [
                    inverse[np.asarray(cell) + starts[volume_id]]
                    for cell in meshes[volume_id]["faces"][domain]
                ]
                boundary_sources.append(
                    {
                        "block": None,
                        "volume": (volume_id, domain),
                        "domain": domain,
                        "cells": cells,
                        "grid": None,
                        "cell_ij": None,
                    }
                )
                if patch is not None:
                    for cell in cells:
                        boundary_cells.append(cell)
                        boundary_patch.append(patch)
    lengths = np.array([len(c) for c in connectivity], dtype=np.int64)
    return MeshAssembly(
        points=merged,
        offsets=np.concatenate([[0], np.cumsum(lengths)]),
        connectivity=np.concatenate(connectivity),
        cell_type=np.concatenate(types),
        cell_domain=np.concatenate(cell_block),
        cell_zone=np.concatenate(cell_zone),
        cell_ij=ijk[:, :2],
        domain_names=[b.name for b in blocks] + [v.name for v in volumes],
        zone_names=zone_names,
        lines=np.zeros((0, 4), dtype=np.int64),
        line_patch=np.asarray(boundary_patch, dtype=np.int32),
        patch_names=patch_names,
        patch_types=patch_types,
        points_before_merge=len(points),
        cell_data=cell_data,
        dimension=3,
        boundary_cells=boundary_cells,
        boundary_patch=np.asarray(boundary_patch, dtype=np.int32),
        boundary_sources=boundary_sources,
    )
