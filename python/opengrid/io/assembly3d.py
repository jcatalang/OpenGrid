"""Unstructured view of a structured 3D multi-block mesh, for the VTKHDF and Gmsh writers.

Cells are hexahedra (VTK_HEXAHEDRON, nodes (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k) and the same
at k+1); points shared by several blocks (interfaces carry identical nodes) are merged.
Boundary faces in patches become quads with outward normals. Blocks that belong to no zone are
collected in the implicit zone "default".
"""

import numpy as np

from opengrid import _core
from opengrid.io.assembly import DEFAULT_ZONE, MeshAssembly
from opengrid.model.blocks import face_array, patch_faces
from opengrid.validation import Issue, ValidationError

VTK_QUAD = 9
VTK_HEXAHEDRON = 12
# Faces whose (first index x second index) direction points into the block.
REVERSED_FACES = ("imin", "kmin", "jmax")


def block_zones(session) -> tuple[list[str], dict[int, int]]:
    zone_names = []
    block_zone = {}
    for zone in session.store.entities("zone"):
        zone_names.append(zone.name)
        for item in zone.params["members"]:
            if session.store.get(item["domain"]).kind == "block":
                block_zone[item["domain"]] = len(zone_names) - 1
    blocks = session.store.entities("block")
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
    store = session.store
    blocks = store.entities("block")
    if not blocks:
        raise ValueError("the model has no blocks")
    zone_names, block_zone = block_zones(session)
    pool = []
    cells = []
    cell_block = []
    cell_zone = []
    cell_ijk = []
    data: dict[str, list] = {}
    offset = 0
    grids = {}
    starts = {}
    for number, block in enumerate(blocks):
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
        cells.append(hexes)
        count = len(hexes)
        cell_block.append(np.full(count, number, dtype=np.int32))
        cell_zone.append(np.full(count, block_zone[block.id], dtype=np.int32))
        kk, jj, ii = np.meshgrid(
            np.arange(nk - 1), np.arange(nj - 1), np.arange(ni - 1), indexing="ij"
        )
        cell_ijk.append(np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1).astype(np.int32))
        if quality:
            for key, values in _core.quality_3d(grid).items():
                data.setdefault(key, []).append(values.reshape(-1))
        offset += ni * nj * nk
    points = np.concatenate(pool)
    merged, inverse = np.unique(points, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    hexes = inverse[np.concatenate(cells)]
    ijk = np.concatenate(cell_ijk)
    cell_data = {key: np.concatenate(v) for key, v in data.items()}
    cell_data["cell_k"] = ijk[:, 2]
    # Boundary quads of patched faces.
    patch_names = [p.name for p in store.entities("patch")]
    patch_types = [p.params["patch_type"] for p in store.entities("patch")]
    quads = np.zeros((0, 4), dtype=np.int64)
    quad_patch = np.zeros(0, dtype=np.int32)
    if topology is not None:
        assigned, _issues = patch_faces(store, topology)
        rows = []
        owners = []
        for (block_id, face), patch in sorted(assigned.items()):
            grid = grids[block_id]
            nk, nj, ni, _ = grid.shape
            index = np.arange(ni * nj * nk).reshape(nk, nj, ni) + starts[block_id]
            ids = face_array(index[..., None], face)[..., 0]
            q = np.stack([ids[:-1, :-1], ids[:-1, 1:], ids[1:, 1:], ids[1:, :-1]], axis=-1)
            q = q.reshape(-1, 4)
            if face in REVERSED_FACES:
                q = q[:, ::-1]
            rows.append(inverse[q])
            owners.append(np.full(len(q), patch_names.index(patch), dtype=np.int32))
        if rows:
            quads = np.concatenate(rows)
            quad_patch = np.concatenate(owners)
    lengths = np.full(len(hexes), 8, dtype=np.int64)
    return MeshAssembly(
        points=merged,
        offsets=np.concatenate([[0], np.cumsum(lengths)]),
        connectivity=hexes.reshape(-1),
        cell_type=np.full(len(hexes), VTK_HEXAHEDRON, dtype=np.uint8),
        cell_domain=np.concatenate(cell_block),
        cell_zone=np.concatenate(cell_zone),
        cell_ij=ijk[:, :2],
        domain_names=[b.name for b in blocks],
        zone_names=zone_names,
        lines=quads,
        line_patch=quad_patch,
        patch_names=patch_names,
        patch_types=patch_types,
        points_before_merge=len(points),
        cell_data=cell_data,
        dimension=3,
    )
