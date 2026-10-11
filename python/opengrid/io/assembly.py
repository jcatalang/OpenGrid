"""Unstructured view of a multi-domain mesh (structured quads and unstructured triangles),
shared by the VTKHDF and Gmsh writers.

Points are merged globally: shared connectors, cuts and the boundaries of unstructured domains
carry bitwise-identical coordinates (the connector nodes), so they merge exactly. Cells are
stored as offsets + connectivity with a VTK cell type per cell (9 quad, 5 triangle), all
counter-clockwise. Boundary line segments carry their patch and run in loop order (domain on
the left). Domains that belong to no zone are collected in an implicit zone called "default".
"""

from dataclasses import dataclass, field

import numpy as np

from opengrid.validation import Issue, ValidationError

DEFAULT_ZONE = "default"
VTK_TRIANGLE = 5
VTK_QUAD = 9


@dataclass
class MeshAssembly:
    points: np.ndarray
    offsets: np.ndarray  # (m + 1,)
    connectivity: np.ndarray  # flat point ids
    cell_type: np.ndarray  # (m,) VTK type
    cell_domain: np.ndarray  # (m,) domain index
    cell_zone: np.ndarray  # (m,) zone index
    cell_ij: np.ndarray  # (m, 2): (i, j) for quads, (triangle index, -1) for triangles
    domain_names: list[str]
    zone_names: list[str]
    lines: np.ndarray  # (k, 2)
    line_patch: np.ndarray  # (k,) patch index
    patch_names: list[str]
    patch_types: list[str]
    points_before_merge: int = 0
    cell_data: dict[str, np.ndarray] = field(default_factory=dict)  # per cell
    # Periodic pairs: (source patch index, target patch index, transform dict,
    # target point ids, matching source point ids)
    periodic: list = field(default_factory=list)

    def cells(self, mask=None) -> list[np.ndarray]:
        """Point ids of each (selected) cell."""
        indices = np.arange(len(self.cell_type))
        if mask is not None:
            indices = indices[mask]
        return [self.connectivity[self.offsets[i] : self.offsets[i + 1]] for i in indices]


def domain_entities(session) -> list:
    return [e for e in session.store if e.kind in ("domain", "unstructured_domain")]


def zone_assignment(session) -> tuple[list[str], dict[str, int]]:
    """(zone names, domain name -> zone index); unzoned domains go to the "default" zone."""
    zone_names = []
    domain_zone = {}
    for zone in session.store.entities("zone"):
        zone_names.append(zone.name)
        for item in zone.params["members"]:
            domain_zone[session.store.get(item["domain"]).name] = len(zone_names) - 1
    unzoned = [e.name for e in domain_entities(session) if e.name not in domain_zone]
    if unzoned:
        if DEFAULT_ZONE in zone_names or DEFAULT_ZONE in session.store:
            raise ValidationError(
                [
                    Issue(
                        DEFAULT_ZONE, f"name is reserved for unzoned domains ({', '.join(unzoned)})"
                    )
                ]
            )
        zone_names.append(DEFAULT_ZONE)
        for name in unzoned:
            domain_zone[name] = len(zone_names) - 1
    return zone_names, domain_zone


def assemble(
    session, topology=None, quality: bool = True, periodic=(), domains=None
) -> MeshAssembly:
    from opengrid import _core

    entities = domain_entities(session)
    if domains is not None:
        entities = [e for e in entities if e.name in set(domains)]
    zone_names, domain_zone = zone_assignment(session)
    pool = []  # point arrays merged at the end
    offset = 0
    cells = []
    types = []
    cell_domain = []
    cell_zone = []
    cell_ij = []
    data: dict[str, list[np.ndarray]] = {}
    for number, entity in enumerate(entities):
        if entity.kind == "domain":
            grid = session.domain_grid(entity.name)
            nj, ni, _ = grid.shape
            pool.append(grid.reshape(-1, 3))
            index = np.arange(ni * nj).reshape(nj, ni) + offset
            q = np.stack(
                [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
            ).reshape(-1, 4)
            cells.extend(q)
            types.append(np.full(len(q), VTK_QUAD, dtype=np.uint8))
            jj, ii = np.meshgrid(np.arange(nj - 1), np.arange(ni - 1), indexing="ij")
            cell_ij.append(np.stack([ii.reshape(-1), jj.reshape(-1)], axis=1).astype(np.int32))
            metrics = None
            if quality:
                metrics = _core.quality_2d(grid)
            offset += ni * nj
            count = len(q)
        else:
            mesh = session.unstructured_mesh(entity.name)
            pool.append(mesh["points"])
            tris = mesh["triangles"] + offset
            quads = mesh["quads"] + offset
            cells.extend(tris)
            cells.extend(quads)
            types.append(np.full(len(tris), VTK_TRIANGLE, dtype=np.uint8))
            types.append(np.full(len(quads), VTK_QUAD, dtype=np.uint8))
            count = len(tris) + len(quads)
            cell_ij.append(
                np.stack([np.arange(count), np.full(count, -1)], axis=1).astype(np.int32)
            )
            metrics = None
            if quality:
                metrics = _core.quality_cells(
                    np.ascontiguousarray(mesh["points"][:, :2]), mesh["triangles"], mesh["quads"]
                )
            offset += len(mesh["points"])
        cell_domain.append(np.full(count, number, dtype=np.int32))
        cell_zone.append(np.full(count, domain_zone[entity.name], dtype=np.int32))
        if metrics is not None:
            for key, values in metrics.items():
                data.setdefault(key, []).append(values.reshape(-1))
    points_before_merge = offset

    # Connector nodes of patch members and periodic pairs: looked up through the same merge.
    connector_slots = {}
    boundary = {}
    if topology is not None:
        boundary = topology.boundary()
    patches = list(session.store.entities("patch"))
    needed = []
    for patch in patches:
        for item in patch.params["members"]:
            if item["connector"] in boundary:
                needed.append(item["connector"])
    for connector_id in sorted(set(needed)):
        nodes = session.connectors.nodes(connector_id).points
        connector_slots[connector_id] = (offset, len(nodes))
        pool.append(np.asarray(nodes))
        offset += len(nodes)

    merged, inverse = np.unique(np.concatenate(pool), axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    lengths = np.array([len(c) for c in cells], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)])
    connectivity = np.zeros(0, dtype=np.int64)
    if cells:
        connectivity = inverse[np.concatenate(cells)]

    def connector_ids(connector_id, reverse=False):
        start, n = connector_slots[connector_id]
        ids = inverse[start : start + n]
        if reverse:
            ids = ids[::-1]
        return ids

    lines = []
    line_patch = []
    for index, patch in enumerate(patches):
        for item in patch.params["members"]:
            use = boundary.get(item["connector"])
            if use is None:
                continue
            ids = connector_ids(item["connector"], use.reversed)
            lines.append(np.stack([ids[:-1], ids[1:]], axis=1))
            line_patch.append(np.full(len(ids) - 1, index, dtype=np.int32))
    periodic_links = []
    names = [p.name for p in patches]
    for source, target, pairs in periodic:
        source_ids = []
        target_ids = []
        for pair in pairs:
            source_ids.append(connector_ids(pair.source))
            target_ids.append(connector_ids(pair.target, not pair.same_direction))
        periodic_links.append(
            (
                names.index(source.name),
                names.index(target.name),
                dict(source.params["periodic"]),
                np.concatenate(target_ids),
                np.concatenate(source_ids),
            )
        )
    if lines:
        line_array = np.concatenate(lines)
        patch_array = np.concatenate(line_patch)
    else:
        line_array = np.zeros((0, 2), dtype=np.int64)
        patch_array = np.zeros(0, dtype=np.int32)
    if not entities:
        raise ValueError("nothing to export: no domains")
    return MeshAssembly(
        points=merged,
        offsets=offsets,
        connectivity=connectivity,
        cell_type=np.concatenate(types),
        cell_domain=np.concatenate(cell_domain),
        cell_zone=np.concatenate(cell_zone),
        cell_ij=np.concatenate(cell_ij),
        domain_names=[e.name for e in entities],
        zone_names=zone_names,
        lines=line_array,
        line_patch=patch_array,
        patch_names=names,
        patch_types=[p.params["patch_type"] for p in patches],
        points_before_merge=points_before_merge,
        cell_data={key: np.concatenate(parts) for key, parts in data.items()},
        periodic=periodic_links,
    )
