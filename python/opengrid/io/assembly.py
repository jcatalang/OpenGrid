"""Unstructured view of a structured multi-domain mesh, shared by the VTKHDF and Gmsh writers.

Points are merged globally (shared connectors and cuts carry bitwise-identical coordinates).
Quads keep their domain and zone; boundary line segments carry their patch and run in loop
order (domain on the left), i.e. boundary normals (dy, -dx) point out of the domain.
Domains that belong to no zone are collected in an implicit zone called "default".
"""

from dataclasses import dataclass, field

import numpy as np

from opengrid.validation import Issue, ValidationError

DEFAULT_ZONE = "default"


@dataclass
class MeshAssembly:
    points: np.ndarray
    quads: np.ndarray  # (m, 4) global point ids, counter-clockwise
    quad_domain: np.ndarray  # (m,) domain index
    quad_zone: np.ndarray  # (m,) zone index
    quad_ij: np.ndarray  # (m, 2) cell indices within the domain
    domain_names: list[str]
    zone_names: list[str]
    lines: np.ndarray  # (k, 2)
    line_patch: np.ndarray  # (k,) patch index
    patch_names: list[str]
    patch_types: list[str]
    cell_data: dict[str, np.ndarray] = field(default_factory=dict)  # per quad
    # Periodic pairs: (source patch index, target patch index, transform dict,
    # target point ids, matching source point ids)
    periodic: list = field(default_factory=list)


def zone_assignment(session) -> tuple[list[str], dict[str, int]]:
    """(zone names, domain name -> zone index); unzoned domains go to the "default" zone."""
    zone_names = []
    domain_zone = {}
    for zone in session.store.entities("zone"):
        zone_names.append(zone.name)
        for item in zone.params["members"]:
            domain_zone[session.store.get(item["domain"]).name] = len(zone_names) - 1
    unzoned = [name for name in session.names("domain") if name not in domain_zone]
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


def assemble(session, topology, quality: bool = True, periodic=()) -> MeshAssembly:
    from opengrid import _core

    domain_names = session.names("domain")
    zone_names, domain_zone = zone_assignment(session)
    all_points = []
    offsets = {}
    quads = []
    quad_domain = []
    quad_zone = []
    quad_ij = []
    data: dict[str, list[np.ndarray]] = {}
    offset = 0
    grids = {}
    for number, name in enumerate(domain_names):
        grid = session.domain_grid(name)
        grids[name] = grid
        nj, ni, _ = grid.shape
        all_points.append(grid.reshape(-1, 3))
        offsets[name] = offset
        index = np.arange(ni * nj).reshape(nj, ni) + offset
        q = np.stack(
            [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
        ).reshape(-1, 4)
        quads.append(q)
        quad_domain.append(np.full(len(q), number, dtype=np.int32))
        quad_zone.append(np.full(len(q), domain_zone[name], dtype=np.int32))
        jj, ii = np.meshgrid(np.arange(nj - 1), np.arange(ni - 1), indexing="ij")
        quad_ij.append(np.stack([ii.reshape(-1), jj.reshape(-1)], axis=1).astype(np.int32))
        if quality:
            for key, values in _core.quality_2d(grid).items():
                data.setdefault(key, []).append(values.reshape(-1))
        offset += ni * nj
    points = np.concatenate(all_points)
    merged, inverse = np.unique(points, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)

    patch_names = []
    patch_types = []
    lines = []
    line_patch = []
    boundary = topology.boundary()
    for patch in session.store.entities("patch"):
        patch_names.append(patch.name)
        patch_types.append(patch.params["patch_type"])
        for item in patch.params["members"]:
            use = boundary.get(item["connector"])
            if use is None:
                continue
            domain = session.store.get(use.domain).name
            nj, ni, _ = grids[domain].shape
            start = use.first
            end = use.last
            if use.reversed:
                start = use.last
                end = use.first
            ids = _node_run(start, end, ni) + offsets[domain]
            segments = np.stack([ids[:-1], ids[1:]], axis=1)
            lines.append(inverse[segments])
            line_patch.append(np.full(len(segments), len(patch_names) - 1, dtype=np.int32))
    if lines:
        line_array = np.concatenate(lines)
        patch_array = np.concatenate(line_patch)
    else:
        line_array = np.zeros((0, 2), dtype=np.int64)
        patch_array = np.zeros(0, dtype=np.int32)
    cell_data = {key: np.concatenate(parts) for key, parts in data.items()}
    periodic_links = []
    for source, target, pairs in periodic:
        source_ids = []
        target_ids = []
        for pair in pairs:
            a = _connector_point_ids(boundary[pair.source], grids, offsets, session, inverse)
            b = _connector_point_ids(boundary[pair.target], grids, offsets, session, inverse)
            if not pair.same_direction:
                b = b[::-1]
            source_ids.append(a)
            target_ids.append(b)
        transform = dict(source.params["periodic"])
        periodic_links.append(
            (
                patch_names.index(source.name),
                patch_names.index(target.name),
                transform,
                np.concatenate(target_ids),
                np.concatenate(source_ids),
            )
        )
    return MeshAssembly(
        points=merged,
        quads=inverse[np.concatenate(quads)],
        quad_domain=np.concatenate(quad_domain),
        quad_zone=np.concatenate(quad_zone),
        quad_ij=np.concatenate(quad_ij),
        domain_names=domain_names,
        zone_names=zone_names,
        lines=line_array,
        line_patch=patch_array,
        patch_names=patch_names,
        patch_types=patch_types,
        cell_data=cell_data,
        periodic=periodic_links,
    )


def _connector_point_ids(use, grids, offsets, session, inverse) -> np.ndarray:
    """Merged point ids of a boundary connector's nodes, in the connector's own direction."""
    domain = session.store.get(use.domain).name
    nj, ni, _ = grids[domain].shape
    return inverse[_node_run(use.first, use.last, ni) + offsets[domain]]


def _node_run(start: tuple[int, int], end: tuple[int, int], ni: int) -> np.ndarray:
    """Flat grid indices (j * ni + i) from start to end along one grid line, inclusive."""
    i0, j0 = start
    i1, j1 = end
    if j0 == j1:
        step = 1
        if i1 < i0:
            step = -1
        i = np.arange(i0, i1 + step, step)
        return j0 * ni + i
    step = 1
    if j1 < j0:
        step = -1
    j = np.arange(j0, j1 + step, step)
    return j * ni + i0
