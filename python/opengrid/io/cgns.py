"""CGNS export of 2D meshes (via the C++ writer on the CGNS library).

One structured zone per structured domain and one unstructured zone per unstructured domain
(zone family = its OpenGrid zone, "default" when unzoned); one Family_t per patch with its BC
type; one BC_t (FamilySpecified, FamilyName = patch) per boundary connector (structured: point
range; unstructured: a BAR_2 element section referenced at EdgeCenter). Interfaces are written in
both directions: GridConnectivity1to1_t between structured zones (including cuts within a
domain), general vertex connectivity (point lists) whenever an unstructured zone is involved;
periodic pairs likewise. CGNS limits names to 32 characters.
"""

import itertools

import numpy as np

from opengrid import _core
from opengrid.io.assembly import zone_assignment
from opengrid.model.patches import CGNS_BC_TYPES
from opengrid.model.topology import one_to_one
from opengrid.validation import Issue, ValidationError

MAX_NAME = 32


def write_cgns(path, session, topology, base_name: str = "OpenGrid", periodic=()) -> dict:
    if not _core.cgns_available():
        raise RuntimeError("OpenGrid was built without CGNS support")
    store = session.store
    domains = store.entities("domain")
    unstructured = store.entities("unstructured_domain")
    # Global zone numbering: structured zones first, then unstructured ones.
    zone_index = {d.id: k for k, d in enumerate(domains)}
    for k, d in enumerate(unstructured):
        zone_index[d.id] = len(domains) + k
    structured_ids = {d.id for d in domains}
    # Zone families: every domain gets one (unzoned domains: the implicit "default" zone, as in
    # the VTKHDF and Gmsh exports), as CGNS recommends.
    zone_names, domain_zone = zone_assignment(session)
    families = [(name, "") for name in zone_names]
    zone_family = {d.id: zone_names[domain_zone[d.name]] for d in [*domains, *unstructured]}
    issues = []
    for patch in store.entities("patch"):
        bc_type = CGNS_BC_TYPES[patch.params["patch_type"]]
        if patch.params["patch_type"] == "periodic":
            bc_type = ""  # periodic faces are connectivity, not boundary conditions
        families.append((patch.name, bc_type))
    for source, _target, _pairs in periodic:
        if source.params["periodic"]["rotation_angle"] != 0.0:
            issues.append(
                Issue(
                    source.name,
                    "rotational periodicity cannot be written to a 2D CGNS base "
                    "(translations only)",
                )
            )
    for name in [d.name for d in [*domains, *unstructured]] + [f[0] for f in families]:
        if len(name) > MAX_NAME:
            issues.append(Issue(name, f"is longer than {MAX_NAME} characters (CGNS limit)"))
    if issues:
        raise ValidationError(issues)

    zones = []
    for d in domains:
        grid = session.domain_grid(d.name)
        nj, ni, _ = grid.shape
        zones.append(
            (
                d.name,
                ni,
                nj,
                grid[..., 0].reshape(-1).tolist(),
                grid[..., 1].reshape(-1).tolist(),
                zone_family[d.id],
            )
        )
    unstructured_zones = []
    meshes = {}
    for d in unstructured:
        mesh = session.unstructured_mesh(d.name)
        meshes[d.id] = _MeshLookup(mesh)
        unstructured_zones.append(
            (
                d.name,
                mesh["points"][:, 0].tolist(),
                mesh["points"][:, 1].tolist(),
                np.ascontiguousarray(mesh["triangles"]),
                np.ascontiguousarray(mesh["quads"]),
                zone_family[d.id],
            )
        )

    def points_of(use, connector_id, reverse=False):
        """1-based point list of the connector's nodes in the use's zone, in connector order."""
        if use.domain in structured_ids:
            n = len(session.connectors.nodes(connector_id).points)
            (i0, j0), (i1, j1) = use.first, use.last
            out = []
            for k in range(n):
                t = k
                if reverse:
                    t = n - 1 - k
                i = i0 + (i1 - i0) * t // max(n - 1, 1)
                j = j0 + (j1 - j0) * t // max(n - 1, 1)
                out.extend([i + 1, j + 1])
            return out
        nodes = session.connectors.nodes(connector_id).points
        if reverse:
            nodes = nodes[::-1]
        return [v + 1 for v in meshes[use.domain].vertices(nodes)]

    boundary = topology.boundary()
    bcs = []
    edge_bcs = []
    names_in_zone: dict[int, set] = {}
    for patch in store.entities("patch"):
        if patch.params["patch_type"] == "periodic":
            continue
        for item in patch.params["members"]:
            if "connector" not in item:
                continue
            use = boundary.get(item["connector"])
            if use is None:
                continue
            zone = zone_index[use.domain]
            name = _unique(patch.name, names_in_zone.setdefault(zone, set()))
            if use.domain in structured_ids:
                imin, jmin, imax, jmax = use.index_range()
                bcs.append((zone, name, patch.name, (imin + 1, jmin + 1, imax + 1, jmax + 1)))
            else:
                lookup = meshes[use.domain]
                nodes = session.connectors.nodes(item["connector"]).points
                edge_bcs.append(
                    (zone - len(domains), name, patch.name, lookup.boundary_edges(nodes))
                )

    connections = []
    vertex_connections = []
    conn_names: dict[int, set] = {}
    for a, b in topology.interfaces():
        label = store.get(a.connector).name
        for this, donor in ((a, b), (b, a)):
            zone = zone_index[this.domain]
            name = _unique(label, conn_names.setdefault(zone, set()))
            if this.domain in structured_ids and donor.domain in structured_ids:
                c = one_to_one(this, donor)
                r = c["range"]
                d = c["donor_range"]
                connections.append(
                    (
                        zone,
                        name,
                        zone_index[donor.domain],
                        (r[0] + 1, r[1] + 1, r[2] + 1, r[3] + 1),
                        (d[0] + 1, d[1] + 1, d[2] + 1, d[3] + 1),
                        c["transform"],
                        None,
                    )
                )
            else:
                vertex_connections.append(
                    (
                        zone,
                        name,
                        zone_index[donor.domain],
                        points_of(this, a.connector),
                        points_of(donor, a.connector),
                        None,
                    )
                )
    for source, _target, pairs in periodic:
        t = source.params["periodic"]["translation"]
        for pair in pairs:
            a = boundary[pair.source]
            b = boundary[pair.target]
            label = store.get(pair.source).name
            for this, donor, shift, this_conn, donor_conn in (
                (a, b, (t[0], t[1]), pair.source, pair.target),
                (b, a, (-t[0], -t[1]), pair.target, pair.source),
            ):
                zone = zone_index[this.domain]
                name = _unique(f"{label}_periodic", conn_names.setdefault(zone, set()))
                if this.domain in structured_ids and donor.domain in structured_ids:
                    c = one_to_one(this, donor, pair.same_direction)
                    r = c["range"]
                    d = c["donor_range"]
                    connections.append(
                        (
                            zone,
                            name,
                            zone_index[donor.domain],
                            (r[0] + 1, r[1] + 1, r[2] + 1, r[3] + 1),
                            (d[0] + 1, d[1] + 1, d[2] + 1, d[3] + 1),
                            c["transform"],
                            shift,
                        )
                    )
                else:
                    # Node k of this connector maps to node k of the donor connector
                    # (same_direction) or to its node n - 1 - k.
                    vertex_connections.append(
                        (
                            zone,
                            name,
                            zone_index[donor.domain],
                            points_of(this, this_conn),
                            points_of(donor, donor_conn, not pair.same_direction),
                            shift,
                        )
                    )
    _core.write_cgns_2d(
        str(path),
        base_name,
        zones,
        unstructured_zones,
        families,
        bcs,
        edge_bcs,
        connections,
        vertex_connections,
    )
    return {
        "zones": len(zones) + len(unstructured_zones),
        "families": len(families),
        "boundaries": len(bcs) + len(edge_bcs),
        "connections": len(connections) + len(vertex_connections),
    }


def _pool_refs(pool, zone_index) -> list:
    """(coordinates, (zone, 1-based place)) of every distinct node of the entries in ``pool``
    (entry, points, faces, places): block faces give (i, j, k) triples, volume faces vertex
    numbers."""
    refs = {}
    for entry, points, _faces, local in pool:
        if entry["block"] is not None:
            zone = zone_index[entry["block"][0]]
            places = [tuple(int(v) + 1 for v in row) for row in local]
        else:
            zone = zone_index[entry["volume"][0]]
            places = [(int(v) + 1,) for v in local]
        for point, place in zip(points, places, strict=True):
            refs.setdefault((zone, place), point)
    return [(point, key) for key, point in refs.items()]


def _periodic_connections(
    session, topology, source, target, grids, zone_index, unique_in, connections, vertex_connections
) -> None:
    """The CGNS connections of one 3D periodic pair (see write_cgns_3d)."""
    from opengrid.model.periodic import apply, cgns_property, inverse
    from opengrid.model.periodic_faces import (
        entry_faces,
        match_points,
        one_to_one_periodic,
        tolerance_of,
    )

    forward = source.params["periodic"]
    backward = inverse(forward)
    label = f"{source.name}_periodic"
    pools = {}
    for patch in (source, target):
        pools[patch.name] = [
            (e, *entry_faces(session, e)) for e in topology["entries"] if e["patch"] == patch.name
        ]
    tolerance = tolerance_of(*[points for pool in pools.values() for _e, points, _f, _l in pool])
    s_pool = pools[source.name]
    t_pool = pools[target.name]
    paired_s = set()
    paired_t = set()
    for a, (ea, *_rest) in enumerate(s_pool):
        if ea["block"] is None:
            continue
        for b, (eb, *_rest) in enumerate(t_pool):
            if eb["block"] is None or b in paired_t:
                continue
            ga = grids[ea["block"][0]]
            gb = grids[eb["block"][0]]
            there = one_to_one_periodic(ga, ea["block"][1], gb, eb["block"][1], forward, tolerance)
            if there is None:
                continue
            back = one_to_one_periodic(gb, eb["block"][1], ga, ea["block"][1], backward, tolerance)
            if back is None:
                continue
            for entry, donor, data, transform in (
                (ea, eb, there, forward),
                (eb, ea, back, backward),
            ):
                zone = zone_index[entry["block"][0]]
                connections.append(
                    (
                        zone,
                        unique_in(zone, label),
                        zone_index[donor["block"][0]],
                        tuple(v + 1 for v in data["range"]),
                        tuple(v + 1 for v in data["donor_range"]),
                        data["transform"],
                        cgns_property(transform),
                    )
                )
            paired_s.add(a)
            paired_t.add(b)
            break
    # The rest as point lists: every remaining source node against the remaining target nodes
    # merged by coordinates (a node shared by two target zones is the image of the same source
    # node in both).
    s_refs = _pool_refs([s_pool[k] for k in range(len(s_pool)) if k not in paired_s], zone_index)
    t_refs = _pool_refs([t_pool[k] for k in range(len(t_pool)) if k not in paired_t], zone_index)
    if not s_refs and not t_refs:
        return
    mismatch = ValidationError(
        [Issue(source.name, f"its nodes do not map onto those of {target.name}")]
    )
    if not s_refs or not t_refs:
        raise mismatch
    s_points = np.array([point for point, _key in s_refs])
    t_points = np.array([point for point, _key in t_refs])
    t_unique, t_inverse = np.unique(t_points, axis=0, return_inverse=True)
    t_inverse = t_inverse.reshape(-1)
    found, _distance = match_points(apply(s_points, forward), t_unique, tolerance)
    if np.any(found < 0) or len(np.unique(found)) != len(t_unique):
        raise mismatch
    members: dict[int, list] = {}
    for k, u in enumerate(t_inverse.tolist()):
        members.setdefault(u, []).append(t_refs[k][1])
    groups: dict[tuple, list] = {}
    for k, u in enumerate(found.tolist()):
        zone_s, place_s = s_refs[k][1]
        for zone_t, place_t in members[u]:
            groups.setdefault((zone_s, zone_t), []).append((place_s, place_t))
    for (zs, zt), pairs in sorted(groups.items()):
        here = [p for p, _q in pairs]
        there = [q for _p, q in pairs]
        for this, donor, mine, theirs, transform in (
            (zs, zt, here, there, forward),
            (zt, zs, there, here, backward),
        ):
            vertex_connections.append(
                (
                    this,
                    unique_in(this, label),
                    donor,
                    [n for place in mine for n in place],
                    [n for place in theirs for n in place],
                    cgns_property(transform),
                )
            )


class _MeshLookup:
    """Vertex ids of connector nodes in an unstructured mesh (exact coordinates: boundary points
    are the connector nodes bitwise) and boundary edges oriented as in the cells."""

    def __init__(self, mesh):
        self.points = mesh["points"]
        self.index = {}
        for k, p in enumerate(self.points):
            self.index.setdefault(p.tobytes(), k)
        self.directed = set()
        for cells in (mesh["triangles"], mesh["quads"]):
            for cell in cells:
                n = len(cell)
                for k in range(n):
                    self.directed.add((int(cell[k]), int(cell[(k + 1) % n])))

    def vertices(self, nodes) -> list[int]:
        return [self.index[np.ascontiguousarray(p).tobytes()] for p in nodes]

    def boundary_edges(self, nodes) -> list[tuple[int, int]]:
        ids = self.vertices(nodes)
        edges = []
        for a, b in itertools.pairwise(ids):
            if (a, b) in self.directed:
                edges.append((a, b))
            else:
                edges.append((b, a))
        return edges


def _unique(name: str, used: set) -> str:
    candidate = name[:MAX_NAME]
    number = 1
    while candidate in used:
        number += 1
        suffix = f"_{number}"
        candidate = name[: MAX_NAME - len(suffix)] + suffix
    used.add(candidate)
    return candidate


CGNS_VOLUME_SECTIONS = {
    "tets": ("Tetrahedra", "TETRA_4"),
    "pyramids": ("Pyramids", "PYRA_5"),
    "prisms": ("Prisms", "PENTA_6"),
    "hexes": ("Hexahedra", "HEXA_8"),
}


def write_cgns_3d(path, session, base_name: str = "OpenGrid") -> dict:
    """3D CGNS of the model's blocks and volumes, in one base.

    * Blocks: one structured zone each (zone family), one BC_t (PointRange) per patched boundary
      face and GridConnectivity1to1_t in both directions for every block interface.
    * Volumes: one unstructured zone each (TETRA_4, PYRA_5, PENTA_6, HEXA_8 sections; zone
      family), one boundary section (TRI_3 / QUAD_4) + BC_t (FaceCenter) per patched face.
    * One Family_t per zone and per patch.
    * Faces shared by a volume and a block or another volume: general connectivity
      (Abutting1to1, Vertex) in both directions, matching the identical nodes of the face.
    * Periodic pairs (connectivity, not BCs, as in 2D): two block faces that map corner to
      corner give 1-to-1 periodic connections, all other nodes general connections with the
      periodic property (RotationCenter, RotationAngle = axis times angle, Translation),
      grouped by zone pair; both in both directions with the inverse transform on the way back.
    """
    from opengrid.io.assembly3d import block_zones
    from opengrid.model.blocks import face_range, one_to_one_3d
    from opengrid.model.volumes import analyze_3d, face_name, volume_cells

    if not _core.cgns_available():
        raise RuntimeError("OpenGrid was built without CGNS support")
    store = session.store
    blocks = store.entities("block")
    volumes = store.entities("volume")
    # Global zone numbering: blocks first, then volumes.
    zone_index = {b.id: k for k, b in enumerate(blocks)}
    for k, v in enumerate(volumes):
        zone_index[v.id] = len(blocks) + k
    zone_names, block_zone = block_zones(session)
    families = [(name, "") for name in zone_names]
    periodic_patches = set()
    for patch in store.entities("patch"):
        bc_type = CGNS_BC_TYPES[patch.params["patch_type"]]
        if patch.params["patch_type"] == "periodic":
            bc_type = ""  # periodic faces are connectivity, not boundary conditions
            periodic_patches.add(patch.name)
        families.append((patch.name, bc_type))
    issues = []
    for name in [e.name for e in [*blocks, *volumes]] + [f[0] for f in families]:
        if len(name) > MAX_NAME:
            issues.append(Issue(name, f"is longer than {MAX_NAME} characters (CGNS limit)"))
    if issues:
        raise ValidationError(issues)
    topology = analyze_3d(session)
    names_in_zone: dict[int, set] = {}

    def unique_in(zone, label):
        return _unique(label, names_in_zone.setdefault(zone, set()))

    # Structured zones and their BCs / 1-to-1 interfaces.
    grids = {b.id: np.ascontiguousarray(session.block_grid(b.name)) for b in blocks}
    block_zones_out = [(b.name, grids[b.id], zone_names[block_zone[b.id]]) for b in blocks]
    bcs = []
    patched = [
        e
        for e in topology["entries"]
        if e["block"] is not None and e["patch"] and e["patch"] not in periodic_patches
    ]
    for entry in sorted(patched, key=lambda e: e["block"]):
        block_id, face = entry["block"]
        zone = zone_index[block_id]
        r = face_range(grids[block_id], face)
        name = unique_in(zone, f"{entry['patch']}_{face}")
        bcs.append((zone, name, entry["patch"], tuple(v + 1 for v in r)))
    connections = []
    for a, b in topology["blocks"].interfaces:
        for this, donor in ((a, b), (b, a)):
            c = one_to_one_3d(session.blocks, this, donor)
            zone = zone_index[this.block]
            name = unique_in(zone, f"{store.get(donor.block).name}_{donor.face}")
            connections.append(
                (
                    zone,
                    name,
                    zone_index[donor.block],
                    tuple(v + 1 for v in c["range"]),
                    tuple(v + 1 for v in c["donor_range"]),
                    c["transform"],
                    None,
                )
            )
    # Unstructured zones: cells, then one boundary section per patched face.
    meshes = {v.id: session.volume_mesh(v.name) for v in volumes}
    zones = []
    bc_count = len(bcs)
    for v in volumes:
        mesh = meshes[v.id]
        zone = zone_index[v.id]
        sections = []
        for key, array, _vtk in volume_cells(mesh):
            if len(array) or key == "tets":
                name, kind = CGNS_VOLUME_SECTIONS[key]
                unique_in(zone, name)
                sections.append((name, kind, array.reshape(-1).tolist(), ""))
        for entry in topology["entries"]:
            if entry["volume"] is None or entry["volume"][0] != v.id or entry["patch"] is None:
                continue
            if entry["patch"] in periodic_patches:
                continue
            cells = mesh["faces"][entry["volume"][1]]
            for size, kind in ((3, "TRI_3"), (4, "QUAD_4")):
                group = [c for c in cells if len(c) == size]
                if not group:
                    continue
                label = f"{entry['patch']}_{face_name(store, entry['volume'][1])}"
                if size == 4:
                    label += "_q"
                name = unique_in(zone, label)
                sections.append((name, kind, np.concatenate(group).tolist(), entry["patch"]))
                bc_count += 1
        points = mesh["points"]
        zones.append(
            (
                v.name,
                points[:, 0].tolist(),
                points[:, 1].tolist(),
                points[:, 2].tolist(),
                sections,
                zone_names[block_zone[v.id]],
            )
        )
    # Faces shared with a volume: nodes matched by their (identical) coordinates. Each side
    # gives {coordinate bytes: 1-based point (vertex number or (i, j, k))}.
    sides: dict[int, list[tuple[int, dict]]] = {}
    for v in volumes:
        pts = meshes[v.id]["points"]
        for d in v.params["faces"] or ():
            ids = np.unique(np.concatenate(meshes[v.id]["faces"][d]))
            sides.setdefault(d, []).append((v.id, {pts[i].tobytes(): (int(i) + 1,) for i in ids}))
    for f in topology["blocks"].boundary:
        if f.domain is None or f.domain not in sides:
            continue
        grid = grids[f.block]
        r = face_range(grid, f.face)
        lookup = {}
        for k in range(r[2], r[5] + 1):
            for j in range(r[1], r[4] + 1):
                for i in range(r[0], r[3] + 1):
                    lookup[grid[k, j, i].tobytes()] = (i + 1, j + 1, k + 1)
        sides[f.domain].append((f.block, lookup))
    vertex_connections = []
    for d, owners in sorted(sides.items()):
        if len(owners) != 2:
            continue
        keys = sorted(owners[0][1])
        if sorted(owners[1][1]) != keys:
            raise ValidationError(
                [Issue(store.get(d).name, "the nodes of its two sides do not match")]
            )
        for (this, here), (donor, there) in (owners, owners[::-1]):
            zone = zone_index[this]
            vertex_connections.append(
                (
                    zone,
                    unique_in(zone, store.get(d).name),
                    zone_index[donor],
                    [n for key in keys for n in here[key]],
                    [n for key in keys for n in there[key]],
                    None,
                )
            )
    for source, target, _match in topology["periodic"]:
        _periodic_connections(
            session,
            topology,
            source,
            target,
            grids,
            zone_index,
            unique_in,
            connections,
            vertex_connections,
        )
    _core.write_cgns_3d(
        str(path), base_name, block_zones_out, zones, families, bcs, connections, vertex_connections
    )
    return {
        "zones": len(block_zones_out) + len(zones),
        "families": len(families),
        "boundaries": bc_count,
        "connections": len(connections) + len(vertex_connections),
    }
