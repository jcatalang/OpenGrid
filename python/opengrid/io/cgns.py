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


def write_cgns_blocks(path, session, topology, base_name: str = "OpenGrid") -> dict:
    """3D structured CGNS of the model's blocks: one zone per block (zone family), one Family_t
    per patch, one BC_t per patched boundary face and GridConnectivity1to1_t in both directions
    for every interface."""
    from opengrid.io.assembly3d import block_zones
    from opengrid.model.blocks import face_range, one_to_one_3d, patch_faces

    if not _core.cgns_available():
        raise RuntimeError("OpenGrid was built without CGNS support")
    store = session.store
    blocks = store.entities("block")
    zone_index = {b.id: k for k, b in enumerate(blocks)}
    zone_names, block_zone = block_zones(session)
    families = [(name, "") for name in zone_names]
    for patch in store.entities("patch"):
        families.append((patch.name, CGNS_BC_TYPES[patch.params["patch_type"]]))
    issues = []
    for name in [b.name for b in blocks] + [f[0] for f in families]:
        if len(name) > MAX_NAME:
            issues.append(Issue(name, f"is longer than {MAX_NAME} characters (CGNS limit)"))
    if issues:
        raise ValidationError(issues)
    zones = []
    for b in blocks:
        zones.append(
            (b.name, np.ascontiguousarray(session.block_grid(b.name)), zone_names[block_zone[b.id]])
        )
    assigned, _issues = patch_faces(store, topology)
    bcs = []
    names_in_zone: dict[int, set] = {}
    for (block_id, face), patch in sorted(assigned.items()):
        zone = zone_index[block_id]
        r = face_range(session.blocks.grid(block_id), face)
        name = _unique(f"{patch}_{face}", names_in_zone.setdefault(zone, set()))
        bcs.append((zone, name, patch, tuple(v + 1 for v in r)))
    connections = []
    for a, b in topology.interfaces:
        for this, donor in ((a, b), (b, a)):
            c = one_to_one_3d(session.blocks, this, donor)
            zone = zone_index[this.block]
            label = f"{store.get(donor.block).name}_{donor.face}"
            name = _unique(label, names_in_zone.setdefault(zone, set()))
            connections.append(
                (
                    zone,
                    name,
                    zone_index[donor.block],
                    tuple(v + 1 for v in c["range"]),
                    tuple(v + 1 for v in c["donor_range"]),
                    c["transform"],
                )
            )
    _core.write_cgns_structured_3d(str(path), base_name, zones, families, bcs, connections)
    return {
        "zones": len(zones),
        "families": len(families),
        "boundaries": len(bcs),
        "connections": len(connections),
    }


def write_cgns_volumes(path, session, base_name: str = "OpenGrid") -> dict:
    """3D unstructured CGNS of the model's volumes: one unstructured zone per volume (TETRA_4,
    PYRA_5 sections; zone family), one Family_t per patch, one boundary section (TRI_3 /
    QUAD_4) + BC_t (FaceCenter) per patched face, vertex connectivity for faces shared by two
    volumes."""
    from opengrid.io.assembly3d import block_zones
    from opengrid.model.volumes import analyze_3d

    if not _core.cgns_available():
        raise RuntimeError("OpenGrid was built without CGNS support")
    store = session.store
    if store.entities("block"):
        raise ValidationError(
            [
                Issue(
                    "CGNS", "models with both structured blocks and volumes are written in Phase 8c"
                )
            ]
        )
    volumes = store.entities("volume")
    zone_index = {v.id: k for k, v in enumerate(volumes)}
    zone_names, block_zone = block_zones(session)
    families = [(name, "") for name in zone_names]
    for patch in store.entities("patch"):
        families.append((patch.name, CGNS_BC_TYPES[patch.params["patch_type"]]))
    issues = []
    for name in [v.name for v in volumes] + [f[0] for f in families]:
        if len(name) > MAX_NAME:
            issues.append(Issue(name, f"is longer than {MAX_NAME} characters (CGNS limit)"))
    if issues:
        raise ValidationError(issues)
    topology = analyze_3d(session)
    meshes = {v.id: session.volume_mesh(v.name) for v in volumes}
    zones = []
    bc_count = 0
    for v in volumes:
        mesh = meshes[v.id]
        sections = [("Tetrahedra", "TETRA_4", mesh["tets"].reshape(-1).tolist(), "")]
        if len(mesh["pyramids"]):
            sections.append(("Pyramids", "PYRA_5", mesh["pyramids"].reshape(-1).tolist(), ""))
        used: set = set()
        for entry in topology["entries"]:
            if entry["volume"] is None or entry["volume"][0] != v.id or entry["patch"] is None:
                continue
            cells = mesh["faces"][entry["volume"][1]]
            for size, kind in ((3, "TRI_3"), (4, "QUAD_4")):
                group = [c for c in cells if len(c) == size]
                if not group:
                    continue
                label = f"{entry['patch']}_{store.get(entry['volume'][1]).name}"
                if size == 4:
                    label += "_q"
                name = _unique(label, used)
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
    # Faces shared by two volumes: matching vertices by coordinates.
    connections = []
    shared: dict[int, list[int]] = {}
    for v in volumes:
        for d in v.params["faces"]:
            shared.setdefault(d, []).append(v.id)
    for d, owners in sorted(shared.items()):
        if len(owners) != 2:
            continue
        a, b = owners
        lookup = {}
        for vid in (a, b):
            pts = meshes[vid]["points"]
            ids = np.unique(np.concatenate(meshes[vid]["faces"][d]))
            lookup[vid] = {pts[i].tobytes(): int(i) for i in ids}
        keys = sorted(lookup[a])
        for this, donor in ((a, b), (b, a)):
            connections.append(
                (
                    zone_index[this],
                    _unique(store.get(d).name, set()),
                    zone_index[donor],
                    [lookup[this][k] + 1 for k in keys],
                    [lookup[donor][k] + 1 for k in keys],
                )
            )
    _core.write_cgns_unstructured_3d(str(path), base_name, zones, families, connections)
    return {
        "zones": len(zones),
        "families": len(families),
        "boundaries": bc_count,
        "connections": len(connections),
    }
