"""CGNS export of structured 2D domains (via the C++ writer on the CGNS library).

One structured zone per domain (zone family = its OpenGrid zone, if any); one Family_t per
patch with its BC type; one BC_t (FamilySpecified, FamilyName = patch) per boundary connector;
GridConnectivity1to1_t in both directions for every interface, including cuts within a domain.
CGNS limits names to 32 characters.
"""

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
    zone_index = {d.id: k for k, d in enumerate(domains)}
    # Zone families: every domain gets one (unzoned domains: the implicit "default" zone, as in
    # the VTKHDF and Gmsh exports), as CGNS recommends.
    zone_names, domain_zone = zone_assignment(session)
    families = [(name, "") for name in zone_names]
    zone_family = {d.id: zone_names[domain_zone[d.name]] for d in domains}
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
    for name in [d.name for d in domains] + [f[0] for f in families]:
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

    boundary = topology.boundary()
    bcs = []
    names_in_zone: dict[int, set] = {}
    for patch in store.entities("patch"):
        if patch.params["patch_type"] == "periodic":
            continue
        for item in patch.params["members"]:
            use = boundary.get(item["connector"])
            if use is None:
                continue
            zone = zone_index[use.domain]
            name = _unique(patch.name, names_in_zone.setdefault(zone, set()))
            imin, jmin, imax, jmax = use.index_range()
            bcs.append((zone, name, patch.name, (imin + 1, jmin + 1, imax + 1, jmax + 1)))

    connections = []
    conn_names: dict[int, set] = {}
    for a, b in topology.interfaces():
        label = store.get(a.connector).name
        for this, donor in ((a, b), (b, a)):
            c = one_to_one(this, donor)
            zone = zone_index[this.domain]
            name = _unique(label, conn_names.setdefault(zone, set()))
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
    for source, _target, pairs in periodic:
        t = source.params["periodic"]["translation"]
        for pair in pairs:
            a = boundary[pair.source]
            b = boundary[pair.target]
            label = store.get(pair.source).name
            for this, donor, shift in ((a, b, (t[0], t[1])), (b, a, (-t[0], -t[1]))):
                c = one_to_one(this, donor, pair.same_direction)
                zone = zone_index[this.domain]
                name = _unique(f"{label}_periodic", conn_names.setdefault(zone, set()))
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
    _core.write_cgns_structured_2d(str(path), base_name, zones, families, bcs, connections)
    return {
        "zones": len(zones),
        "families": len(families),
        "boundaries": len(bcs),
        "connections": len(connections),
    }


def _unique(name: str, used: set) -> str:
    candidate = name[:MAX_NAME]
    number = 1
    while candidate in used:
        number += 1
        suffix = f"_{number}"
        candidate = name[: MAX_NAME - len(suffix)] + suffix
    used.add(candidate)
    return candidate
