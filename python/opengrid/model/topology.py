"""Boundary topology of structured 2D domains: which connectors lie on the boundary, which are
interfaces between (or within) domains, 1-to-1 connectivity and the patch validator.

Every use of a connector in a domain edge is an :class:`EdgeUse` with its node index range in
that domain. A connector used once is a boundary connector (it must belong to exactly one
patch before export); used twice it is an interface (two domains, or a cut within one domain,
detected automatically, no patch needed); more than twice is non-manifold (an error).
Index ranges are 0-based here and converted to 1-based for CGNS by the exporters.
"""

from dataclasses import dataclass, field

from opengrid.validation import Issue

# Face of each edge: (axis normal to the face, True if it is the max face, running axis).
# Axis 0 = i, 1 = j.
EDGE_FACES = {0: (1, False, 0), 1: (0, True, 1), 2: (1, True, 0), 3: (0, False, 1)}


@dataclass(frozen=True)
class EdgeUse:
    domain: int  # entity id
    edge: int
    item: int
    connector: int  # entity id
    reversed: bool
    first: tuple[int, int]  # (i, j) of the connector's first node in this domain
    last: tuple[int, int]  # (i, j) of the connector's last node

    def face(self) -> tuple[int, bool, int]:
        return EDGE_FACES[self.edge]

    def index_range(self) -> tuple[int, int, int, int]:
        """(imin, jmin, imax, jmax), 0-based, sorted."""
        i0, j0 = self.first
        i1, j1 = self.last
        return (min(i0, i1), min(j0, j1), max(i0, i1), max(j0, j1))


@dataclass
class Topology:
    uses: dict[int, list[EdgeUse]] = field(default_factory=dict)  # connector id -> uses
    dimensions: dict[int, tuple[int, int]] = field(default_factory=dict)  # domain id -> (ni, nj)

    def boundary(self) -> dict[int, EdgeUse]:
        return {c: u[0] for c, u in self.uses.items() if len(u) == 1}

    def interfaces(self) -> list[tuple[EdgeUse, EdgeUse]]:
        return [(u[0], u[1]) for c, u in sorted(self.uses.items()) if len(u) == 2]

    def non_manifold(self) -> dict[int, list[EdgeUse]]:
        return {c: u for c, u in self.uses.items() if len(u) > 2}


def _position(edge: int, p: int, ni: int, nj: int) -> tuple[int, int]:
    """(i, j) of the p-th node along an edge in loop order."""
    if edge == 0:
        return (p, 0)
    if edge == 1:
        return (ni - 1, p)
    if edge == 2:
        return (ni - 1 - p, nj - 1)
    return (0, nj - 1 - p)


def analyze(store, session) -> Topology:
    topo = Topology()
    for entity in store.entities("connector"):
        topo.uses[entity.id] = []
    # All domains in id order (both kinds), so uses and interface labels follow creation order.
    for domain in store:
        if domain.kind == "unstructured_domain":
            for loop in domain.params["loops"]:
                for m, item in enumerate(loop):
                    # No structured index range: edge -1 marks a loop use.
                    topo.uses[item["connector"]].append(
                        EdgeUse(
                            domain.id,
                            -1,
                            m,
                            item["connector"],
                            item["reversed"],
                            (-1, -1),
                            (-1, -1),
                        )
                    )
            continue
        if domain.kind != "domain":
            continue
        ni, nj = session.domain_dimensions(domain.name)
        topo.dimensions[domain.id] = (ni, nj)
        for k, edge in enumerate(domain.params["edges"]):
            p = 0
            for m, item in enumerate(edge):
                n = store.get(item["connector"]).params["dimension"]
                start = _position(k, p, ni, nj)
                end = _position(k, p + n - 1, ni, nj)
                first = start
                last = end
                if item["reversed"]:
                    first = end
                    last = start
                topo.uses[item["connector"]].append(
                    EdgeUse(domain.id, k, m, item["connector"], item["reversed"], first, last)
                )
                p += n - 1
    return topo


def one_to_one(a: EdgeUse, b: EdgeUse, same_direction: bool = True) -> dict:
    """CGNS 1-to-1 connection from use ``a`` to use ``b`` (0-based ranges).

    The connector's first/last nodes are matched in both domains; PointRange (in ``a``) is
    sorted, PointRangeDonor lists the donor images of its begin and end, and the transform maps
    a's index axes to b's (1-based, signed).
    """
    a_norm, a_max, a_run = a.face()
    b_norm, b_max, b_run = b.face()
    begin_a, end_a = sorted([a.first, a.last])
    # Donor image of a node of a: the same connector node in b (interfaces), or its periodic
    # image (``same_direction``: a's first node maps to b's first node).
    b_first = b.first
    b_last = b.last
    if not same_direction:
        b_first = b.last
        b_last = b.first
    if begin_a == a.first:
        begin_b, end_b = b_first, b_last
    else:
        begin_b, end_b = b_last, b_first
    running_a = end_a[a_run] - begin_a[a_run]
    running_b = end_b[b_run] - begin_b[b_run]
    sign_run = 1
    if (running_a > 0) != (running_b > 0):
        sign_run = -1
    sign_normal = -1
    if a_max != b_max:
        sign_normal = 1
    transform = [0, 0]
    transform[a_run] = sign_run * (b_run + 1)
    transform[a_norm] = sign_normal * (b_norm + 1)
    return {
        "range": (begin_a[0], begin_a[1], end_a[0], end_a[1]),
        "donor_range": (begin_b[0], begin_b[1], end_b[0], end_b[1]),
        "transform": tuple(transform),
    }


@dataclass
class BoundaryReport:
    boundary: list[str]
    interfaces: list[str]
    unused: list[str]
    unassigned: list[str]
    issues: list[Issue]
    # (source patch, target patch, [ConnectorPair]) for every consistent periodic pair
    periodic: list = field(default_factory=list)
    assignment: dict = field(default_factory=dict)  # boundary label -> patch name or None
    ids: dict = field(default_factory=dict)  # boundary label -> entity id to select

    @property
    def ok(self) -> bool:
        return not self.issues


def check_periodic(store, session, issues: list[Issue], matcher=None) -> list:
    """Every periodic patch is paired exactly once (from one side) and its boundary maps onto
    the partner's node by node: ``matcher(session, source, target)`` returns (the match, its
    issues); by default the 2D connector match (opengrid.model.periodic.match), in 3D the face
    match (opengrid.model.periodic_faces.match_faces). Returns [(source, target, match)]."""
    from opengrid.model.periodic import match

    if matcher is None:
        matcher = match
    periodic = [p for p in store.entities("patch") if p.params["patch_type"] == "periodic"]
    targeted: dict[int, list[str]] = {}
    for patch in store.entities("patch"):
        pairing = patch.params.get("periodic")
        if pairing is not None:
            targeted.setdefault(pairing["partner"], []).append(patch.name)
    pairs = []
    for patch in periodic:
        pairing = patch.params.get("periodic")
        sources = targeted.get(patch.id, [])
        if pairing is None:
            if not sources:
                issues.append(
                    Issue(patch.name, "periodic patch has no partner (use pair_periodic)")
                )
            elif len(sources) > 1:
                issues.append(
                    Issue(
                        patch.name,
                        f"is the periodic partner of several patches: {', '.join(sources)}",
                    )
                )
            continue
        partner = store.get(pairing["partner"])
        if partner.id == patch.id:
            issues.append(Issue(patch.name, "is paired with itself"))
            continue
        if partner.params["patch_type"] != "periodic":
            issues.append(
                Issue(
                    patch.name,
                    f"partner {partner.name} is a {partner.params['patch_type']} patch, "
                    "not periodic",
                )
            )
            continue
        if partner.params.get("periodic") is not None or sources:
            issues.append(
                Issue(patch.name, f"pair {patch.name} and {partner.name} from one side only")
            )
            continue
        if session is None:
            continue
        matched, match_issues = matcher(session, patch, partner)
        issues.extend(match_issues)
        if not match_issues:
            pairs.append((patch, partner, matched))
    return pairs


def validate_patches(store, topology: Topology, session=None) -> BoundaryReport:
    """Every boundary connector in exactly one patch; patches contain only boundary
    connectors; no non-manifold connectors; each domain in at most one zone; periodic patches
    paired consistently (node-by-node check when ``session`` is given)."""
    issues: list[Issue] = []
    boundary = topology.boundary()
    membership: dict[int, list[str]] = {}
    for patch in store.entities("patch"):
        for item in patch.params["members"]:
            if "connector" not in item:
                from opengrid.model.patches import member_label

                issues.append(
                    Issue(
                        member_label(item, store),
                        f"is a face member of patch {patch.name}, but the model has no blocks "
                        "(2D patches contain connectors)",
                    )
                )
                continue
            membership.setdefault(item["connector"], []).append(patch.name)
    for connector_id, uses in sorted(topology.non_manifold().items()):
        name = store.get(connector_id).name
        domains = sorted({store.get(u.domain).name for u in uses})
        issues.append(
            Issue(
                name, f"is used {len(uses)} times (domains {', '.join(domains)}); at most 2 allowed"
            )
        )
    unassigned = []
    for connector_id in sorted(boundary):
        name = store.get(connector_id).name
        patches = membership.get(connector_id, [])
        domain = store.get(boundary[connector_id].domain).name
        if not patches:
            unassigned.append(name)
            issues.append(Issue(name, f"lies on the boundary of {domain} but belongs to no patch"))
        elif len(patches) > 1:
            issues.append(Issue(name, f"belongs to several patches: {', '.join(patches)}"))
    for connector_id, patches in sorted(membership.items()):
        uses = topology.uses.get(connector_id, [])
        name = store.get(connector_id).name
        if len(uses) == 2:
            domains = " and ".join(store.get(u.domain).name for u in uses)
            issues.append(
                Issue(
                    name,
                    f"is an interface ({domains}), not a boundary; "
                    f"remove it from {', '.join(patches)}",
                )
            )
        elif not uses:
            issues.append(Issue(name, f"is in patch {', '.join(patches)} but used by no domain"))
    zones_of: dict[int, list[str]] = {}  # domains of either kind
    for zone in store.entities("zone"):
        for item in zone.params["members"]:
            zones_of.setdefault(item["domain"], []).append(zone.name)
    for domain_id, zones in sorted(zones_of.items()):
        if len(zones) > 1:
            issues.append(
                Issue(store.get(domain_id).name, f"belongs to several zones: {', '.join(zones)}")
            )
    interfaces = []
    for a, b in topology.interfaces():
        label = store.get(a.connector).name
        if a.domain == b.domain:
            interfaces.append(f"{label} (cut in {store.get(a.domain).name})")
        else:
            interfaces.append(f"{label} ({store.get(a.domain).name} | {store.get(b.domain).name})")
    unused = [store.get(c).name for c, u in sorted(topology.uses.items()) if not u]
    boundary_names = [store.get(c).name for c in sorted(boundary)]
    periodic = check_periodic(store, session, issues)
    assignment = {}
    ids = {}
    for connector_id in sorted(boundary):
        name = store.get(connector_id).name
        patches = membership.get(connector_id, [])
        assignment[name] = None
        if patches:
            assignment[name] = patches[0]
        ids[name] = connector_id
    return BoundaryReport(
        boundary_names, interfaces, unused, unassigned, issues, periodic, assignment, ids
    )
