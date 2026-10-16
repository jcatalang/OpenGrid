"""CGNS export of straight-sided high-order elements (Phase 11f).

One base (cell and physical dimension 2 or 3) with one unstructured zone per OpenGrid zone
("Elements", or "Elements_<zone>" with several; zone family = the OpenGrid zone). Each zone
has:

* one element section per element kind, in the SIDS complete Lagrange types (TRI_6 ..
  HEXA_125, orders up to 4) and the SIDS node order (``io/lagrange.py``);
* one boundary section per patch and boundary kind (BAR_n in 2D; TRI_n, QUAD_n in 3D), each
  referenced by a BC_t (FamilySpecified, FamilyName = the patch; EdgeCenter / FaceCenter);
  periodic patches have no BC (their family has no BC type), as in the linear export;
* general vertex connectivity (Abutting1to1, point lists) to every other zone it shares nodes
  with, in both directions, and periodic pairs likewise with their translation.

One Family_t per zone and per patch (with its BC type). CGNS limits names to 32 characters.
"""

import numpy as np

from opengrid import _core
from opengrid.io import lagrange
from opengrid.io.cgns import MAX_NAME, _unique
from opengrid.model.patches import CGNS_BC_TYPES
from opengrid.validation import Issue, ValidationError

SECTION_NAMES = {
    "line": "Bars",
    "triangle": "Triangles",
    "quad": "Quadrilaterals",
    "tetra": "Tetrahedra",
    "pyramid": "Pyramids",
    "prism": "Prisms",
    "hexa": "Hexahedra",
}


def _rows(a: np.ndarray) -> np.ndarray:
    """Rows of an int64 array as one opaque value each (for sorting and searching)."""
    a = np.ascontiguousarray(a, dtype=np.int64)
    return a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1]))).reshape(-1)


def _face_owner(mesh) -> dict:
    """Zone of the cell owning each boundary element (kind -> (count,) zone index)."""
    nv = {kind: len(lagrange.VERTICES[kind]) for kind in lagrange.KINDS}
    keys = []
    zones = []
    for kind, conn in mesh.cells.items():
        corners = conn[:, : nv[kind]]
        if mesh.dim == 3:
            sides = lagrange.GMSH[kind]["faces"]
        else:
            sides = lagrange.GMSH[kind]["edges"]
        for side in sides:
            ids = np.sort(corners[:, list(side)], axis=1)
            if ids.shape[1] < 4 and mesh.dim == 3:
                ids = np.column_stack([ids, np.full(len(ids), -1, dtype=np.int64)])
            keys.append(ids)
            zones.append(mesh.cell_zone[kind])
    keys = _rows(np.vstack(keys))
    zones = np.concatenate(zones)
    order = np.argsort(keys, kind="stable")
    keys = keys[order]
    zones = zones[order]
    out = {}
    for kind, conn in mesh.faces.items():
        ids = np.sort(conn[:, : nv[kind]], axis=1)
        if ids.shape[1] < 4 and mesh.dim == 3:
            ids = np.column_stack([ids, np.full(len(ids), -1, dtype=np.int64)])
        wanted = _rows(ids)
        at = np.searchsorted(keys, wanted)
        at = np.minimum(at, len(keys) - 1)
        if np.any(keys[at] != wanted):
            raise ValidationError([Issue("mesh", "a boundary element is not a face of a cell")])
        out[kind] = zones[at]
    return out


def write_cgns_elements(path, mesh, base_name: str = "OpenGrid") -> dict:
    if not _core.cgns_available():
        raise RuntimeError("OpenGrid was built without CGNS support")
    p = mesh.order
    issues = []
    if p > lagrange.CGNS_MAX_ORDER:
        issues.append(
            Issue("high_order", f"CGNS defines elements up to order {lagrange.CGNS_MAX_ORDER}")
        )
    used_zones = sorted({int(z) for tags in mesh.cell_zone.values() for z in np.unique(tags)})
    used_patches = sorted({int(q) for tags in mesh.face_patch.values() for q in np.unique(tags)})
    names = [mesh.zone_names[z] for z in used_zones] + [mesh.patch_names[q] for q in used_patches]
    for name in [base_name, *names]:
        if len(name) > MAX_NAME:
            issues.append(Issue(name, f"is longer than {MAX_NAME} characters (CGNS limit)"))
    for _source, _target, transform, _t, _s in mesh.periodic:
        if transform.get("rotation_angle", 0.0) != 0.0:
            issues.append(
                Issue(
                    "periodic", "rotational periodicity is not written to CGNS (translation only)"
                )
            )
    if issues:
        raise ValidationError(issues)
    families = [(mesh.zone_names[z], "") for z in used_zones]
    periodic_patches = set()
    for q in used_patches:
        bc_type = CGNS_BC_TYPES.get(mesh.patch_types[q], "BCTypeUserDefined")
        if mesh.patch_types[q] == "periodic":
            bc_type = ""  # periodic faces are connectivity, not boundary conditions
            periodic_patches.add(q)
        families.append((mesh.patch_names[q], bc_type))
    owner = {}
    if mesh.faces:
        owner = _face_owner(mesh)
    zones = []
    zone_nodes = []  # per CGNS zone: sorted global node ids
    names_in_zone: list[set] = []
    # Zone names must differ from the families (both are children of the base).
    base_names = {name for name, _ in families}
    for z in used_zones:
        parts = [conn[mesh.cell_zone[kind] == z] for kind, conn in mesh.cells.items()]
        nodes = np.unique(np.concatenate([c.reshape(-1) for c in parts]))
        zone_nodes.append(nodes)
        used = set()
        names_in_zone.append(used)
        sections = []
        for kind in mesh.cells:
            rows = mesh.connectivity(kind, "cgns")[mesh.cell_zone[kind] == z]
            if len(rows):
                local = np.searchsorted(nodes, rows)
                name = _unique(SECTION_NAMES[kind], used)
                sections.append((name, lagrange.cgns_type(kind, p), local.reshape(-1).tolist(), ""))
        for q in used_patches:
            if q in periodic_patches:
                continue
            kinds = [k for k in mesh.faces if np.any((mesh.face_patch[k] == q) & (owner[k] == z))]
            for kind in kinds:
                mask = (mesh.face_patch[kind] == q) & (owner[kind] == z)
                rows = mesh.connectivity(kind, "cgns", boundary=True)[mask]
                local = np.searchsorted(nodes, rows)
                label = mesh.patch_names[q]
                if len(kinds) > 1:
                    label += "_" + SECTION_NAMES[kind][:4].lower()
                name = _unique(label, used)
                sections.append(
                    (
                        name,
                        lagrange.cgns_type(kind, p),
                        local.reshape(-1).tolist(),
                        mesh.patch_names[q],
                    )
                )
        xyz = mesh.points[nodes]
        label = "Elements"
        if len(used_zones) > 1:
            label = f"Elements_{mesh.zone_names[z]}"
        zones.append(
            (
                _unique(label, base_names),
                xyz[:, 0].tolist(),
                xyz[:, 1].tolist(),
                xyz[:, 2].tolist(),
                sections,
                mesh.zone_names[z],
            )
        )
    connections = []
    for a in range(len(zones)):
        for b in range(len(zones)):
            if a == b:
                continue
            shared = np.intersect1d(zone_nodes[a], zone_nodes[b])
            if len(shared) == 0:
                continue
            connections.append(
                (
                    a,
                    _unique(zones[b][0], names_in_zone[a]),
                    b,
                    (np.searchsorted(zone_nodes[a], shared) + 1).tolist(),
                    (np.searchsorted(zone_nodes[b], shared) + 1).tolist(),
                    None,
                )
            )
    for source, _target, transform, target_ids, source_ids in mesh.periodic:
        t = transform["translation"]
        label = mesh.patch_names[source]
        groups: dict = {}
        for tid, sid in zip(target_ids.tolist(), source_ids.tolist(), strict=True):
            zt = _zone_of(zone_nodes, tid)
            zs = _zone_of(zone_nodes, sid)
            groups.setdefault((zs, zt), []).append((sid, tid))
        for (zs, zt), pairs in sorted(groups.items()):
            s_ids = np.array([x for x, _ in pairs])
            t_ids = np.array([y for _, y in pairs])
            for this, donor, here, there, shift in (
                (zs, zt, s_ids, t_ids, (t[0], t[1])),
                (zt, zs, t_ids, s_ids, (-t[0], -t[1])),
            ):
                connections.append(
                    (
                        this,
                        _unique(f"{label}_periodic", names_in_zone[this]),
                        donor,
                        (np.searchsorted(zone_nodes[this], here) + 1).tolist(),
                        (np.searchsorted(zone_nodes[donor], there) + 1).tolist(),
                        shift,
                    )
                )
    _core.write_cgns_elements(str(path), base_name, mesh.dim, zones, families, connections)
    return {
        "zones": len(zones),
        "families": len(families),
        "sections": sum(len(z[4]) for z in zones),
        "connections": len(connections),
        "nodes": int(sum(len(n) for n in zone_nodes)),
        "elements": mesh.elements,
        "order": p,
    }


def _zone_of(zone_nodes, node: int) -> int:
    for k, nodes in enumerate(zone_nodes):
        at = np.searchsorted(nodes, node)
        if at < len(nodes) and nodes[at] == node:
            return k
    raise ValueError(f"node {node} is in no zone")
