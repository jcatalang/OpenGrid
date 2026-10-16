"""Gmsh MSH 4.1 (ASCII) export of a mesh assembly with physical groups (2D: surfaces and boundary
curves; 3D: volumes and boundary surfaces), linear or straight-sided high-order elements.

Physical groups: one 2D group per zone (quads and triangles), one 1D group per patch
(boundary lines). Each
zone is a surface entity and each patch a curve entity carrying its physical tag; all nodes are
listed in one block attached to the first surface. Paired periodic patches are written to the
$Periodic section (slave = target curve, master = source curve, affine transform and node
correspondences). Reference: Gmsh reference manual, section
"MSH file format" (version 4.1).
"""

from pathlib import Path

import numpy as np

GMSH_LINE = 1
GMSH_TRIANGLE = 2
GMSH_QUAD = 3
GMSH_HEXAHEDRON = 5
GMSH_TETRAHEDRON = 4
GMSH_PYRAMID = 7
GMSH_PRISM = 6
VTK_TO_GMSH = {
    5: GMSH_TRIANGLE,
    9: GMSH_QUAD,
    12: GMSH_HEXAHEDRON,
    10: GMSH_TETRAHEDRON,
    14: GMSH_PYRAMID,
    13: GMSH_PRISM,
}


def write_gmsh(path, assembly) -> dict:
    """Linear elements of a mesh assembly."""
    dim = assembly.dimension
    zones = [z for z in range(len(assembly.zone_names)) if np.any(assembly.cell_zone == z)]
    patches = [p for p in range(len(assembly.patch_names)) if np.any(assembly.line_patch == p)]

    def boundary_of(p):
        if dim == 3:
            return [assembly.boundary_cells[k] for k in np.flatnonzero(assembly.line_patch == p)]
        return list(assembly.lines[assembly.line_patch == p])

    blocks = []
    for p in patches:
        cells = boundary_of(p)
        if dim == 3:
            for size, kind in ((3, GMSH_TRIANGLE), (4, GMSH_QUAD)):
                group = [c for c in cells if len(c) == size]
                if group:
                    blocks.append((dim - 1, p + 1, kind, np.asarray(group)))
        else:
            blocks.append((dim - 1, p + 1, GMSH_LINE, np.asarray(cells)))
    for z in zones:
        for vtk_type, gmsh_type in sorted(VTK_TO_GMSH.items()):
            mask = (assembly.cell_zone == z) & (assembly.cell_type == vtk_type)
            if np.any(mask):
                blocks.append((dim, z + 1, gmsh_type, np.array(assembly.cells(mask))))
    return _write(
        path,
        assembly.points,
        dim,
        [(z, assembly.zone_names[z]) for z in zones],
        [(p, assembly.patch_names[p]) for p in patches],
        blocks,
        assembly.periodic,
    )


def write_gmsh_elements(path, mesh) -> dict:
    """Straight-sided order-p elements (``model/elements.py``): Gmsh's high-order element
    types (e.g. 9 = 6-node triangle, 12 = 27-node hexahedron, 93 = 125-node hexahedron) with
    the same physical groups, entities and periodic links as the linear export."""
    from opengrid.io import lagrange

    dim = mesh.dim
    zones = sorted({int(z) for tags in mesh.cell_zone.values() for z in np.unique(tags)})
    patches = sorted({int(p) for tags in mesh.face_patch.values() for p in np.unique(tags)})
    blocks = []
    for p in patches:
        for kind, conn in mesh.faces.items():
            rows = conn[mesh.face_patch[kind] == p]
            if len(rows):
                blocks.append((dim - 1, p + 1, lagrange.gmsh_type(kind, mesh.order), rows))
    for z in zones:
        for kind, conn in mesh.cells.items():
            rows = conn[mesh.cell_zone[kind] == z]
            if len(rows):
                blocks.append((dim, z + 1, lagrange.gmsh_type(kind, mesh.order), rows))
    info = _write(
        path,
        mesh.points,
        dim,
        [(z, mesh.zone_names[z]) for z in zones],
        [(p, mesh.patch_names[p]) for p in patches],
        blocks,
        mesh.periodic,
    )
    info["order"] = mesh.order
    return info


def _write(path, points, dim, zones, patches, blocks, periodic) -> dict:
    """zones / patches: (index, name); blocks: (dim, entity tag, Gmsh type, (count, nodes)
    0-based node ids); periodic: (source, target, transform, target ids, source ids)."""
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames"]
    lines.append(str(len(zones) + len(patches)))
    for p, name in patches:
        lines.append(f'{dim - 1} {p + 1} "{name}"')
    for z, name in zones:
        lines.append(f'{dim} {z + 1} "{name}"')
    lines.append("$EndPhysicalNames")

    def bbox(ids):
        coords = points[np.unique(ids)]
        low = coords.min(axis=0)
        high = coords.max(axis=0)
        return " ".join(repr(float(v)) for v in (*low, *high))

    lines.append("$Entities")
    if dim == 3:
        lines.append(f"0 0 {len(patches)} {len(zones)}")
    else:
        lines.append(f"0 {len(patches)} {len(zones)} 0")
    for entity_dim, group in ((dim - 1, patches), (dim, zones)):
        for index, _ in group:
            ids = np.concatenate(
                [
                    np.asarray(c).reshape(-1)
                    for d, e, _, c in blocks
                    if d == entity_dim and e == index + 1
                ]
            )
            lines.append(f"{index + 1} {bbox(ids)} 1 {index + 1} 0")
    lines.append("$EndEntities")

    n = len(points)
    lines.append("$Nodes")
    lines.append(f"1 {n} 1 {n}")
    lines.append(f"{dim} {zones[0][0] + 1} 0 {n}")
    lines.extend(str(k + 1) for k in range(n))
    lines.extend(" ".join(repr(float(c)) for c in point) for point in points)
    lines.append("$EndNodes")

    total = sum(len(b[3]) for b in blocks)
    lines.append("$Elements")
    lines.append(f"{len(blocks)} {total} 1 {total}")
    tag = 1
    for block_dim, entity, kind, cells in blocks:
        lines.append(f"{block_dim} {entity} {kind} {len(cells)}")
        for row in cells:
            lines.append(f"{tag} " + " ".join(str(int(v) + 1) for v in row))
            tag += 1
    lines.append("$EndElements")
    if periodic:
        # Master = source patch entity, slave = target entity; the affine matrix maps master
        # coordinates onto slave coordinates (row-major 4 x 4).
        from opengrid.model.periodic import affine_matrix

        lines.append("$Periodic")
        lines.append(str(len(periodic)))
        for source, target, transform, target_ids, source_ids in periodic:
            lines.append(f"{dim - 1} {target + 1} {source + 1}")
            matrix = affine_matrix(transform).reshape(-1)
            lines.append("16 " + " ".join(repr(float(v)) for v in matrix))
            lines.append(str(len(target_ids)))
            for slave, master in zip(target_ids, source_ids, strict=True):
                lines.append(f"{int(slave) + 1} {int(master) + 1}")
        lines.append("$EndPeriodic")
    Path(path).write_text("\n".join(lines) + "\n")
    return {
        "nodes": n,
        "elements": total,
        "physical_groups": len(zones) + len(patches),
        "periodic_links": len(periodic),
    }
