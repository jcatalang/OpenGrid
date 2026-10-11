"""Gmsh MSH 4.1 (ASCII) export of a mesh assembly with physical groups (2D: surfaces and boundary
curves; 3D blocks: volumes of hexahedra and boundary surfaces of quads).

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
VTK_TO_GMSH = {5: GMSH_TRIANGLE, 9: GMSH_QUAD, 12: GMSH_HEXAHEDRON}


def write_gmsh(path, assembly) -> dict:
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames"]
    zones = [z for z in range(len(assembly.zone_names)) if np.any(assembly.cell_zone == z)]
    patches = [p for p in range(len(assembly.patch_names)) if np.any(assembly.line_patch == p)]
    dim = assembly.dimension
    boundary_kind = GMSH_LINE
    if dim == 3:
        boundary_kind = GMSH_QUAD
    lines.append(str(len(zones) + len(patches)))
    for p in patches:
        lines.append(f'{dim - 1} {p + 1} "{assembly.patch_names[p]}"')
    for z in zones:
        lines.append(f'{dim} {z + 1} "{assembly.zone_names[z]}"')
    lines.append("$EndPhysicalNames")

    def bbox(ids):
        coords = assembly.points[np.unique(ids)]
        low = coords.min(axis=0)
        high = coords.max(axis=0)
        return " ".join(repr(float(v)) for v in (*low, *high))

    lines.append("$Entities")
    if dim == 3:
        lines.append(f"0 0 {len(patches)} {len(zones)}")
    else:
        lines.append(f"0 {len(patches)} {len(zones)} 0")
    for p in patches:
        ids = assembly.lines[assembly.line_patch == p]
        lines.append(f"{p + 1} {bbox(ids)} 1 {p + 1} 0")
    for z in zones:
        ids = np.concatenate(assembly.cells(assembly.cell_zone == z))
        lines.append(f"{z + 1} {bbox(ids)} 1 {z + 1} 0")
    lines.append("$EndEntities")

    n = len(assembly.points)
    lines.append("$Nodes")
    lines.append(f"1 {n} 1 {n}")
    lines.append(f"{dim} {zones[0] + 1} 0 {n}")
    lines.extend(str(k + 1) for k in range(n))
    lines.extend(" ".join(repr(float(c)) for c in point) for point in assembly.points)
    lines.append("$EndNodes")

    blocks = []
    tag = 1
    for p in patches:
        cells = assembly.lines[assembly.line_patch == p]
        blocks.append((dim - 1, p + 1, boundary_kind, cells))
    for z in zones:
        for vtk_type, gmsh_type in sorted(VTK_TO_GMSH.items()):
            mask = (assembly.cell_zone == z) & (assembly.cell_type == vtk_type)
            if np.any(mask):
                blocks.append((dim, z + 1, gmsh_type, np.array(assembly.cells(mask))))
    total = sum(len(b[3]) for b in blocks)
    lines.append("$Elements")
    lines.append(f"{len(blocks)} {total} 1 {total}")
    for dim, entity, kind, cells in blocks:
        lines.append(f"{dim} {entity} {kind} {len(cells)}")
        for row in cells:
            lines.append(f"{tag} " + " ".join(str(int(v) + 1) for v in row))
            tag += 1
    lines.append("$EndElements")
    if assembly.periodic:
        # Master = source patch curve, slave = target curve; the affine matrix maps master
        # coordinates onto slave coordinates (row-major 4 x 4).
        from opengrid.model.periodic import affine_matrix

        lines.append("$Periodic")
        lines.append(str(len(assembly.periodic)))
        for source, target, transform, target_ids, source_ids in assembly.periodic:
            lines.append(f"1 {target + 1} {source + 1}")
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
        "periodic_links": len(assembly.periodic),
    }
