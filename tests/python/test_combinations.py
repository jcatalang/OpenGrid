"""Phase 12b: end-to-end checks of feature combinations (templates, wall layers, smoothing
groups, sketches, hybrid meshes, periodicity, zones, 3D layers) through every export.

One checker for every model: the merged mesh is consistent (every point used, every boundary
element an edge / face of a cell, the patches cover exactly the mesh boundary), every export
is read back by its reference reader with the right counts, the curved high-order meshes are
valid and their boundary matches, the journal replays to the same model and the project
round-trips. The connector-row bug of Phase 11g (patch lines pointing at stray points after
group smoothing) is the class of fault these catch.
"""

import shutil
import subprocess
from collections import Counter

import h5py
import numpy as np
import pytest
from test_periodic import channel
from test_templates import airfoil
from test_unstructured import hybrid

import opengrid
from opengrid import _core
from opengrid.io import lagrange
from opengrid.model.shape import jacobians, jacobians_at_centre
from opengrid.project import to_json

VTK_NAMES = {5: "triangle", 9: "quad", 10: "tetra", 12: "hexahedron", 13: "wedge", 14: "pyramid"}


# ------------------------------------------------------------------ the checker
def _sides(cells, types, dim):
    """Edges (2D) or faces (3D) of every cell as sorted vertex tuples."""
    out = []
    for cell, t in zip(cells, types, strict=True):
        cell = [int(v) for v in cell]
        if dim == 2:
            n = len(cell)
            out.extend(tuple(sorted((cell[k], cell[(k + 1) % n]))) for k in range(n))
        else:
            from opengrid.model.volumes import CELL_FACES

            key = {10: "tets", 14: "pyramids", 13: "prisms", 12: "hexes"}[int(t)]
            for face in CELL_FACES[key]:
                out.append(tuple(sorted(cell[k] for k in face)))
    return out


def check_assembly(assembly) -> dict:
    cells = assembly.cells()
    dim = assembly.dimension
    used = np.zeros(len(assembly.points), dtype=bool)
    used[assembly.connectivity] = True
    assert used.all(), f"{int((~used).sum())} points used by no cell"
    assert len(np.unique(np.round(assembly.points, 12), axis=0)) == len(assembly.points)
    counts = Counter(_sides(cells, assembly.cell_type, dim))
    boundary = {side for side, n in counts.items() if n == 1}
    assert all(n <= 2 for n in counts.values()), "a side shared by more than two cells"
    if dim == 2:
        patched = {tuple(sorted(int(v) for v in line)) for line in assembly.lines}
    else:
        patched = {tuple(sorted(int(v) for v in face)) for face in assembly.boundary_cells}
    assert patched <= set(counts), "a boundary element that is not a side of a cell"
    assert patched == boundary, (
        f"patches cover {len(patched & boundary)} of {len(boundary)} boundary sides, "
        f"{len(patched - boundary)} interior"
    )
    return {"cells": len(cells), "points": len(assembly.points), "boundary": len(boundary)}


def _cgnscheck(path) -> None:
    checker = shutil.which("cgnscheck")
    if checker is None:
        return
    result = subprocess.run(
        [checker, str(path)], capture_output=True, text=True, check=False, timeout=600
    )
    log = result.stdout + result.stderr
    assert result.returncode == 0 and "ERROR" not in log.upper(), log


def check_exports(og, tmp_path, counts) -> None:
    vtk = pytest.importorskip("vtk")
    meshio = pytest.importorskip("meshio")
    # VTKHDF: every cell once, across the collection's blocks.
    path = tmp_path / "m.vtkhdf"
    og.export_vtkhdf(path)
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    output = reader.GetOutput()
    total = 0
    boundary_cells = 0
    if output.IsA("vtkPartitionedDataSetCollection"):
        iterator = output.NewIterator()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            block = iterator.GetCurrentDataObject()
            if block is not None:
                volume = sum(
                    1
                    for c in range(block.GetNumberOfCells())
                    if block.GetCell(c).GetCellDimension() == og_dim(og)
                )
                total += volume
                boundary_cells += block.GetNumberOfCells() - volume
            iterator.GoToNextItem()
    else:
        total = output.GetNumberOfCells()
    assert total == counts["cells"]
    if boundary_cells:
        assert boundary_cells == counts["boundary"]
    # Gmsh: cells by type and the physical groups.
    path = tmp_path / "m.msh"
    summary = og.export_gmsh(path)
    mesh = meshio.read(path)
    by_type = Counter()
    for block in mesh.cells:
        by_type[block.type] += len(block.data)
    volume_types = {"tetra", "pyramid", "wedge", "hexahedron"}
    boundary_types = {"triangle", "quad"}
    if og_dim(og) == 2:
        volume_types = {"triangle", "quad"}
        boundary_types = {"line"}
    assert sum(n for t, n in by_type.items() if t in volume_types) == counts["cells"]
    assert sum(n for t, n in by_type.items() if t in boundary_types) == counts["boundary"]
    assert summary["elements"] == counts["cells"] + counts["boundary"]
    # CGNS: the library's checker, and the zones exist.
    if _core.cgns_available():
        path = tmp_path / "m.cgns"
        og.export_cgns(path)
        _cgnscheck(path)
        with h5py.File(path, "r") as f:
            base = next(k for k in f if f[k].attrs.get("label") == b"CGNSBase_t")
            zones = [k for k in f[base] if f[base][k].attrs.get("label") == b"Zone_t"]
            assert zones


def og_dim(og) -> int:
    if og.is_3d():
        return 3
    return 2


def check_high_order(og, tmp_path, counts, orders=(2, 3)) -> None:
    for p in orders:
        mesh = og.high_order(p)
        assert mesh.elements == counts["cells"]
        assert mesh.boundary_elements == counts["boundary"]
        for kind, conn in mesh.cells.items():
            det = jacobians(kind, p, mesh.points[conn], mesh.dim)
            sign = np.sign(jacobians_at_centre(kind, p, mesh.points[conn], mesh.dim))
            assert (det * sign[:, None]).min() > 0, f"{kind} folded at order {p}"
        # Shared nodes: the boundary elements use the cells' nodes.
        used = np.zeros(len(mesh.points), dtype=bool)
        for conn in mesh.cells.values():
            used[conn.reshape(-1)] = True
        assert used.all()
        for conn in mesh.faces.values():
            assert used[conn.reshape(-1)].all()
        if p == 2:
            path = tmp_path / "ho.msh"
            og.export_high_order(path, order=2)
            meshio = pytest.importorskip("meshio")
            mesh2 = meshio.read(path)
            second = {"tetra10", "pyramid14", "wedge18", "hexahedron27"}
            second_boundary = {"triangle6", "quad9"}
            if mesh.dim == 2:
                second = {"triangle6", "quad9"}
                second_boundary = {"line3"}
            assert sum(len(b.data) for b in mesh2.cells if b.type in second) == counts["cells"]
            boundary = sum(len(b.data) for b in mesh2.cells if b.type in second_boundary)
            assert boundary == counts["boundary"]
            if _core.cgns_available():
                og.export_high_order(tmp_path / "ho.cgns", order=2)
                _cgnscheck(tmp_path / "ho.cgns")
        if len(lagrange.KINDS) and p == 3 and _core.cgns_available():
            og.export_high_order(tmp_path / "ho3.cgns", order=3)
            _cgnscheck(tmp_path / "ho3.cgns")


def check_replay_and_project(og, tmp_path) -> None:
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)
    path = tmp_path / "m.ogp"
    og.save(path)
    again = opengrid.Session()
    again.load(path)
    assert to_json(again.store) == to_json(og.store)
    a = og._linear_assembly()
    b = again._linear_assembly()
    assert np.array_equal(a.points, b.points) and np.array_equal(a.connectivity, b.connectivity)


def check_model(og, tmp_path, orders=(2, 3)) -> dict:
    assert og.boundary_report().ok, [str(i) for i in og.boundary_report().issues]
    counts = check_assembly(og._linear_assembly())
    check_exports(og, tmp_path, counts)
    check_high_order(og, tmp_path, counts, orders)
    check_replay_and_project(og, tmp_path)
    return counts


# ------------------------------------------------------------------ models
def boundary_connectors(og) -> list[str]:
    return sorted(og.store.get(c).name for c in og.topology().boundary())


def patch_everything(og, wall_curves=(), wall_name="walls", other="far") -> None:
    """Patches for every boundary connector (2D) or face (3D): the walls and the rest."""
    if og.is_3d():
        from opengrid.model.volumes import analyze_3d

        store = og.store
        members = []
        seen = set()
        for entry in analyze_3d(og)["entries"]:
            if entry["patch"] is not None:
                continue
            domains = [store.get(key[1]).name for key in entry["keys"] if key[0] == "domain"]
            if domains:
                if domains[0] not in seen:
                    seen.add(domains[0])
                    members.append(domains[0])
            else:  # a block face without a domain
                _kind, block, face = next(key for key in entry["keys"] if key[0] == "block")
                members.append({"block": store.get(block).name, "face": face})
        og.create_patch(wall_name, "wall", members)
        return
    ids = {og.store.get(c).id for c in wall_curves}
    walls, rest = [], []
    for name in boundary_connectors(og):
        segments = og.store.get(name).params["segments"]
        if ids and all(s["curve"] in ids for s in segments):
            walls.append(name)
        else:
            rest.append(name)
    if walls:
        og.create_patch(wall_name, "wall", walls)
    if rest:
        og.create_patch(other, "farfield", rest)


def cgrid_layers_group():
    og = opengrid.Session()
    body = airfoil(og)
    og.create_template("cgrid", body=body, wall=1e-5, layers=12, radial_points=(41, 25), name="c")
    og.smooth_group("c-group")
    patch_everything(og, body, "airfoil")
    return og


def channel_periodic_zones():
    og = channel()
    og.pair_periodic("left", "right")
    og.create_zone("bottom", ["lower"])
    og.create_zone("top", ["upper"])
    return og


def hybrid_ring():
    return hybrid()


def cylinder_layers():
    from test_inflation import cylinder

    return cylinder(
        inflation={"connectors": ["c-wall"], "first_height": 1e-3, "growth_rate": 1.2, "layers": 8}
    )


def ogrid_group():
    from test_groups import circle_ogrid

    og = opengrid.Session()
    names = circle_ogrid(og)
    og.smooth_group(og.create_smoothing_group(names, tolerance=1e-10))
    patch_everything(og, [f"a{k}" for k in range(4)], "wall")
    return og


def sketch_square():
    from test_sketch import square

    og = square()
    patch_everything(og)
    return og


def pipe_template():
    og = opengrid.Session()
    og.create_template(
        "pipe", boundary=[og.create_circle((0, 0), 0.5, name="c")], length=2.0, span_points=5
    )
    og.smooth_group(og.names("smoothing_group")[0])
    patch_everything(og)
    return og


def bump_channel_layers():
    from test_templates import channel

    og = opengrid.Session()
    og.create_template(
        "bump_channel",
        **channel(og, bump=True),
        blocks=(3, 1),
        wall=5e-3,
        length=1.0,
        three_d_blocks=True,
        layers=6,
    )
    og.smooth_group(og.names("smoothing_group")[0])
    patch_everything(og)
    return og


def sphere_box_layers():
    pytest.importorskip("netgen.meshing")
    from test_curved_elements import sphere_in_box

    og, _radius = sphere_in_box("structured")
    return og


CASES = {
    "cgrid_layers_group": (cgrid_layers_group, (2, 3)),
    "channel_periodic_zones": (channel_periodic_zones, (2, 3)),
    "hybrid_ring": (hybrid_ring, (2, 3)),
    "cylinder_layers": (cylinder_layers, (2, 3)),
    "ogrid_group": (ogrid_group, (2, 3)),
    "sketch_square": (sketch_square, (2, 3)),
    "pipe_template": (pipe_template, (2,)),
    "bump_channel_layers": (bump_channel_layers, (2,)),
    "sphere_box_layers": (sphere_box_layers, (2,)),
}


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_combination(case, tmp_path):
    build, orders = CASES[case]
    og = build()
    counts = check_model(og, tmp_path, orders)
    assert counts["cells"] > 0 and counts["boundary"] > 0
