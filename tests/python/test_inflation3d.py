"""Phase 8b: 3D inflation layers (prisms / hexahedra) in volumes, sliding on planar faces."""

import itertools
import math
import shutil
import subprocess
from collections import Counter

import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.model.volumes import cell_volumes, total_volume, volume_cells
from opengrid.project import to_json

pytest.importorskip("netgen.meshing")

from test_volumes import box


def body_in_box(kind="unstructured", n=11, layers=6):
    og = opengrid.Session()
    outer = box(og, (0, 0, 0), (4, 3, 3), 9, "o", "structured")
    inner = box(og, (1.5, 1, 1), (2.5, 2, 2), n, "i", kind)
    inflation = {"faces": inner, "first_height": 0.01, "growth_rate": 1.25, "layers": layers}
    og.create_volume(
        outer + inner, generator={"max_size": 0.4, "inflation": inflation}, name="fluid"
    )
    og.create_patch("far", "farfield", outer)
    og.create_patch("body", "wall", inner)
    return og, outer, inner


def cube_on_floor():
    """A unit cube standing on the floor of a channel; the floor (unstructured, planar, with
    the cube's footprint as a hole) is the sliding plane."""
    og = opengrid.Session()
    lo, hi = (0.0, 0.0, 0.0), (5.0, 4.0, 3.0)
    corners = {c: tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3)) for c in _cube()}
    edges = {}

    def edge(a, b, pts, tag, n):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"{tag}{len(edges)}"
            og.create_line(pts[key[0]], pts[key[1]], name=name)
            edges[key] = og.create_connector(name, n, name=f"c{name}")
        return edges[key], key != (a, b)

    loops = {
        "xmin": [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
        "xmax": [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
        "ymin": [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)],
        "ymax": [(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
        "zmax": [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
        "zmin": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
    }
    outer = []
    floor_loop = []
    for name, cs in loops.items():
        ed = [edge(cs[k], cs[(k + 1) % 4], corners, "b", 6) for k in range(4)]
        if name == "zmin":
            floor_loop = [e[0] for e in ed]
            continue
        og.create_domain([[e] for e in ed], name=f"far_{name}")
        outer.append(f"far_{name}")
    edges.clear()
    clo, chi = (2.0, 1.5, 0.0), (3.0, 2.5, 1.0)
    cube = {c: tuple(clo[k] + c[k] * (chi[k] - clo[k]) for k in range(3)) for c in _cube()}
    walls = []
    for name, cs in loops.items():
        ed = [edge(cs[k], cs[(k + 1) % 4], cube, "k", 9)[0] for k in range(4)]
        if name == "zmin":
            foot = ed
            continue
        og.create_unstructured_domain(connectors=ed, name=f"cube_{name}")
        walls.append(f"cube_{name}")
    og.create_unstructured_domain(connectors=floor_loop + foot, name="floor")
    og.create_volume(
        [*outer, "floor", *walls],
        generator={
            "max_size": 0.5,
            "inflation": {"faces": walls, "first_height": 0.01, "growth_rate": 1.25, "layers": 8},
        },
        name="fluid",
    )
    og.create_patch("far", "farfield", outer)
    og.create_patch("walls", "wall", [*walls, "floor"])
    return og, walls


def _cube():
    return list(itertools.product((0, 1), repeat=3))


def test_prism_layers_on_a_body():
    og, _outer, inner = body_in_box()
    mesh = og.volume_mesh("fluid")
    walls = {og.entity(f).id for f in inner}
    n_wall = sum(len(mesh["faces"][w]) for w in walls)
    info = mesh["inflation"]
    assert info["max_layers"] == 6
    assert info["stops"] == {}
    assert len(mesh["prisms"]) == 6 * n_wall
    assert len(mesh["hexes"]) == 0
    assert total_volume(mesh) == pytest.approx(35.0, rel=1e-12)
    for key, cells, _ in volume_cells(mesh):
        assert np.all(cell_volumes(mesh["points"], key, cells) > 0), key
    q = og.volume_quality("fluid")
    assert q["scaled_jacobian"].min() > 0.1
    # Achieved heights: every wall node got 6 layers, nominal height first * (g^6 - 1) / (g - 1).
    wall_nodes = np.unique(np.concatenate([np.concatenate(mesh["faces"][w]) for w in walls]))
    assert np.all(info["layers"][wall_nodes] == 6)
    expected = 0.01 * (1.25**6 - 1.0) / 0.25
    assert np.allclose(info["height"][wall_nodes], expected, rtol=1e-12)
    report = og.boundary_report()
    assert report.ok, report.issues
    assert "prisms" in og.volume_status("fluid")


def test_hexahedral_layers_on_structured_walls():
    og, _outer, _inner = body_in_box(kind="structured", n=9, layers=4)
    mesh = og.volume_mesh("fluid")
    assert len(mesh["hexes"]) == 4 * 6 * 64
    assert len(mesh["prisms"]) == 0
    # The front is quads: capped by pyramids towards the tetrahedra.
    assert len(mesh["pyramids"]) >= 6 * 64
    assert total_volume(mesh) == pytest.approx(35.0, rel=1e-12)
    q = og.volume_quality("fluid")
    assert q["volume"].min() > 0


def test_layers_slide_on_a_floor(tmp_path):
    og, _walls = cube_on_floor()
    mesh = og.volume_mesh("fluid")
    info = mesh["inflation"]
    floor = og.entity("floor").id
    assert info["sliding"] == (floor,)
    assert total_volume(mesh) == pytest.approx(5 * 4 * 3 - 1, rel=1e-12)
    for key, cells, _ in volume_cells(mesh):
        assert np.all(cell_volumes(mesh["points"], key, cells) > 0), key
    # The re-meshed floor: in the plane, its cells are the layer sides (quads / triangles) and
    # triangles; it covers the floor area minus the footprint.
    cells = mesh["faces"][floor]
    pts = mesh["points"]
    assert np.all(pts[np.unique(np.concatenate(cells))][:, 2] == 0.0)
    area = 0.0
    for c in cells:
        p = pts[c]
        area += 0.5 * np.cross(p[1] - p[0], p[2] - p[0])[2]
        if len(c) == 4:
            area += 0.5 * np.cross(p[2] - p[0], p[3] - p[0])[2]
    assert abs(area) == pytest.approx(20.0 - 1.0, rel=1e-12)
    assert Counter(len(c) for c in cells)[4] > 0
    # Floor nodes along the cube rose along the floor, as high as their columns.
    report = og.boundary_report()
    assert report.ok, report.issues
    # Exports.
    meshio = pytest.importorskip("meshio")
    summary = og.export_gmsh(tmp_path / "f.msh")
    read = meshio.read(tmp_path / "f.msh")
    counts = Counter()
    for block in read.cells:
        counts[block.type] += len(block.data)
    assert counts["wedge"] == len(mesh["prisms"])
    assert counts["tetra"] == len(mesh["tets"])
    assert counts["pyramid"] == len(mesh["pyramids"])
    assert summary["physical_groups"] == 3
    vtk = pytest.importorskip("vtk")
    og.export_vtkhdf(tmp_path / "f.vtkhdf")
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(tmp_path / "f.vtkhdf"))
    reader.Update()
    grid = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.Update()
    from vtk.util.numpy_support import vtk_to_numpy

    volumes = vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Volume"))
    assert volumes.min() > 0
    assert volumes.sum() == pytest.approx(59.0, rel=1e-9)
    if _core.cgns_available():
        path = tmp_path / "f.cgns"
        og.export_cgns(path)
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
            assert "error" not in out.lower() and "warning" not in out.lower(), out
    # Journal and project round trips give the same mesh.
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)
    og.save(tmp_path / "f.ogp")
    loaded = opengrid.Session()
    loaded.load(tmp_path / "f.ogp")
    again = loaded.volume_mesh("fluid")
    assert np.array_equal(again["points"], mesh["points"])
    assert np.array_equal(again["prisms"], mesh["prisms"])


def test_without_sliding_the_layers_taper_at_the_floor():
    og, _walls = cube_on_floor()
    generator = dict(og.entity("fluid").params["generator"])
    inflation = dict(generator["inflation"])
    inflation["slide"] = False
    og.set_volume_generator("fluid", inflation=inflation)
    mesh = og.volume_mesh("fluid")
    assert mesh["inflation"]["sliding"] == ()
    assert mesh["inflation"]["stops"]["fixed"] > 0
    assert total_volume(mesh) == pytest.approx(59.0, rel=1e-12)
    # The floor keeps its own triangulation.
    floor = og.entity("floor")
    own = og.unstructured_mesh("floor")
    assert len(mesh["faces"][floor.id]) == len(own["triangles"])


def test_inflation_settings_are_validated():
    og, _outer, inner = body_in_box()
    with pytest.raises(opengrid.ValidationError, match="not a face of the volume"):
        og.set_volume_generator("fluid", inflation={"faces": ["o-e0"], "first_height": 0.1})
    with pytest.raises(opengrid.ValidationError, match="first_height"):
        og.set_volume_generator("fluid", inflation={"faces": inner})
    with pytest.raises(opengrid.ValidationError, match="unexpected key"):
        og.set_volume_generator(
            "fluid", inflation={"faces": inner, "first_height": 0.1, "fans": True}
        )
    with pytest.raises(opengrid.ValidationError, match="growth_rate"):
        og.set_volume_generator(
            "fluid", inflation={"faces": inner, "first_height": 0.1, "growth_rate": 0.9}
        )
    params = og.entity("fluid").params["generator"]["inflation"]
    assert params["faces"] == tuple(og.entity(f).id for f in inner)
    assert params["slide"] is True
    # The journal names the wall faces.
    assert "'i-f00'" in og.journal_text()


def test_total_thickness_limits_the_layers():
    og, _outer, inner = body_in_box(layers=30)
    og.set_volume_generator(
        "fluid",
        inflation={
            "faces": inner,
            "first_height": 0.01,
            "growth_rate": 1.25,
            "total_thickness": 0.05,
            "transition": False,
        },
    )
    info = og.volume_mesh("fluid")["inflation"]
    # Smallest N with 0.01 (1.25^N - 1) / 0.25 >= 0.05: N = 4.
    assert info["max_layers"] == 4
    assert math.isclose(info["height"].max(), 0.01 * (1.25**4 - 1) / 0.25)


def test_gui_volume_inflation_editor(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    outer = box(og, (0, 0, 0), (4, 3, 3), 7, "o", "structured")
    inner = box(og, (1.5, 1, 1), (2.5, 2, 2), 7, "i", "unstructured")
    og.create_volume(outer + inner, generator={"max_size": 0.5}, name="fluid")
    window.set_selection([og.entity("fluid").id])
    _kind, editor = window.property_panel._editors["generator"]
    editor.walls.setText(", ".join(inner))
    editor.first_height.setText("0.02")
    editor.layers.setValue(3)
    assert window.property_panel.apply()
    inflation = og.entity("fluid").params["generator"]["inflation"]
    assert inflation["layers"] == 3 and inflation["first_height"] == 0.02
    window.set_selection([og.entity("fluid").id])
    assert window.property_panel.read_changes() == (None, {})
    assert "fluid" in window.model_tree.entity_names("Blocks")
    assert "prisms" in og.volume_status("fluid")
    window.quality_action.setChecked(True)
    assert window.show_worst_cell() is not None
    window.close()
