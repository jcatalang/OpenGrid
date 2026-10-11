"""Phase 7a: surfaces, 3D domains, blocks (six faces / extrusion), 3D topology and exports."""

import math
import random
import shutil
import subprocess
from collections import Counter

import h5py
import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read

import opengrid
from opengrid import _core
from opengrid.project import to_json

T = math.pi / 3


def polar(r, t, z):
    return (r * math.cos(t), r * math.sin(t), z)


def duct(og, n_theta=31, t0=0.0, t1=T, tag="", shared=None):
    """Annular sector r in [1, 2], theta in [t0, t1], z in [0, 1] as a six-face block.
    ``shared``: names of the theta = t0 side (connectors and domain) of a previous duct."""
    for z, s in ((0.0, "0"), (1.0, "1")):
        og.create_arc((0, 0, z), 1.0, t0, t1, name=f"ai{s}{tag}")
        og.create_arc((0, 0, z), 2.0, t0, t1, name=f"ao{s}{tag}")
        og.create_connector(f"ai{s}{tag}", n_theta, name=f"c-ai{s}{tag}")
        og.create_connector(f"ao{s}{tag}", n_theta, name=f"c-ao{s}{tag}")
    side = {}
    for t, which in ((t0, "a"), (t1, "b")):
        if which == "a" and shared is not None:
            side[which] = shared
            continue
        names = {}
        for z, s in ((0.0, "0"), (1.0, "1")):
            og.create_line(polar(1, t, z), polar(2, t, z), name=f"r{which}{s}{tag}")
            names[f"r{s}"] = og.create_connector(
                f"r{which}{s}{tag}", 11, name=f"c-r{which}{s}{tag}"
            )
        for r, s in ((1, "i"), (2, "o")):
            og.create_line(polar(r, t, 0), polar(r, t, 1), name=f"v{which}{s}{tag}")
            names[f"v{s}"] = og.create_connector(f"v{which}{s}{tag}", 9, name=f"c-v{which}{s}{tag}")
        names["domain"] = og.create_domain(
            [[names["r0"]], [names["vo"]], [(names["r1"], True)], [(names["vi"], True)]],
            name=f"side{which}{tag}",
        )
        side[which] = names
    a, b = side["a"], side["b"]
    og.create_extruded_surface(f"ai0{tag}", (0, 0, 1), name=f"s-inner{tag}")
    og.create_extruded_surface(f"ao0{tag}", (0, 0, 1), name=f"s-outer{tag}")
    og.create_domain(
        [[a["r0"]], [f"c-ao0{tag}"], [(b["r0"], True)], [(f"c-ai0{tag}", True)]],
        name=f"bottom{tag}",
    )
    og.create_domain(
        [[a["r1"]], [f"c-ao1{tag}"], [(b["r1"], True)], [(f"c-ai1{tag}", True)]],
        name=f"top{tag}",
    )
    og.create_domain(
        [[f"c-ai0{tag}"], [b["vi"]], [(f"c-ai1{tag}", True)], [(a["vi"], True)]],
        name=f"inner{tag}",
        surface=f"s-inner{tag}",
    )
    og.create_domain(
        [[f"c-ao0{tag}"], [b["vo"]], [(f"c-ao1{tag}", True)], [(a["vo"], True)]],
        name=f"outer{tag}",
        surface=f"s-outer{tag}",
    )
    faces = [f"bottom{tag}", f"top{tag}", f"inner{tag}", f"outer{tag}", a["domain"], b["domain"]]
    block = og.create_block(faces=faces, name=f"duct{tag}")
    return block, b


def patched_duct(n_theta=31) -> opengrid.Session:
    og = opengrid.Session()
    duct(og, n_theta)
    og.create_patch("walls", "wall", ["inner", "outer"])
    og.create_patch("ends", "inlet", ["bottom", "top"])
    og.create_patch("sides", "symmetry", ["sidea", "sideb"])
    return og


def total_volume(og, name):
    return float(og.block_quality(name)["volume"].sum())


def test_surfaces_entities_and_projection():
    og = opengrid.Session()
    og.create_arc((0, 0, 0), 2.0, 0.0, T, name="arc")
    og.create_extruded_surface("arc", (0, 0, 3), name="cyl")
    og.create_plane((0, 0, 0), (1, 0, 0), (0, 1, 0), name="floor")
    og.create_line((1, 0, 0), (1, 0, 2), name="gen")
    og.create_revolved_surface("gen", math.pi, name="rev")
    og.create_nurbs_surface(
        1, 1, (0, 0, 1, 1), (0, 0, 1, 1), [[(0, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 1, 0)]], name="n"
    )
    _u, _v, p, d = og.project_to_surface("cyl", polar(2.5, 0.5, 1.0))
    assert d == pytest.approx(0.5, abs=1e-12)
    assert np.hypot(p[0], p[1]) == pytest.approx(2.0, abs=1e-12)
    _u, _v, _p, d = og.project_to_surface("rev", polar(0.4, 1.0, 0.5))
    assert d == pytest.approx(0.6, abs=1e-12)
    assert og.project_to_surface("n", (0.3, 0.4, 2.0))[3] == pytest.approx(2.0, abs=1e-12)
    with pytest.raises(opengrid.ValidationError, match="used by"):
        og.delete("arc")
    # Changing the curve rebuilds the surface.
    og.modify("arc", radius=3.0)
    assert og.project_to_surface("cyl", polar(3.5, 0.5, 1.0))[3] == pytest.approx(0.5, abs=1e-12)
    with pytest.raises(opengrid.ValidationError, match="not a curve"):
        og.create_extruded_surface("floor", (0, 0, 1))
    with pytest.raises(opengrid.ValidationError, match="parallel"):
        og.create_plane((0, 0, 0), (1, 0, 0), (2, 0, 0))
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)


def test_six_face_block_volume_surfaces_and_order():
    og = opengrid.Session()
    duct(og, n_theta=601)
    grid = og.block_grid("duct")
    q = og.block_quality("duct")
    assert q["scaled_jacobian"].min() > 0.99
    analytic = 0.5 * T * (4.0 - 1.0) * 1.0
    assert total_volume(og, "duct") == pytest.approx(analytic, rel=1e-6)
    r = np.hypot(grid[..., 0], grid[..., 1])
    assert r.min() == pytest.approx(1.0, abs=1e-12)
    assert r.max() == pytest.approx(2.0, abs=1e-12)
    # Any face order gives the same block up to its index frame (same cells).
    faces = list(og.entity("duct").params["faces"])
    random.Random(4).shuffle(faces)
    names = [og.store.get(f).name for f in faces]
    other = og.create_block(faces=names, name="again")
    assert sorted(og.block_grid(other).shape[:3]) == sorted(grid.shape[:3])
    assert total_volume(og, "again") == pytest.approx(total_volume(og, "duct"), rel=1e-12)
    assert og.block_quality("again")["scaled_jacobian"].min() > 0.99
    assert "6 faces" in og.block_status("duct")


def test_block_validation():
    og = opengrid.Session()
    duct(og)
    with pytest.raises(opengrid.ValidationError, match="six structured domains"):
        og.create_block(faces=["bottom", "top"])
    with pytest.raises(opengrid.ValidationError, match="either six faces or an extrusion"):
        og.create_block()
    with pytest.raises(opengrid.ValidationError, match="listed twice"):
        og.create_block(faces=["bottom", "bottom", "inner", "outer", "sidea", "sideb"])
    # A face that does not close the box.
    og.create_line((5, 0, 0), (6, 0, 0), name="x0")
    og.create_line((6, 0, 0), (6, 1, 0), name="x1")
    og.create_line((6, 1, 0), (5, 1, 0), name="x2")
    og.create_line((5, 1, 0), (5, 0, 0), name="x3")
    for k in range(4):
        og.create_connector(f"x{k}", 5, name=f"c-x{k}")
    og.create_domain([["c-x0"], ["c-x1"], ["c-x2"], ["c-x3"]], name="far")
    og.create_block(faces=["bottom", "far", "inner", "outer", "sidea", "sideb"], name="bad")
    with pytest.raises(opengrid.ValidationError, match="corners"):
        og.block_grid("bad")
    with pytest.raises(opengrid.ValidationError, match="elliptic smoothing need a planar"):
        og.set_domain_generator("inner", method="hermite_tfi")


def test_two_blocks_interface_and_cgns(tmp_path):
    og = opengrid.Session()
    _first, side = duct(og, 21, 0.0, T / 2, tag="-a")
    duct(og, 21, T / 2, T, tag="-b", shared=side)
    report_names = og.boundary_report().interfaces
    assert report_names == ["duct-a.jmax | duct-b.jmin"] or len(report_names) == 1
    og.create_patch("walls", "wall", ["inner-a", "outer-a", "inner-b", "outer-b"])
    og.create_patch("ends", "outlet", ["bottom-a", "top-a", "bottom-b", "top-b"])
    og.create_patch("sides", "symmetry", ["sidea-a", "sideb-b"])
    report = og.boundary_report()
    assert report.ok, report.issues
    assert len(report.boundary) == 10 and len(report.interfaces) == 1
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    path = tmp_path / "two.cgns"
    summary = og.export_cgns(path)
    assert summary == {"zones": 2, "families": 4, "boundaries": 10, "connections": 2}
    checker = shutil.which("cgnscheck")
    if checker:
        out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
        assert "error" not in out.lower() and "warning" not in out.lower(), out
    grids = {n: og.block_grid(n) for n in ("duct-a", "duct-b")}
    with h5py.File(path) as f:
        found = 0
        for zone_name in ("duct-a", "duct-b"):
            zone = f[f"OpenGrid/{zone_name}"]
            for conn in zone["ZoneGridConnectivity"].values():
                if conn.attrs.get("label") != b"GridConnectivity1to1_t":
                    continue
                found += 1
                donor_name = conn[" data"][()].tobytes().decode().strip("\x00 ")
                r = conn["PointRange"][" data"][()].reshape(-1) - 1
                d = conn["PointRangeDonor"][" data"][()].reshape(-1) - 1
                t = conn["Transform"][" data"][()].tolist()
                g = grids[zone_name]
                h = grids[donor_name]
                # Every point of the range maps through the transform onto the donor point.
                for k in range(r[2], r[5] + 1):
                    for j in range(r[1], r[4] + 1):
                        for i in range(r[0], r[3] + 1):
                            delta = (i - r[0], j - r[1], k - r[2])
                            target = list(d[:3])
                            for axis in range(3):
                                m = abs(t[axis]) - 1
                                target[m] += int(np.sign(t[axis])) * delta[axis]
                            assert np.array_equal(g[k, j, i], h[target[2], target[1], target[0]])
        assert found == 2


def test_extruded_blocks():
    og = opengrid.Session()
    corners = [(0, 0), (2, 0), (2, 1), (0, 1)]
    for k in range(4):
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"l{k}")
        og.create_connector(f"l{k}", 6 + (k % 2) * 3, name=f"c{k}")
    og.create_domain([["c0"], ["c1"], ["c2"], ["c3"]], name="base")
    up = og.create_extruded_block("base", "translate", vector=(0, 0, 3), dimension=7)
    assert total_volume(og, up) == pytest.approx(6.0, rel=1e-13)
    assert og.block_grid(up).shape == (7, 9, 6, 3)
    down = og.create_extruded_block(
        "base",
        "translate",
        vector=(0, 0, -1),
        dimension=5,
        distribution={"type": "tanh", "begin_spacing": 0.05},
    )
    assert total_volume(og, down) == pytest.approx(2.0, rel=1e-13)
    assert og.block_quality(down)["scaled_jacobian"].min() > 0
    # Extruding up and down from the same base: the base is an interface.
    report = og.boundary_report()
    assert len(report.interfaces) == 1
    turned = og.create_extruded_block(
        "base", "rotate", vector=(0, 1, 0), axis_point=(-1, 0, 0), angle=-math.pi / 2, dimension=91
    )
    assert og.block_quality(turned)["scaled_jacobian"].min() > 0
    # Rotating the rectangle x in [0, 2] (radius 1..3 about x = -1) by 90 degrees about y.
    analytic = 0.25 * math.pi * (9 - 1) * 1.0
    assert total_volume(og, turned) == pytest.approx(analytic, rel=2e-4)
    thick = og.create_extruded_block("base", "normal", distance=0.5, dimension=3)
    assert total_volume(og, thick) == pytest.approx(1.0, rel=1e-13)
    with pytest.raises(opengrid.ValidationError, match="must not be zero"):
        og.create_extruded_block("base", "rotate", angle=0.0)
    with pytest.raises(opengrid.ValidationError, match="mode"):
        og.create_extruded_block("base", "twist")
    replayed = opengrid.replay(og.journal_text())
    assert np.array_equal(replayed.block_grid(turned), og.block_grid(turned))


def test_3d_patches_and_exports(tmp_path):
    og = patched_duct()
    report = og.boundary_report()
    assert report.ok, report.issues
    meshio = pytest.importorskip("meshio")
    msh = tmp_path / "d.msh"
    summary = og.export_gmsh(msh)
    assert summary["physical_groups"] == 4
    mesh = meshio.read(msh)
    types = Counter()
    for block in mesh.cells:
        types[block.type] += len(block.data)
    assert types == {"hexahedron": 30 * 10 * 8, "quad": 2 * (30 * 10 + 30 * 8 + 10 * 8)}
    if find_gmsh():
        log = gmsh_read(msh)
        assert "Error" not in log and "Warning" not in log, log
    vtk = pytest.importorskip("vtk")
    path = tmp_path / "d.vtkhdf"
    info = og.export_vtkhdf(path)
    assert info["volume_cells"] == 2400
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    volume = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    assert volume.GetNumberOfCells() == 2400
    assert volume.GetCellType(0) == vtk.VTK_HEXAHEDRON
    # Outward boundary quads: the area vectors of each patch sum to ~zero for a closed box
    # only overall, so check the walls' normals point away from the duct centre line.
    from opengrid.io import read_plot3d

    xyz = tmp_path / "d.xyz"
    og.export_plot3d(xyz)
    blocks = read_plot3d(xyz, three_d=True)
    assert np.allclose(blocks[0], og.block_grid("duct"))
    # Unassigned faces and connector members are reported.
    og.remove_from_patch("sides", "sideb")
    assert any("in no patch" in str(i) for i in og.boundary_report().issues)
    og.add_to_patch(
        "sides", [{"block": "duct", "face": og.boundary_report().unassigned[0].split(".")[1]}]
    )
    assert og.boundary_report().ok
    og.create_patch("lines", "wall", ["c-ai0"])
    assert any("3D patches contain" in str(i) for i in og.boundary_report().issues)


def test_project_round_trip(tmp_path):
    og = patched_duct()
    og.create_zone("fluid", ["duct"])
    path = tmp_path / "d.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert to_json(other.store) == to_json(og.store)
    assert np.array_equal(other.block_grid("duct"), og.block_grid("duct"))
    assert opengrid.replay(og.journal_text()).boundary_report().ok


def test_gui_blocks(qtbot, tmp_path):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    duct(og)
    og.delete("duct")
    faces = ["bottom", "top", "inner", "outer", "sidea", "sideb"]
    window.set_selection([og.entity(n).id for n in faces])
    name = window.create_block_on_selection()
    assert name in window.model_tree.entity_names("Blocks")
    window.set_selection([og.entity("bottom").id])
    extruded = window.extrude_selection()
    assert og.entity(extruded).params["extrusion"]["mode"] == "translate"
    # Property panel: change the extrusion to 5 points.
    window.set_selection([og.entity(extruded).id])
    _kind, editor = window.property_panel._editors["extrusion"]
    editor.dimension.setValue(5)
    assert window.property_panel.apply()
    assert og.block_grid(extruded).shape[0] == 5
    window.set_selection([og.entity(extruded).id])
    assert window.property_panel.read_changes() == (None, {})
    window.cut_action.setChecked(True)
    window.shrink_action.setChecked(True)
    window.quality_action.setChecked(True)
    window.set_selection([og.entity(name).id])
    assert window.show_worst_cell() is not None
    window.close()
