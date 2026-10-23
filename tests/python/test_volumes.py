"""Phase 8a: unstructured domains in 3D, closed shells and tetrahedral volumes (Netgen)."""

import itertools
import math
import shutil
import subprocess
from collections import Counter

import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.model.volumes import total_volume
from opengrid.project import to_json

pytest.importorskip("netgen.meshing")


def box(og, lo, hi, n, tag, kind, shared=None):
    """Six faces of an axis-aligned box (12 connectors, 6 domains); returns the face names.
    ``shared``: faces to reuse instead of creating them, by (axis, side)."""
    corners = {}
    for c in itertools.product((0, 1), repeat=3):
        corners[c] = tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3))
    edges = {}

    def connector(c, d):
        key = (corners[c], corners[d])
        if key not in edges:
            name = f"{tag}-e{len(edges)}"
            existing = [
                n
                for n in og.names("line")
                if og.entity(n).params["p0"] == corners[c]
                and og.entity(n).params["p1"] == corners[d]
            ]
            if existing:
                edges[key] = "c-" + existing[0]
            else:
                og.create_line(corners[c], corners[d], name=name)
                edges[key] = og.create_connector(name, n, name=f"c-{name}")
        return edges[key]

    faces = []
    for axis in range(3):
        for side in (0, 1):
            if shared is not None and (axis, side) in shared:
                faces.append(shared[(axis, side)])
                continue
            a, b = [x for x in range(3) if x != axis]

            def v(i, j, axis=axis, side=side, a=a, b=b):
                c = [0, 0, 0]
                c[axis] = side
                c[a] = i
                c[b] = j
                return tuple(c)

            names = [
                connector(v(0, 0), v(1, 0)),
                connector(v(1, 0), v(1, 1)),
                connector(v(0, 1), v(1, 1)),
                connector(v(0, 0), v(0, 1)),
            ]
            face = f"{tag}-f{axis}{side}"
            if kind == "structured":
                og.create_domain(
                    [[names[0]], [names[1]], [(names[2], True)], [(names[3], True)]], name=face
                )
            else:
                og.create_unstructured_domain(connectors=names, name=face)
            faces.append(face)
    return faces


def box_with_body():
    og = opengrid.Session()
    outer = box(og, (0, 0, 0), (4, 3, 3), 13, "o", "structured")
    inner = box(og, (1.5, 1, 1), (2.5, 2, 2), 9, "i", "unstructured")
    og.create_volume(outer + inner, generator={"max_size": 0.4}, name="fluid")
    og.create_patch("walls", "farfield", outer)
    og.create_patch("body", "wall", inner)
    return og


def test_unstructured_domains_in_3d():
    og = opengrid.Session()
    pts = [(1, 0, 0), (1, 2, 0), (1, 2, 2), (1, 0, 2)]
    for k in range(4):
        og.create_line(pts[k], pts[(k + 1) % 4], name=f"e{k}")
        og.create_connector(f"e{k}", 21, name=f"c{k}")
    og.create_circle((1, 1, 1), 0.4, normal=(1, 0, 0), name="hole")
    og.create_connector("hole", 41, name="ch")
    og.create_unstructured_domain(connectors=["c0", "c1", "c2", "c3", "ch"], name="plane")
    mesh = og.unstructured_mesh("plane")
    assert np.ptp(mesh["points"][:, 0]) == 0.0
    q = og.domain_quality("plane")
    polygon = 0.5 * 40 * 0.16 * math.sin(2 * math.pi / 40)
    assert q["area"].sum() == pytest.approx(4.0 - polygon, rel=1e-12)
    assert q["min_angle"].min() > 20.0
    # On a cylinder.
    og.create_arc((0, 0, 0), 2.0, 0.0, math.pi / 2, name="a0")
    og.create_arc((0, 0, 3), 2.0, 0.0, math.pi / 2, name="a1")
    og.create_line((2, 0, 0), (2, 0, 3), name="l0")
    og.create_line((0, 2, 0), (0, 2, 3), name="l1")
    for name, n in (("a0", 31), ("a1", 31), ("l0", 21), ("l1", 21)):
        og.create_connector(name, n, name=f"c-{name}")
    og.create_extruded_surface("a0", (0, 0, 3), name="cyl")
    og.create_unstructured_domain(
        connectors=["c-a0", "c-a1", "c-l0", "c-l1"], surface="cyl", name="patch"
    )
    mesh = og.unstructured_mesh("patch")
    assert np.allclose(np.hypot(mesh["points"][:, 0], mesh["points"][:, 1]), 2.0, atol=1e-12)
    assert og.domain_quality("patch")["area"].sum() == pytest.approx(3 * math.pi, rel=2e-3)
    # A non-planar boundary without a surface is rejected.
    og.create_line((0, 2, 3), (2, 0, 0), name="skew")
    og.create_connector("skew", 5, name="c-skew")
    with pytest.raises(opengrid.ValidationError, match="not planar"):
        og.create_unstructured_domain(connectors=["c-a0", "c-l1", "c-skew"], name="bad")
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)


def test_volume_with_inner_body():
    og = box_with_body()
    mesh = og.volume_mesh("fluid")
    assert mesh["shells"] == 2
    assert total_volume(mesh) == pytest.approx(36.0 - 1.0, rel=1e-12)
    q = og.volume_quality("fluid")
    assert q["volume"].min() > 0
    assert q["scaled_jacobian"].min() > 0.05
    # Boundary nodes are the face nodes, bitwise.
    face_nodes = np.concatenate([og.domain_grid(f).reshape(-1, 3) for f in og.names("domain")])
    surface = mesh["points"][: mesh["n_surface"]]
    assert {p.tobytes() for p in np.unique(face_nodes, axis=0)} <= {p.tobytes() for p in surface}
    # Conformity: every tet face appears twice, except the boundary faces (once).
    faces = Counter()
    for t in mesh["tets"]:
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            faces[tuple(sorted(int(t[k]) for k in f))] += 1
    boundary = {tuple(sorted(int(v) for v in c)) for cells in mesh["faces"].values() for c in cells}
    assert {f for f, n in faces.items() if n == 1} == boundary
    assert max(faces.values()) == 2
    report = og.boundary_report()
    assert report.ok, report.issues
    assert len(report.boundary) == 12
    assert "tetrahedra" in og.volume_status("fluid")


def test_shell_errors():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (1, 1, 1), 5, "b", "structured")
    og.create_volume(faces[:-1], name="open")
    with pytest.raises(opengrid.ValidationError, match="closed shells"):
        og.volume_mesh("open")
    with pytest.raises(opengrid.ValidationError, match="at least two"):
        og.create_volume(faces[:1])
    with pytest.raises(opengrid.ValidationError, match="pyramids"):
        og.create_volume(faces, generator={"pyramids": "maybe"})


def test_pyramids_always_and_never():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (2, 1, 1), 6, "b", "structured")
    og.create_volume(faces, generator={"max_size": 0.3, "pyramids": "always"}, name="p")
    mesh = og.volume_mesh("p")
    quads = 6 * 5 * 5  # 6 points on every edge
    assert len(mesh["pyramids"]) == quads
    assert total_volume(mesh) == pytest.approx(2.0, rel=1e-12)
    assert og.volume_quality("p")["volume"].min() > 0
    og.set_volume_generator("p", pyramids="never")
    mesh = og.volume_mesh("p")
    assert len(mesh["pyramids"]) == 0
    assert total_volume(mesh) == pytest.approx(2.0, rel=1e-12)


def test_two_volumes_share_a_face(tmp_path):
    og = opengrid.Session()
    left = box(og, (0, 0, 0), (1, 1, 1), 7, "l", "unstructured")
    right = box(og, (1, 0, 0), (2, 1, 1), 7, "r", "unstructured", shared={(0, 0): left[1]})
    og.create_volume(left, generator={"max_size": 0.25}, name="left")
    og.create_volume(right, generator={"max_size": 0.25}, name="right")
    outer = [f for f in left + right if f != left[1]]
    og.create_patch("walls", "wall", outer)
    report = og.boundary_report()
    assert report.ok, report.issues
    assert len(report.interfaces) == 1
    a = og.volume_mesh("left")
    b = og.volume_mesh("right")
    shared = og.entity(left[1]).id
    na = {a["points"][i].tobytes() for c in a["faces"][shared] for i in c}
    nb = {b["points"][i].tobytes() for c in b["faces"][shared] for i in c}
    assert na == nb
    assert total_volume(a) + total_volume(b) == pytest.approx(2.0, rel=1e-12)
    if _core.cgns_available():
        path = tmp_path / "two.cgns"
        summary = og.export_cgns(path)
        assert summary["zones"] == 2 and summary["connections"] == 2
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
            assert "error" not in out.lower() and "warning" not in out.lower(), out


def test_volume_exports(tmp_path):
    og = box_with_body()
    meshio = pytest.importorskip("meshio")
    msh = tmp_path / "v.msh"
    summary = og.export_gmsh(msh)
    mesh = meshio.read(msh)
    counts = Counter()
    for block in mesh.cells:
        counts[block.type] += len(block.data)
    volume = og.volume_mesh("fluid")
    assert counts["tetra"] == len(volume["tets"])
    assert counts["triangle"] == sum(len(c) for c in volume["faces"].values())
    assert summary["physical_groups"] == 3
    vtk = pytest.importorskip("vtk")
    path = tmp_path / "v.vtkhdf"
    og.export_vtkhdf(path)
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    cells = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    assert cells.GetCellType(0) == vtk.VTK_TETRA
    with pytest.raises(opengrid.ValidationError, match="structured blocks only"):
        og.export_plot3d(tmp_path / "v.xyz")
    og.remove_from_patch("body", "i-f00")
    assert any("in no patch" in str(i) for i in og.boundary_report().issues)


def test_gui_volume(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    faces = box(og, (0, 0, 0), (1, 1, 1), 5, "b", "structured")
    window.set_selection([og.entity(f).id for f in faces])
    name = window.create_volume_on_selection()
    assert name in window.model_tree.entity_names("Blocks")
    window.set_selection([og.entity(name).id])
    _kind, editor = window.property_panel._editors["generator"]
    editor.max_size.setText("0.2")
    assert window.property_panel.apply()
    assert og.entity(name).params["generator"]["max_size"] == 0.2
    window.set_selection([og.entity(name).id])
    assert window.property_panel.read_changes() == (None, {})
    window.cut_action.setChecked(True)
    window.quality_action.setChecked(True)
    assert window.show_worst_cell() is not None
    window.close()


def disc(og, n=61, generator=None):
    og.create_circle((0, 0), 1.0, name="c")
    og.create_connector("c", n, name="cc")
    og.create_unstructured_domain(connectors=["cc"], generator=generator, name="disc")
    mesh = og.unstructured_mesh("disc")
    points = np.asarray(mesh["points"])
    area = 0.0
    for cells in (mesh["triangles"], mesh["quads"]):
        for cell in cells:
            p = points[list(cell)]
            for k in range(1, len(cell) - 1):
                area += 0.5 * abs(np.cross(p[k] - p[0], p[k + 1] - p[0])[2])
    return area


def test_extruded_volume_prisms_and_hexahedra_with_patched_ends_and_sides(tmp_path):
    og = opengrid.Session()
    area = disc(og, generator={"element_type": "quad_dominant"})
    mesh2d = og.unstructured_mesh("disc")
    name = og.create_extruded_volume("disc", "translate", vector=(0, 0, 0.5), dimension=6)
    mesh = og.volume_mesh(name)
    assert len(mesh["prisms"]) == 5 * len(mesh2d["triangles"])
    assert len(mesh["hexes"]) == 5 * len(mesh2d["quads"]) > 0
    assert len(mesh["tets"]) == 0 and len(mesh["points"]) == 6 * len(mesh2d["points"])
    assert abs(total_volume(mesh) - 0.5 * area) < 1e-12 * area
    quality = og.volume_quality(name)
    assert quality["scaled_jacobian"].min() > 0.3 and quality["volume"].min() > 0
    assert "translate extrusion of disc, 5 layers" in og.volume_status(name)
    faces = mesh["faces"]
    domain = og.entity("disc").id
    connector = og.entity("cc").id
    assert set(faces) == {domain, ("top", domain), ("side", connector)}
    assert len(faces[("side", connector)]) == 5 * 60
    # The side strip stands on the connector's nodes, in order.
    nodes = og.connector_nodes("cc").points
    quad = faces[("side", connector)][0]
    assert np.allclose(mesh["points"][quad[0]], nodes[0])
    assert np.allclose(mesh["points"][quad[1]], nodes[1])
    # Outward faces: bottom down, top up, sides away from the axis.
    p = mesh["points"]
    bottom = faces[domain][0]
    assert np.cross(p[bottom[1]] - p[bottom[0]], p[bottom[2]] - p[bottom[0]])[2] < 0
    top = faces[("top", domain)][0]
    assert np.cross(p[top[1]] - p[top[0]], p[top[2]] - p[top[0]])[2] > 0
    outward = np.cross(p[quad[1]] - p[quad[0]], p[quad[3]] - p[quad[0]])
    assert outward @ p[quad[0]] > 0
    # Patches: the domain (bottom), the volume's top and the connector (side wall).
    og.create_patch("back", "symmetry", ["disc"])
    og.create_patch("front", "symmetry", [{"volume": name, "face": "top"}])
    og.create_patch("wall", "wall", ["cc"])
    assert og.boundary_report().issues == []
    summary = og.export_cgns(tmp_path / "slab.cgns")
    # One section per patch and face type: bottom and top (triangles + quads), side (quads).
    assert summary["boundaries"] == 5 and summary["zones"] == 1
    meshio = pytest.importorskip("meshio")
    og.export_gmsh(tmp_path / "slab.msh")
    counts = Counter()
    for block in meshio.read(tmp_path / "slab.msh").cells:
        counts[block.type] += len(block.data)
    assert counts["wedge"] == len(mesh["prisms"]) and counts["hexahedron"] == len(mesh["hexes"])
    assert counts["quad"] == 5 * 60 + 2 * len(mesh2d["quads"])
    assert counts["triangle"] == 2 * len(mesh2d["triangles"])
    vtk = pytest.importorskip("vtk")
    og.export_vtkhdf(tmp_path / "slab.vtkhdf")
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(tmp_path / "slab.vtkhdf"))
    reader.Update()
    cells = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    assert cells.GetNumberOfCells() == len(mesh["prisms"]) + len(mesh["hexes"])
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.save(tmp_path / "slab.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "slab.ogp")
    assert to_json(other.store) == to_json(og.store)
    # One less undo step: the volume goes, its patches first.
    og.delete(["front"])
    assert og.boundary_report().issues != []


def test_extruded_volume_modes_mirroring_and_refusals():
    og = opengrid.Session()
    area = disc(og)
    name = og.create_extruded_volume("disc", "translate", vector=(0, 0, -0.5), dimension=4)
    mesh = og.volume_mesh(name)
    assert mesh["extrusion"]["mirrored"] and "base mirrored" in og.volume_status(name)
    assert (
        abs(total_volume(mesh) - 0.5 * area) < 1e-12 and og.volume_quality(name)["volume"].min() > 0
    )
    # Rotation: a quarter turn about the x axis two units away; Pappus with straight sides.
    extrusion = dict(og.entity(name).params["extrusion"])
    extrusion.update(mode="rotate", vector=(1, 0, 0), axis_point=(0, 2, 0), angle=math.pi / 2)
    og.modify(name, extrusion=extrusion)
    step = math.pi / 2 / 3
    expected = area * 2.0 * (math.pi / 2) * math.sin(step) / step
    assert abs(total_volume(og.volume_mesh(name)) - expected) < 1e-9 * expected
    assert og.volume_quality(name)["volume"].min() > 0
    extrusion.update(mode="normal", distance=0.3)
    og.modify(name, extrusion=extrusion)
    assert abs(total_volume(og.volume_mesh(name)) - 0.3 * area) < 1e-12
    # An axis through the domain folds it: reported, not hidden.
    extrusion.update(mode="rotate", axis_point=(0, 0, 0), angle=math.pi / 2)
    og.modify(name, extrusion=extrusion)
    assert "flat or inverted" in og.volume_status(name)
    with pytest.raises(opengrid.ValidationError, match="translate, rotate, normal"):
        og.create_extruded_volume("disc", "hyperbolic", distance=0.1)
    og.create_spline([(0, 0), (1, 0.1), (2, 0)], name="b")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_line((2, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    for curve in ("b", "r", "t", "l"):
        og.create_connector(curve, 11, name=f"c{curve}")
    og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="d")
    with pytest.raises(opengrid.ValidationError, match="not an unstructured domain"):
        og.create_extruded_volume("d", "translate", vector=(0, 0, 1))
    with pytest.raises(opengrid.ValidationError, match=r"faces .* or an extrusion"):
        og.create_volume(faces=["disc", "d"], extrusion={"domain": "disc"})
    with pytest.raises(opengrid.ValidationError, match="no inflation"):
        og.create_extruded_volume(
            "disc", vector=(0, 0, 1), generator={"inflation": {"faces": [], "first_height": 0.01}}
        )


def test_gui_extrudes_an_unstructured_domain_to_a_volume(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    disc(og, n=31)
    window.set_selection([og.entity("disc").id])
    window.dimension_box.setValue(4)
    name = window.extrude_selection()
    assert name is not None and og.entity(name).kind == "volume"
    mesh = og.volume_mesh(name)
    assert mesh["extrusion"]["layers"] == 3 and len(mesh["prisms"]) > 0
    assert abs(og.entity(name).params["extrusion"]["vector"][2]) > 0  # along the normal
    window.set_selection([og.entity(name).id])
    _kind, editor = window.property_panel._editors["extrusion"]
    editor.dimension.setValue(6)
    assert window.property_panel.apply()
    assert og.volume_mesh(name)["extrusion"]["layers"] == 5
    # Patches from the selection: the volume gives both ends, the connector the side wall.
    window.set_selection([og.entity(name).id])
    assert window.create_patch_from_selection("ends", "symmetry") == "ends"
    window.set_selection([og.entity("cc").id])
    assert window.create_patch_from_selection("rim", "wall") == "rim"
    assert og.boundary_report().issues == []
    window.close()


def test_one_volume_meshed_from_two_threads_at_once():
    """The window meshes in a worker thread while its main thread may ask for the same volume
    (the boundary check): Netgen's one task manager per process deadlocked when both ran it.
    The second request waits for the first and takes its result."""
    import threading

    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (1, 1, 1), 5, "b", "unstructured_domain")
    og.create_volume(faces, generator={"max_size": 0.2}, name="v")
    results = {}
    worker = threading.Thread(target=lambda: results.update(worker=og.volume_mesh("v")))
    worker.start()
    results["main"] = og.volume_mesh("v")
    worker.join(timeout=120)
    assert not worker.is_alive()
    assert results["worker"] is results["main"]


def test_the_boundary_check_waits_for_a_long_operation(qtbot):
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    calls = []
    window.validate_boundary = lambda *a, **k: calls.append(1)
    window._long_running = True
    window._refresh_boundary()
    assert calls == [] and window._boundary_timer.isActive()
    window._long_running = False
    window._refresh_boundary()
    assert calls == [1]
    window.close()
