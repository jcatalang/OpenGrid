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
