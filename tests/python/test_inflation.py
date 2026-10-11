"""Phase 6: boundary-layer inflation in unstructured domains."""

import shutil
import subprocess
from collections import Counter

import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read

import opengrid
from opengrid import _core
from opengrid.project import to_json

INFLATION = {"connectors": ["c-wall"], "first_height": 1e-4, "growth_rate": 1.2}


def cylinder(og=None, **generator) -> opengrid.Session:
    if og is None:
        og = opengrid.Session()
    og.create_circle((0, 0), 0.5, name="wall")
    og.create_circle((0, 0), 10.0, name="far")
    og.create_connector("wall", 201, name="c-wall")
    og.create_connector("far", 81, name="c-far")
    og.create_unstructured_domain(connectors=["c-wall", "c-far"], generator=generator, name="flow")
    og.create_patch("p-wall", "wall", ["c-wall"])
    og.create_patch("p-far", "farfield", ["c-far"])
    return og


def test_layers_first_height_and_validity():
    og = cylinder(inflation=INFLATION)
    mesh = og.unstructured_mesh("flow")
    assert mesh["layers"] > 10
    assert mesh["warnings"] == []
    q = og.domain_quality("flow")
    assert q["scaled_jacobian"].min() > 0
    plain = cylinder().domain_quality("flow")
    assert q["area"].sum() == pytest.approx(plain["area"].sum(), rel=1e-12)
    # First layer: 200 nodes at 1e-4 from the wall polygon (curved wall: within 2 %).
    wall = og.connector_nodes("c-wall").points[:, :2]
    first = mesh["points"][280:480, :2]
    a = wall[:-1]
    b = wall[1:]
    ab = b - a
    distance = np.full(len(first), np.inf)
    for k in range(len(a)):
        t = np.clip(((first - a[k]) @ ab[k]) / (ab[k] @ ab[k]), 0.0, 1.0)
        distance = np.minimum(distance, np.linalg.norm(first - (a[k] + t[:, None] * ab[k]), axis=1))
    assert np.allclose(distance / 1e-4, 1.0, atol=0.02)
    assert "inflation layers" in og.domain_status("flow")


def test_validation_and_round_trip(tmp_path):
    og = cylinder()
    with pytest.raises(opengrid.ValidationError, match="first_height"):
        og.set_unstructured_generator("flow", inflation={"connectors": ["c-wall"]})
    og.create_line((20, 0), (21, 0), name="elsewhere")
    og.create_connector("elsewhere", 3, name="c-elsewhere")
    with pytest.raises(opengrid.ValidationError, match="not part of the domain's loops"):
        og.set_unstructured_generator(
            "flow", inflation={"connectors": ["c-elsewhere"], "first_height": 1e-4}
        )
    with pytest.raises(opengrid.ValidationError, match="layers"):
        og.set_unstructured_generator("flow", inflation={**INFLATION, "layers": 0})
    og.set_unstructured_generator("flow", inflation=INFLATION)
    stored = og.entity("flow").params["generator"]["inflation"]
    assert stored["connectors"] == (og.entity("c-wall").id,)
    path = tmp_path / "i.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert to_json(other.store) == to_json(og.store)
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)
    # Renaming the wall connector keeps the reference (ids are stored).
    og.rename("c-wall", "c-body")
    assert og.unstructured_mesh("flow")["layers"] > 10
    og.undo()
    og.undo()
    assert og.entity("flow").params["generator"]["inflation"] is None


@pytest.mark.parametrize("element_type", ["triangles", "quads"])
def test_exports_conform(tmp_path, element_type):
    og = cylinder(inflation=INFLATION, element_type=element_type)
    mesh = og.unstructured_mesh("flow")
    if element_type == "quads":
        assert len(mesh["triangles"]) == 0
    meshio = pytest.importorskip("meshio")
    path = tmp_path / "c.msh"
    og.export_gmsh(path)
    read = meshio.read(path)
    edges = Counter()
    for block in read.cells:
        if block.type == "line":
            continue
        for cell in block.data:
            for k in range(len(cell)):
                edges[tuple(sorted((int(cell[k]), int(cell[(k + 1) % len(cell)]))))] += 1
    assert max(edges.values()) == 2
    assert sum(1 for c in edges.values() if c == 1) == 200 + 80  # wall and farfield only
    if find_gmsh():
        log = gmsh_read(path)
        assert "Error" not in log and "Warning" not in log, log
    if _core.cgns_available():
        cgns = tmp_path / "c.cgns"
        og.export_cgns(cgns)
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run(
                [checker, str(cgns)], capture_output=True, text=True, check=False
            ).stdout
            assert "error" not in out.lower() and "warning" not in out.lower(), out


def test_gui_inflation_fields(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    session = cylinder(window.session)
    window.set_selection([session.entity("flow").id])
    panel = window.property_panel
    _kind, editor = panel._editors["generator"]
    assert editor.walls.text() == ""
    editor.walls.setText("c-wall")
    editor.first_height.setText("1e-4")
    assert panel.apply()
    stored = session.entity("flow").params["generator"]["inflation"]
    assert stored["connectors"] == (session.entity("c-wall").id,)
    assert session.unstructured_mesh("flow")["layers"] > 10
    window.set_selection([session.entity("flow").id])
    assert window.property_panel.read_changes() == (None, {})
    _kind, editor = window.property_panel._editors["generator"]
    editor.first_height.setText("")
    assert not window.property_panel.apply()
    assert "first height" in window.property_panel._error.text()
    window.close()
