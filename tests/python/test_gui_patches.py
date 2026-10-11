import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.boundary_dock import INTERFACE_COLOR, PATCH_PALETTE, UNASSIGNED_COLOR
from opengrid.app.main_window import MainWindow

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    spec = importlib.util.spec_from_file_location("step", EXAMPLES / "03_backward_step.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.build(win.session)
    yield win
    win.close()


def _ids(window, *names):
    return [window.session.entity(n).id for n in names]


def test_patch_workflow_and_validation(window, tmp_path):
    report = window.validate_boundary()
    assert len(report.unassigned) == 8
    assert "8 unassigned" in window.boundary_dock.summary.text()
    window.set_selection(_ids(window, "c-in-bottom", "c-in-top", "c-step", "c-lo-bottom"))
    assert window.create_patch_from_selection("walls", "wall") == "walls"
    window.set_selection(_ids(window, "c-up-top"))
    assert window.add_selection_to_patch("walls")
    window.set_selection(_ids(window, "c-inlet"))
    window.create_patch_from_selection("inflow", "inlet")
    window.set_selection(_ids(window, "c-out-lo", "c-out-up"))
    window.create_patch_from_selection("outflow", "outlet")
    report = window.validate_boundary()
    assert report.ok
    assert "OK" in window.boundary_dock.summary.text()
    assert window.model_tree.entity_names("Patches") == ["walls", "inflow", "outflow"]
    window.set_selection(_ids(window, "lower-block", "upper-block"))
    assert window.create_zone_from_selection("expansion") == "expansion"
    assert window.model_tree.entity_names("Zones") == ["expansion"]
    assert window.export_cgns(str(tmp_path / "s.cgns"))["zones"] == 3
    assert window.export_gmsh(str(tmp_path / "s.msh"))["physical_groups"] == 5
    assert window.export_vtkhdf(str(tmp_path / "s.vtkhdf"))["boundary_cells"] > 0


def test_invalid_export_is_reported(window, tmp_path):
    assert window.export_cgns(str(tmp_path / "x.cgns")) is None
    assert "belongs to no patch" in window.log_panel.toPlainText()


def test_patch_colouring(window):
    window.set_selection(_ids(window, "c-step"))
    window.create_patch_from_selection("walls", "wall")
    window.patch_color_action.setChecked(True)
    colors = window.scene.patch_colors
    assert colors[window.session.entity("c-step").id] == PATCH_PALETTE[0]
    assert colors[window.session.entity("c-mid").id] == INTERFACE_COLOR
    assert colors[window.session.entity("c-inlet").id] == UNASSIGNED_COLOR
    window.patch_color_action.setChecked(False)
    assert window.scene.patch_colors is None


def test_property_panel_edits_patch(window):
    window.session.create_patch("walls", "wall", ["c-step"])
    window.set_selection(_ids(window, "walls"))
    panel = window.property_panel
    _kind, combo = panel._editors["patch_type"]
    combo.setCurrentText("symmetry")
    _kind, members = panel._editors["members"]
    members.setPlainText("c-step\nc-lo-bottom\n")
    assert panel.apply()
    entity = window.session.entity("walls")
    assert entity.params["patch_type"] == "symmetry"
    assert len(entity.params["members"]) == 2
    window.set_selection(_ids(window, "walls"))
    assert window.property_panel.read_changes() == (None, {})
    members_editor = window.property_panel._editors["members"][1]
    members_editor.setPlainText("nope")
    assert not window.property_panel.apply()
    assert "no entity 'nope'" in window.property_panel._error.text()
