import importlib.util
from pathlib import Path

import pytest

import opengrid

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.boundary_dock import INTERFACE_COLOR, PATCH_PALETTE, UNASSIGNED_COLOR
from opengrid.app.main_window import MainWindow
from opengrid.validation import ValidationError

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
    assert "8 of 8 boundaries have no type" in window.boundary_dock.status.text()
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
    assert "ready to export" in window.boundary_dock.status.text()
    # No Patches category: each curve shows its boundary type.
    curves = window.model_tree.category_item("Curves")
    rows = {curves.child(i).text(0): curves.child(i).text(1) for i in range(curves.childCount())}
    assert rows["c-inlet"].endswith("\u00b7 inlet") and rows["c-step"].endswith("\u00b7 wall")
    assert "\u00b7" not in rows["c-mid"]  # an interface: no type
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


def test_boundary_type_from_the_selection_and_the_menus(window):
    og = window.session
    window.set_selection(_ids(window, "c-in-bottom", "c-in-top"))
    assert window.set_boundary_type_on_selection("wall") == "wall"
    assert window.patch_color_action.isChecked()  # the colouring shows the result
    walls = og.entity("wall").params
    assert walls["patch_type"] == "wall" and len(walls["members"]) == 2
    window.set_selection(_ids(window, "c-step"))
    assert window.set_boundary_type_on_selection("wall") == "wall"  # joins the existing one
    assert len(og.entity("wall").params["members"]) == 3
    window.set_selection(_ids(window, "c-inlet"))
    inlet = window.set_boundary_type_on_selection("inlet")
    assert inlet == "inlet-2"  # "inlet" is a curve of this model: names are unique
    assert og.entity(inlet).params["patch_type"] == "inlet"
    # A member moves between patches: one patch per boundary.
    window.set_selection(_ids(window, "c-step"))
    assert window.set_boundary_type_on_selection("inlet") == inlet
    assert len(og.entity("wall").params["members"]) == 2
    assert len(og.entity(inlet).params["members"]) == 2
    # Named patches of the type are used or made; the type is checked.
    assert og.set_boundary_type(["c-up-top"], "outlet", name="exit") == "exit"
    assert og.entity("exit").params["patch_type"] == "outlet"
    with pytest.raises(ValidationError, match="not wall"):
        og.set_boundary_type(["c-up-top"], "wall", name="exit")
    assert og.set_boundary_type(["c-up-top"], "wall") == "wall"
    assert og.entity("exit").params["members"] == ()
    # The menus carry it: the context menu (with a selection) and the Create menu.
    window.set_selection(_ids(window, "c-lo-bottom"))
    menu = window.context_menu()
    names = [a.text() for a in menu.actions() if a.menu() is not None]
    assert "Boundary type" in names
    submenu = next(a.menu() for a in menu.actions() if a.text() == "Boundary type")
    assert [a.text() for a in submenu.actions()][:3] == ["wall", "inlet", "outlet"]
    submenu.actions()[3].trigger()  # symmetry
    assert og.entity("symmetry").params["patch_type"] == "symmetry"
    window.set_selection([])
    assert "Boundary type" not in [a.text() for a in window.context_menu().actions()]
    assert window.set_boundary_type_on_selection("wall") is None  # a hint, no patch
    from opengrid.project.ogp import to_json

    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_pair_as_periodic_and_extrude_model_from_the_gui(window):
    og = window.session
    window.set_selection(_ids(window, "c-inlet"))
    assert window.pair_selection_as_periodic() is None  # one side: a hint
    window.set_selection(_ids(window, "c-in-top", "c-lo-bottom"))
    names = window.pair_selection_as_periodic()
    assert names == ["periodic-1a", "periodic-1b"]
    assert window.selected_names() == ["periodic-1a"]  # the panel shows the pairing
    assert window.patch_color_action.isChecked()
    assert "periodic" in window.property_panel._editors
    assert "Pair as periodic" not in [a.text() for a in window.context_menu().actions()]
    window.set_selection(_ids(window, "periodic-1a", "periodic-1b"))
    assert "Pair as periodic" in [a.text() for a in window.context_menu().actions()]
    og.undo()
    # The whole 2D model, one cell thick, with periodic front and back.
    for name in og.unassigned_boundary_connectors():
        og.set_boundary_type([name], "wall")
    names = window.extrude_model(thickness=0.05, ends="periodic")
    assert names and names[-2:] == ["front", "back"]
    assert og.is_3d() and window.navigation_mode == "3d"
    assert og.boundary_report().ok
    assert window.extrude_model(thickness=0.05, ends="periodic") is None  # 3D already
