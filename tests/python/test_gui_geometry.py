"""GUI behaviour driven through world-space handlers (no real mouse events needed)."""

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from helpers import naca0012_selig, write_selig

from opengrid.app.main_window import MainWindow
from opengrid.app.scene import pick_entity
from opengrid.app.tools import snap, snap_targets


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    win.resize(1000, 800)
    yield win
    win.close()


def _wpp(window, value=0.01):
    # Pin the world-per-pixel scale so snapping/picking tolerances are deterministic.
    window.scene.world_per_pixel = lambda: value


def test_line_tool_with_snapping(window):
    _wpp(window)
    window.session.create_point((1.0, 1.0), name="anchor")
    window.set_tool("line")
    window.handle_click((0.0, 0.0, 0.0))
    assert window.tool.points == [(0.0, 0.0, 0.0)]
    # Within 10 px (0.1 world units) of the point: snaps exactly onto it.
    window.handle_click((1.05, 0.97, 0.0))
    line = window.session.entity("line-1")
    assert line.params["p1"] == (1.0, 1.0, 0.0)
    assert window.tool.points == []
    assert window.model_tree.entity_names() == ["anchor", "line-1"]


def test_snap_to_curve_end(window):
    _wpp(window)
    window.session.create_line((0, 0), (2, 0), name="l")
    targets = snap_targets(window.session)
    assert [t.label for t in targets] == ["l start", "l end"]
    snapped, target = snap((2.03, 0.04, 0.0), targets, 0.1)
    assert snapped == (2.0, 0.0, 0.0)
    assert target.label == "l end"
    snapped, target = snap((2.5, 0.0, 0.0), targets, 0.1)
    assert target is None
    assert snapped == (2.5, 0.0, 0.0)


def test_arc_polyline_spline_tools(window):
    _wpp(window, 1e-4)
    window.set_tool("arc")
    for xyz in [(1, 0, 0), (0, 1, 0), (-1, 0, 0)]:
        window.handle_click(xyz)
    assert window.session.entity("arc-1").params["radius"] == pytest.approx(1.0)

    window.set_tool("polyline")
    for xyz in [(0, -1, 0), (1, -2, 0), (2, -1, 0)]:
        window.handle_click(xyz)
    assert "polyline-1" not in window.session.store
    assert window.finish_action.isEnabled()
    window.finish_tool()
    assert len(window.session.entity("polyline-1").params["points"]) == 3

    window.set_tool("spline")
    window.handle_click((0, 3, 0))
    window.handle_move((1, 4, 0))
    window.handle_click((1, 4, 0))
    window.handle_move((2, 3.5, 0))  # spline preview with 3 points
    window.handle_click((2, 3, 0))
    window.cancel_tool()  # Esc drops the clicks
    assert "spline-1" not in window.session.store
    for xyz in [(0, 3, 0), (1, 4, 0), (2, 3, 0)]:
        window.handle_click(xyz)
    window.finish_tool()
    assert window.session.entity("spline-1").kind == "spline"


def test_invalid_tool_input_is_reported_not_raised(window):
    _wpp(window, 1e-6)
    window.set_tool("arc")
    for xyz in [(0, 0, 0), (1, 1, 0), (2, 2, 0)]:  # collinear
        window.handle_click(xyz)
    assert len(window.session.store) == 0
    assert "collinear" in window.log_panel.toPlainText()


def test_pick_select_and_delete(window):
    _wpp(window)
    window.session.create_line((0, 0), (2, 0), name="bottom")
    window.session.create_circle((1, 0), 0.5, name="ring")
    window.session.create_point((2, 0), name="corner")
    window.set_tool("select")
    window.handle_click((1.0, 0.05, 0.0))
    assert window.selected_names() == ["bottom"]
    window.handle_click((1.0, 0.52, 0.0), additive=True)
    assert window.selected_names() == ["bottom", "ring"]
    # A point coinciding with a curve end wins the pick.
    assert pick_entity(window.session, (2.0, 0.01, 0.0), 0.1) == "corner"
    window.handle_click((5.0, 5.0, 0.0))
    assert window.selected_names() == []
    window.set_selection([window.session.entity("ring").id])
    window.delete_selected()
    assert "ring" not in window.session.store
    window.undo()
    assert "ring" in window.session.store


def test_tree_selection_drives_property_panel(window):
    window.session.create_line((0, 0), (1, 0), name="chord")
    entity_id = window.session.entity("chord").id
    window.model_tree.select_ids([entity_id])
    window.model_tree._emit_selection()
    assert window.selected_ids == [entity_id]
    assert window.property_panel.entity_id == entity_id


def test_property_panel_edit_and_undo(window):
    window.session.create_line((0, 0), (1, 0), name="chord")
    window.set_selection([window.session.entity("chord").id])
    panel = window.property_panel
    _kind, editor = panel._editors["p1"]
    editor.setText("2, 0.5")
    _name_kind, name_editor = panel._editors["name"]
    name_editor.setText("wing-chord")
    assert panel.apply()
    entity = window.session.entity("wing-chord")
    assert entity.params["p1"] == (2.0, 0.5, 0.0)
    window.undo()  # rename
    window.undo()  # modify
    assert window.session.entity("chord").params["p1"] == (1.0, 0.0, 0.0)


def test_property_panel_reports_errors(window):
    window.session.create_spline([(0, 0), (1, 1), (2, 0)], name="s")
    window.set_selection([window.session.entity("s").id])
    panel = window.property_panel
    _kind, editor = panel._editors["points"]
    editor.setPlainText("0 0\n1 x\n")
    assert not panel.apply()
    assert "not a number" in panel._error.text()
    editor.setPlainText("0 0\n0 0\n")
    assert not panel.apply()
    assert "coincides" in panel._error.text()
    assert window.session.entity("s").params["points"][1] == (1.0, 1.0, 0.0)


def test_property_panel_every_kind_round_trips_without_changes(window):
    og = window.session
    og.create_point((1, 2, 3))
    og.create_line((0, 0), (1, 0))
    og.create_arc((0, 0), 1.0, 0.0, 2.0)
    og.create_polyline([(0, 0), (1, 1), (2, 0)])
    og.create_spline([(0, 0), (1, 1), (2, 0)], method="clamped", start_tangent=(0, 1))
    og.create_nurbs(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)])
    for entity in list(og.store):
        window.set_selection([entity.id])
        name, changes = window.property_panel.read_changes()
        assert name is None, entity.kind
        assert changes == {}, entity.kind


def test_file_round_trip_and_import(window, tmp_path):
    dat = tmp_path / "naca0012.dat"
    write_selig(dat, naca0012_selig())
    assert window.import_airfoil(str(dat), mode="split") == ["naca0012-upper", "naca0012-lower"]
    assert window.is_modified()
    project = tmp_path / "p.ogp"
    assert window.save_project_as(str(project))
    assert not window.is_modified()
    assert "p.ogp" in window.windowTitle()
    window.new_project(confirm=False)
    assert len(window.session.store) == 0
    assert window.open_project(str(project))
    assert window.model_tree.entity_names() == ["naca0012-upper", "naca0012-lower"]
    # The console's `og` keeps pointing at the live session after a load.
    assert window.console.run_line("print(len(og.store))").strip() == "2"


def test_open_invalid_project_is_reported(window, tmp_path):
    bad = tmp_path / "bad.ogp"
    bad.write_text("{}")
    assert not window.open_project(str(bad))
    assert "not an OpenGrid project" in window.log_panel.toPlainText()


def test_save_journal_replays(window, tmp_path):
    import opengrid

    window.session.create_line((0, 0), (1, 0))
    window.session.create_circle((0, 0), 3.0)
    window.undo()
    journal = tmp_path / "j.py"
    window.save_journal(str(journal))
    replayed = opengrid.replay(journal.read_text())
    assert replayed.names() == ["line-1"]
