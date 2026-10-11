import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.main_window import MainWindow
from opengrid.app.scene import pick_entity


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    win.scene.world_per_pixel = lambda: 0.01
    og = win.session
    og.create_line((0, 0), (1, 0), name="a")
    og.create_line((1, 0), (1, 1), name="b")
    og.create_arc_3pt((1, 1), (0.5, 1.5), (0, 1), name="c")
    yield win
    win.close()


def _ids(window, *names):
    return [window.session.entity(n).id for n in names]


def test_connectors_on_selected_curves(window):
    window.set_selection(_ids(window, "a", "b"))
    window.dimension_box.setValue(11)
    names = window.create_connectors_on_selection()
    assert names == ["connector-1", "connector-2"]
    assert window.model_tree.entity_names("Connectors") == names
    assert window.selected_names() == names
    assert len(window.session.connector_nodes("connector-1").points) == 11
    window.undo()
    assert window.model_tree.entity_names("Connectors") == []


def test_chain_connector_from_selection(window):
    window.set_selection(_ids(window, "c", "a", "b"))
    name = window.create_chain_connector_on_selection()
    nodes = window.session.connector_nodes(name)
    # Starts at the free end of the first selected curve that has one ("c" ends at (0, 1)).
    assert nodes.points[0] == pytest.approx([0.0, 1.0, 0.0])
    assert nodes.points[-1] == pytest.approx([0.0, 0.0, 0.0])
    assert len(window.session.entity(name).params["segments"]) == 3


def test_no_selection_is_reported(window):
    window.set_selection([])
    assert window.create_connectors_on_selection() == []
    assert "Select" in window.statusBar().currentMessage()


def test_pick_prefers_connector_over_its_curve(window):
    window.session.create_connector("a", 5, name="ca")
    assert pick_entity(window.session, (0.5, 0.01, 0.0), 0.1) == "ca"
    window.set_tool("select")
    window.handle_click((0.5, 0.01, 0.0))
    assert window.selected_names() == ["ca"]


def test_property_panel_edits_connector(window):
    window.session.create_connector("a", 11, name="ca")
    window.set_selection(_ids(window, "ca"))
    panel = window.property_panel
    _kind, dimension = panel._editors["dimension"]
    dimension.setValue(31)
    _kind, dist = panel._editors["distribution"]
    dist.type.setCurrentText("tanh")
    dist.begin.setText("0.005")
    dist.end.setText("0.02")
    _kind, breakpoints = panel._editors["breakpoints"]
    breakpoints.setPlainText("0.5 15 geometric 0.02 - -")
    assert panel.apply()
    nodes = window.session.connector_nodes("ca")
    assert len(nodes.points) == 31
    assert abs(nodes.begin_spacing - 0.005) / 0.005 < 1e-8
    params = window.session.entity("ca").params
    assert params["breakpoints"][0]["index"] == 15
    # Re-reading an unchanged panel yields no changes.
    window.set_selection(_ids(window, "ca"))
    assert window.property_panel.read_changes() == (None, {})
    window.undo()
    assert window.session.entity("ca").params["dimension"] == 11


def test_property_panel_connector_errors(window):
    window.session.create_connector("a", 11, name="ca")
    window.set_selection(_ids(window, "ca"))
    panel = window.property_panel
    _kind, breakpoints = panel._editors["breakpoints"]
    breakpoints.setPlainText("0.5 x uniform")
    assert not panel.apply()
    assert "not an integer" in panel._error.text()
    breakpoints.setPlainText("")
    _kind, dist = panel._editors["distribution"]
    dist.type.setCurrentText("geometric")
    assert not panel.apply()
    assert "exactly one" in panel._error.text()


def test_stretching_warning_and_threshold(window):
    window.session.create_connector("a", 21, {"type": "geometric", "growth_ratio": 1.4}, name="ca")
    assert "max stretching ratio 1.4 exceeds 1.3" in window.log_panel.toPlainText()
    item_text = window.model_tree.category_item("Connectors").child(0).text(1)
    assert "⚠" in item_text
    window.set_threshold(1.5)
    assert "⚠" not in window.model_tree.category_item("Connectors").child(0).text(1)


def test_colouring_and_nodes_toggle(window):
    window.session.create_connector("a", 21, {"type": "tanh", "begin_spacing": 0.01}, name="ca")
    entity_id = window.session.entity("ca").id
    assert window.scene._actors[entity_id] == [f"entity-{entity_id}", f"entity-{entity_id}-nodes"]
    window.nodes_action.setChecked(False)
    assert window.scene._actors[entity_id] == [f"entity-{entity_id}"]
    window.color_action.setChecked(False)
    assert not window.scene.color_by_stretching
    window.render_check = window.grab()
    assert not window.render_check.isNull()


def test_modifying_curve_redraws_connector(window):
    window.session.create_connector("b", 5, name="cb")
    window.session.modify("b", p1=(1.0, 3.0))
    nodes = window.session.connector_nodes("cb")
    assert np.allclose(nodes.points[-1], [1.0, 3.0, 0.0])
