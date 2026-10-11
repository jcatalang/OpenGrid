import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.main_window import MainWindow
from opengrid.app.scene import pick_entity


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    og = win.session
    og.create_circle((0, 0), 1.0, name="hole")
    og.create_circle((0, 0), 4.0, name="outer")
    og.create_connector("hole", 41, name="c-hole")
    og.create_connector("outer", 61, name="c-outer")
    win.scene.world_per_pixel = lambda: 0.01
    yield win
    win.close()


def _ids(window, *names):
    return [window.session.entity(n).id for n in names]


def test_create_from_selection_tree_and_pick(window):
    window.set_selection(_ids(window, "c-hole", "c-outer"))
    name = window.create_unstructured_on_selection()
    assert name == "unstructured_domain-1"
    count = len(window.session.unstructured_mesh(name)["triangles"])
    item = window.model_tree.category_item("Domains").child(0)
    assert (item.text(0), item.text(1)) == (name, f"{count} tris")
    assert pick_entity(window.session, (2.5, 0.3, 0.0), 0.05) == name
    assert pick_entity(window.session, (0.2, 0.1, 0.0), 0.05) is None  # inside the hole
    assert str(count) in window.quality_dock.stats.text()


def test_property_panel_generator_edit(window):
    name = window.session.create_unstructured_domain(connectors=["c-hole", "c-outer"])
    before = len(window.session.unstructured_mesh(name)["triangles"])
    window.set_selection(_ids(window, name))
    panel = window.property_panel
    _kind, editor = panel._editors["generator"]
    editor.max_size.setText("0.15")
    editor.sources.setPlainText("2 0 0.02 0.05\n")
    assert panel.apply()
    generator = window.session.entity(name).params["generator"]
    assert generator["max_size"] == 0.15
    assert generator["sources"][0]["size"] == 0.02
    assert len(window.session.unstructured_mesh(name)["triangles"]) > before
    window.set_selection(_ids(window, name))
    assert window.property_panel.read_changes() == (None, {})
    _kind, editor = window.property_panel._editors["generator"]
    editor.sources.setPlainText("1 2 3")
    assert not window.property_panel.apply()
    assert "sources line 1" in window.property_panel._error.text()
    window.undo()
    assert len(window.session.unstructured_mesh(name)["triangles"]) == before


def test_quality_colouring_and_worst_triangle(window):
    name = window.session.create_unstructured_domain(connectors=["c-hole", "c-outer"])
    window.quality_dock.metric.setCurrentIndex(3)  # min angle
    window.quality_action.setChecked(True)
    result = window.show_worst_cell()
    assert result[0] == name
    q = window.session.domain_quality(name)["min_angle"].reshape(-1)
    assert q[result[1]] == q.min()
    assert not window.grab().isNull()


def test_algorithm_combo(window):
    name = window.session.create_unstructured_domain(connectors=["c-hole", "c-outer"])
    window.set_selection(_ids(window, name))
    _kind, editor = window.property_panel._editors["generator"]
    assert editor.min_angle.isEnabled()
    editor.algorithm.setCurrentText("advancing_front")
    assert not editor.min_angle.isEnabled()
    assert not editor.refine.isEnabled()
    assert window.property_panel.apply()
    assert window.session.entity(name).params["generator"]["algorithm"] == "advancing_front"
    assert window.session.domain_quality(name)["min_angle"].min() >= 25.0
