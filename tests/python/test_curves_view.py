"""One curve in the window: a database curve and the connector of its points shown, selected,
edited and deleted as one item (model/curves.py), and structured domains straight from
curves."""

import numpy as np
import pytest

import opengrid
from opengrid.model.curves import CurveIndex
from opengrid.validation import ValidationError


def test_pairing_rules():
    og = opengrid.Session()
    og.create_line((0, 0), (1, 0), name="a")  # with its points: one curve
    og.create_line((1, 0), (1, 1), name="b")  # split in two: two curves along it
    og.create_line((1, 1), (0, 1), name="c")  # no points yet
    og.create_line((0, 1), (0, 0), name="d")
    og.create_line((0, 0), (-1, 0), name="e")  # d and e: one chain
    og.create_connector("a", 11, name="ca")
    og.create_connector([{"curve": "b", "t0": 0.0, "t1": 0.5}], 6, name="cb1")
    og.create_connector([{"curve": "b", "t0": 0.5, "t1": 1.0}], 6, name="cb2")
    og.create_chain_connector(["d", "e"], 21, name="cde")
    index = CurveIndex(og.store)

    def ids(*names):
        return [og.entity(n).id for n in names]

    assert index.connector_of == {og.entity("a").id: og.entity("ca").id}
    assert index.representative(og.entity("a").id) == og.entity("ca").id
    assert index.members(og.entity("ca").id) == ids("a", "ca")
    assert index.members(og.entity("a").id) == ids("a", "ca")
    assert index.members(og.entity("cb1").id) == ids("cb1")
    listed = [e.name for e in og.store if index.is_listed(e) and e.kind != "point"]
    assert listed == ["c", "ca", "cb1", "cb2", "cde"]  # a, b, d, e stand behind their points
    texts = {name: index.describe(og.entity(name)) for name in listed}
    assert texts == {
        "c": "line, no points",
        "ca": "line a, 11 pts",
        "cb1": "6 pts, part of a line",
        "cb2": "6 pts, part of a line",
        "cde": "21 pts, along 2 curves",
    }
    # Automatic names say nothing: not repeated.
    og.create_line((5, 0), (6, 0))
    og.create_connectors_on_curves(["line-1"], 5)
    assert CurveIndex(og.store).describe(og.entity("curve-1")) == "line, 5 pts"


def test_automatic_names_of_points_are_curves():
    og = opengrid.Session()
    lines = og.create_rectangle(width=2, height=1)
    assert og.create_connectors_on_curves(lines[:2], 5) == ["curve-1", "curve-2"]
    og.create_connector(lines[2], 5)
    assert og.names("connector")[-1] == "curve-3"
    # A replay ends with the same counters (explicit names advance them).
    again = opengrid.replay(og.journal_text())
    again.create_connector(lines[3], 5)
    assert again.names("connector")[-1] == "curve-4"


def test_structured_domain_from_curves_is_one_step():
    og = opengrid.Session()
    lines = og.create_rectangle(width=2, height=1)
    og.create_connector(lines[0], 21, name="bottom")  # one side has points already
    domain = og.structured_domain_from_curves(lines[1:], 21, connectors=["bottom"])
    assert og.domain_dimensions(domain) == (21, 21)
    assert sorted(og.names("connector")) == ["bottom", "curve-1", "curve-2", "curve-3"]
    og.undo()
    assert og.names("connector") == ["bottom"] and og.names("domain") == []
    og.redo()
    replayed = opengrid.replay(og.journal_text())
    assert replayed.domain_dimensions(domain) == (21, 21)
    # Curves that bound no domain (an open loop): refused, nothing left behind.
    with pytest.raises(ValidationError):
        og.structured_domain_from_curves(og.create_rectangle(center=(5, 0))[:3], 5)
    assert sorted(og.names("connector")) == ["bottom", "curve-1", "curve-2", "curve-3"]


# ------------------------------------------------------------------------ window
@pytest.fixture
def window(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    yield win
    win.close()


def test_a_curve_is_one_row_one_selection_one_panel(window):
    og = window.session
    og.create_line((0, 0), (1, 0), name="a")
    og.create_connector("a", 11, name="ca")
    tree = window.model_tree
    assert tree.entity_names() == ["ca"]
    window.set_selection([og.entity("a").id])  # the shape (the viewport picks it)
    assert window.selected_names() == ["ca"]
    assert [item.text(0) for item in tree.selectedItems()] == ["ca"]
    panel = window.property_panel
    assert "curve" in panel._title.text() and "line, 11 pts" in panel._title.text()
    assert {"shape:p0", "shape:p1", "dimension", "distribution"} <= set(panel._editors)
    # The shape's end moves from the panel; the points follow.
    panel._editors["shape:p1"][1].setText("2, 0, 0")
    assert panel.apply()
    assert og.entity("a").params["p1"] == (2.0, 0.0, 0.0)
    assert og.connector_nodes("ca").points[-1][0] == pytest.approx(2.0)
    # Both edits in one Apply: the shape, then the points.
    panel._editors["shape:p0"][1].setText("-1, 0, 0")
    panel._editors["dimension"][1].setValue(21)
    assert panel.apply()
    assert og.entity("a").params["p0"] == (-1.0, 0.0, 0.0)
    assert len(og.connector_nodes("ca").points) == 21


def test_curve_without_points_gets_them_from_its_panel(window):
    og = window.session
    og.create_line((0, 0), (1, 0), name="a")
    window.dimension_box.setValue(7)
    window.set_selection([og.entity("a").id])
    assert "line, no points" in window.property_panel._title.text()
    button = window.property_panel._form.itemAt(
        window.property_panel._form.rowCount() - 1, window.property_panel._form.ItemRole.FieldRole
    ).widget()
    assert button.text() == "Give it points"
    button.click()
    assert og.names("connector") == ["curve-1"]
    assert window.selected_names() == ["curve-1"]
    assert len(og.connector_nodes("curve-1").points) == 7


def test_deleting_a_curve_takes_its_shape(window, monkeypatch):
    og = window.session
    og.create_line((0, 0), (1, 0), name="a")
    og.create_connector("a", 11, name="ca")
    og.create_line((0, 0), (0, 1), name="kept")
    window.set_selection([og.entity("ca").id])
    window.delete_selected()
    assert "a" not in og.store and "ca" not in og.store and "kept" in og.store
    assert "Deleted ca" in window.statusBar().currentMessage()
    og.undo()
    assert "a" in og.store and "ca" in og.store


def test_structured_domain_from_selected_curves(window):
    og = window.session
    lines = og.create_rectangle(width=2, height=1)
    window.dimension_box.setValue(9)
    window.set_selection([og.entity(n).id for n in lines])
    name = window.create_domain_on_selection()
    assert og.domain_dimensions(name) == (9, 9)
    assert window.model_tree.entity_names() == ["curve-1", "curve-2", "curve-3", "curve-4"]
    assert window.model_tree.entity_names("Domains") == [name]
    window.undo()  # one step
    assert og.names("connector") == [] and og.names("domain") == []


def test_points_are_previewed_before_meshing(window):
    og = window.session
    og.create_line((0, 0), (1, 0), name="a")
    og.create_line((1, 0), (1, 1), name="b")
    window.dimension_box.setValue(6)
    window.set_selection([og.entity("a").id, og.entity("b").id])
    preview = window.scene._overlays["points-preview"].mapper.dataset
    assert preview.n_points == 12  # 6 points on each curve, nothing created
    assert og.names("connector") == []
    window.dimension_box.setValue(9)  # follows the 'points' box
    assert window.scene._overlays["points-preview"].mapper.dataset.n_points == 18
    # A curve with points: its panel's edits are previewed before Apply.
    window.set_selection([])
    assert window.scene._overlays.get("points-preview") is None
    og.create_connector("a", 11, name="ca")
    window.set_selection([og.entity("ca").id])
    panel = window.property_panel
    panel._editors["dimension"][1].setValue(4)
    points = window.scene._overlays["points-preview"].mapper.dataset.points
    assert len(points) == 4 and len(og.connector_nodes("ca").points) == 11
    assert np.allclose(points[:, 0], [0, 1 / 3, 2 / 3, 1])
    assert panel.apply()  # applied: the preview goes, the points are there
    assert window.scene._overlays.get("points-preview") is None
    assert len(og.connector_nodes("ca").points) == 4
    # The API: where the points would go, changing nothing.
    assert og.preview_points("b", 3).tolist() == [[1, 0, 0], [1, 0.5, 0], [1, 1, 0]]


def test_selected_curves_are_bold(window):
    from opengrid.app.scene import SELECTED_CURVE_WIDTH

    og = window.session
    og.create_line((0, 0), (1, 0), name="a")
    og.create_connector("a", 11, name="ca")
    og.create_line((0, 0), (0, 1), name="plain")
    og.create_line((2, 0), (3, 0), name="other")
    og.create_connector("other", 5, name="co")  # not selected: drawn normally
    window.set_selection([og.entity("ca").id, og.entity("plain").id])
    for layer in ("connectors-selected", "curves-selected"):
        assert window.scene.layer(layer).prop.line_width == SELECTED_CURVE_WIDTH
    assert window.scene.layer("connectors").prop.line_width < SELECTED_CURVE_WIDTH
    # Coloured by patch: the selected curve keeps its patch's colour, in bold.
    og.set_boundary_type(["ca"], "inlet")
    window.patch_color_action.setChecked(True)
    window.set_selection([og.entity("ca").id])
    selected = window.scene.layer("connectors-selected")
    colour = window.scene.layer_data("connectors-selected").cell_data["colors"][0]
    expected = window.scene.patch_colors[og.entity("ca").id]
    import pyvista as pv

    assert tuple(colour) == pv.Color(expected).int_rgb
    assert selected.prop.line_width == SELECTED_CURVE_WIDTH
