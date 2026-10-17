"""GUI interactions driven through the real event chain: Qt mouse and key events on the
viewport widget -> VTK's interactor -> the window's observers and the camera style.

The off-screen window has no interactor of its own (``viewport.iren`` is None); the fixture
attaches one to the widget's generic interactor, exactly as the on-screen window has, so a
``qtbot.mousePress`` on the widget reaches the same observers a real click does. The earlier
GUI tests call ``handle_click`` / ``handle_press`` in world space and never reach this chain,
which is where the zoom-window bug lived.
"""

import math
from unittest import mock

import numpy as np
import pytest
from PySide6.QtCore import QPoint, Qt
from PySide6.QtTest import QTest

from opengrid.app import main_window as mw
from opengrid.app.main_window import MainWindow


class Live:
    """A window whose viewport takes real Qt events."""

    def __init__(self, qtbot):
        from pyvista.plotting.render_window_interactor import RenderWindowInteractor

        self.qtbot = qtbot
        self.window = MainWindow(off_screen=True)
        qtbot.addWidget(self.window)
        self.window.resize(1200, 800)
        self.window.show()
        viewport = self.window.viewport
        iren = RenderWindowInteractor(viewport, interactor=viewport._Iren)
        viewport.iren = iren
        iren.initialize()  # enables VTK's event delivery (a disabled interactor drops events)
        mw.set_navigation_style(iren)
        self.window._connect_viewport_events()
        self.widget = viewport
        # Off-screen the render window has no size until something resizes it: give it the
        # widget's (pixel ratio 1), so display and widget coordinates agree.
        viewport.ren_win.SetSize(max(viewport.width(), 1), max(viewport.height(), 1))
        self.og = self.window.session
        # QTest does not update QApplication.mouseButtons(): the window asks it whether the
        # button is still down during a drag.
        self._buttons = Qt.MouseButton.NoButton
        self._patch = mock.patch.object(
            mw, "_left_button_down", lambda: bool(self._buttons & Qt.MouseButton.LeftButton)
        )
        self._patch.start()
        qtbot.addWidget(self.window)

    def close(self):
        self._patch.stop()
        self.window.close()

    def fit(self):
        self.window.scene.fit()
        self.widget.render()

    def display(self, xyz) -> QPoint:
        """Widget coordinates (Qt, y down) of a world point."""
        renderer = self.widget.renderer
        renderer.SetWorldPoint(*xyz, 1.0)
        renderer.WorldToDisplay()
        x, y, _z = renderer.GetDisplayPoint()
        scale = self.widget.ren_win.GetSize()[0] / max(self.widget.width(), 1)  # pixel ratio
        return QPoint(round(x / scale), round(self.widget.height() - 1 - y / scale))

    def press(self, xyz, modifier=Qt.KeyboardModifier.NoModifier):
        self._buttons = Qt.MouseButton.LeftButton
        QTest.mousePress(self.widget, Qt.MouseButton.LeftButton, modifier, self.display(xyz))

    def move(self, xyz, modifier=Qt.KeyboardModifier.NoModifier):
        QTest.mouseMove(self.widget, self.display(xyz))
        # QTest.mouseMove alone does not deliver a move while a button is held off-screen:
        # send the event the widget would get.
        from PySide6.QtCore import QEvent, QPointF
        from PySide6.QtGui import QMouseEvent

        point = QPointF(self.display(xyz))
        event = QMouseEvent(
            QEvent.Type.MouseMove,
            point,
            point,
            Qt.MouseButton.NoButton,
            self._buttons,
            modifier,
        )
        self.widget.mouseMoveEvent(event)

    def release(self, xyz, modifier=Qt.KeyboardModifier.NoModifier):
        QTest.mouseRelease(self.widget, Qt.MouseButton.LeftButton, modifier, self.display(xyz))
        self._buttons = Qt.MouseButton.NoButton

    def double_click(self, xyz):
        """Press, release, double-click, release: what Qt sends for a real double click (the
        last release matters: QTest's button state is process-wide and later tests ask it)."""
        point = self.display(xyz)
        button = Qt.MouseButton.LeftButton
        modifier = Qt.KeyboardModifier.NoModifier
        QTest.mouseClick(self.widget, button, modifier, point)
        QTest.mouseDClick(self.widget, button, modifier, point)
        QTest.mouseRelease(self.widget, button, modifier, point)
        self._buttons = Qt.MouseButton.NoButton

    def click(self, xyz, shift: bool = False):
        modifier = Qt.KeyboardModifier.NoModifier
        if shift:
            modifier = Qt.KeyboardModifier.ShiftModifier
        self.press(xyz, modifier)
        self.release(xyz, modifier)

    def drag(self, start, end, steps: int = 8, modifier=Qt.KeyboardModifier.NoModifier):
        self.press(start, modifier)
        a = np.asarray(start, dtype=float)
        b = np.asarray(end, dtype=float)
        for k in range(1, steps + 1):
            self.move(a + (b - a) * k / steps, modifier)
        self.release(end, modifier)

    def key(self, key, widget=None):
        if widget is None:
            widget = self.widget
        widget.setFocus()
        QTest.keyClick(widget, key)

    def actor(self, name: str):
        return self.widget.renderer.actors.get(name)


@pytest.fixture
def live(qtbot):
    pytest.importorskip("pyvistaqt")
    live = Live(qtbot)
    og = live.og
    og.create_spline([(0, 0), (1, 0.1), (2, 0)], name="b")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_line((2, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    og.create_connector("b", 41, name="cb")
    og.create_connector("t", 41, name="ct")
    og.create_connector("r", 21, name="cr")
    og.create_connector("l", 21, name="cl")
    og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="d")
    live.fit()
    yield live
    live.close()


def test_the_event_chain_reaches_the_window(live):
    seen = []
    live.window.handle_click = lambda xyz, additive=False: seen.append((tuple(xyz), additive))
    live.click((2.0, 0.5, 0.0))
    assert len(seen) == 1
    xyz, additive = seen[0]
    assert math.dist(xyz[:2], (2.0, 0.5)) < 0.02 and not additive  # within a pixel or two
    live.click((2.0, 0.5, 0.0), shift=True)
    assert seen[1][1] is True


def test_clicks_select_and_shift_adds(live):
    live.click((2.0, 0.5, 0.0))
    assert live.window.selected_names() == ["cr"]
    live.click((1.0, 1.0, 0.0), shift=True)
    assert sorted(live.window.selected_names()) == ["cr", "ct"]
    live.click((1.0, 1.0, 0.0), shift=True)  # toggles off
    assert live.window.selected_names() == ["cr"]
    live.click((1.0, 0.5, 0.0))  # inside the domain
    assert live.window.selected_names() == ["d"]
    live.click((5.0, 5.0, 0.0))  # nothing there
    assert live.window.selected_names() == []


def test_line_tool_two_clicks_with_a_live_preview(live):
    live.window.set_tool("line")
    before = set(live.og.names("line"))
    live.click((3.0, 0.0, 0.0))
    assert len(live.window.tool.points) == 1
    live.move((3.5, 0.5, 0.0))
    assert live.actor("tool-preview") is not None  # the rubber band follows the pointer
    live.click((4.0, 1.0, 0.0))
    (made,) = set(live.og.names("line")) - before
    p0, p1 = live.og.entity(made).params["p0"], live.og.entity(made).params["p1"]
    assert math.dist(p0[:2], (3.0, 0.0)) < 0.02 and math.dist(p1[:2], (4.0, 1.0)) < 0.02
    assert live.actor("tool-preview") is None
    assert not live.window.zoom_box.active


def test_drag_in_a_drawing_tool_zooms_and_adds_no_point(live):
    live.window.set_tool("line")
    camera = live.widget.renderer.GetActiveCamera()
    scale = camera.GetParallelScale()
    live.drag((0.2, 0.2, 0.0), (1.8, 0.9, 0.0))
    assert camera.GetParallelScale() < scale
    assert live.window.tool.points == []
    assert not live.window.zoom_box.active


def test_escape_in_the_viewport_cancels_then_clears(live):
    live.window.set_tool("line")
    live.click((3.0, 0.0, 0.0))
    live.key(Qt.Key.Key_Escape)
    assert live.window.tool.points == []
    assert live.window.tool.name == "line"
    live.key(Qt.Key.Key_Escape)
    assert live.window.tool.name == "select"
    live.click((2.0, 0.5, 0.0))
    assert live.window.selected_names() == ["cr"]
    live.key(Qt.Key.Key_Escape)
    assert live.window.selected_names() == []


def test_return_finishes_a_polyline(live):
    live.window.set_tool("polyline")
    for point in ((3, 0), (3.5, 0.5), (4, 0)):
        live.click((*point, 0.0))
    before = set(live.og.names("polyline"))
    live.key(Qt.Key.Key_Return)
    assert len(set(live.og.names("polyline")) - before) == 1


def test_handle_drag_moves_a_control_point(live):
    live.click((1.0, 0.1, 0.0))
    assert live.window.selected_names() == ["cb"]  # the connector lies on the curve: it wins
    # ... and offers the curve's handles, so the curve can still be dragged from the viewport.
    assert [h.entity for h in live.window.visible_handles()] == ["b", "b", "b"]
    journal = len(live.og.journal_text().splitlines())
    live.drag((1.0, 0.1, 0.0), (1.0, 0.3, 0.0))
    moved = live.og.entity("b").params["points"][1]
    assert math.dist(moved[:2], (1.0, 0.3)) < 0.02
    assert len(live.og.journal_text().splitlines()) == journal + 1
    assert not live.window.zoom_box.active  # the drag took the button before the zoom box
    assert live.window.selected_names() == ["cb"]
    grid = live.og.domain_grid("d")
    assert abs(grid[0, 20, 1] - 0.3) < 0.05  # the connector and the domain followed


def test_freehand_stroke_becomes_curves(live):
    live.window.set_tool("freehand")
    before = set(live.og.names("nurbs"))
    live.press((3.0, 0.0, 0.0))
    for k in range(1, 30):
        live.move((3 + 0.05 * k, 0.3 * math.sin(k / 5), 0.0))
    live.release((4.5, 0.3 * math.sin(6.0), 0.0))
    made = set(live.og.names("nurbs")) - before
    assert len(made) >= 1
    assert not live.window.zoom_box.active and live.window._stroke is None


def test_sketch_tool_four_clicks(live):
    live.window.set_tool("sketch")
    for point in ((3, 0), (4, 0), (4, 1), (3, 1)):
        live.click((*point, 0.0))
    assert len(live.og.names("sketch")) == 1
    sketch = live.og.entity(live.og.names("sketch")[0])
    assert len(sketch.params["vertices"]) == 4


def test_tool_switch_clears_points_and_preview(live):
    live.window.set_tool("line")
    live.click((3.0, 0.0, 0.0))
    live.move((3.5, 0.5, 0.0))
    live.window.tool_actions["select"].trigger()
    assert live.window.tool.name == "select" and live.window.tool.points == []
    assert live.actor("tool-preview") is None


def test_boundary_dock_buttons(live):
    dock = live.window.boundary_dock
    live.click((2.0, 0.5, 0.0))
    dock.name.setText("walls")
    dock.patch_type.setCurrentText("wall")
    live.qtbot.mouseClick(dock.new_patch, Qt.MouseButton.LeftButton)
    cr, ct = live.og.entity("cr").id, live.og.entity("ct").id
    assert [m["connector"] for m in live.og.entity("walls").params["members"]] == [cr]
    live.click((1.0, 1.0, 0.0))
    live.qtbot.mouseClick(dock.add, Qt.MouseButton.LeftButton)
    assert [m["connector"] for m in live.og.entity("walls").params["members"]] == [cr, ct]
    live.qtbot.mouseClick(dock.validate, Qt.MouseButton.LeftButton)
    assert "unassigned" in dock.summary.text()


def test_delete_key_in_the_viewport(live):
    live.window.set_tool("point")
    live.click((5.0, 5.0, 0.0))
    (name,) = live.og.names("point")
    live.window.set_tool("select")
    live.click((5.0, 5.0, 0.0))
    assert live.window.selected_names() == [name]
    live.key(Qt.Key.Key_Backspace)
    assert live.og.names("point") == []


def test_model_tree_click_selects_in_the_viewport(live):
    tree = live.window.model_tree
    item = next(
        tree.topLevelItem(c).child(k)
        for c in range(tree.topLevelItemCount())
        for k in range(tree.topLevelItem(c).childCount())
        if tree.topLevelItem(c).child(k).text(0) == "cr"
    )
    tree.scrollToItem(item)  # the docks are short: the item may be scrolled out of view
    rect = tree.visualItemRect(item)
    live.qtbot.mouseClick(tree.viewport(), Qt.MouseButton.LeftButton, pos=rect.center())
    assert live.window.selected_names() == ["cr"]
    assert live.window.property_panel.entity_id == live.og.entity("cr").id


def test_wheel_zooms_and_f_fits(live):
    from PySide6.QtCore import QPointF
    from PySide6.QtGui import QWheelEvent

    camera = live.widget.renderer.GetActiveCamera()
    scale = camera.GetParallelScale()
    point = QPointF(live.display((1.0, 0.5, 0.0)))
    event = QWheelEvent(
        point,
        live.widget.mapToGlobal(point.toPoint()).toPointF(),
        QPoint(0, 0),
        QPoint(0, 120),
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
        Qt.ScrollPhase.NoScrollPhase,
        False,
    )
    live.widget.wheelEvent(event)
    assert camera.GetParallelScale() != scale  # two-finger scroll / wheel zooms
    live.key(Qt.Key.Key_F)
    assert abs(camera.GetParallelScale() - scale) < 1e-9 * scale  # F fits the view again


def test_shift_drag_pans_in_a_drawing_tool(live):
    live.window.set_tool("line")
    camera = live.widget.renderer.GetActiveCamera()
    position = camera.GetPosition()
    live.drag((0.5, 0.5, 0.0), (1.5, 0.5, 0.0), modifier=Qt.KeyboardModifier.ShiftModifier)
    assert camera.GetPosition() != position
    assert live.window.tool.points == [] and not live.window.zoom_box.active


def test_property_panel_edit_by_keyboard(live):
    from PySide6.QtWidgets import QSpinBox

    panel = live.window.property_panel
    live.click((2.0, 0.5, 0.0))  # a side of the domain: its opposite side must match
    assert panel.entity_id == live.og.entity("cr").id
    spin = next(s for s in panel.findChildren(QSpinBox) if s.value() == 21)  # the point count
    spin.lineEdit().selectAll()
    live.qtbot.keyClicks(spin, "31")
    live.qtbot.mouseClick(panel.apply_button, Qt.MouseButton.LeftButton)
    assert live.og.entity("cr").params["dimension"] == 21  # refused, with the reason shown
    assert "cl" in panel._error.text() or "opposite" in panel._error.text().lower()
    live.og.create_line((3, 0), (4, 0), name="free")
    live.og.create_connector("free", 11, name="cf")
    live.click((3.5, 0.0, 0.0))
    assert panel.entity_id == live.og.entity("cf").id
    spin = next(s for s in panel.findChildren(QSpinBox) if s.value() == 11)
    spin.lineEdit().selectAll()
    live.qtbot.keyClicks(spin, "31")
    live.qtbot.mouseClick(panel.apply_button, Qt.MouseButton.LeftButton)
    assert live.og.entity("cf").params["dimension"] == 31
    live.window.undo()
    assert live.og.entity("cf").params["dimension"] == 11


def test_freehand_loop_released_near_its_start_meshes_with_one_click(live):
    window = live.window
    window.set_tool("freehand")
    before = set(live.og.names("nurbs"))
    centre = np.array([1.0, 0.5, 0.0])
    radius = 0.3
    pixel = window.scene.world_per_pixel()
    start = centre + np.array([radius, 0.0, 0.0])
    live.press(start)
    for angle in np.linspace(0.0, 2 * math.pi, 60)[1:-1]:
        live.move(centre + np.array([radius * math.cos(angle), radius * math.sin(angle), 0.0]))
    live.release(start - np.array([0.0, 3 * pixel, 0.0]))  # three pixels short of the start
    (curve,) = set(live.og.names("nurbs")) - before
    geometry = live.og.geometry(curve)
    assert np.linalg.norm(np.subtract(geometry.eval(1.0), geometry.eval(0.0))) == 0.0
    assert window.selected_names() == [curve]  # ready for the next action
    window.dimension_box.setValue(41)
    window.unstructured_action.trigger()
    (domain,) = live.og.names("unstructured_domain")
    assert window.selected_names() == [domain]
    assert "triangles" in live.og.domain_status(domain)
    (loop,) = live.og.entity(domain).params["loops"]
    assert len(loop) == 1
    assert live.og.entity(loop[0]["connector"]).params["dimension"] == 41
    window.undo_action.trigger()  # one step: the connector and the domain
    assert live.og.names("unstructured_domain") == [] and live.og.names("connector") == [
        "cb",
        "ct",
        "cr",
        "cl",
    ]


def test_actions_explain_a_wrong_selection(live):
    window = live.window
    window.set_selection([live.og.entity("d").id])
    window.unstructured_action.trigger()
    assert (
        window.statusBar().currentMessage()
        == "d is a domain: select the connectors or curves of the boundary loops"
    )
    window.set_selection([])
    window.unstructured_action.trigger()
    assert (
        window.statusBar().currentMessage()
        == "Select the connectors or curves of the boundary loops first"
    )
    assert live.og.names("unstructured_domain") == []


def test_double_click_adds_and_removes_a_spline_point(live):
    window = live.window
    live.click((1.0, 0.1, 0.0))  # the connector on the spline b
    assert window.selected_names() == ["cb"]
    before = len(live.og.entity("b").params["points"])
    journal = len(live.og.journal_text().splitlines())
    live.double_click((0.5, 0.07, 0.0))  # on the spline, between its first two points
    points = live.og.entity("b").params["points"]
    assert len(points) == before + 1
    assert window.statusBar().currentMessage().startswith("Added a point to b")
    assert len(live.og.journal_text().splitlines()) == journal + 1
    assert live.og.domain_quality("d")["scaled_jacobian"].min() > 0.5  # d follows b
    live.double_click(points[1])  # the new handle
    assert len(live.og.entity("b").params["points"]) == before
    assert window.statusBar().currentMessage().startswith("Removed point 1 of b")
    assert window.selected_names() == ["cb"]
    window.set_selection([])
    live.double_click((1.0, 0.5, 0.0))  # inside the domain, on nothing editable
    assert window.statusBar().currentMessage().startswith("Double-click a polyline")
    assert len(live.og.entity("b").params["points"]) == before


def test_hover_highlights_the_entity_under_the_pointer(live):
    """A resting pointer highlights what a click would pick; the status bar shows the 2D
    coordinates; the highlight goes when the entity is selected or the pointer leaves."""
    window = live.window
    scene = window.scene
    live.move((2.0, 0.5, 0.0))  # on connector cr
    live.qtbot.waitUntil(lambda: scene.hovered == live.og.entity("cr").id, timeout=2000)
    assert scene._overlays.get("hover") is not None
    shown = window.coordinates_label.text().split()  # "x 1.99968   y 0.5" (pixel rounding)
    assert abs(float(shown[1]) - 2.0) < 0.01 and abs(float(shown[3]) - 0.5) < 0.01
    live.move((5.0, 5.0, 0.0))  # nothing there
    live.qtbot.waitUntil(lambda: scene.hovered is None, timeout=2000)
    assert scene._overlays.get("hover") is None
    live.move((1.0, 0.5, 0.0))  # inside the domain
    live.qtbot.waitUntil(lambda: scene.hovered == live.og.entity("d").id, timeout=2000)
    live.click((1.0, 0.5, 0.0))  # selected: the hover highlight makes way
    assert window.selected_names() == ["d"] and scene.hovered is None
    window.set_tool("line")  # drawing tools do not hover
    live.move((2.0, 0.5, 0.0))
    live.qtbot.wait(80)
    assert scene.hovered is None
    window.set_tool("select")
    window._on_vtk_leave()
    assert window.coordinates_label.text() == ""


def test_context_menu_offers_what_applies_to_the_selection(live):
    window = live.window
    og = live.og

    def labels(menu):
        return [a.text().replace("&", "") for a in menu.actions() if not a.isSeparator()]

    window.set_selection([])
    names = labels(window.context_menu())
    assert "Delete" not in names and "Reset camera" in names and "Topology" in names
    window.set_selection([og.entity("b").id])  # a curve
    names = labels(window.context_menu())
    assert "Connectors on curves" in names and "Structured domain" not in names
    assert "Zoom to selection" in names and "Delete" in names
    window.set_selection([og.entity("cr").id])  # a connector
    names = labels(window.context_menu())
    assert "Structured domain" in names and "Unstructured domain" in names
    assert "Connectors on curves" not in names
    window.set_selection([og.entity("d").id])  # a domain
    names = labels(window.context_menu())
    assert "Extrude domain" in names and "Smoothing group" in names
    # A right click (no drag) opens it in the viewport: exercised through the handlers.
    opened = []
    window.show_context_menu = lambda position: opened.append(position)
    window._right_press_position = (100.0, 100.0)
    window.viewport.iren.interactor.SetEventInformation(101, 100, 0, 0, "\0", 0, None)
    window._on_vtk_right_release()
    assert len(opened) == 1
    window._right_press_position = (100.0, 100.0)
    window.viewport.iren.interactor.SetEventInformation(160, 140, 0, 0, "\0", 0, None)
    window._on_vtk_right_release()  # a drag: no menu
    assert len(opened) == 1
