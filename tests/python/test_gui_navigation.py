"""GUI navigation: zoom to a window (app/zoom_box.py) and docks whose dividers move freely."""

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QScrollArea

from opengrid.app.main_window import MainWindow
from opengrid.app.zoom_box import zoom_to_box

PROJECT = Path(__file__).resolve().parents[2] / "examples" / "output" / "04_airfoil_cgrid.ogp"


def _window(qtbot) -> MainWindow:
    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    window.resize(1200, 800)
    window.open_project(str(PROJECT))
    window.show()
    window.viewport.ren_win.SetSize(1000, 700)  # off-screen: give the view a size
    window.viewport.render()
    return window


def test_zoom_to_a_window(qtbot):
    window = _window(qtbot)
    renderer = window.viewport.renderer
    camera = renderer.GetActiveCamera()
    assert camera.GetParallelProjection()
    width, height = renderer.GetSize()
    scale = camera.GetParallelScale()
    # The upper right quarter of the view, by its world corners.
    corners = []
    for x, y in ((0.5 * width, 0.5 * height), (width, height)):
        renderer.SetDisplayPoint(x, y, 0.5)
        renderer.DisplayToWorld()
        w = renderer.GetWorldPoint()
        corners.append(np.array(w[:3]) / w[3])
    factor = zoom_to_box(renderer, 0.5 * width, 0.5 * height, width, height)
    assert abs(factor - 2.0) < 1e-12
    assert abs(camera.GetParallelScale() - scale / 2) < 1e-9 * scale
    centre = 0.5 * (corners[0] + corners[1])
    assert np.allclose(camera.GetFocalPoint()[:2], centre[:2], atol=1e-9 * scale)


def test_zoom_box_outline_and_clicks(qtbot):
    window = _window(qtbot)
    box = window.zoom_box
    camera = window.viewport.renderer.GetActiveCamera()
    scale = camera.GetParallelScale()
    box.start(100, 100)
    assert not box.update(102, 101)  # too small to be a window yet
    assert not box.actor.GetVisibility()
    assert box.update(300, 250)
    assert box.actor.GetVisibility()
    assert box.finish(300, 250)
    assert not box.active and not box.actor.GetVisibility()
    assert camera.GetParallelScale() < scale
    zoomed = camera.GetParallelScale()
    box.start(10, 10)
    assert not box.finish(12, 11)  # a click: no zoom
    assert camera.GetParallelScale() == zoomed
    assert window.zoom_window_action.isChecked()


def test_forms_scroll_so_docks_resize(qtbot):
    window = _window(qtbot)
    window.set_selection([window.session.entity("c-c-lower").id])
    panel = window.property_panel
    qtbot.waitUntil(lambda: panel.minimumSizeHint().height() > 300, timeout=3000)
    dock = window.property_panel.parentWidget()
    while dock is not None and not isinstance(dock, QScrollArea):
        dock = dock.parentWidget()
    assert isinstance(dock, QScrollArea)
    assert window.property_panel.minimumSizeHint().height() > 300  # a long form...
    assert dock.minimumSizeHint().height() < 150  # ...in a dock that can shrink


def _live_interactor(window):
    """A VTK interactor on the off-screen view with the window's navigation and observers, so
    events run through the real chain (camera style, focus, priorities)."""
    from unittest import mock

    from pyvista.plotting.render_window_interactor import RenderWindowInteractor

    from opengrid.app import main_window

    iren = RenderWindowInteractor(window.viewport)
    iren.interactor.SetRenderWindow(window.viewport.ren_win)
    window.viewport.iren = iren
    main_window.set_navigation_style(iren)
    window._connect_viewport_events()
    held = mock.patch.object(main_window, "_left_button_down", return_value=True)
    return iren.interactor, held


def test_zoom_window_through_the_event_chain(qtbot):
    """The box must finish on release: the camera style grabs the interactor's focus when it
    sees the press, after which nothing else got the moves and the release (the box stayed)."""
    window = _window(qtbot)
    interactor, held = _live_interactor(window)
    box = window.zoom_box
    camera = window.viewport.renderer.GetActiveCamera()
    clicks = []
    window.handle_click = lambda xyz, additive=False: clicks.append(xyz)

    def event(name, x, y, shift=0):
        interactor.SetEventInformation(x, y, 0, shift, "\0", 0, None)
        interactor.InvokeEvent(name)

    scale = camera.GetParallelScale()
    with held:
        event("LeftButtonPressEvent", 100, 100)
        event("MouseMoveEvent", 300, 260)
        assert box.active and box.actor.GetVisibility()
        event("LeftButtonReleaseEvent", 300, 260)
    assert not box.active and not box.actor.GetVisibility()
    assert camera.GetParallelScale() < scale and clicks == []
    event("LeftButtonPressEvent", 400, 400)
    event("LeftButtonReleaseEvent", 401, 400)
    assert len(clicks) == 1 and not box.active  # a click selects
    with held:
        event("LeftButtonPressEvent", 100, 100)
        event("MouseMoveEvent", 250, 220)
    window.cancel_tool()  # Esc
    assert not box.active and not box.actor.GetVisibility()
    with held:
        event("LeftButtonPressEvent", 100, 100)
        event("MouseMoveEvent", 250, 220)
    event("MouseMoveEvent", 260, 230)  # the release was missed: no button any more
    assert not box.active and not box.actor.GetVisibility()
    position = camera.GetPosition()
    with held:
        event("LeftButtonPressEvent", 300, 300, shift=1)
        event("MouseMoveEvent", 350, 330, shift=1)
        event("LeftButtonReleaseEvent", 350, 330, shift=1)
    assert camera.GetPosition() != position and not box.active  # Shift + drag pans


def test_high_order_overlay(qtbot):
    """View -> Show high-order elements: the curved edges follow the geometry; straight-sided
    ones are chords; the overlay follows model changes and goes away when hidden."""
    from test_curved_elements import ring

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = ring(16, og=window.session)  # in the window's own session (it listens to it)
    from opengrid.app.high_order_view import SAMPLES, element_edges

    def wall_samples(mesh):
        """Sampled points of the drawn edges whose ends are on the wall circle."""
        edges = element_edges(mesh)
        actor = window.viewport.renderer.actors["high-order-edges"]
        points = np.asarray(actor.mapper.dataset.points).reshape(len(edges), SAMPLES, 3)
        ends = np.linalg.norm(mesh.points[edges[:, :2], :2], axis=2)
        wall = np.all(np.abs(ends - 1.0) < 1e-9, axis=1)
        assert wall.sum() == 16
        return np.linalg.norm(points[wall][..., :2], axis=-1)

    mesh = window.show_high_order(order=4, curved=True)
    assert mesh is not None and mesh.order == 4
    assert np.abs(wall_samples(mesh) - 1.0).max() < 1e-6  # on the circle between the nodes
    assert "high-order-nodes" in window.viewport.renderer.actors
    assert "order 4" in window.statusBar().currentMessage()
    flat = window.show_high_order(order=4, curved=False)
    assert flat is not None
    assert np.abs(wall_samples(flat) - 1.0).max() > 1e-3  # chords
    before = len(element_edges(flat))
    og.modify("c-cut", dimension=9)  # a model change refreshes the overlay
    actor = window.viewport.renderer.actors["high-order-edges"]
    assert len(np.asarray(actor.mapper.dataset.points)) > before * SAMPLES
    window.hide_high_order()
    assert "high-order-edges" not in window.viewport.renderer.actors
