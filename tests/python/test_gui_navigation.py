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
