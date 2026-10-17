"""The camera stays where it is while the model is edited: a handle drag must not refit or
flip the view (the user: every control-point move reset the camera and mirrored the view)."""

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.main_window import MainWindow


def camera_state(window):
    camera = window.scene.plotter.camera
    return (tuple(camera.position), tuple(camera.focal_point), camera.parallel_scale)


def test_handle_drag_keeps_the_camera(qtbot):
    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    og.create_nurbs(
        3, [0, 0, 0, 0, 0.5, 1, 1, 1, 1], [(0, 0), (1, 1), (2, -1), (3, 1), (4, 0)], name="n"
    )
    og.create_connector("n", 21, name="cn")
    window.scene.world_per_pixel = lambda: 0.002
    camera = window.scene.plotter.camera
    camera.position = (2.0, 0.0, 10.0)
    camera.focal_point = (2.0, 0.0, 0.0)
    camera.parallel_scale = 1.5
    before = camera_state(window)
    window.set_selection([og.entity("cn").id])
    assert window.handle_press((2.0, -1.0, 0.0))
    window.handle_drag((2.2, -0.4, 0.0))
    assert window.scene.has_geometry()  # the whole model is being previewed: still geometry
    assert window.handle_release((2.5, 0.5, 0.0)) == ["n"]
    assert og.entity("n").params["control_points"][2] == (2.5, 0.5, 0.0)
    assert camera_state(window) == before
    # A fit never flips to the far side, however often it runs.
    window.scene.fit()
    first = camera_state(window)
    assert first[0][2] > first[1][2]  # looking down from +z
    window.scene.fit()
    assert camera_state(window) == first
    # The view key does toggle, on purpose.
    window.scene.view_preset("xy")
    flipped = camera_state(window)
    assert flipped[0][2] < flipped[1][2]
    window.close()
