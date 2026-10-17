"""The reference grid (app/reference_grid.py): a step that follows the zoom, axes, labels in
views along z, none in other views, persisted with the view options."""

from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings

from opengrid.app.main_window import MainWindow
from opengrid.app.reference_grid import grid_step

PROJECT = Path(__file__).resolve().parents[2] / "examples" / "output" / "04_airfoil_cgrid.ogp"


def test_grid_step_is_a_round_fraction_of_the_view():
    assert grid_step(10.0) == (2.0, 4)
    assert grid_step(1.0) == (0.2, 4)
    assert grid_step(0.37) == (0.05, 5)
    assert grid_step(50.0) == (10.0, 10)
    assert grid_step(0.0) == (1.0, 10)


def test_reference_grid_follows_the_camera_and_the_view(qtbot, tmp_path):
    settings = QSettings(str(tmp_path / "grid.ini"), QSettings.Format.IniFormat)
    window = MainWindow(off_screen=True, settings=settings)
    qtbot.addWidget(window)
    window.resize(1200, 800)
    window.open_project(str(PROJECT))
    window.show()
    window.viewport.ren_win.SetSize(1000, 700)
    window.scene.fit()
    grid = window.scene.grid
    assert not grid.enabled and grid._actors == []
    window.grid_action.setChecked(True)
    assert grid.enabled and len(grid._actors) >= 2  # minor and major lines (and axes)
    assert len(grid._labels) > 4  # a view along z: coordinates along the edges
    names = {a.name for a in grid._actors}
    assert {"reference-grid-xy-minor", "reference-grid-xy-major"} <= names
    assert not any("-xz-" in n or "-yz-" in n for n in names)  # one plane in an axis view
    centre, extent = grid.visible_extent()
    step = grid.step
    assert grid_step(extent)[0] == step
    camera = window.viewport.renderer.GetActiveCamera()
    camera.SetParallelScale(camera.GetParallelScale() / 10.0)
    grid.update()
    assert abs(grid.step - step / 10.0) < 1e-12 * step
    # The lines cover the view around its centre.
    major = next(a for a in grid._actors if a.name == "reference-grid-xy-major")
    bounds = np.array(major.GetBounds())
    assert bounds[0] <= centre[0] <= bounds[1] and bounds[2] <= centre[1] <= bounds[3]
    # A view along y shows the XZ plane with its edge labels; along x the YZ plane.
    window.view_actions["xz"].trigger()
    names = {a.name for a in grid._actors}
    assert "reference-grid-xz-minor" in names and not any("-xy-" in n for n in names)
    assert len(grid._labels) > 4
    window.view_actions["yz"].trigger()
    assert any(a.name == "reference-grid-yz-minor" for a in grid._actors)
    # A tilted view: the three planes (the facing one with minor lines) and labels along
    # the axes that face the camera.
    window.view_actions["iso"].trigger()
    names = {a.name for a in grid._actors}
    assert {
        "reference-grid-xy-major",
        "reference-grid-xz-major",
        "reference-grid-yz-major",
    } <= names
    assert sum(1 for n in names if n.endswith("-minor")) == 1
    assert grid._labels and all(hasattr(label, "GetPosition") for label in grid._labels[:-1])
    assert any("axis-x" in n for n in names)  # the planes meet on the world axes here
    window.view_actions["xy"].trigger()
    assert grid._labels
    window.grid_action.setChecked(False)
    assert grid._actors == [] and grid._labels == []
    window.grid_action.setChecked(True)
    window.close()
    again = MainWindow(
        off_screen=True, settings=QSettings(str(tmp_path / "grid.ini"), QSettings.Format.IniFormat)
    )
    qtbot.addWidget(again)
    assert again.grid_action.isChecked() and again.scene.grid.enabled
    again.close()
