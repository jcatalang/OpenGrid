"""Cell quality on the mesh (the Quality panel): every cell coloured by the metric, or only the
cells worse than a threshold filled, the rest as lines; in 3D also the cells inside blocks and
volumes."""

import numpy as np
import pytest

from opengrid.app.quality_dock import THRESHOLDS, worse_than


@pytest.fixture
def window(qtbot):
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    yield window
    window.close()


def worse_cells(window, name):
    return window.scene.worse_cells(window.session.entity(name).id)


def test_2d_cells_worse_than_a_threshold(window):
    og = window.session
    domain = og.unstructured_domain_from_curves(og.create_rectangle(width=2, height=1), 21)
    dock = window.quality_dock
    dock.metric.setCurrentIndex(dock.metric.findData("min_angle"))
    assert dock.threshold.value() == THRESHOLDS["min_angle"]
    # The panel's switch and the View menu entry are one setting.
    dock.colour.setChecked(True)
    assert window.quality_action.isChecked() and window.scene.quality_metric == "min_angle"
    window.quality_action.setChecked(False)
    assert not dock.colour.isChecked() and window.scene.quality_metric is None
    # Only the cells worse than 45 degrees: they are filled, the domain is drawn as lines.
    dock.threshold.setValue(45.0)
    dock.only_worse.setChecked(True)  # (turns the colouring on)
    assert window.quality_action.isChecked() and window.scene.quality_threshold == 45.0
    values = np.asarray(og.domain_quality(domain)["min_angle"]).ravel()
    expected = int(np.count_nonzero(worse_than(values, "min_angle", 45.0)))
    assert 0 < expected < len(values)
    assert worse_cells(window, domain).n_cells == expected == window.scene.worse_cell_count()
    assert f"{expected} worse than 45" in dock.stats.text()
    # Selecting it recolours its lines only.
    window.set_selection([og.entity(domain).id])
    # A metric change offers that metric's threshold.
    dock.metric.setCurrentIndex(dock.metric.findData("scaled_jacobian"))
    assert window.scene.quality_threshold == THRESHOLDS["scaled_jacobian"]
    dock.only_worse.setChecked(False)
    assert window.scene.quality_threshold is None and worse_cells(window, domain) is None


def test_3d_cells_inside_a_volume(window, qtbot):
    og = window.session
    base = og.unstructured_domain_from_curves(og.create_rectangle(width=2, height=1), 11)
    volume = og.create_extruded_volume(base, "translate", vector=(0, 0, 1), dimension=6)
    dock = window.quality_dock
    dock.metric.setCurrentIndex(dock.metric.findData("min_angle"))
    dock.threshold.setValue(50.0)
    dock.only_worse.setChecked(True)
    values = np.asarray(og.volume_quality(volume)["min_angle"]).ravel()
    expected = int(np.count_nonzero(worse_than(values, "min_angle", 50.0)))
    assert 0 < expected < len(values)
    cells = worse_cells(window, volume)
    assert cells.n_cells == expected
    # Cells inside, not only those on its faces.
    centres = np.asarray(cells.cell_centers().points)
    inside = (centres[:, 2] > 0.25) & (centres[:, 2] < 0.75)
    assert np.any(inside)
    # The block's faces are still picked (its outline carries them).
    assert window.scene._face_members[og.entity(volume).id]
