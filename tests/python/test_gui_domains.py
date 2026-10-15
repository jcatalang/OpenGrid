import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.main_window import MainWindow


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    og = win.session
    og.create_circle((0, 0), 1.0, name="inner")
    og.create_circle((0, 0), 6.0, name="outer")
    og.create_line((1, 0), (6, 0), name="cut")
    og.create_connector("inner", 41, name="c-inner")
    og.create_connector("outer", 41, name="c-outer")
    og.create_connector("cut", 21, {"type": "tanh", "begin_spacing": 0.02}, name="c-cut")
    yield win
    win.close()


def _ids(window, *names):
    return [window.session.entity(n).id for n in names]


def test_domain_from_selected_connectors(window):
    window.set_selection(_ids(window, "c-inner", "c-outer", "c-cut"))
    name = window.create_domain_on_selection()
    assert name == "domain-1"
    assert window.model_tree.entity_names("Domains") == ["domain-1"]
    assert window.model_tree.category_item("Domains").child(0).text(1) == "41 x 21"
    assert window.selected_names() == ["domain-1"]
    assert "800 cells" in window.quality_dock.stats.text()


def test_quality_colouring_histogram_and_worst_cell(window):
    window.session.create_structured_domain(connectors=["c-inner", "c-outer", "c-cut"], name="ring")
    window.set_selection([])
    window.quality_dock.metric.setCurrentIndex(
        list(window.quality_dock.metric.itemData(i) for i in range(7)).index("aspect_ratio")
    )
    assert "all domains" in window.quality_dock.stats.text()
    window.quality_action.setChecked(True)
    assert window.scene.quality_metric == "aspect_ratio"
    name, i, j = window.show_worst_cell()
    assert name == "ring"
    q = window.session.domain_quality("ring")["aspect_ratio"]
    assert q[j, i] == q.max()
    window.quality_action.setChecked(False)
    assert window.scene.quality_metric is None


def test_generator_edit_and_threaded_smoothing(window, qtbot):
    window.session.create_structured_domain(connectors=["c-inner", "c-outer", "c-cut"], name="ring")
    window.set_selection(_ids(window, "ring"))
    panel = window.property_panel
    _kind, editor = panel._editors["generator"]
    editor.elliptic_enabled.setChecked(True)
    editor.control.setCurrentText("laplace")
    editor.wall_boxes[0].setChecked(True)
    editor.wall_spacing[0].setText("0.01")
    assert panel.apply()
    generator = window.session.entity("ring").params["generator"]
    assert generator["elliptic"]["walls"][0]["spacing"] == 0.01
    assert "⚠" in window.model_tree.category_item("Domains").child(0).text(1)  # pending
    # Unchanged panel -> no changes.
    window.set_selection(_ids(window, "ring"))
    assert window.property_panel.read_changes() == (None, {})
    assert window.run_smoother("ring")
    window.smoother_dock.wait()
    qtbot.waitUntil(lambda: not window.smoother_dock.running, timeout=30000)
    assert "converged" in window.session.domain_status("ring")
    assert "⚠" not in window.model_tree.category_item("Domains").child(0).text(1)
    grid = window.session.domain_grid("ring", compute=False)
    spacing = np.linalg.norm(grid[1] - grid[0], axis=1)
    assert np.allclose(spacing[1:-1], 0.01, rtol=0.05)


def test_cancel_smoothing(window, qtbot):
    window.session.create_structured_domain(
        connectors=["c-inner", "c-outer", "c-cut"],
        generator={
            # Plain relaxation: the accelerated one reaches round-off before it can be cancelled.
            "elliptic": {
                "control": "laplace",
                "max_iterations": 10**8,
                "tolerance": 1e-300,
                "acceleration": 0,
            }
        },
        name="ring",
    )
    assert window.run_smoother("ring")
    qtbot.waitUntil(lambda: window.smoother_dock._progress.iteration > 3, timeout=10000)
    window.smoother_dock.cancel()
    qtbot.waitUntil(lambda: not window.smoother_dock.running, timeout=10000)
    assert "cancelled" in window.smoother_dock.status.text()
    assert window.session.domain_status("ring").startswith("smoothing pending")


def test_rotate_button_and_exports(window, tmp_path):
    name = window.session.create_structured_domain(connectors=["c-inner", "c-outer", "c-cut"])
    window.set_selection(_ids(window, name))
    window.property_panel.rotate_button.click()
    assert window.session.domain_dimensions(name) == (21, 41)
    summary = window.export_vtkhdf(str(tmp_path / "ring.vtkhdf"))
    assert summary["cells"] == 800
    assert summary["points"] == 41 * 21 - 21  # the cut column is merged
    assert window.export_plot3d(str(tmp_path / "ring.xyz"))
    assert (tmp_path / "ring.xyz").stat().st_size > 0
