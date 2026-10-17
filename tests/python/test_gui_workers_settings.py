"""Long operations in a worker thread (app/worker.py) and the settings that persist between
runs (window layout, view options, recent files, the last directory)."""

import threading
import time
from pathlib import Path

import pytest
from PySide6.QtCore import QSettings, QTimer

from opengrid.app.main_window import MainWindow
from opengrid.app.worker import run_blocking

PROJECT = Path(__file__).resolve().parents[2] / "examples" / "output" / "04_airfoil_cgrid.ogp"


def test_run_blocking_keeps_the_gui_thread_alive(qtbot):
    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    ticks = []
    timer = QTimer(window)
    timer.setInterval(10)
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start()
    result = run_blocking(window, "Sleeping…", lambda: time.sleep(0.3) or 42, grace_ms=50)
    timer.stop()
    assert result == 42
    assert len(ticks) >= 5  # events were processed meanwhile
    with pytest.raises(ValueError, match="boom"):
        run_blocking(window, "Failing…", lambda: (_ for _ in ()).throw(ValueError("boom")))
    # Cancel calls back; the function sees it and returns.
    stop = threading.Event()

    def wait_for_cancel():
        stop.wait(5.0)
        return "stopped"

    QTimer.singleShot(
        150,
        lambda: [
            d.cancel()
            for d in window.findChildren(__import__("PySide6.QtWidgets").QtWidgets.QProgressDialog)
        ],
    )
    assert (
        run_blocking(window, "Cancellable…", wait_for_cancel, cancel=stop.set, grace_ms=20)
        == "stopped"
    )
    assert stop.is_set()
    window.close()


def test_long_operations_refresh_once_and_report_errors(qtbot):
    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    refreshes = []
    original = window._on_model_changed
    window._on_model_changed = lambda change: (refreshes.append(change), original(change))

    def several_commands():
        og.create_line((0, 0), (1, 0))
        og.create_line((1, 0), (1, 1))
        return "done"

    assert window._run_long("Drawing…", several_commands) == "done"
    assert len(refreshes) == 1 and sorted(refreshes[0].ids()) == [1, 2]
    assert window.scene.layer("curves") is not None
    assert window._run_long("Failing…", lambda: og.create_line((0, 0), (0, 0))) is None
    assert "coincide" in window.statusBar().currentMessage()
    window.close()


def test_settings_persist_between_windows(qtbot, tmp_path):
    settings = QSettings(str(tmp_path / "opengrid.ini"), QSettings.Format.IniFormat)
    window = MainWindow(off_screen=True, settings=settings)
    qtbot.addWidget(window)
    assert window.recent_menu.actions() == [] and not window.recent_menu.isEnabled()
    window.dimension_box.setValue(33)
    window.nodes_action.setChecked(False)
    window.zoom_window_action.setChecked(False)
    assert window.open_project(str(PROJECT))
    assert window._last_dir == "" or window._last_dir  # set by the dialogs only
    window._last_dir = str(PROJECT.parent)
    window.close()
    again = MainWindow(
        off_screen=True,
        settings=QSettings(str(tmp_path / "opengrid.ini"), QSettings.Format.IniFormat),
    )
    qtbot.addWidget(again)
    assert again.dimension_box.value() == 33
    assert not again.nodes_action.isChecked() and not again.zoom_window_action.isChecked()
    assert again._last_dir == str(PROJECT.parent)
    recent = [a.text() for a in again.recent_menu.actions()]
    assert recent == [PROJECT.name] and again.recent_menu.isEnabled()
    again.recent_menu.actions()[0].trigger()
    assert again.session.project_path == PROJECT
    again.close()
