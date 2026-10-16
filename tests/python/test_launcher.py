"""The installed macOS launcher (packaging/macos/install_launcher.sh) must match its scripts:
skipped where no launcher of this checkout is installed."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "check_launcher", ROOT / "tools" / "check_launcher.py"
)
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


def test_installed_launcher_is_current():
    issues = check.problems(Path("/Applications/OpenGrid.app"))
    if issues is None:
        pytest.skip("no OpenGrid.app launcher of this checkout installed")
    assert not issues, "; ".join(issues) + ": rerun packaging/macos/install_launcher.sh"


def test_file_open_events_open_projects(qtbot):
    from PySide6.QtCore import QCoreApplication
    from PySide6.QtGui import QFileOpenEvent

    from opengrid.__main__ import _file_open_filter
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    opener = _file_open_filter(window)
    app = QCoreApplication.instance()
    app.installEventFilter(opener)
    try:
        project = ROOT / "examples" / "output" / "04_airfoil_cgrid.ogp"
        QCoreApplication.sendEvent(app, QFileOpenEvent(str(project)))
        assert window.session.project_path == project
        assert "04_airfoil_cgrid.ogp" in window.windowTitle()
    finally:
        app.removeEventFilter(opener)


def test_update_notice(qtbot):
    from PySide6.QtWidgets import QLabel

    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    window._check_for_update()  # unchanged code: nothing shown
    assert all("restart" not in label.text() for label in window.statusBar().findChildren(QLabel))
    window._code_stamp = 0.0  # as if the code changed after the start
    window._check_for_update()
    texts = [label.text() for label in window.statusBar().findChildren(QLabel)]
    assert any("restart it to use the new version" in text for text in texts)
