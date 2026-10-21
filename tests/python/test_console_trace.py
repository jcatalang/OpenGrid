"""The window's Python trace (app/trace.py, the Script tab), the IPython console and their
place in the lower right corner."""

import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from PySide6.QtCore import Qt

import opengrid
from opengrid.app.console import PythonConsole, ipython_available, make_console
from opengrid.app.main_window import MainWindow
from opengrid.app.trace import SessionTrace, new_trace_path


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    yield win
    win.close()


def test_trace_follows_the_session_across_new_and_open(tmp_path):
    og = opengrid.Session()
    trace = SessionTrace(og, tmp_path / "trace.py")
    added = []
    trace.on_append = added.append
    og.create_line((0, 0), (1, 0))
    og.undo()
    og.create_line((0, 0), (2, 0), name="long")
    og.save(tmp_path / "p.ogp")
    trace.note(f"og.save({str(tmp_path / 'p.ogp')!r})")
    og.new()
    og.create_line((0, 0), (0, 3), name="other")
    og.load(tmp_path / "p.ogp")
    og.create_connector("long", 5, name="pts")
    text = trace.text()
    calls = [line.split("(")[0] for line in text.splitlines() if line.startswith("og.")]
    assert calls == [
        "og.create_line",
        "og.undo",
        "og.create_line",
        "og.save",
        "og.new",
        "og.create_line",
        "og.load",
        "og.create_connector",
    ]
    assert "".join(added) in text and len(added) == 8
    # The trace replays: what it ends with is the session.
    replayed = opengrid.replay(text)
    assert replayed.names() == og.names() == ["long", "pts"]


def test_old_traces_are_removed(tmp_path):
    old = tmp_path / "session-20000101-000000-1.py"
    old.write_text("# old\n")
    os.utime(old, (0, 0))
    recent = tmp_path / "session-20990101-000000-2.py"
    recent.write_text("# recent\n")
    path = new_trace_path(tmp_path)
    assert not old.exists() and recent.exists()
    assert path.parent == tmp_path and path.name.startswith("session-")


def test_script_tab_shows_what_the_window_does(window, tmp_path):
    og = window.session
    panel = window.script_panel
    assert panel.text().startswith("# OpenGrid session trace")
    window.scene.world_per_pixel = lambda: 0.001  # no snapping between the clicks
    window.set_tool("line")
    window.handle_click((0.0, 0.0, 0.0))
    window.handle_click((1.0, 0.0, 0.0))
    assert og.names() == ["line-1"]
    assert "og.create_line(p0=(0.0, 0.0, 0.0), p1=(1.0, 0.0, 0.0), name='line-1')" in panel.text()
    window.save_project_as(str(tmp_path / "drawn.ogp"))
    assert f"og.save({str(tmp_path / 'drawn.ogp')!r})" in panel.text()
    saved = window.save_trace(str(tmp_path / "script.py"))
    assert Path(saved).read_text() == window.trace.text()
    # Console commands are traced too, and the console knows where the trace is.
    assert window.console.namespace["trace_path"] == str(window.trace.path)
    window.console.run_line("og.create_line((0, 0), (0, 1), name='typed')")
    assert "name='typed')" in panel.text().splitlines()[-1]
    trace_dir = window.trace.path.parent
    window.close()
    assert not trace_dir.exists()  # a test window's trace goes with it


def test_console_lives_in_the_lower_right_corner(window):
    dock = window.console_dock_widget
    assert window.dockWidgetArea(dock) == Qt.DockWidgetArea.RightDockWidgetArea
    assert (
        window.dockWidgetArea(window.quality_dock_widget) == Qt.DockWidgetArea.RightDockWidgetArea
    )
    tabs = [window.console_tabs.tabText(i) for i in range(window.console_tabs.count())]
    assert tabs == ["Console", "Script", "Log"]
    assert window.console_tabs.currentWidget() is window.console
    # Nothing in the bottom area: the viewport keeps the window's height.
    from PySide6.QtWidgets import QDockWidget

    bottom = [
        d.windowTitle()
        for d in window.findChildren(QDockWidget)
        if window.dockWidgetArea(d) == Qt.DockWidgetArea.BottomDockWidgetArea
    ]
    assert bottom == []


@pytest.mark.skipif(not ipython_available(), reason="qtconsole / ipykernel not installed")
def test_the_console_is_ipython(window):
    console = window.console
    assert console.kind == "ipython"
    assert console.run_line("21 * 2").strip() == "42"
    assert console.run_line("print(type(og).__name__)").strip() == "Session"
    console.run_line("og.create_line((0, 0), (1, 0))")
    assert window.model_tree.entity_names() == ["line-1"]
    assert "ZeroDivisionError" in console.run_line("1 / 0")
    assert console.run_line("get_ipython().__class__.__name__").strip().endswith("Shell'")


def test_plain_console_without_ipython(monkeypatch, qtbot):
    monkeypatch.setenv("OPENGRID_PLAIN_CONSOLE", "1")
    console = make_console({"x": 20}, banner="hello")
    qtbot.addWidget(console)
    assert isinstance(console, PythonConsole) and console.kind == "python"
    assert console.run_line("print(x + 1)").strip() == "21"
    assert "hello" in console.output.toPlainText()
