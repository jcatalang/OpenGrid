import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from opengrid.app.main_window import MainWindow
from opengrid.app.model_tree import CATEGORIES


@pytest.fixture
def window(qtbot):
    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    yield win
    win.close()


def test_model_tree_categories(window):
    assert window.model_tree.category_names() == list(CATEGORIES)
    assert CATEGORIES == ("Database", "Connectors", "Domains", "Blocks", "Patches", "Zones")


def test_console_evaluates(window):
    output = window.console.run_line("print(1 + 1)")
    assert output.strip() == "2"
    window.console.run_line("x = 21")
    assert window.console.run_line("print(x * 2)").strip() == "42"


def test_console_has_opengrid(window):
    assert "opengrid" in window.console.namespace
    assert window.console.run_line("print(opengrid.core_version())").strip() == "0.0.1"


def test_window_renders(window, tmp_path):
    window.viewport.render()
    image = window.grab()
    assert not image.isNull()
    assert image.width() > 0
