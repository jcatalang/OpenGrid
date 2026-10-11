"""Main application window: viewport, model tree, property panel, log and Python console."""

import logging

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import QDockWidget, QFormLayout, QLabel, QMainWindow, QWidget

import opengrid
from opengrid.app.console import PythonConsole
from opengrid.app.log_panel import LogPanel
from opengrid.app.model_tree import ModelTree

log = logging.getLogger("opengrid")


class PropertyPanel(QWidget):
    """Shows editable properties of the selected entity (filled from Phase 1)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QFormLayout(self)
        self.placeholder = QLabel("No selection")
        self.placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addRow(self.placeholder)


class MainWindow(QMainWindow):
    def __init__(self, parent=None, off_screen: bool = False):
        super().__init__(parent)
        self.setWindowTitle(f"OpenGrid {opengrid.__version__}")
        self.resize(1400, 900)

        self.viewport = self._create_viewport(off_screen)
        self.setCentralWidget(self.viewport)

        self.model_tree = ModelTree()
        self.property_panel = PropertyPanel()
        self.log_panel = LogPanel()
        self.console = PythonConsole(namespace={"opengrid": opengrid, "window": self})

        self._add_dock("Model", self.model_tree, Qt.DockWidgetArea.LeftDockWidgetArea)
        self._add_dock("Properties", self.property_panel, Qt.DockWidgetArea.LeftDockWidgetArea)
        console_dock = self._add_dock(
            "Console", self.console, Qt.DockWidgetArea.BottomDockWidgetArea
        )
        log_dock = self._add_dock("Log", self.log_panel, Qt.DockWidgetArea.BottomDockWidgetArea)
        self.tabifyDockWidget(log_dock, console_dock)
        console_dock.raise_()

        self._create_menus()

        logging.getLogger("opengrid").addHandler(self.log_panel.handler)
        logging.getLogger("opengrid").setLevel(logging.INFO)
        log.info("OpenGrid %s (core %s) started", opengrid.__version__, opengrid.core_version())

    def _create_viewport(self, off_screen: bool):
        from pyvistaqt import QtInteractor

        viewport = QtInteractor(self, off_screen=off_screen)
        viewport.set_background("#2b3036", top="#59636e")
        viewport.add_axes()
        viewport.enable_parallel_projection()
        return viewport

    def _add_dock(self, title: str, widget: QWidget, area: Qt.DockWidgetArea) -> QDockWidget:
        dock = QDockWidget(title, self)
        dock.setObjectName(f"dock_{title.lower()}")
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        return dock

    def _create_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        quit_action = QAction("&Quit", self)
        quit_action.setShortcut(QKeySequence.StandardKey.Quit)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        edit_menu = self.menuBar().addMenu("&Edit")
        # Undo/redo are wired to the CommandStack in Phase 1.
        self.undo_action = QAction("&Undo", self)
        self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.setEnabled(False)
        self.redo_action = QAction("&Redo", self)
        self.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.redo_action.setEnabled(False)
        edit_menu.addAction(self.undo_action)
        edit_menu.addAction(self.redo_action)

        view_menu = self.menuBar().addMenu("&View")
        reset_action = QAction("&Reset camera", self)
        reset_action.triggered.connect(self.viewport.reset_camera)
        view_menu.addAction(reset_action)

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt override)
        logging.getLogger("opengrid").removeHandler(self.log_panel.handler)
        self.viewport.close()
        super().closeEvent(event)
