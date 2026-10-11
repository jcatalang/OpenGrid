"""Main application window: viewport, model tree, property panel, log and Python console."""

import logging
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QMainWindow,
    QMessageBox,
    QToolBar,
    QWidget,
)

import opengrid
from opengrid.api import Session
from opengrid.app.console import PythonConsole
from opengrid.app.log_panel import LogPanel
from opengrid.app.model_tree import ModelTree
from opengrid.app.property_panel import PropertyPanel
from opengrid.app.scene import Scene, pick_entity
from opengrid.app.tools import TOOLS, Tool, snap, snap_targets
from opengrid.validation import ValidationError

log = logging.getLogger("opengrid")

SNAP_PIXELS = 10.0
PICK_PIXELS = 8.0
CLICK_MOVE_PIXELS = 4.0
PROJECT_FILTER = "OpenGrid project (*.ogp)"


class MainWindow(QMainWindow):
    def __init__(self, parent=None, off_screen: bool = False, session: Session | None = None):
        super().__init__(parent)
        if session is None:
            session = Session()
        self.session = session
        self._off_screen = off_screen
        self.selected_ids: list[int] = []
        self.tool: Tool = TOOLS["select"](session)
        self._saved_journal_length = 0
        self._press_position: tuple[int, int] | None = None
        self.resize(1400, 900)

        self.viewport = self._create_viewport(off_screen)
        self.setCentralWidget(self.viewport)
        self.scene = Scene(self.viewport, session)

        self.model_tree = ModelTree()
        self.property_panel = PropertyPanel(session)
        self.log_panel = LogPanel()
        self.console = PythonConsole(
            namespace={"opengrid": opengrid, "og": session, "window": self}
        )

        self._add_dock("Model", self.model_tree, Qt.DockWidgetArea.LeftDockWidgetArea)
        self._add_dock("Properties", self.property_panel, Qt.DockWidgetArea.LeftDockWidgetArea)
        console_dock = self._add_dock(
            "Console", self.console, Qt.DockWidgetArea.BottomDockWidgetArea
        )
        log_dock = self._add_dock("Log", self.log_panel, Qt.DockWidgetArea.BottomDockWidgetArea)
        self.tabifyDockWidget(log_dock, console_dock)
        console_dock.raise_()

        self._create_actions()
        self._create_menus()
        self._create_toolbars()
        self._connect_viewport_events()

        self.model_tree.selection_changed.connect(self.set_selection)
        self.property_panel.message.connect(self._on_panel_message)
        session.commands.add_listener(self._on_model_changed)

        logging.getLogger("opengrid").addHandler(self.log_panel.handler)
        logging.getLogger("opengrid").setLevel(logging.INFO)
        self._refresh_all()
        log.info("OpenGrid %s (core %s) started", opengrid.__version__, opengrid.core_version())

    # ------------------------------------------------------------------ setup
    def _create_viewport(self, off_screen: bool):
        from pyvistaqt import QtInteractor

        viewport = QtInteractor(self, off_screen=off_screen)
        viewport.set_background("#2b3036", top="#59636e")
        viewport.add_axes()
        viewport.enable_parallel_projection()
        viewport.view_xy()
        if viewport.iren is not None:
            # 2D navigation: left-drag pans, wheel / two-finger scroll zooms, right-drag zooms.
            viewport.enable_2d_style()
        return viewport

    def _add_dock(self, title: str, widget: QWidget, area: Qt.DockWidgetArea) -> QDockWidget:
        dock = QDockWidget(title, self)
        dock.setObjectName(f"dock_{title.lower()}")
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        return dock

    def _action(self, text: str, slot, shortcut=None, tip: str = "") -> QAction:
        action = QAction(text, self)
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip:
            action.setStatusTip(tip)
        action.triggered.connect(slot)
        return action

    def _create_actions(self) -> None:
        self.new_action = self._action("&New", self.new_project, QKeySequence.StandardKey.New)
        self.open_action = self._action("&Open…", self.open_project, QKeySequence.StandardKey.Open)
        self.save_action = self._action("&Save", self.save_project, QKeySequence.StandardKey.Save)
        self.save_as_action = self._action(
            "Save &As…", self.save_project_as, QKeySequence.StandardKey.SaveAs
        )
        self.import_dat_action = self._action("Import &airfoil (.dat)…", self.import_airfoil)
        self.import_csv_action = self._action("Import &points (CSV)…", self.import_csv)
        self.save_journal_action = self._action("Save &journal…", self.save_journal)
        self.quit_action = self._action("&Quit", self.close, QKeySequence.StandardKey.Quit)

        self.undo_action = self._action("&Undo", self.undo, QKeySequence.StandardKey.Undo)
        self.redo_action = self._action("&Redo", self.redo, QKeySequence.StandardKey.Redo)
        self.delete_action = self._action(
            "&Delete", self.delete_selected, QKeySequence.StandardKey.Delete
        )
        self.delete_action.setShortcuts(
            [QKeySequence(QKeySequence.StandardKey.Delete), QKeySequence(Qt.Key.Key_Backspace)]
        )
        self.delete_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.viewport.addAction(self.delete_action)
        self.model_tree.addAction(self.delete_action)

        self.finish_action = self._action(
            "&Finish", self.finish_tool, QKeySequence(Qt.Key.Key_Return), "Finish the current tool"
        )
        self.finish_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.viewport.addAction(self.finish_action)
        self.cancel_action = self._action(
            "&Cancel", self.cancel_tool, QKeySequence(Qt.Key.Key_Escape), "Cancel / clear selection"
        )
        self.cancel_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.viewport.addAction(self.cancel_action)

        self.fit_action = self._action("&Fit view", self.scene.fit, QKeySequence("F"))
        self.fit_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.viewport.addAction(self.fit_action)

        self.tool_group = QActionGroup(self)
        self.tool_actions: dict[str, QAction] = {}
        labels = {
            "select": "Select",
            "point": "Point",
            "line": "Line",
            "arc": "Arc (3 points)",
            "polyline": "Polyline",
            "spline": "Spline",
        }
        for name, label in labels.items():
            action = QAction(label, self)
            action.setCheckable(True)
            action.triggered.connect(lambda _checked=False, n=name: self.set_tool(n))
            self.tool_group.addAction(action)
            self.tool_actions[name] = action
        self.tool_actions["select"].setChecked(True)

    def _create_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        for action in (self.new_action, self.open_action, self.save_action, self.save_as_action):
            file_menu.addAction(action)
        file_menu.addSeparator()
        file_menu.addAction(self.import_dat_action)
        file_menu.addAction(self.import_csv_action)
        file_menu.addSeparator()
        file_menu.addAction(self.save_journal_action)
        file_menu.addSeparator()
        file_menu.addAction(self.quit_action)

        edit_menu = self.menuBar().addMenu("&Edit")
        edit_menu.addAction(self.undo_action)
        edit_menu.addAction(self.redo_action)
        edit_menu.addSeparator()
        edit_menu.addAction(self.delete_action)

        create_menu = self.menuBar().addMenu("&Create")
        for name in ("point", "line", "arc", "polyline", "spline"):
            create_menu.addAction(self.tool_actions[name])
        create_menu.addSeparator()
        create_menu.addAction(self.finish_action)
        create_menu.addAction(self.cancel_action)

        view_menu = self.menuBar().addMenu("&View")
        view_menu.addAction(self.fit_action)

    def _create_toolbars(self) -> None:
        toolbar = QToolBar("Tools", self)
        toolbar.setObjectName("toolbar_tools")
        for action in self.tool_actions.values():
            toolbar.addAction(action)
        toolbar.addSeparator()
        toolbar.addAction(self.finish_action)
        toolbar.addAction(self.cancel_action)
        toolbar.addSeparator()
        toolbar.addAction(self.undo_action)
        toolbar.addAction(self.redo_action)
        toolbar.addAction(self.fit_action)
        self.addToolBar(toolbar)

    def _connect_viewport_events(self) -> None:
        iren = self.viewport.iren
        if iren is None:
            return  # off-screen (tests): handlers are driven directly in world space
        iren.add_observer("LeftButtonPressEvent", self._on_vtk_press)
        iren.add_observer("LeftButtonReleaseEvent", self._on_vtk_release)
        iren.add_observer("MouseMoveEvent", self._on_vtk_move)

    # ------------------------------------------------------- viewport events
    def _event_position(self) -> tuple[int, int]:
        x, y = self.viewport.iren.get_event_position()
        return int(x), int(y)

    def _on_vtk_press(self, *_args) -> None:
        self._press_position = self._event_position()

    def _on_vtk_release(self, *_args) -> None:
        if self._press_position is None:
            return
        x, y = self._event_position()
        x0, y0 = self._press_position
        self._press_position = None
        if abs(x - x0) + abs(y - y0) > CLICK_MOVE_PIXELS:
            return  # a drag (pan), not a click
        interactor = self.viewport.iren.interactor
        additive = bool(interactor.GetShiftKey() or interactor.GetControlKey())
        self.handle_click(self.scene.display_to_world(x, y), additive)

    def _on_vtk_move(self, *_args) -> None:
        x, y = self._event_position()
        self.handle_move(self.scene.display_to_world(x, y))

    # ---------------------------------------------------- world-space handlers
    def snap_point(self, xyz):
        tolerance = SNAP_PIXELS * self.scene.world_per_pixel()
        return snap(xyz, snap_targets(self.session), tolerance)

    def handle_click(self, xyz, additive: bool = False) -> None:
        """A click at world position ``xyz`` (from the viewport, or from tests)."""
        if self.tool.name == "select":
            tolerance = PICK_PIXELS * self.scene.world_per_pixel()
            name = pick_entity(self.session, xyz, tolerance)
            self._select_picked(name, additive)
            return
        snapped, target = self.snap_point(xyz)
        if target is not None:
            log.info("snapped to %s", target.label)
        created = self._run(lambda: self.tool.click(snapped))
        if created:
            log.info("created %s", created)
        self._update_tool_feedback(snapped)

    def handle_move(self, xyz) -> None:
        if self.tool.name == "select":
            return
        snapped, target = self.snap_point(xyz)
        if target is None:
            self.scene.show_snap(None)
        else:
            self.scene.show_snap(snapped)
        self.scene.show_preview(self.tool.preview(snapped))

    def _select_picked(self, name: str | None, additive: bool) -> None:
        if name is None:
            if not additive:
                self.set_selection([])
            return
        entity_id = self.session.store.get(name).id
        if additive:
            ids = list(self.selected_ids)
            if entity_id in ids:
                ids.remove(entity_id)
            else:
                ids.append(entity_id)
            self.set_selection(sorted(ids))
        else:
            self.set_selection([entity_id])

    # ------------------------------------------------------------------ tools
    def set_tool(self, name: str) -> None:
        self.tool.cancel()
        self.tool = TOOLS[name](self.session)
        self.tool_actions[name].setChecked(True)
        self.scene.show_preview([])
        self.scene.show_snap(None)
        self._update_tool_feedback(None)

    def finish_tool(self) -> None:
        created = self._run(self.tool.finish)
        if created:
            log.info("created %s", created)
        self._update_tool_feedback(None)

    def cancel_tool(self) -> None:
        if self.tool.points:
            self.tool.cancel()
            self.scene.show_preview([])
        elif self.tool.name != "select":
            self.set_tool("select")
        else:
            self.set_selection([])
        self._update_tool_feedback(None)

    def _update_tool_feedback(self, cursor) -> None:
        self.finish_action.setEnabled(self.tool.can_finish())
        prompt = self.tool.prompt()
        if prompt:
            self.statusBar().showMessage(prompt)
        else:
            self.statusBar().clearMessage()
        if cursor is not None:
            self.scene.show_preview(self.tool.preview(cursor))
        elif not self.tool.points:
            self.scene.show_preview([])

    # ------------------------------------------------------- model / selection
    def set_selection(self, ids) -> None:
        ids = [i for i in ids if i in self.session.store]
        self.selected_ids = list(ids)
        self.scene.set_selection(ids)
        self.model_tree.select_ids(ids)
        if len(ids) == 1:
            self.property_panel.show_entity(ids[0])
        else:
            self.property_panel.show_entity(None)
        self.delete_action.setEnabled(bool(ids))

    def selected_names(self) -> list[str]:
        return [self.session.store.get(i).name for i in self.selected_ids]

    def _on_model_changed(self, change) -> None:
        had_geometry = len(self.scene._actors) > 0
        if change is None:
            self.scene.rebuild()
        else:
            self.scene.update_entities(change.ids())
        self.model_tree.populate(self.session.store)
        self.set_selection(self.selected_ids)
        self._update_edit_actions()
        if not had_geometry and len(self.session.store) > 0:
            self.scene.fit()

    def _refresh_all(self) -> None:
        self.scene.rebuild()
        self.model_tree.populate(self.session.store)
        self.set_selection([])
        self._update_edit_actions()
        self._update_tool_feedback(None)

    def _update_edit_actions(self) -> None:
        commands = self.session.commands
        self.undo_action.setEnabled(commands.can_undo)
        self.redo_action.setEnabled(commands.can_redo)
        if commands.can_undo:
            self.undo_action.setText(f"&Undo {commands.undo_label()}")
        else:
            self.undo_action.setText("&Undo")
        if commands.can_redo:
            self.redo_action.setText(f"&Redo {commands.redo_label()}")
        else:
            self.redo_action.setText("&Redo")
        self._update_title()

    def _update_title(self) -> None:
        name = "untitled"
        if self.session.project_path is not None:
            name = self.session.project_path.name
        dirty = ""
        if self.is_modified():
            dirty = " *"
        self.setWindowTitle(f"{name}{dirty} — OpenGrid {opengrid.__version__}")

    def is_modified(self) -> bool:
        return len(self.session.commands.journal.lines) != self._saved_journal_length

    def _on_panel_message(self, text: str, is_error: bool) -> None:
        if is_error:
            log.warning("%s", text)
        else:
            log.info("%s", text)
        self.statusBar().showMessage(text.splitlines()[0], 5000)

    # ------------------------------------------------------------ edit actions
    def undo(self) -> None:
        if self.session.commands.can_undo:
            self.session.undo()

    def redo(self) -> None:
        if self.session.commands.can_redo:
            self.session.redo()

    def delete_selected(self) -> None:
        names = self.selected_names()
        if names:
            self._run(lambda: self.session.delete(names))

    def _run(self, function):
        """Run a model operation; report validation errors instead of raising."""
        try:
            return function()
        except (ValidationError, KeyError, ValueError) as error:
            text = str(error)
            log.warning("%s", text)
            self.statusBar().showMessage(text.splitlines()[0], 8000)
            return None

    # ------------------------------------------------------------ file actions
    def _confirm_discard(self) -> bool:
        if not self.is_modified() or len(self.session.store) == 0:
            return True
        answer = QMessageBox.question(
            self,
            "Unsaved changes",
            "Discard unsaved changes?",
            QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
        )
        return answer == QMessageBox.StandardButton.Discard

    def new_project(self, confirm: bool = True) -> None:
        if confirm and not self._confirm_discard():
            return
        self.session.new()
        self._saved_journal_length = len(self.session.commands.journal.lines)
        self._refresh_all()

    def open_project(self, path: str | None = None) -> bool:
        if path is None:
            if not self._confirm_discard():
                return False
            path, _filter = QFileDialog.getOpenFileName(self, "Open project", "", PROJECT_FILTER)
            if not path:
                return False
        try:
            self.session.load(path)
        except (ValidationError, OSError) as error:
            self._report_error("Cannot open project", error)
            return False
        self._saved_journal_length = len(self.session.commands.journal.lines)
        self._refresh_all()
        self.scene.fit()
        log.info("opened %s", path)
        return True

    def save_project(self) -> bool:
        if self.session.project_path is None:
            return self.save_project_as()
        return self._save_to(self.session.project_path)

    def save_project_as(self, path: str | None = None) -> bool:
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(self, "Save project", "", PROJECT_FILTER)
            if not path:
                return False
            if not path.endswith(".ogp"):
                path += ".ogp"
        return self._save_to(Path(path))

    def _save_to(self, path: Path) -> bool:
        try:
            self.session.save(path)
        except OSError as error:
            self._report_error("Cannot save project", error)
            return False
        self._saved_journal_length = len(self.session.commands.journal.lines)
        self._update_title()
        log.info("saved %s", path)
        return True

    def import_airfoil(self, path: str | None = None, mode: str | None = None) -> list[str]:
        if path is None:
            path, _filter = QFileDialog.getOpenFileName(
                self, "Import airfoil", "", "Airfoil coordinates (*.dat *.txt);;All files (*)"
            )
            if not path:
                return []
        if mode is None:
            mode, ok = QInputDialog.getItem(
                self, "Import airfoil", "Create as:", ["spline", "split", "polyline"], 0, False
            )
            if not ok:
                return []
        try:
            names = self.session.import_airfoil_dat(path, mode=mode)
        except ValidationError as error:
            self._report_error("Cannot import airfoil", error)
            return []
        log.info("imported %s as %s", Path(path).name, ", ".join(names))
        self.scene.fit()
        return names

    def import_csv(self, path: str | None = None, mode: str | None = None) -> list[str]:
        if path is None:
            path, _filter = QFileDialog.getOpenFileName(
                self, "Import points", "", "Point lists (*.csv *.txt *.dat);;All files (*)"
            )
            if not path:
                return []
        if mode is None:
            mode, ok = QInputDialog.getItem(
                self, "Import points", "Create as:", ["spline", "polyline"], 0, False
            )
            if not ok:
                return []
        try:
            names = self.session.import_csv(path, mode=mode)
        except ValidationError as error:
            self._report_error("Cannot import points", error)
            return []
        log.info("imported %s as %s", Path(path).name, ", ".join(names))
        self.scene.fit()
        return names

    def save_journal(self, path: str | None = None) -> None:
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(
                self, "Save journal", "journal.py", "Python script (*.py)"
            )
            if not path:
                return
        self.session.save_journal(path)
        log.info("journal saved to %s", path)

    def _report_error(self, title: str, error: Exception) -> None:
        log.error("%s: %s", title, error)
        if self.isVisible() and not self._off_screen:
            QMessageBox.critical(self, title, str(error))

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt override)
        logging.getLogger("opengrid").removeHandler(self.log_panel.handler)
        if self._on_model_changed in self.session.commands._listeners:
            self.session.commands.remove_listener(self._on_model_changed)
        self.viewport.close()
        super().closeEvent(event)
