"""Main application window: viewport, model tree, property panel, log and Python console."""

import logging
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QSpinBox,
    QToolBar,
    QWidget,
)

import opengrid
from opengrid.api import Session
from opengrid.app.boundary_dock import (
    INTERFACE_COLOR,
    PATCH_PALETTE,
    UNASSIGNED_COLOR,
    BoundaryDock,
)
from opengrid.app.console import PythonConsole
from opengrid.app.log_panel import LogPanel
from opengrid.app.model_tree import ModelTree
from opengrid.app.property_panel import PropertyPanel
from opengrid.app.quality_dock import METRICS, QualityDock, worst_cell
from opengrid.app.scene import Scene, pick_entity
from opengrid.app.smoother_dock import SmootherDock
from opengrid.app.tools import TOOLS, Tool, snap, snap_targets
from opengrid.validation import ValidationError

log = logging.getLogger("opengrid")

SNAP_PIXELS = 10.0
PICK_PIXELS = 8.0
CLICK_MOVE_PIXELS = 4.0
PROJECT_FILTER = "OpenGrid project (*.ogp)"
DEFAULT_DIMENSION = 21
CURVE_KINDS = ("line", "arc", "polyline", "spline", "nurbs")


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

        model_dock = self._add_dock("Model", self.model_tree, Qt.DockWidgetArea.LeftDockWidgetArea)
        properties_dock = self._add_dock(
            "Properties", self.property_panel, Qt.DockWidgetArea.LeftDockWidgetArea
        )
        console_dock = self._add_dock(
            "Console", self.console, Qt.DockWidgetArea.BottomDockWidgetArea
        )
        log_dock = self._add_dock("Log", self.log_panel, Qt.DockWidgetArea.BottomDockWidgetArea)
        self.quality_dock = QualityDock()
        self.smoother_dock = SmootherDock(session)
        quality_dock = self._add_dock(
            "Quality", self.quality_dock, Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.smoother_dock_widget = self._add_dock(
            "Smoother", self.smoother_dock, Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.boundary_dock = BoundaryDock()
        self.boundary_dock_widget = self._add_dock(
            "Boundary", self.boundary_dock, Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.tabifyDockWidget(quality_dock, self.smoother_dock_widget)
        self.tabifyDockWidget(quality_dock, self.boundary_dock_widget)
        quality_dock.raise_()
        self.quality_dock_widget = quality_dock
        self.tabifyDockWidget(log_dock, console_dock)
        console_dock.raise_()
        self.resizeDocks([model_dock, properties_dock], [320, 480], Qt.Orientation.Vertical)

        self._create_actions()
        self._create_menus()
        self._create_toolbars()
        self._connect_viewport_events()

        self.model_tree.selection_changed.connect(self.set_selection)
        self.property_panel.message.connect(self._on_panel_message)
        self.property_panel.smooth_requested.connect(self.run_smoother)
        self.smoother_dock.finished.connect(self._on_smoothing_finished)
        self.quality_dock.metric_changed.connect(self._on_metric_changed)
        self.quality_dock.locate_requested.connect(self.show_worst_cell)
        self.boundary_dock.patch_requested.connect(self.create_patch_from_selection)
        self.boundary_dock.add_requested.connect(self.add_selection_to_patch)
        self.boundary_dock.remove_requested.connect(self.remove_selection_from_patch)
        self.boundary_dock.zone_requested.connect(self.create_zone_from_selection)
        self.boundary_dock.validate_requested.connect(self.validate_boundary)
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

        self.connectors_action = self._action(
            "Connectors on curves",
            self.create_connectors_on_selection,
            tip="One connector on each selected curve",
        )
        self.chain_action = self._action(
            "Chain connector",
            self.create_chain_connector_on_selection,
            tip="One connector along the selected curves (ordered end to end)",
        )
        self.color_action = QAction("Colour connectors by stretching ratio", self)
        self.color_action.setCheckable(True)
        self.color_action.setChecked(True)
        self.color_action.toggled.connect(self.set_color_by_stretching)
        self.nodes_action = QAction("Show connector nodes", self)
        self.nodes_action.setCheckable(True)
        self.nodes_action.setChecked(True)
        self.nodes_action.toggled.connect(self.set_show_nodes)
        self.threshold_action = self._action("Stretching ratio threshold…", self.ask_threshold)
        self.unstructured_action = self._action(
            "Unstructured domain",
            self.create_unstructured_on_selection,
            tip="Triangulated domain bounded by the selected connectors (closed loops)",
        )
        self.domain_action = self._action(
            "Structured domain",
            self.create_domain_on_selection,
            tip="Structured domain from the selected connectors (edges assembled automatically)",
        )
        self.quality_action = QAction("Colour domains by quality", self)
        self.quality_action.setCheckable(True)
        self.quality_action.toggled.connect(self.set_quality_colouring)
        self.patch_color_action = QAction("Colour connectors by patch", self)
        self.patch_color_action.setCheckable(True)
        self.patch_color_action.toggled.connect(self.set_patch_colouring)
        self.export_cgns_action = self._action("Export CGNS…", self.export_cgns)
        self.export_gmsh_action = self._action("Export Gmsh…", self.export_gmsh)
        self.export_vtkhdf_action = self._action("Export VTKHDF…", self.export_vtkhdf)
        self.export_plot3d_action = self._action("Export Plot3D…", self.export_plot3d)

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
        file_menu.addAction(self.export_vtkhdf_action)
        file_menu.addAction(self.export_cgns_action)
        file_menu.addAction(self.export_gmsh_action)
        file_menu.addAction(self.export_plot3d_action)
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
        create_menu.addAction(self.connectors_action)
        create_menu.addAction(self.chain_action)
        create_menu.addAction(self.domain_action)
        create_menu.addAction(self.unstructured_action)
        create_menu.addSeparator()
        create_menu.addAction(self.finish_action)
        create_menu.addAction(self.cancel_action)

        view_menu = self.menuBar().addMenu("&View")
        view_menu.addAction(self.fit_action)
        view_menu.addSeparator()
        view_menu.addAction(self.color_action)
        view_menu.addAction(self.nodes_action)
        view_menu.addAction(self.threshold_action)
        view_menu.addAction(self.quality_action)
        view_menu.addAction(self.patch_color_action)

    def _create_toolbars(self) -> None:
        toolbar = QToolBar("Tools", self)
        toolbar.setObjectName("toolbar_tools")
        for action in self.tool_actions.values():
            toolbar.addAction(action)
        toolbar.addSeparator()
        toolbar.addAction(self.finish_action)
        toolbar.addAction(self.cancel_action)
        toolbar.addSeparator()
        toolbar.addAction(self.connectors_action)
        toolbar.addAction(self.chain_action)
        toolbar.addAction(self.domain_action)
        toolbar.addAction(self.unstructured_action)
        toolbar.addWidget(QLabel(" dim "))
        self.dimension_box = QSpinBox()
        self.dimension_box.setRange(2, 100000)
        self.dimension_box.setValue(DEFAULT_DIMENSION)
        self.dimension_box.setToolTip("Number of points for new connectors")
        toolbar.addWidget(self.dimension_box)
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
        domains = self._selected_domains()
        if len(domains) == 1:
            self.smoother_dock.set_domain(domains[0])
        else:
            self.smoother_dock.set_domain(None)
        self.refresh_quality()

    def selected_names(self) -> list[str]:
        return [self.session.store.get(i).name for i in self.selected_ids]

    def _on_model_changed(self, change) -> None:
        had_geometry = len(self.scene._actors) > 0
        if change is None:
            self.scene.rebuild()
        else:
            self.scene.update_entities(change.ids())
            self._report_connector_warnings(change.ids())
        self.model_tree.populate(self.session)
        if self.scene.patch_colors is not None:
            self.scene.patch_colors = self._patch_colors()
            self.scene.update_entities([])
        self.set_selection(self.selected_ids)
        self._update_edit_actions()
        if not had_geometry and len(self.session.store) > 0:
            self.scene.fit()

    def _report_connector_warnings(self, ids) -> None:
        for entity_id in ids:
            if entity_id not in self.session.store:
                continue
            entity = self.session.store.get(entity_id)
            if entity.kind != "connector":
                continue
            nodes = self.session.connector_nodes(entity.name)
            for warning in nodes.warnings:
                log.warning("%s: %s", entity.name, warning)
            if nodes.max_stretching_ratio > self.session.stretching_threshold:
                log.warning(
                    "%s: max stretching ratio %.4g exceeds %g",
                    entity.name,
                    nodes.max_stretching_ratio,
                    self.session.stretching_threshold,
                )

    def _refresh_all(self) -> None:
        self.scene.rebuild()
        self.model_tree.populate(self.session)
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

    # ---------------------------------------------------------------- domains
    def _selected_domains(self) -> list[str]:
        names = []
        for entity_id in self.selected_ids:
            if entity_id not in self.session.store:
                continue
            entity = self.session.store.get(entity_id)
            if entity.kind in ("domain", "unstructured_domain"):
                names.append(entity.name)
        return names

    def _quality_scope(self) -> tuple[list[str], str]:
        names = self._selected_domains()
        if names:
            return names, ", ".join(names)
        return self.session.names("domain") + self.session.names(
            "unstructured_domain"
        ), "all domains"

    def create_domain_on_selection(self) -> str | None:
        connectors = []
        for name in self.selected_names():
            if self.session.store.get(name).kind == "connector":
                connectors.append(name)
        if not connectors:
            self.statusBar().showMessage("Select the connectors of the domain first", 5000)
            return None
        name = self._run(lambda: self.session.create_structured_domain(connectors=connectors))
        if name:
            ni, nj = self.session.domain_dimensions(name)
            log.info("created %s (%d x %d)", name, ni, nj)
            self.set_selection([self.session.store.get(name).id])
        return name

    def create_unstructured_on_selection(self) -> str | None:
        connectors = self._selected_of_kind("connector")
        if not connectors:
            self.statusBar().showMessage("Select the connectors of the boundary loops first", 5000)
            return None
        name = self._run(lambda: self.session.create_unstructured_domain(connectors=connectors))
        if name:
            log.info("created %s: %s", name, self.session.domain_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def run_smoother(self, name: str) -> bool:
        self.smoother_dock_widget.raise_()
        started = self.smoother_dock.start(name)
        if started:
            log.info("smoothing %s…", name)
        return started

    def _on_smoothing_finished(self, name: str, result) -> None:
        if isinstance(result, Exception):
            log.error("smoothing %s failed: %s", name, result)
            return
        for warning in result.warnings:
            log.warning("%s: %s", name, warning)
        if result.cancelled:
            log.info("%s: smoothing cancelled after %d iterations", name, result.iterations)
        else:
            log.info("%s: %s", name, self.session.domain_status(name))
        if name in self.session.store:
            self.scene.update_entities([self.session.store.get(name).id])
        self.model_tree.populate(self.session)
        self.property_panel.refresh()
        self.refresh_quality()

    def set_quality_colouring(self, enabled: bool) -> None:
        if enabled:
            self.scene.quality_metric = self.quality_dock.metric_key()
        else:
            self.scene.quality_metric = None
        self.scene.rebuild()

    def _on_metric_changed(self, key: str) -> None:
        if self.scene.quality_metric is not None:
            self.scene.quality_metric = key
            self.scene.rebuild()
        self.refresh_quality()

    def refresh_quality(self) -> None:
        names, scope = self._quality_scope()
        if not names:
            self.quality_dock.show_values(None, scope)
            return
        key = self.quality_dock.metric_key()
        values = []
        for name in names:
            values.append(
                self.scene._display_quality(self.session.store.get(name))[key].reshape(-1)
            )
        self.quality_dock.show_values(np.concatenate(values), scope)

    def _cell_corners(self, name: str, j: int, i: int) -> np.ndarray:
        if self.session.store.get(name).kind == "unstructured_domain":
            mesh = self.session.unstructured_mesh(name)
            n = len(mesh["triangles"])
            if i < n:
                return mesh["points"][mesh["triangles"][i]]
            return mesh["points"][mesh["quads"][i - n]]
        grid = self.session.domain_grid(name, compute=False)
        return np.array([grid[j, i], grid[j, i + 1], grid[j + 1, i + 1], grid[j + 1, i]])

    def show_worst_cell(self) -> tuple[str, int, int] | None:
        names, _scope = self._quality_scope()
        if not names:
            return None
        key = self.quality_dock.metric_key()
        larger_is_worse = METRICS[key][1]
        best = None
        for name in names:
            values = self.scene._display_quality(self.session.store.get(name))[key]
            j, i = worst_cell(values, larger_is_worse)
            value = float(values[j, i])
            better = best is None
            if not better:
                if larger_is_worse:
                    better = value > best[0]
                else:
                    better = value < best[0]
            if better:
                best = (value, name, j, i)
        value, name, j, i = best
        self.scene.show_worst_cell(self._cell_corners(name, j, i))
        log.info("worst %s: %.4g in %s at cell (i=%d, j=%d)", key, value, name, i, j)
        return name, i, j

    def export_vtkhdf(self, path: str | None = None) -> dict | None:
        if not self.session.names("domain"):
            self.statusBar().showMessage("No domains to export", 5000)
            return None
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(
                self, "Export VTKHDF", "", "VTKHDF (*.vtkhdf *.hdf)"
            )
            if not path:
                return None
        summary = self._run(lambda: self.session.export_vtkhdf(path))
        if summary:
            log.info("exported %s: %s", path, summary)
        return summary

    def export_plot3d(self, path: str | None = None, binary: bool = True) -> bool:
        if not self.session.names("domain"):
            self.statusBar().showMessage("No domains to export", 5000)
            return False
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(
                self, "Export Plot3D", "", "Plot3D (*.xyz *.x *.p3d)"
            )
            if not path:
                return False
        self._run(lambda: self.session.export_plot3d(path, binary=binary))
        log.info("exported %s", path)
        return True

    # ------------------------------------------------------- patches and zones
    def _selected_of_kind(self, kind: str) -> list[str]:
        return [n for n in self.selected_names() if self.session.store.get(n).kind == kind]

    def create_patch_from_selection(self, name: str, patch_type: str = "wall") -> str | None:
        connectors = self._selected_of_kind("connector")
        if not name:
            self.statusBar().showMessage("Enter a patch name", 5000)
            return None
        created = self._run(lambda: self.session.create_patch(name, patch_type, connectors))
        if created:
            log.info(
                "created patch %s (%s) with %d connector(s)", created, patch_type, len(connectors)
            )
            self.validate_boundary()
        return created

    def add_selection_to_patch(self, name: str) -> bool:
        connectors = self._selected_of_kind("connector")
        if not connectors:
            self.statusBar().showMessage("Select connectors first", 5000)
            return False
        done = self._run(lambda: self.session.add_to_patch(name, connectors))
        self.validate_boundary()
        return bool(done)

    def remove_selection_from_patch(self, name: str) -> bool:
        connectors = self._selected_of_kind("connector")
        done = self._run(lambda: self.session.remove_from_patch(name, connectors))
        self.validate_boundary()
        return bool(done)

    def create_zone_from_selection(self, name: str) -> str | None:
        domains = self._selected_of_kind("domain")
        if not name:
            self.statusBar().showMessage("Enter a zone name", 5000)
            return None
        created = self._run(lambda: self.session.create_zone(name, domains))
        if created:
            log.info("created zone %s with %d domain(s)", created, len(domains))
            self.validate_boundary()
        return created

    def validate_boundary(self):
        report = self.session.boundary_report()
        self.boundary_dock.show_report(report)
        if self.scene.patch_colors is not None:
            self.scene.patch_colors = self._patch_colors(report)
            self.scene.rebuild()
        return report

    def _patch_colors(self, report=None) -> dict[int, str]:
        topology = self.session.topology()
        colors = {}
        for connector_id, uses in topology.uses.items():
            if len(uses) == 2:
                colors[connector_id] = INTERFACE_COLOR
            elif len(uses) == 1:
                colors[connector_id] = UNASSIGNED_COLOR
        for index, patch in enumerate(self.session.store.entities("patch")):
            for item in patch.params["members"]:
                colors[item["connector"]] = PATCH_PALETTE[index % len(PATCH_PALETTE)]
        return colors

    def set_patch_colouring(self, enabled: bool) -> None:
        if enabled:
            self.scene.patch_colors = self._patch_colors()
        else:
            self.scene.patch_colors = None
        self.scene.rebuild()

    def export_cgns(self, path: str | None = None) -> dict | None:
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(self, "Export CGNS", "", "CGNS (*.cgns)")
            if not path:
                return None
        summary = self._run(lambda: self.session.export_cgns(path))
        if summary:
            log.info("exported %s: %s", path, summary)
        self.validate_boundary()
        return summary

    def export_gmsh(self, path: str | None = None) -> dict | None:
        if path is None:
            path, _filter = QFileDialog.getSaveFileName(self, "Export Gmsh", "", "Gmsh (*.msh)")
            if not path:
                return None
        summary = self._run(lambda: self.session.export_gmsh(path))
        if summary:
            log.info("exported %s: %s", path, summary)
        self.validate_boundary()
        return summary

    # -------------------------------------------------------------- connectors
    def _selected_curves(self) -> list[str]:
        names = []
        for name in self.selected_names():
            if self.session.store.get(name).kind in CURVE_KINDS:
                names.append(name)
        return names

    def create_connectors_on_selection(self) -> list[str]:
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select one or more curves first", 5000)
            return []
        names = self._run(
            lambda: self.session.create_connectors_on_curves(curves, self.dimension_box.value())
        )
        if names:
            log.info("created %s", ", ".join(names))
            self.set_selection([self.session.store.get(n).id for n in names])
        return names or []

    def create_chain_connector_on_selection(self) -> str | None:
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select the curves of the chain first", 5000)
            return None
        name = self._run(
            lambda: self.session.create_chain_connector(curves, self.dimension_box.value())
        )
        if name:
            log.info("created %s", name)
            self.set_selection([self.session.store.get(name).id])
        return name

    def set_color_by_stretching(self, enabled: bool) -> None:
        self.scene.color_by_stretching = enabled
        self.scene.rebuild()

    def set_show_nodes(self, enabled: bool) -> None:
        self.scene.show_nodes = enabled
        self.scene.rebuild()

    def set_threshold(self, value: float) -> None:
        self.session.stretching_threshold = value
        self.scene.rebuild()
        self.model_tree.populate(self.session)
        self.property_panel.refresh()
        log.info("stretching ratio threshold set to %g", value)

    def ask_threshold(self) -> None:
        value, ok = QInputDialog.getDouble(
            self,
            "Stretching ratio threshold",
            "Warn when consecutive spacings differ by more than:",
            self.session.stretching_threshold,
            1.0,
            100.0,
            3,
        )
        if ok:
            self.set_threshold(value)

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
