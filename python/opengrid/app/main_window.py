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


def compute_summary() -> str:
    """One line on the compute devices: CPU threads and the GPU (OpenCL) if one passed the
    self-test. Rendering is VTK's OpenGL on the system's GPU driver either way."""
    text = f"Compute: {opengrid.threads()} CPU threads"
    gpu = opengrid.gpu_info()
    if gpu["available"]:
        text += (
            f"; GPU {gpu['device']} ({gpu['version'].strip()}, {gpu['compute_units']} compute "
            f"units), experimental GPU smoothing {opengrid.gpu_mode()} (opengrid.set_gpu)"
        )
    else:
        text += f"; no GPU compute ({gpu['reason']})"
    return text


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
        self._drag = None  # handle drag in progress (app.handles.Drag)
        self._stroke: list | None = None  # freehand stroke in progress
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
        log.info("%s", compute_summary())

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
        self.import_cad_action = self._action("Import &CAD (STEP / IGES)…", self.import_cad)
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
        self.block_action = self._action(
            "Block from 6 faces",
            self.create_block_on_selection,
            tip="Structured block (3D TFI) from six selected structured domains",
        )
        self.extrude_action = self._action(
            "Extrude domain",
            self.extrude_selection,
            tip="Block by extruding the selected structured domain (edit the extrusion afterwards)",
        )
        self.network_action = self._action(
            "Domains from curves",
            self.domains_from_selected_curves,
            tip="Structured domains from the selected curves drawn across each other (regions "
            "with four corners); point counts from the 'dim' box",
        )
        self.group_action = self._action(
            "Smoothing group",
            self.create_group_on_selection,
            tip="Smooth the selected structured domains together as one grid (interfaces move); "
            "selected connectors become walls",
        )
        self.volume_action = self._action(
            "Volume from faces",
            self.create_volume_on_selection,
            tip="Unstructured volume (tetrahedra) enclosed by the selected domains",
        )
        self.cut_action = QAction("Cut blocks (plane x)", self)
        self.cut_action.setCheckable(True)
        self.cut_action.toggled.connect(self.set_cut_plane)
        self.shrink_action = QAction("Shrink block cells", self)
        self.shrink_action.setCheckable(True)
        self.shrink_action.toggled.connect(self.set_shrink)
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
        self.export_curved_high_order_action = self._action(
            "Export curved high-order elements (structured)…", self.export_curved_high_order
        )
        self.export_high_order_action = self._action(
            "Export high-order elements…", self.export_high_order
        )

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
            "nurbs": "NURBS",
            "freehand": "Freehand",
            "sketch": "Sketch block",
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
        file_menu.addAction(self.import_cad_action)
        file_menu.addSeparator()
        file_menu.addAction(self.export_vtkhdf_action)
        file_menu.addAction(self.export_cgns_action)
        file_menu.addAction(self.export_gmsh_action)
        file_menu.addAction(self.export_plot3d_action)
        file_menu.addAction(self.export_high_order_action)
        file_menu.addAction(self.export_curved_high_order_action)
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
        for name in ("point", "line", "arc", "polyline", "spline", "nurbs", "freehand", "sketch"):
            create_menu.addAction(self.tool_actions[name])
        create_menu.addSeparator()
        create_menu.addAction(self.connectors_action)
        create_menu.addAction(self.chain_action)
        create_menu.addAction(self.domain_action)
        create_menu.addAction(self.unstructured_action)
        create_menu.addAction(self.network_action)
        create_menu.addAction(self.group_action)
        create_menu.addAction(self.block_action)
        create_menu.addAction(self.extrude_action)
        create_menu.addAction(self.volume_action)
        create_menu.addSeparator()
        create_menu.addAction(self.finish_action)
        create_menu.addAction(self.cancel_action)

        topology_menu = self.menuBar().addMenu("&Topology")
        self.topology_actions = {
            "split_connector": self._action("Split connector at node…", self.split_connector),
            "join_connectors": self._action("Join two connectors", self.join_connectors),
            "split_domain": self._action("Split structured domain…", self.split_domain),
            "join_domains": self._action("Join two structured domains", self.join_domains),
            "split_block": self._action("Split block…", self.split_block),
            "join_blocks": self._action("Join two blocks", self.join_blocks),
            "ogrid": self._action("O-grid (butterfly) domain or block…", self.ogrid),
            "sketch_assign": self._action(
                "Assign selected geometry to sketch",
                self.assign_selected_geometry,
                tip="Assign the selected curves (surfaces) to the sketch edges (faces) near them",
            ),
            "sketch_regenerate": self._action("Regenerate sketch", self.regenerate_sketch),
            "sketch_report": self._action("Sketch report", self.sketch_report),
            "template_ogrid": self._action(
                "Template: O-grid around the selected body",
                lambda: self.template_from_selection("ogrid"),
            ),
            "template_cgrid": self._action(
                "Template: C-grid around the selected body (sharp trailing edge)",
                lambda: self.template_from_selection("cgrid"),
            ),
            "template_butterfly": self._action(
                "Template: butterfly inside the selected closed curve",
                lambda: self.template_from_selection("butterfly"),
            ),
            "auto_topology": self._action(
                "Automatic topology (2D) from the selected curves",
                self.auto_topology_from_selection,
                tip="Cross-field block decomposition of the region the curves bound",
            ),
        }
        for key, action in self.topology_actions.items():
            if key in ("split_domain", "split_block", "ogrid", "sketch_assign", "template_ogrid"):
                topology_menu.addSeparator()
            topology_menu.addAction(action)

        view_menu = self.menuBar().addMenu("&View")
        view_menu.addAction(self.fit_action)
        view_menu.addSeparator()
        view_menu.addAction(self.color_action)
        view_menu.addAction(self.nodes_action)
        view_menu.addAction(self.threshold_action)
        view_menu.addAction(self.quality_action)
        view_menu.addAction(self.patch_color_action)
        view_menu.addSeparator()
        view_menu.addAction(self.cut_action)
        view_menu.addAction(self.shrink_action)

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
        toolbar.addAction(self.block_action)
        toolbar.addAction(self.extrude_action)
        toolbar.addAction(self.volume_action)
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
        # Drags of handles and freehand strokes take the left button before the camera does:
        # high-priority observers that stop the event when they use it.
        interactor = iren.interactor
        self._grab_tags = {}
        for event, handler in (
            ("LeftButtonPressEvent", self._grab_press),
            ("MouseMoveEvent", self._grab_move),
            ("LeftButtonReleaseEvent", self._grab_release),
        ):
            self._grab_tags[event] = interactor.AddObserver(event, handler, 10.0)

    def _abort(self, interactor, event: str) -> None:
        interactor.GetCommand(self._grab_tags[event]).SetAbortFlag(1)

    def _grab_press(self, interactor, event) -> None:
        x, y = interactor.GetEventPosition()
        if self.handle_press(self.scene.display_to_world(x, y)):
            self._abort(interactor, event)

    def _grab_move(self, interactor, event) -> None:
        if self._drag is None and self._stroke is None:
            return
        x, y = interactor.GetEventPosition()
        self.handle_drag(self.scene.display_to_world(x, y))
        self._abort(interactor, event)

    def _grab_release(self, interactor, event) -> None:
        if self._drag is None and self._stroke is None:
            return
        x, y = interactor.GetEventPosition()
        self.handle_release(self.scene.display_to_world(x, y))
        self._press_position = None
        self._abort(interactor, event)

    # ------------------------------------------------- drags and freehand strokes
    def visible_handles(self):
        from opengrid.app.handles import handles_for

        if self.tool.name != "select":
            return []
        return handles_for(self.session, self.selected_names())

    def _show_handles(self) -> None:
        self.scene.show_handles([h.xyz for h in self.visible_handles()])

    def handle_press(self, xyz) -> bool:
        """Start a handle drag (Select tool, near a handle of a selected entity) or a freehand
        stroke (Freehand tool). Returns True when the press is used."""
        from opengrid.app.handles import Drag, pick_handle

        if getattr(self.tool, "stroke", False):
            self._stroke = [tuple(xyz)]
            return True
        if self.tool.name != "select":
            return False
        tolerance = PICK_PIXELS * self.scene.world_per_pixel()
        handle = pick_handle(self.visible_handles(), xyz, tolerance)
        if handle is None:
            return False
        self._drag = Drag(self.session, handle)
        self.statusBar().showMessage(f"Dragging {handle.entity}.{handle.field}")
        return True

    def handle_drag(self, xyz) -> None:
        if self._stroke is not None:
            self._stroke.append(tuple(xyz))
            self.scene.show_preview(self._stroke)
            return
        if self._drag is None:
            return
        self._drag.preview(xyz, self.scene.preview_entities)

    def handle_release(self, xyz) -> list[str]:
        """Finish a stroke (returns the curves) or a drag (returns [entity])."""
        if self._stroke is not None:
            stroke = [*self._stroke, tuple(xyz)]
            self._stroke = None
            self.scene.show_preview([])
            created = self._run(lambda: self.tool.create_stroke(stroke))
            if created:
                log.info("freehand: %s", ", ".join(created))
            return created or []
        if self._drag is None:
            return []
        drag = self._drag
        self._drag = None
        name = self._run(lambda: drag.finish(xyz))
        self.scene.end_preview()
        self.scene.update_entities(drag.affected)
        self._show_handles()
        if not name:
            return []
        self._smooth_affected(drag.affected)
        return [name]

    def _smooth_affected(self, ids) -> None:
        """After a drag: rerun the smoothing of an affected group or domain (the preview
        showed TFI)."""
        store = self.session.store
        for entity_id in ids:
            entity = store.get(entity_id)
            if entity.kind == "smoothing_group":
                self.run_smoother(entity.name)
                return
        for entity_id in ids:
            entity = store.get(entity_id)
            if entity.kind == "domain" and entity.params["generator"]["elliptic"] is not None:
                self.run_smoother(entity.name)
                return

    def domains_from_selected_curves(self) -> list[str]:
        from opengrid.model.connectors import CURVE_KINDS

        curves = [n for n in self.selected_names() if self.session.store.get(n).kind in CURVE_KINDS]
        if not curves:
            self.statusBar().showMessage("Select the curves of the network first", 5000)
            return []
        total = sum(self.session.length(c) for c in curves)
        spacing = total / len(curves) / max(self.dimension_box.value() - 1, 1)
        created = self._run(lambda: self.session.domains_from_curves(curves, spacing=spacing))
        if created:
            log.info("created %s", ", ".join(created))
            for warning in self.session.last_warnings:
                log.warning(warning)
            self.set_selection([self.session.store.get(n).id for n in created])
        return created or []

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
            if self.session.store.get(created).kind == "sketch":
                self._after_sketch(created)
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
        if name == "sketch":
            sketches = self._selected_of_kind("sketch")
            if sketches:
                self.tool.sketch = sketches[0]
            self.tool.points_per_edge = self.dimension_box.value()
        self._stroke = None
        self._drag = None
        self._show_handles()
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
        self._show_handles()
        self.model_tree.select_ids(ids)
        if len(ids) == 1:
            self.property_panel.show_entity(ids[0])
        else:
            self.property_panel.show_entity(None)
        self.delete_action.setEnabled(bool(ids))
        domains = self._selected_domains()
        groups = self._selected_of_kind("smoothing_group")
        blocks = self._selected_of_kind("block")
        if len(groups) == 1 and not domains:
            self.smoother_dock.set_domain(groups[0])
        elif len(blocks) == 1 and not domains:
            self.smoother_dock.set_domain(blocks[0])
        elif len(domains) == 1:
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
            if entity.kind in ("domain", "unstructured_domain", "block", "volume"):
                names.append(entity.name)
        return names

    def create_block_on_selection(self) -> str | None:
        faces = self._selected_of_kind("domain")
        if len(faces) != 6:
            self.statusBar().showMessage("Select the six structured domains of the block", 5000)
            return None
        name = self._run(lambda: self.session.create_block(faces=faces))
        if name:
            log.info("created %s: %s", name, self.session.block_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def create_volume_on_selection(self) -> str | None:
        faces = self._selected_of_kind("domain") + self._selected_of_kind("unstructured_domain")
        if len(faces) < 2:
            self.statusBar().showMessage("Select the domains enclosing the volume", 5000)
            return None
        name = self._run(lambda: self.session.create_volume(faces))
        if name:
            log.info("created %s: %s", name, self.session.volume_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def extrude_selection(self) -> str | None:
        domains = self._selected_of_kind("domain")
        if len(domains) != 1:
            self.statusBar().showMessage("Select one structured domain to extrude", 5000)
            return None
        grid = self.session.domain_grid(domains[0])
        size = float(np.linalg.norm(np.ptp(grid.reshape(-1, 3), axis=0)))
        normal = np.cross(grid[0, -1] - grid[0, 0], grid[-1, 0] - grid[0, 0])
        if np.linalg.norm(normal) == 0:
            normal = np.array([0.0, 0.0, 1.0])
        vector = tuple(float(v) for v in 0.25 * size * normal / np.linalg.norm(normal))
        name = self._run(
            lambda: self.session.create_extruded_block(
                domains[0], "translate", vector=vector, dimension=self.dimension_box.value()
            )
        )
        if name:
            log.info("created %s: %s", name, self.session.block_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def set_cut_plane(self, enabled: bool) -> None:
        """Cut the blocks with the plane x = centre of the model (normal +x)."""
        self.scene.cut_plane = None
        if enabled and self.session.is_3d():
            parts = [self.session.block_grid(n).reshape(-1, 3) for n in self.session.names("block")]
            parts += [self.session.volume_mesh(n)["points"] for n in self.session.names("volume")]
            points = np.concatenate(parts)
            centre = tuple(float(v) for v in 0.5 * (points.min(axis=0) + points.max(axis=0)))
            self.scene.cut_plane = (centre, (1.0, 0.0, 0.0))
        self.scene.rebuild()

    def set_shrink(self, enabled: bool) -> None:
        self.scene.shrink = enabled
        self.scene.rebuild()

    def _quality_scope(self) -> tuple[list[str], str]:
        names = self._selected_domains()
        if names:
            return names, ", ".join(names)
        if self.session.is_3d():
            return self.session.names("block") + self.session.names("volume"), "all blocks"
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

    def create_group_on_selection(self) -> str | None:
        blocks = self._selected_of_kind("block")
        domains = self._selected_of_kind("domain")
        if blocks:
            # 3D: selected domains (block faces) become walls.
            name = self._run(
                lambda: self.session.create_smoothing_group(blocks=blocks, walls=domains)
            )
            if name:
                log.info("created %s from %d blocks", name, len(blocks))
                self.set_selection([self.session.store.get(name).id])
                self.run_smoother(name)
            return name
        if not domains:
            self.statusBar().showMessage(
                "Select the structured domains (2D) or blocks (3D) of the group", 5000
            )
            return None
        walls = self._selected_of_kind("connector")
        name = self._run(lambda: self.session.create_smoothing_group(domains, walls=walls))
        if name:
            log.info("created %s from %d domains", name, len(domains))
            self.set_selection([self.session.store.get(name).id])
            self.run_smoother(name)
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
        elif self.session.store.get(name).kind == "smoothing_group":
            log.info("%s: %s", name, self.session.group_status(name))
        else:
            log.info("%s: %s", name, self.session.domain_status(name))
        if name in self.session.store:
            entity = self.session.store.get(name)
            ids = [entity.id]
            if entity.kind == "smoothing_group":
                ids = list(entity.params["domains"]) + list(entity.params["blocks"])
                cached = self.session.groups.result(entity.id)
                if cached is not None:
                    ids += sorted(cached.roles)  # moved interfaces, sliding connectors / faces
                    for d in getattr(cached, "faces", {}):
                        for edge in self.session.store.get(d).params["edges"]:
                            ids += [item["connector"] for item in edge]
                ids = sorted(set(ids))
            self.scene.update_entities(ids)
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
        if self.session.store.get(name).kind == "volume":
            from opengrid.model.volumes import CELL_FACES, volume_cells

            mesh = self.session.volume_mesh(name)
            for key, array, _vtk in volume_cells(mesh):
                if i < len(array):
                    # A polyline over the cell's faces (outline of the worst cell).
                    order = [v for face in CELL_FACES[key] for v in (*face, face[0])]
                    return mesh["points"][array[i]][order]
                i -= len(array)
        if self.session.store.get(name).kind == "block":
            grid = self.session.block_grid(name)
            nk, nj, ni, _ = grid.shape
            k, jj, ii = np.unravel_index(i, (nk - 1, nj - 1, ni - 1))
            c = grid[k : k + 2, jj : jj + 2, ii : ii + 2].reshape(-1, 3)
            # Bottom face loop then top face loop (outlined as one polyline).
            return c[[0, 1, 3, 2, 0, 4, 5, 7, 6, 4]]
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
            if values.ndim == 3:
                # Blocks: (k, j, i); the cell is located by its flat index.
                if larger_is_worse:
                    flat = int(np.argmax(values))
                else:
                    flat = int(np.argmin(values))
                j, i = 0, flat
                value = float(values.reshape(-1)[flat])
            else:
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

    def export_high_order(
        self, path: str | None = None, order: int | None = None, nodes: str = "equispaced"
    ) -> dict | None:
        """Straight-sided order-p elements of every mesh cell: Gmsh, CGNS, VTKHDF Lagrange
        cells, or .h5 with the chosen node set."""
        if order is None:
            from PySide6.QtWidgets import QInputDialog

            order, ok = QInputDialog.getInt(
                self, "High-order elements", "Order (CGNS: up to 4):", 2, 1, 6
            )
            if not ok:
                return None
        if path is None:
            path, chosen = QFileDialog.getSaveFileName(
                self,
                "Export high-order elements",
                "",
                "Gmsh (*.msh);;CGNS (*.cgns);;VTKHDF Lagrange cells (*.vtkhdf);;"
                "HDF5, equispaced nodes (*.h5);;HDF5, Gauss-Lobatto nodes (*.h5);;"
                "HDF5, Gauss nodes (*.h5)",
            )
            if not path:
                return None
            if "Lobatto" in chosen:
                nodes = "gll"
            elif "Gauss nodes" in chosen:
                nodes = "gauss"
        summary = self._run(lambda: self.session.export_high_order(path, order=order, nodes=nodes))
        if summary:
            log.info("exported %s: %s", path, summary)
        return summary

    def export_curved_high_order(
        self, path: str | None = None, order: int | None = None, nodes: str = "equispaced"
    ) -> dict | None:
        """Curved order-p elements of the structured domains (blocks in 3D): VTKHDF Lagrange
        cells, or .h5 with the chosen node set."""
        if order is None:
            from PySide6.QtWidgets import QInputDialog

            order, ok = QInputDialog.getInt(
                self,
                "Curved high-order elements",
                "Order (interval counts must be multiples):",
                2,
                1,
                10,
            )
            if not ok:
                return None
        if path is None:
            path, chosen = QFileDialog.getSaveFileName(
                self,
                "Export curved high-order elements",
                "",
                "VTKHDF Lagrange cells (*.vtkhdf);;HDF5, Gauss-Lobatto nodes (*.h5);;"
                "HDF5, Gauss nodes (*.h5)",
            )
            if not path:
                return None
            if "Lobatto" in chosen:
                nodes = "gll"
            elif "Gauss nodes" in chosen:
                nodes = "gauss"
        summary = self._run(
            lambda: self.session.export_curved_high_order(path, order=order, nodes=nodes)
        )
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
        connectors = self._selected_of_kind("connector") + self._selected_of_kind("domain")
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
                if "connector" not in item:
                    continue
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

    # --------------------------------------------------------------- sketches
    def _after_sketch(self, name: str) -> None:
        """After a sketch changed (its warnings are logged by the session): redraw it and smooth
        its group (live progress)."""
        entity = self.session.store.get(name)
        ids = [entity.id] + [i for _role, i in entity.params["outputs"]]
        self.scene.update_entities(ids)
        self.model_tree.populate(self.session)
        for role, entity_id in entity.params["outputs"]:
            if role == "group":
                self.run_smoother(self.session.store.get(entity_id).name)

    def assign_selected_geometry(self) -> str | None:
        from opengrid.model.connectors import CURVE_KINDS
        from opengrid.model.entities import SURFACE_KINDS

        sketch = self._one("sketch", "sketch (and the curves or surfaces to assign)")
        if sketch is None:
            return None
        geometry = [
            n
            for n in self.selected_names()
            if self.session.store.get(n).kind in (*CURVE_KINDS, *SURFACE_KINDS)
        ]
        if not geometry:
            self.statusBar().showMessage("Select curves or surfaces with the sketch", 5000)
            return None
        name = self._run(lambda: self.session.assign_sketch_geometry(sketch, geometry))
        if name:
            log.info("%s: assigned %s", name, ", ".join(geometry))
            self._after_sketch(name)
        return name

    def _selected_curves(self) -> list[str]:
        from opengrid.model.connectors import CURVE_KINDS

        return [n for n in self.selected_names() if self.session.store.get(n).kind in CURVE_KINDS]

    def template_from_selection(self, kind: str) -> str | None:
        """A 2D template from the selected curves (the body, or the boundary)."""
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select the curves of the body first", 5000)
            return None
        key = "body"
        if kind == "butterfly":
            key = "boundary"
        name = self._run(lambda: self.session.create_template(kind, **{key: curves}))
        if name:
            log.info("template %s: %s", kind, name)
            self._after_sketch(name)
        return name

    def auto_topology_from_selection(self) -> str | None:
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select the curves bounding the region first", 5000)
            return None
        name = self._run(lambda: self.session.auto_topology(curves))
        if name:
            log.info("automatic topology: %s", name)
            self.scene.rebuild()
            self._after_sketch(name)
        return name

    def regenerate_sketch(self) -> str | None:
        sketch = self._one("sketch", "sketch")
        if sketch is None:
            return None
        name = self._run(lambda: self.session.regenerate_sketch(sketch))
        if name:
            self._after_sketch(name)
        return name

    def sketch_report(self) -> str | None:
        sketch = self._one("sketch", "sketch")
        if sketch is None:
            return None
        text = self.session.sketch_report(sketch)
        for line in text.splitlines():
            log.info("%s", line)
        return text

    # ------------------------------------------------------- topology editing
    def _one(self, kind: str, label: str) -> str | None:
        names = self._selected_of_kind(kind)
        if len(names) != 1:
            self.statusBar().showMessage(f"Select one {label}", 5000)
            return None
        return names[0]

    def _two(self, kind: str, label: str) -> list[str] | None:
        names = self._selected_of_kind(kind)
        if len(names) != 2:
            self.statusBar().showMessage(f"Select two {label}", 5000)
            return None
        return names

    def _ask_int(self, title: str, label: str, low: int, high: int) -> int | None:
        value, ok = QInputDialog.getInt(self, title, label, (low + high) // 2, low, high)
        if not ok:
            return None
        return value

    def _after_topology(self, names) -> None:
        names = [n for n in names or [] if n]
        if names:
            log.info("topology edit: %s", ", ".join(names))
            ids = [self.session.store.get(n).id for n in names if n in self.session.store]
            self.set_selection(ids)

    def split_connector(self, index: int | None = None) -> list[str] | None:
        name = self._one("connector", "connector")
        if name is None:
            return None
        n = len(self.session.connector_nodes(name).points)
        if index is None:
            index = self._ask_int("Split connector", f"Node (1 … {n - 2}):", 1, n - 2)
            if index is None:
                return None
        names = self._run(lambda: self.session.split_connector(name, index))
        self._after_topology(names)
        return names

    def join_connectors(self) -> str | None:
        names = self._two("connector", "connectors sharing an end")
        if names is None:
            return None
        joined = self._run(lambda: self.session.join_connectors(*names))
        self._after_topology([joined])
        return joined

    def split_domain(self, direction: str | None = None, index: int | None = None):
        name = self._one("domain", "structured domain")
        if name is None:
            return None
        ni, nj = self.session.domain_dimensions(name)
        if direction is None:
            direction, ok = QInputDialog.getItem(
                self, "Split domain", "Along grid line:", ["i", "j"], 0, False
            )
            if not ok:
                return None
        size = ni
        if direction == "j":
            size = nj
        if index is None:
            index = self._ask_int("Split domain", f"{direction} (1 … {size - 2}):", 1, size - 2)
            if index is None:
                return None
        names = self._run(lambda: self.session.split_domain(name, **{direction: index}))
        self._after_topology(names)
        return names

    def join_domains(self) -> str | None:
        names = self._two("domain", "structured domains sharing an edge")
        if names is None:
            return None
        joined = self._run(lambda: self.session.join_domains(*names))
        self._after_topology([joined])
        return joined

    def split_block(self, axis: str | None = None, index: int | None = None):
        name = self._one("block", "block")
        if name is None:
            return None
        nk, nj, ni = self.session.block_grid(name).shape[:3]
        if axis is None:
            axis, ok = QInputDialog.getItem(self, "Split block", "Axis:", ["i", "j", "k"], 0, False)
            if not ok:
                return None
        size = {"i": ni, "j": nj, "k": nk}[axis]
        if index is None:
            index = self._ask_int("Split block", f"{axis} (1 … {size - 2}):", 1, size - 2)
            if index is None:
                return None
        names = self._run(lambda: self.session.split_block(name, axis, index))
        self._after_topology(names)
        return names

    def join_blocks(self) -> str | None:
        names = self._two("block", "blocks sharing a face")
        if names is None:
            return None
        joined = self._run(lambda: self.session.join_blocks(*names))
        self._after_topology([joined])
        return joined

    def ogrid(self, fraction: float | None = None) -> list[str] | None:
        domains = self._selected_of_kind("domain")
        blocks = self._selected_of_kind("block")
        if len(domains) + len(blocks) != 1:
            self.statusBar().showMessage("Select one structured domain or block", 5000)
            return None
        if fraction is None:
            fraction, ok = QInputDialog.getDouble(
                self, "O-grid", "Core size (fraction of the outer size):", 0.5, 0.05, 0.95, 2
            )
            if not ok:
                return None
        if domains:
            names = self._run(lambda: self.session.ogrid_domain(domains[0], fraction=fraction))
        else:
            names = self._run(lambda: self.session.ogrid_block(blocks[0], fraction=fraction))
        self._after_topology(names)
        return names

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

    def import_cad(self, path: str | None = None, domains: str | None = None):
        """Import a STEP / IGES file (Session.import_cad); the report goes to the log."""
        if path is None:
            path, _filter = QFileDialog.getOpenFileName(
                self,
                "Import CAD",
                "",
                "CAD files (*.step *.stp *.iges *.igs);;All files (*)",
            )
            if not path:
                return None
        if domains is None:
            domains, ok = QInputDialog.getItem(
                self,
                "Import CAD",
                "Create faces as domains:",
                ["unstructured", "auto (structured if four-sided)", "none"],
                0,
                False,
            )
            if not ok:
                return None
            domains = domains.split(" ")[0]
        try:
            result = self.session.import_cad(path, domains=domains)
        except (ValidationError, RuntimeError) as error:
            self._report_error("Cannot import CAD", error)
            return None
        for line in result.report().splitlines():
            log.info("%s", line)
        self.scene.fit()
        return result

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
