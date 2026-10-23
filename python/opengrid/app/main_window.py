"""Main application window: viewport, model tree, property panel, log and Python console."""

import logging
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from PySide6.QtCore import QEvent, QSettings, QSize, Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QFrame,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMenu,
    QMessageBox,
    QScrollArea,
    QSpinBox,
    QTabWidget,
    QToolBar,
    QVBoxLayout,
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
from opengrid.app.console import make_console
from opengrid.app.icons import icon
from opengrid.app.log_panel import LogPanel
from opengrid.app.message_bar import MessageBar
from opengrid.app.model_tree import ModelTree
from opengrid.app.property_panel import PropertyPanel
from opengrid.app.quality_dock import METRICS, QualityDock, worst_cell
from opengrid.app.scene import Scene, pick_entity
from opengrid.app.script_panel import ScriptPanel
from opengrid.app.smoother_dock import SmootherDock
from opengrid.app.tools import TOOLS, Tool, snap, snap_targets
from opengrid.app.trace import SessionTrace, new_trace_path
from opengrid.app.zoom_box import ZoomBox
from opengrid.model.curves import CurveIndex
from opengrid.model.entities import SURFACE_KINDS
from opengrid.model.patches import PATCH_TYPES
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
HOVER_DELAY_MS = 30  # the highlight follows the pointer after it rests this long
UPDATE_CHECK_MS = 30000  # how often to look for updated code
SETTINGS_ORGANIZATION = "OpenGrid"
SETTINGS_APPLICATION = "OpenGrid"
RECENT_FILES = 8
TOOLBAR_ICON_SIZE = 20
WINDOW_STATE_KEY = "window/state_v3"  # bumped when the default toolbar or dock layout changes
PROJECT_FILTER = "OpenGrid project (*.ogp)"
DEFAULT_DIMENSION = 21
# A 3D boundary report slower than this is not rerun after each change (the Boundary panel
# then waits for Refresh).
SLOW_REPORT_SECONDS = 1.0
CONSOLE_BANNER = (
    "OpenGrid console. og: this window's session (every change shows in the window), "
    "window: the window, opengrid: the package, trace_path: the Python of everything done "
    "in this window (the Script tab). Example: og.create_line((0, 0), (1, 0))"
)
CURVE_KINDS = ("line", "arc", "polyline", "spline", "nurbs")
MESHED_KINDS = ("domain", "unstructured_domain", "block", "volume")
INSERTABLE_FIELD = {"polyline": "points", "spline": "points", "nurbs": "control_points"}


# Mouse bindings of the two navigation modes (pyvista's custom trackball style). 2D: left-drag
# zooms to a window (app/zoom_box.py; View menu: or pans), Shift + left-drag and middle-drag
# pan, wheel / two-finger scroll and right-drag zoom, Ctrl + left / middle-drag spin, Ctrl +
# right-drag rotates. 3D: left-drag rotates, Shift + left-drag and middle-drag pan, Alt + left-
# drag zooms to a window, wheel and right-drag zoom, Ctrl + left-drag spins.
NAVIGATION = {
    "2d": {
        "left": "pan",
        "shift_left": "pan",
        "control_left": "spin",
        "middle": "pan",
        "shift_middle": "dolly",
        "control_middle": "spin",
        "right": "dolly",
        "shift_right": "dolly",
        "control_right": "rotate",
    },
    "3d": {
        "left": "rotate",
        "shift_left": "pan",
        "control_left": "spin",
        "middle": "pan",
        "shift_middle": "dolly",
        "control_middle": "spin",
        "right": "dolly",
        "shift_right": "pan",
        "control_right": "rotate",
    },
}
NAVIGATION_HINTS = {
    "2d": "2D navigation: left-drag zooms to a window or pans, Shift + drag pans, wheel zooms; "
    "F re-centres",
    "3d": "3D navigation: left-drag rotates, Shift + drag pans, wheel zooms, Alt + drag zooms "
    "to a window; 1 / 2 / 3: the XY / XZ / YZ planes, 7: isometric, F: re-centre",
}
VIEW_NAMES = {
    "xy": "XY plane from +Z (top)",
    "-xy": "XY plane from -Z (bottom)",
    "xz": "XZ plane from -Y (front)",
    "-xz": "XZ plane from +Y (back)",
    "yz": "YZ plane from +X (side)",
    "-yz": "YZ plane from -X",
    "iso": "isometric",
}


def set_navigation_style(iren, mode: str = "2d") -> None:
    """Install the mouse bindings of a navigation mode (NAVIGATION) on an interactor."""
    iren.enable_custom_trackball_style(**NAVIGATION[mode])


def _setting_bool(value, default: bool) -> bool:
    """A QSettings value as a bool (ini files hand strings back)."""
    if value is None:
        return default
    if isinstance(value, str):
        return value.lower() in ("true", "1", "yes")
    return bool(value)


def _tooltip(action: QAction) -> None:
    """Tooltip from the action's text, shortcut and status tip."""
    text = action.text().replace("&", "")
    shortcut = action.shortcut().toString(QKeySequence.SequenceFormat.NativeText)
    if shortcut:
        text += f"  ({shortcut})"
    tip = action.statusTip()
    if tip:
        text += f"\n{tip}"
    action.setToolTip(text)


class _Batched:
    """The change of a batch of commands: the union of their entity ids."""

    def __init__(self, ids):
        self._ids = list(ids)

    def ids(self) -> list[int]:
        return self._ids


def _left_button_down() -> bool:
    """Whether the left mouse button is held now (Qt's view of the pointer)."""
    from PySide6.QtWidgets import QApplication

    return bool(QApplication.mouseButtons() & Qt.MouseButton.LeftButton)


def code_stamp() -> float:
    """Latest modification time of OpenGrid's code: its Python files and its extension modules
    (an editable install keeps those in site-packages)."""
    import sys

    package = Path(opengrid.__file__).resolve().parent
    paths = list(package.rglob("*.py"))
    for name in ("opengrid._core", "opengrid._occt"):
        module = sys.modules.get(name)
        if module is not None and getattr(module, "__file__", None):
            paths.append(Path(module.__file__))
    latest = 0.0
    for path in paths:
        try:
            latest = max(latest, path.stat().st_mtime)
        except OSError:
            continue  # removed meanwhile
    return latest


def version_label() -> str:
    """The version, and the source commit when running from a checkout ("+": with changes)."""
    from opengrid.revision import source_revision

    revision = source_revision()
    if revision is None:
        return opengrid.__version__
    return f"{opengrid.__version__} ({revision})"


class MainWindow(QMainWindow):
    def __init__(
        self,
        parent=None,
        off_screen: bool = False,
        session: Session | None = None,
        settings: QSettings | None = None,
    ):
        super().__init__(parent)
        if session is None:
            session = Session()
        self.session = session
        self._off_screen = off_screen
        # Window layout, view options, recent files and the last directory persist between
        # runs (an off-screen window, i.e. the tests, writes a throw-away file).
        if settings is None:
            if off_screen:
                import os
                import tempfile

                descriptor, name = tempfile.mkstemp(suffix=".ini", prefix="opengrid-")
                os.close(descriptor)
                settings = QSettings(name, QSettings.Format.IniFormat)
            else:
                settings = QSettings(SETTINGS_ORGANIZATION, SETTINGS_APPLICATION)
        self.settings = settings
        self._last_dir = ""
        self._recent: list[str] = []
        self.selected_ids: list[int] = []
        # Selected faces of blocks and volumes (patch members with ids: {"block": id, "face":
        # "imin"}, {"volume": id, "face": "top"}, a side wall {"connector": id}, ...).
        self.selected_faces: list[dict] = []
        self.tool: Tool = TOOLS["select"](session)
        self._saved_journal_length = 0
        self._press_position: tuple[int, int] | None = None
        self._pick_cycle: dict | None = None  # what is under the last click (front to back)
        self._drag = None  # handle drag in progress (app.handles.Drag)
        self._stroke: list | None = None  # freehand stroke in progress
        # Mouse bindings: "2d" or "3d" (NAVIGATION); the mode follows the model until the user
        # picks one (View menu). Alt is tracked here: the Qt -> VTK adapter forwards Ctrl and
        # Shift only.
        self.navigation_mode = "2d"
        self._navigation_chosen = False
        self._alt_pressed = False
        # Refreshes of the commands run inside batch_updates() are coalesced into one.
        self._batch_depth = 0
        self._batched: tuple[set[int], bool] | None = None  # (ids, whole model)
        self.resize(1400, 900)

        self.viewport = self._create_viewport(off_screen)
        # Warnings and errors are shown above the view (app/message_bar.py), not only logged.
        self.message_bar = MessageBar()
        central = QWidget()
        central_layout = QVBoxLayout(central)
        central_layout.setContentsMargins(0, 0, 0, 0)
        central_layout.setSpacing(0)
        central_layout.addWidget(self.message_bar)
        central_layout.addWidget(self.viewport, 1)
        self.setCentralWidget(central)
        self.zoom_box = ZoomBox(self.viewport.renderer, min_pixels=int(CLICK_MOVE_PIXELS) + 1)
        self.scene = Scene(self.viewport, session)

        self.model_tree = ModelTree()
        self.property_panel = PropertyPanel(session)
        self.log_panel = LogPanel()
        # Everything done in the window, as Python (app/trace.py); a window that is never
        # shown (the tests) writes it to a folder of its own, removed when it closes.
        self._trace_dir = None
        if off_screen:
            import tempfile

            self._trace_dir = Path(tempfile.mkdtemp(prefix="opengrid-trace-"))
            trace_path = self._trace_dir / "session.py"
        else:
            trace_path = new_trace_path()
        self.trace = SessionTrace(session, trace_path)
        self.script_panel = ScriptPanel(self.trace)
        self.script_panel.save_requested.connect(self.save_trace)
        self.script_panel.runner = self.run_script
        self.console = make_console(
            namespace={
                "opengrid": opengrid,
                "og": session,
                "window": self,
                "trace_path": str(trace_path),
            },
            banner=CONSOLE_BANNER,
        )

        model_dock = self._add_dock("Model", self.model_tree, Qt.DockWidgetArea.LeftDockWidgetArea)
        # Forms scroll inside their docks, so the docks' dividers can be moved freely.
        properties_dock = self._add_dock(
            "Properties", self.property_panel, Qt.DockWidgetArea.LeftDockWidgetArea, scroll=True
        )
        self.quality_dock = QualityDock()
        self.smoother_dock = SmootherDock(session)
        quality_dock = self._add_dock(
            "Quality", self.quality_dock, Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.smoother_dock_widget = self._add_dock(
            "Smoother", self.smoother_dock, Qt.DockWidgetArea.RightDockWidgetArea, scroll=True
        )
        self.boundary_dock = BoundaryDock()
        self.boundary_dock_widget = self._add_dock(
            "Boundary", self.boundary_dock, Qt.DockWidgetArea.RightDockWidgetArea, scroll=True
        )
        self.tabifyDockWidget(quality_dock, self.smoother_dock_widget)
        self.tabifyDockWidget(quality_dock, self.boundary_dock_widget)
        quality_dock.raise_()
        self.quality_dock_widget = quality_dock
        # The console, the script and the log share the lower right corner, under the panels
        # (the viewport keeps the whole height): tabs in one dock (a second group of tabified
        # docks leaves a stray tab bar over the viewport in Qt 6.11).
        self.console_tabs = QTabWidget()
        self.console_tabs.setDocumentMode(True)
        self.console_tabs.addTab(self.console, "Console")
        self.console_tabs.addTab(self.script_panel, "Script")
        self.console_tabs.addTab(self.log_panel, "Log")
        console_dock = self._add_dock(
            "Console", self.console_tabs, Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.console_dock_widget = console_dock
        self.resizeDocks([model_dock, properties_dock], [320, 480], Qt.Orientation.Vertical)
        self.resizeDocks([quality_dock, console_dock], [480, 360], Qt.Orientation.Vertical)
        self.resizeDocks([model_dock, quality_dock], [360, 460], Qt.Orientation.Horizontal)

        self._create_actions()
        self._create_menus()
        self._create_toolbars()
        self._connect_viewport_events()
        self._restore_settings()

        self.model_tree.selection_changed.connect(self.set_selection)
        self.model_tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.model_tree.customContextMenuRequested.connect(
            lambda position: self.show_context_menu(
                self.model_tree.viewport().mapToGlobal(position)
            )
        )
        self.coordinates_label = QLabel("")
        self.coordinates_label.setMinimumWidth(180)
        self.statusBar().addPermanentWidget(self.coordinates_label)
        self._last_pointer: tuple[float, float] | None = None
        self._right_press_position: tuple[float, float] | None = None
        self._report_seconds = 0.0  # how long the last 3D boundary report took
        self._long_running = False  # a worker thread runs a long operation (_run_long)
        self._boundary_timer = QTimer(self)  # the Boundary panel follows the model (coalesced)
        self._boundary_timer.setSingleShot(True)
        self._boundary_timer.setInterval(100)
        self._boundary_timer.timeout.connect(self._refresh_boundary)
        self._hover_timer = QTimer(self)
        self._hover_timer.setSingleShot(True)
        self._hover_timer.setInterval(HOVER_DELAY_MS)
        self._hover_timer.timeout.connect(self._update_hover)
        self.property_panel.message.connect(self._on_panel_message)
        self.property_panel.smooth_requested.connect(self.run_smoother)
        self.property_panel.points_requested.connect(self.give_points)
        self.property_panel.inflation_requested.connect(lambda name: self.inflation_on_selection())
        self.property_panel.points_preview.connect(self._show_edit_preview)
        self.dimension_box.valueChanged.connect(lambda _n: self._update_points_preview())
        self.smoother_dock.finished.connect(self._on_smoothing_finished)
        self.quality_dock.metric_changed.connect(self._on_metric_changed)
        self.quality_dock.colouring_toggled.connect(self.quality_action.setChecked)
        self.quality_dock.threshold_changed.connect(self.set_quality_threshold)
        self.quality_dock.locate_requested.connect(self.show_worst_cell)
        self.boundary_dock.type_requested.connect(self._set_type_from_dock)
        self.boundary_dock.pair_requested.connect(self.pair_selection_as_periodic)
        self.boundary_dock.clear_requested.connect(self.clear_boundary_type_on_selection)
        self.boundary_dock.zone_requested.connect(self.create_zone_from_selection)
        self.boundary_dock.refresh_requested.connect(self.validate_boundary)
        self.boundary_dock.entity_chosen.connect(lambda entity_id: self.set_selection([entity_id]))
        self.boundary_dock.face_chosen.connect(lambda face: self.set_selection([], faces=[face]))
        self.property_panel.entity_requested.connect(lambda i: self.set_selection([i]))
        session.commands.add_listener(self._on_model_changed)

        # The Log tab shows INFO and above, the bar WARNING and above; details (DEBUG) go to the
        # launcher's log file only.
        logging.getLogger("opengrid").addHandler(self.log_panel.handler)
        logging.getLogger("opengrid").addHandler(self.message_bar.handler)
        logging.getLogger("opengrid").setLevel(logging.DEBUG)
        self.message_bar.log_requested.connect(self.show_log)
        self._refresh_all()
        log.info("OpenGrid %s", version_label())
        log.debug("core %s", opengrid.core_version())
        self._code_stamp = code_stamp()
        self._update_timer = QTimer(self)
        self._update_timer.timeout.connect(self._check_for_update)
        self._update_timer.start(UPDATE_CHECK_MS)
        log.debug("%s", compute_summary())

    def show_log(self) -> None:
        """The Log tab, in front."""
        self.console_dock_widget.show()
        self.console_dock_widget.raise_()
        self.console_tabs.setCurrentWidget(self.log_panel)

    # ------------------------------------------------------------------ setup
    def _create_viewport(self, off_screen: bool):
        from pyvistaqt import QtInteractor

        viewport = QtInteractor(self, off_screen=off_screen)
        viewport.set_background("#2b3036", top="#59636e")
        viewport.add_axes()
        viewport.enable_parallel_projection()
        viewport.view_xy()
        if viewport.iren is not None:
            set_navigation_style(viewport.iren)
        return viewport

    def _add_dock(
        self, title: str, widget: QWidget, area: Qt.DockWidgetArea, scroll: bool = False
    ) -> QDockWidget:
        dock = QDockWidget(title, self)
        dock.setObjectName(f"dock_{title.lower()}")
        if scroll:
            area_widget = QScrollArea(dock)
            area_widget.setWidgetResizable(True)
            area_widget.setFrameShape(QFrame.Shape.NoFrame)
            area_widget.setWidget(widget)
            widget = area_widget
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        return dock

    def _action(self, text: str, slot, shortcut=None, tip: str = "") -> QAction:
        action = QAction(text, self)
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip:
            action.setStatusTip(tip)
        # QAction.triggered carries a "checked" bool: it must not land in the slot's first
        # optional parameter (a spacing, a thickness, a confirmation flag).
        action.triggered.connect(lambda _checked=False: slot())
        _tooltip(action)
        return action

    def _decorate(self) -> None:
        """Icons and tooltips (text, shortcut, what it does) of the toolbar actions."""
        icons = {
            "select": self.tool_actions["select"],
            "point": self.tool_actions["point"],
            "line": self.tool_actions["line"],
            "arc": self.tool_actions["arc"],
            "polyline": self.tool_actions["polyline"],
            "spline": self.tool_actions["spline"],
            "nurbs": self.tool_actions["nurbs"],
            "freehand": self.tool_actions["freehand"],
            "sketch": self.tool_actions["sketch"],
            "sketch-unstructured": self.sketch_unstructured_action,
            "finish": self.finish_action,
            "cancel": self.cancel_action,
            "connectors": self.connectors_action,
            "chain": self.chain_action,
            "structured": self.domain_action,
            "unstructured": self.unstructured_action,
            "network": self.network_action,
            "group": self.group_action,
            "block": self.block_action,
            "extrude": self.extrude_action,
            "volume": self.volume_action,
            "undo": self.undo_action,
            "redo": self.redo_action,
            "reset-camera": self.fit_action,
            "zoom-selection": self.zoom_selection_action,
            "view-xy": self.view_actions["xy"],
            "view-xz": self.view_actions["xz"],
            "view-yz": self.view_actions["yz"],
            "view-iso": self.view_actions["iso"],
            "navigation-3d": self.navigation_3d_action,
            "perspective": self.perspective_action,
            "zoom-window": self.zoom_window_action,
            "grid": self.grid_action,
        }
        for name, action in icons.items():
            action.setIcon(icon(name))
            _tooltip(action)

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
            "Points on curves",
            self.create_connectors_on_selection,
            tip="Give each selected curve without points the 'points' box's points",
        )
        self.chain_action = self._action(
            "Chain of curves",
            self.create_chain_connector_on_selection,
            tip="One curve of points along the selected curves (ordered end to end)",
        )
        self.color_action = QAction("Colour curves by stretching ratio", self)
        self.color_action.setCheckable(True)
        self.color_action.setChecked(True)
        self.color_action.toggled.connect(self.set_color_by_stretching)
        self.nodes_action = QAction("Show curve points", self)
        self.nodes_action.setCheckable(True)
        self.nodes_action.setChecked(True)
        self.nodes_action.toggled.connect(self.set_show_nodes)
        self.threshold_action = self._action("Stretching ratio threshold…", self.ask_threshold)
        self.unstructured_action = self._action(
            "Unstructured domain",
            self.create_unstructured_on_selection,
            tip="Triangulated domain bounded by the selected curves (closed loops; curves without "
            "points get the 'points' box's)",
        )
        self.inflation_action = self._action(
            "Inflation layers…",
            self.inflation_on_selection,
            tip="Boundary layers grown from the selected curves into the unstructured domains "
            "they bound (first height, growth, layers; running into the neighbouring curves)",
        )
        self.domain_action = self._action(
            "Structured domain",
            self.create_domain_on_selection,
            tip="Structured domain bounded by the selected curves (sides found automatically; "
            "curves without points get the 'points' box's)",
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
            "with four corners); point counts from the 'points' box",
        )
        self.group_action = self._action(
            "Smoothing group",
            self.create_group_on_selection,
            tip="Smooth the selected structured domains together as one grid (interfaces move); "
            "selected curves become walls",
        )
        self.volume_action = self._action(
            "Volume from faces",
            self.create_volume_on_selection,
            tip="Unstructured volume (tetrahedra) enclosed by the selected domains",
        )
        self.mesh_surfaces_action = self._action(
            "&Mesh surfaces…",
            self.mesh_selected_surfaces,
            tip="Surface meshes (unstructured domains) on the selected surfaces, their edges "
            "and curves shared, ready for Volume from faces",
        )
        self.periodic_pair_action = self._action(
            "Pair as periodic",
            self.pair_selection_as_periodic,
            tip="The two selected sides (two curves; 3D: two faces, or one extruded volume "
            "for its ends) become a periodic pair; two selected periodic patches are paired",
        )
        self.extrude_model_action = self._action(
            "Extrude 2D model (one cell thick)…",
            self.extrude_model,
            tip="Every domain extruded along z: blocks and volumes with the 2D patches on "
            "their sides and the front / back as a periodic pair (or symmetry planes)",
        )
        self.cut_action = QAction("Cut blocks (plane x)", self)
        self.cut_action.setCheckable(True)
        self.cut_action.toggled.connect(self.set_cut_plane)
        self.shrink_action = QAction("Shrink block cells", self)
        self.shrink_action.setCheckable(True)
        self.shrink_action.toggled.connect(self.set_shrink)
        self.quality_action = QAction("Colour the mesh by quality", self)
        self.quality_action.setStatusTip(
            "Colour the cells by the Quality panel's metric (there: only the cells worse than a "
            "value, also inside 3D blocks and volumes)"
        )
        self.quality_action.setCheckable(True)
        self.quality_action.toggled.connect(self.set_quality_colouring)
        self.patch_color_action = QAction("Colour curves by patch", self)
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

        self.fit_action = self._action(
            "&Reset camera",
            self.scene.reset_camera,
            QKeySequence("F"),
            "Centre the view on the whole model, keeping the viewing direction",
        )
        self.zoom_selection_action = self._action(
            "&Zoom to selection",
            self.zoom_to_selection,
            QKeySequence("Z"),
            "Centre the view on the selected entities (the whole model when nothing is selected)",
        )
        self.view_actions: dict[str, QAction] = {}
        for key, label, shortcut, tip in (
            ("xy", "XY plane (top, from +Z)", "1", "Look down the z axis; again: from -Z"),
            ("xz", "XZ plane (front, from -Y)", "2", "Look along the y axis; again: from +Y"),
            ("yz", "YZ plane (side, from +X)", "3", "Look along the x axis; again: from -X"),
            ("iso", "Isometric view", "7", "Look from the (1, 1, 1) corner"),
            ("-xy", "XY plane from -Z (bottom)", None, "Look up the z axis"),
            ("-xz", "XZ plane from +Y (back)", None, "Look along -y"),
            ("-yz", "YZ plane from -X (other side)", None, "Look along +x"),
        ):
            sequence = None
            if shortcut is not None:
                sequence = QKeySequence(shortcut)
            action = self._action(
                label, lambda _checked=False, k=key: self.view_from(k), sequence, tip
            )
            self.view_actions[key] = action
        self.perspective_action = QAction("&Perspective projection", self)
        self.perspective_action.setCheckable(True)
        self.perspective_action.setShortcut(QKeySequence("V"))
        self.perspective_action.setStatusTip(
            "Perspective instead of parallel projection (the 2D tools need the parallel view)"
        )
        self.perspective_action.toggled.connect(self.set_perspective)
        self.grid_action = QAction("Reference &grid", self)
        self.grid_action.setCheckable(True)
        self.grid_action.setShortcut(QKeySequence("G"))
        self.grid_action.setStatusTip(
            "Grid of the z = 0 plane at a step that follows the zoom, with the axes and, in a "
            "view along z, the coordinates along the edges"
        )
        self.grid_action.toggled.connect(self.set_grid)
        self.hide_action = self._action(
            "&Hide selected",
            self.hide_selected,
            QKeySequence("H"),
            tip="Take the selected entities out of the view (a far-field box around a body); "
            "Shift+H shows everything again",
        )
        self.show_all_action = self._action(
            "Show &all",
            self.show_all,
            QKeySequence("Shift+H"),
            tip="Show the hidden entities again",
        )
        self.navigation_3d_action = QAction("&3D navigation (left-drag rotates)", self)
        self.navigation_3d_action.setCheckable(True)
        self.navigation_3d_action.setStatusTip(
            "Left-drag rotates the view (Shift: pan, Alt: zoom to a window); off: the 2D "
            "bindings. Follows the model (blocks, volumes, surfaces) until set here."
        )
        self.navigation_3d_action.toggled.connect(self._on_navigation_toggled)
        # Single-key view and camera shortcuts work wherever the focus is in the window,
        # except in text fields (which take their own keys first).
        for action in (
            self.fit_action,
            self.zoom_selection_action,
            *self.view_actions.values(),
            self.perspective_action,
            self.grid_action,
        ):
            action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            self.addAction(action)
        self.high_order_action = QAction("Show high-order elements…", self)
        self.high_order_action.setCheckable(True)
        self.high_order_action.setStatusTip(
            "Overlay the curved (or straight-sided) high-order elements of the mesh"
        )
        self.high_order_action.toggled.connect(self._toggle_high_order)
        self._high_order_settings = (4, True)  # order, curved
        self.zoom_window_action = QAction("Left-drag zooms to a window", self)
        self.zoom_window_action.setCheckable(True)
        self.zoom_window_action.setChecked(True)
        self.zoom_window_action.setStatusTip(
            "Left-drag (three-finger drag) draws a window and zooms onto it; off: it pans. "
            "Shift + drag always pans."
        )

        self.tool_group = QActionGroup(self)
        self.tool_actions: dict[str, QAction] = {}
        tools = {
            "select": (
                "Select",
                "S",
                "Click to select, Shift + click to add; drag a handle to move it",
            ),
            "point": ("Point", "P", "Click to place a point"),
            "line": ("Line", "L", "Click the two end points"),
            "arc": ("Arc (3 points)", "A", "Click the start, a point on the arc and the end"),
            "polyline": ("Polyline", "W", "Click the points; Enter finishes, Esc cancels"),
            "spline": ("Spline", "B", "Click the points of a cubic spline; Enter finishes"),
            "nurbs": ("NURBS", "N", "Click the control points; Enter finishes"),
            "freehand": (
                "Freehand",
                "D",
                "Press and draw; release to fit curves (near the start: a loop)",
            ),
            "sketch": (
                "Sketch block",
                "K",
                "Click the four corners of a block of the sketch topology",
            ),
        }
        for name, (label, key, tip) in tools.items():
            action = QAction(label, self)
            action.setCheckable(True)
            action.setShortcut(QKeySequence(key))
            action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            action.setStatusTip(tip)
            action.triggered.connect(lambda _checked=False, n=name: self.set_tool(n))
            self.tool_group.addAction(action)
            self.tool_actions[name] = action
            self.addAction(action)
        self.tool_actions["select"].setChecked(True)
        self.sketch_unstructured_action = QAction("Unstructured sketch blocks", self)
        self.sketch_unstructured_action.setCheckable(True)
        self.sketch_unstructured_action.setStatusTip(
            "New sketch blocks are meshed unstructured on their curves (the sketch's walls "
            "on them become inflation layers) instead of structured"
        )
        self.sketch_unstructured_action.toggled.connect(self._on_sketch_unstructured)
        self._decorate()
        for action in self.actions_with_single_keys():
            self.model_tree.passthrough_keys.add(int(action.shortcut()[0].key()))

    def _create_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        for action in (self.new_action, self.open_action):
            file_menu.addAction(action)
        self.recent_menu = file_menu.addMenu("Open &recent")
        self._update_recent_menu()
        for action in (self.save_action, self.save_as_action):
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
        create_menu.addAction(self.sketch_unstructured_action)
        create_menu.addMenu(self.shapes_menu(create_menu))
        create_menu.addSeparator()
        create_menu.addAction(self.connectors_action)
        create_menu.addAction(self.chain_action)
        create_menu.addAction(self.domain_action)
        create_menu.addAction(self.unstructured_action)
        create_menu.addAction(self.inflation_action)
        create_menu.addAction(self.network_action)
        create_menu.addAction(self.group_action)
        create_menu.addAction(self.block_action)
        create_menu.addAction(self.extrude_action)
        create_menu.addAction(self.mesh_surfaces_action)
        create_menu.addAction(self.volume_action)
        create_menu.addSeparator()
        create_menu.addMenu(self.boundary_type_menu(create_menu))
        create_menu.addAction(self.periodic_pair_action)
        create_menu.addAction(self.extrude_model_action)
        create_menu.addSeparator()
        create_menu.addAction(self.finish_action)
        create_menu.addAction(self.cancel_action)

        topology_menu = self.menuBar().addMenu("&Topology")
        self.topology_actions = {
            "split_connector": self._action("Split curve at point…", self.split_connector),
            "join_connectors": self._action("Join two curves", self.join_connectors),
            "split_domain": self._action("Split structured domain…", self.split_domain),
            "join_domains": self._action("Join two structured domains", self.join_domains),
            "convert_domain": self._action(
                "Convert domain (structured <-> unstructured)",
                self.convert_selected_domain,
                tip="The selected domain changes kind on the same curves (an unstructured "
                "one needs a single loop of four edges)",
            ),
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
        view_menu.addAction(self.zoom_selection_action)
        view_menu.addSeparator()
        for key in ("xy", "xz", "yz", "iso"):
            view_menu.addAction(self.view_actions[key])
        opposite = view_menu.addMenu("From the other side")
        for key in ("-xy", "-xz", "-yz"):
            opposite.addAction(self.view_actions[key])
        view_menu.addSeparator()
        view_menu.addAction(self.navigation_3d_action)
        view_menu.addAction(self.zoom_window_action)
        view_menu.addAction(self.perspective_action)
        view_menu.addAction(self.grid_action)
        view_menu.addSeparator()
        view_menu.addAction(self.hide_action)
        view_menu.addAction(self.show_all_action)
        view_menu.addSeparator()
        view_menu.addAction(self.color_action)
        view_menu.addAction(self.nodes_action)
        view_menu.addAction(self.threshold_action)
        view_menu.addAction(self.quality_action)
        view_menu.addAction(self.patch_color_action)
        view_menu.addAction(self.high_order_action)
        view_menu.addSeparator()
        view_menu.addAction(self.cut_action)
        view_menu.addAction(self.shrink_action)

    def _create_toolbars(self) -> None:
        toolbar = QToolBar("Tools", self)
        toolbar.setObjectName("toolbar_tools")
        toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        toolbar.setIconSize(QSize(TOOLBAR_ICON_SIZE, TOOLBAR_ICON_SIZE))
        for action in self.tool_actions.values():
            toolbar.addAction(action)
        toolbar.addAction(self.sketch_unstructured_action)
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
        toolbar.addWidget(QLabel(" points "))
        self.dimension_box = QSpinBox()
        self.dimension_box.setRange(2, 100000)
        self.dimension_box.setValue(DEFAULT_DIMENSION)
        self.dimension_box.setToolTip(
            "Points given to curves that have none: Points on curves and Chain of curves (per "
            "curve), Structured / Unstructured domain and Domains from curves (per curve), "
            "Sketch block (per edge); the layers of Extrude domain. Change a curve's points "
            "later in its property panel."
        )
        toolbar.addWidget(self.dimension_box)
        toolbar.addSeparator()
        toolbar.addAction(self.undo_action)
        toolbar.addAction(self.redo_action)
        self.addToolBar(toolbar)
        view_toolbar = QToolBar("View", self)
        view_toolbar.setObjectName("toolbar_view")
        view_toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        view_toolbar.setIconSize(QSize(TOOLBAR_ICON_SIZE, TOOLBAR_ICON_SIZE))
        view_toolbar.addAction(self.fit_action)
        view_toolbar.addAction(self.zoom_selection_action)
        for key in ("xy", "xz", "yz", "iso"):
            view_toolbar.addAction(self.view_actions[key])
        view_toolbar.addSeparator()
        view_toolbar.addAction(self.navigation_3d_action)
        view_toolbar.addAction(self.zoom_window_action)
        view_toolbar.addAction(self.perspective_action)
        view_toolbar.addAction(self.grid_action)
        view_toolbar.addSeparator()
        view_toolbar.addAction(self.network_action)
        view_toolbar.addAction(self.group_action)
        self.addToolBarBreak()  # the view toolbar on its own row: the tools fill the first
        self.addToolBar(view_toolbar)
        self.view_toolbar = view_toolbar

    def _connect_viewport_events(self) -> None:
        # Double clicks come from Qt (VTK has no such event): add / remove curve points.
        self.viewport.installEventFilter(self)
        iren = self.viewport.iren
        if iren is None:
            return  # off-screen (tests): handlers are driven directly in world space
        set_navigation_style(iren, self.navigation_mode)
        interactor = iren.interactor
        # Clicks, moves and right clicks: passive observers, invoked even while the camera
        # style holds the interactor's focus (it grabs it on every press it sees).
        for event, handler in (
            ("LeftButtonPressEvent", self._on_vtk_press),
            ("LeftButtonReleaseEvent", self._on_vtk_release),
            ("MouseMoveEvent", self._on_vtk_move),
            ("RightButtonPressEvent", self._on_vtk_right_press),
            ("RightButtonReleaseEvent", self._on_vtk_right_release),
            ("LeaveEvent", self._on_vtk_leave),
        ):
            tag = interactor.AddObserver(event, handler)
            interactor.GetCommand(tag).SetPassiveObserver(True)
        # The reference grid follows the camera: after the style's moves (low priority) and
        # at the end of every interaction.
        for event in (
            "InteractionEvent",
            "EndInteractionEvent",
            "MouseWheelForwardEvent",
            "MouseWheelBackwardEvent",
        ):
            tag = interactor.AddObserver(event, self._on_camera_moved, -1.0)
            interactor.GetCommand(tag).SetPassiveObserver(True)
        # Drags of handles and freehand strokes take the left button before the camera does:
        # high-priority observers that stop the event when they use it.
        self._grab_tags = {}
        for event, handler in (
            ("LeftButtonPressEvent", self._grab_press),
            ("MouseMoveEvent", self._grab_move),
            ("LeftButtonReleaseEvent", self._grab_release),
        ):
            self._grab_tags[event] = interactor.AddObserver(event, handler, 10.0)
        # Zoom window: after handle drags and strokes (10), before the camera style (0).
        self._zoom_tags = {}
        for event, handler in (
            ("LeftButtonPressEvent", self._zoom_press),
            ("MouseMoveEvent", self._zoom_move),
            ("LeftButtonReleaseEvent", self._zoom_release),
        ):
            self._zoom_tags[event] = interactor.AddObserver(event, handler, 5.0)

    def _zoom_press(self, interactor, event) -> None:
        if self.navigation_mode == "3d":
            if not self._alt_pressed:
                return  # left-drag rotates; Alt + left-drag is the zoom window
        else:
            if not self.zoom_window_action.isChecked():
                return
            if interactor.GetShiftKey() or interactor.GetControlKey():
                return  # Shift: pan, Ctrl: spin (the camera style)
        x, y = interactor.GetEventPosition()
        self.zoom_box.start(x, y)
        # The press is used here: the camera style must not see it, or it starts a pan and
        # grabs the interactor's focus, after which only it gets the moves and the release
        # (the box would never finish). The click logic still gets the press position.
        self._press_position = self._event_position()
        interactor.GetCommand(self._zoom_tags[event]).SetAbortFlag(1)

    def _zoom_move(self, interactor, event) -> None:
        if not self.zoom_box.active:
            return
        x, y = interactor.GetEventPosition()
        if not _left_button_down():
            # The release never came (e.g. outside the window): end the box as a release.
            self._end_zoom(x, y)
            return
        self.zoom_box.update(x, y)
        self.viewport.render()
        interactor.GetCommand(self._zoom_tags[event]).SetAbortFlag(1)  # no hover meanwhile

    def _zoom_release(self, interactor, event) -> None:
        if not self.zoom_box.active:
            return
        x, y = interactor.GetEventPosition()
        self._end_zoom(x, y)
        # Not aborted: a short drag is still a click (_on_vtk_release decides).

    def _end_zoom(self, x: float, y: float) -> None:
        self.zoom_box.finish(x, y)
        self.scene.camera_changed()  # zoomed, or the outline removed; the grid follows

    def _on_camera_moved(self, *_args) -> None:
        if self.scene.grid.enabled:
            self.scene.grid.update()

    def hide_selected(self) -> None:
        """View > Hide selected (H): the selection leaves the view until Show all (Shift+H);
        hidden entities are neither drawn nor picked, so what they enclosed can be clicked."""
        ids = [i for i in self.selected_ids if i in self.session.store]
        if not ids:
            self.statusBar().showMessage(self._selection_hint("the entities to hide"), 5000)
            return
        self.set_selection([])
        self.scene.set_hidden(self.scene.hidden | set(ids))
        self.statusBar().showMessage(
            f"{len(ids)} hidden ({len(self.scene.hidden)} in all): Shift+H shows everything", 6000
        )

    def show_all(self) -> None:
        count = len(self.scene.hidden)
        if count:
            self.scene.set_hidden(set())
        self.statusBar().showMessage(f"{count} entities shown again", 4000)

    def set_grid(self, enabled: bool) -> None:
        self.scene.grid.set_enabled(enabled)
        if enabled:
            self.statusBar().showMessage(
                f"Reference grid: major lines every {self.scene.grid.step:g}; in 2D, clicks "
                "snap to its minor lines (G toggles it)",
                6000,
            )

    def cancel_zoom_box(self) -> None:
        """Drop a zoom window in progress (Esc, tool change)."""
        if self.zoom_box.active or self.zoom_box.actor.GetVisibility():
            self.zoom_box.cancel()
            self._press_position = None
            self.viewport.render()

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
            # Released within the pick tolerance of the start: a closed loop.
            close_tolerance = PICK_PIXELS * self.scene.world_per_pixel()
            created = self._run(
                lambda: self.tool.create_stroke(stroke, close_tolerance=close_tolerance)
            )
            if created:
                log.info("freehand: %s", ", ".join(created))
                # The new curves are the selection: "Points on curves" or
                # "Unstructured domain" apply to them directly.
                self.set_selection([self.session.store.get(n).id for n in created])
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
        curves = self._selected_unmeshed_curves()  # (its points come from the network)
        if not curves:
            if self._selected_curves():
                self._has_points_message()
                return []
            self.statusBar().showMessage(self._selection_hint("the curves of the network"), 5000)
            return []
        total = sum(self.session.length(c) for c in curves)
        spacing = total / len(curves) / max(self.dimension_box.value() - 1, 1)
        created = self._run(lambda: self.session.domains_from_curves(curves, spacing=spacing))
        if created:
            log.info("created %s", ", ".join(created))  # (the session logs its warnings)
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
        self.handle_click_at(x, y, additive)

    def _on_vtk_move(self, *_args) -> None:
        x, y = self._event_position()
        self._last_pointer = (x, y)
        self.handle_move(self.scene.display_to_world(x, y))
        self._update_coordinates(x, y)
        if self._hover_wanted():
            self._hover_timer.start()
        elif self.scene.hovered is not None:
            self.scene.set_hover(None)

    def _on_vtk_leave(self, *_args) -> None:
        self._hover_timer.stop()
        self.scene.set_hover(None)
        self.coordinates_label.setText("")

    def _on_vtk_right_press(self, *_args) -> None:
        self._right_press_position = self._event_position()

    def _on_vtk_right_release(self, *_args) -> None:
        """A right click (no drag: right-drag zooms) opens the context menu."""
        if self._right_press_position is None:
            return
        x, y = self._event_position()
        x0, y0 = self._right_press_position
        self._right_press_position = None
        if abs(x - x0) + abs(y - y0) > CLICK_MOVE_PIXELS:
            return
        from PySide6.QtGui import QCursor

        self.show_context_menu(QCursor.pos())

    # ------------------------------------------------------------ hover, coordinates
    def _hover_wanted(self) -> bool:
        return (
            self.tool.name == "select"
            and self._drag is None
            and self._stroke is None
            and not self.zoom_box.active
            and not _left_button_down()
        )

    def _update_hover(self) -> None:
        """Highlight the entity under the pointer (run shortly after the last move)."""
        if not self._hover_wanted() or self._last_pointer is None:
            self.scene.set_hover(None)
            return
        x, y = self._last_pointer
        entity_id = self.pick_at(x, y)
        self.scene.set_hover(entity_id)

    def pick_at(self, x: float, y: float) -> int | None:
        """The entity under display pixel (x, y): along the view ray in a 3D view, in the
        z = 0 plane otherwise (Scene.pick_display / pick_entity)."""
        if self._ray_picking():
            entity_id, _point = self.scene.pick_display(x, y)
            if entity_id is None or entity_id not in self.session.store:
                return None
        else:
            tolerance = PICK_PIXELS * self.scene.world_per_pixel()
            name = pick_entity(self.session, self.scene.display_to_world(x, y), tolerance)
            if name is None:
                return None
            entity_id = self.session.store.get(name).id
        # A curve's shape and its points are one curve: the points stand for it.
        return CurveIndex(self.session.store).representative(entity_id)

    def _update_coordinates(self, x: float, y: float) -> None:
        """World coordinates of the pointer in the status bar (2D views: the z = 0 plane)."""
        if self._ray_picking():
            self.coordinates_label.setText("")
            return
        wx, wy, _wz = self.scene.display_to_world(x, y)
        self.coordinates_label.setText(f"x {wx:.6g}   y {wy:.6g}")

    # ------------------------------------------------------------- context menu
    def context_menu(self) -> QMenu:
        """The right-click menu: the actions that apply to the selection, topology edits,
        zoom and delete, then the views."""
        from PySide6.QtWidgets import QMenu

        menu = QMenu(self)
        kinds = {self.session.store.get(i).kind for i in self.selected_ids}
        curves = bool(kinds & set(CURVE_KINDS))
        connectors = "connector" in kinds
        domains = bool(kinds & {"domain", "unstructured_domain"})
        structured = "domain" in kinds
        groups = [
            (self.connectors_action, curves),
            (self.chain_action, curves),
            (self.network_action, curves),
            (self.unstructured_action, curves or connectors),
            (self.domain_action, connectors or curves),
            (self.inflation_action, connectors or (domains and self.session.is_3d())),
            (self.extrude_action, structured or "unstructured_domain" in kinds),
            (self.group_action, structured),
            (self.block_action, structured),
            (self.volume_action, domains),
            (self.topology_actions["convert_domain"], domains and len(self.selected_ids) == 1),
        ]
        for action, applies in groups:
            if applies:
                menu.addAction(action)
        if any(applies for _action, applies in groups):
            menu.addSeparator()
        members = self._patch_members_from_selection()
        if members:
            menu.addMenu(self.boundary_type_menu(menu))
            if len(members) == 2 or len(self._selected_of_kind("patch")) == 2:
                menu.addAction(self.periodic_pair_action)
            menu.addSeparator()
        elif len(self._selected_of_kind("patch")) == 2:
            menu.addAction(self.periodic_pair_action)
            menu.addSeparator()
        topology = menu.addMenu("Topology")
        for action in self.topology_actions.values():
            topology.addAction(action)
        menu.addSeparator()
        if self.selected_ids:
            menu.addAction(self.zoom_selection_action)
            menu.addAction(self.delete_action)
            menu.addSeparator()
        menu.addAction(self.fit_action)
        views = menu.addMenu("View from")
        for key in ("xy", "xz", "yz", "iso", "-xy", "-xz", "-yz"):
            views.addAction(self.view_actions[key])
        menu.addAction(self.navigation_3d_action)
        return menu

    def show_context_menu(self, global_position) -> None:
        menu = self.context_menu()
        menu.exec(global_position)

    # ---------------------------------------------------- world-space handlers
    def snap_point(self, xyz):
        """Nearest point, curve end or sketch vertex within a few pixels; otherwise, with the
        reference grid shown in 2D, the nearest grid point (minor lines), so drawn blocks
        have exact, equal sides (periodic sides must be images of each other)."""
        from opengrid.app.tools import SnapTarget

        tolerance = SNAP_PIXELS * self.scene.world_per_pixel()
        snapped, target = snap(xyz, snap_targets(self.session), tolerance)
        grid = self.scene.grid
        if target is None and grid.enabled and self.navigation_mode == "2d":
            step = grid.step / max(grid.divisions, 1)
            point = (
                float(round(xyz[0] / step) * step),
                float(round(xyz[1] / step) * step),
                float(xyz[2]),
            )
            return point, SnapTarget(point, f"grid ({step:g})")
        return snapped, target

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt virtual)
        if watched is self.viewport and event.type() in (
            QEvent.Type.MouseButtonPress,
            QEvent.Type.MouseButtonRelease,
            QEvent.Type.MouseMove,
        ):
            self._alt_pressed = bool(event.modifiers() & Qt.KeyboardModifier.AltModifier)
        if (
            watched is self.viewport
            and event.type() == QEvent.Type.MouseButtonDblClick
            and event.button() == Qt.MouseButton.LeftButton
        ):
            position = event.position()
            scale = self.viewport.ren_win.GetSize()[0] / max(self.viewport.width(), 1)
            x = position.x() * scale
            y = (self.viewport.height() - 1 - position.y()) * scale
            self.handle_double_click(self.scene.display_to_world(x, y))
            return True  # not a press for VTK (the click before it did the selecting)
        return super().eventFilter(watched, event)

    def handle_double_click(self, xyz) -> str | None:
        """Select tool: a double-click on a handle removes that point (undo restores it); on a
        polyline, spline or NURBS curve - or a connector lying on one - it adds a point there
        (a control point by knot insertion: the curve keeps its shape). Returns the curve."""
        from opengrid.app.handles import pick_handle

        if self.tool.name != "select":
            return None
        tolerance = PICK_PIXELS * self.scene.world_per_pixel()
        handle = pick_handle(self.visible_handles(), xyz, tolerance)
        if handle is not None and handle.field in ("points", "control_points"):
            name = self._run(lambda: self.session.remove_handle(handle.entity, handle.index))
            if name:
                self._after_curve_edit(name, f"Removed point {handle.index} of {name}")
            return name
        picked = pick_entity(self.session, xyz, tolerance)
        curve = self._curve_under(picked, xyz)
        if curve is None:
            self.statusBar().showMessage(
                "Double-click a polyline, spline or NURBS curve to add a point there, or one "
                "of its handles to remove it",
                5000,
            )
            return None
        name = self._run(lambda: self.session.insert_handle(curve, xyz))
        if name:
            if curve not in self.selected_names() and picked not in self.selected_names():
                self.set_selection([self.session.store.get(curve).id])
            entity = self.session.store.get(curve)
            count = len(entity.params[INSERTABLE_FIELD[entity.kind]])
            self._after_curve_edit(name, f"Added a point to {name} ({count} now)")
        return name

    def _curve_under(self, picked: str | None, xyz) -> str | None:
        """The editable curve a double-click lands on: the picked curve, or the closest of the
        curves a picked connector runs on."""
        from opengrid.model.topology_edit import INSERTABLE

        if picked is None:
            return None
        entity = self.session.store.get(picked)
        if entity.kind in INSERTABLE:
            return picked
        if entity.kind != "connector":
            return None
        best = None
        best_distance = None
        for segment in entity.params["segments"]:
            curve = self.session.store.get(segment["curve"])
            if curve.kind not in INSERTABLE:
                continue
            distance = float(self.session.geometry(curve.name).project(tuple(xyz))[2])
            if best_distance is None or distance < best_distance:
                best = curve.name
                best_distance = distance
        return best

    def _after_curve_edit(self, name: str, message: str) -> None:
        from opengrid.app.handles import dependents_closure

        log.info("%s", message)  # (the session logs its warnings)
        self.scene.update_entities(dependents_closure(self.session.store, name))
        self._show_handles()
        self.statusBar().showMessage(f"{message}; undo restores the curve", 5000)

    def handle_click_at(self, x: float, y: float, additive: bool = False) -> None:
        """A click at display pixel (x, y). Selection in a 3D view (3D navigation, or a view
        not looking along z) picks what is drawn under the pointer (Scene.pick_display); the
        2D tools and the 2D selection work in the z = 0 plane."""
        if self.tool.name == "select" and self._ray_picking():
            # Clicking again at the same place (same view) takes the next thing behind: the
            # faces hidden by others, then the blocks and volumes (pick_stack).
            camera = self.viewport.renderer.GetActiveCamera()
            view = (tuple(camera.GetPosition()), tuple(camera.GetFocalPoint()))
            cycle = self._pick_cycle
            same = (
                cycle is not None
                and not additive
                and cycle["view"] == view
                and abs(cycle["x"] - x) <= CLICK_MOVE_PIXELS
                and abs(cycle["y"] - y) <= CLICK_MOVE_PIXELS
            )
            if same:
                stack = cycle["stack"]
                index = (cycle["index"] + 1) % max(len(stack), 1)
            else:
                stack = self.scene.pick_stack(x, y)
                index = 0
            self._pick_cycle = {"x": x, "y": y, "view": view, "stack": stack, "index": index}
            if stack:
                entity_id, face = stack[index]
                if entity_id not in self.session.store:
                    self._pick_cycle = None
                    return
                if face is not None:
                    self._select_face(face, entity_id, additive, (index, len(stack)))
                    return
                self._select_picked(self.session.store.get(entity_id).name, additive)
                if len(stack) > 1:
                    self.statusBar().showMessage(
                        f"{index + 1} of {len(stack)} under the pointer: click again for the "
                        "next one behind",
                        6000,
                    )
                return
            self._select_picked(None, additive)
            return
        self.handle_click(self.scene.display_to_world(x, y), additive)

    def _ray_picking(self) -> bool:
        """Whether clicks pick along the view ray: in 3D navigation, or whenever the view does
        not look along z (the 2D picker works in the z = 0 plane)."""
        if self.navigation_mode == "3d":
            return True
        direction = self.viewport.renderer.GetActiveCamera().GetDirectionOfProjection()
        return abs(direction[2]) < 0.999

    def handle_click(self, xyz, additive: bool = False) -> None:
        """A click at world position ``xyz`` in the z = 0 plane (from the viewport in a 2D
        view, or from tests)."""
        if self.tool.name == "select":
            tolerance = PICK_PIXELS * self.scene.world_per_pixel()
            name = pick_entity(self.session, xyz, tolerance)
            self._select_picked(name, additive)
            return
        snapped, target = self.snap_point(xyz)
        if target is not None:
            log.debug("snapped to %s", target.label)
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

    def _select_face(self, face: dict, owner: int, additive: bool, place=(0, 1)) -> None:
        """A click on a face of a block or volume: that face (a boundary to give a type);
        Shift adds or removes it. ``place``: (index, count) among what is under the pointer."""
        if additive:
            faces = list(self.selected_faces)
            if face in faces:
                faces.remove(face)
            else:
                faces.append(face)
            self.set_selection(self.selected_ids, faces=faces)
            return
        self.set_selection([], faces=[face])
        from opengrid.model.patches import member_label

        index, count = place
        self.statusBar().showMessage(
            f"Face {member_label(face, self.session.store)} ({index + 1} of {count} under the "
            "pointer: click again for the next one behind): give it a boundary type in its "
            "panel, the Boundary panel or the right-click menu",
            8000,
        )

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
        self.cancel_zoom_box()
        self.tool.cancel()
        self.tool = TOOLS[name](self.session)
        if name == "sketch":
            sketches = self._selected_of_kind("sketch")
            if not sketches:
                sketches = self.session.names("sketch")  # one sketch: blocks go to it
            if len(sketches) == 1:
                self.tool.sketch = sketches[0]
            self.tool.points_per_edge = self.dimension_box.value()
            self.tool.unstructured = self.sketch_unstructured_action.isChecked()
        self._stroke = None
        self._drag = None
        self._show_handles()
        self.tool_actions[name].setChecked(True)
        self.scene.show_preview([])
        self.scene.show_snap(None)
        self._update_tool_feedback(None)

    def _on_sketch_unstructured(self, checked: bool) -> None:
        if self.tool.name == "sketch":
            self.tool.unstructured = checked

    def finish_tool(self) -> None:
        created = self._run(self.tool.finish)
        if created:
            log.info("created %s", created)
        self._update_tool_feedback(None)

    def cancel_tool(self) -> None:
        if self.zoom_box.active or self.zoom_box.actor.GetVisibility():
            self.cancel_zoom_box()  # Esc drops a zoom window first
            return
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
    def set_selection(self, ids, faces=None) -> None:
        """Select entities (ids) and faces of blocks and volumes (``faces``: patch members
        with ids; none by default)."""
        # A curve's shape and its points are one curve, selected as its points (the panel
        # shows both; model/curves.py).
        index = CurveIndex(self.session.store)
        chosen = []
        for entity_id in ids:
            if entity_id in self.session.store:
                representative = index.representative(entity_id)
                if representative not in chosen:
                    chosen.append(representative)
        ids = chosen
        self.selected_ids = list(ids)
        store = self.session.store
        self.selected_faces = [
            dict(f) for f in faces or [] if all(v in store for k, v in f.items() if k != "face")
        ]
        self.scene.set_selection(ids)
        self.scene.set_face_selection(self.selected_faces)
        if self.scene.hovered in ids:
            self.scene.set_hover(None)
        self._show_handles()
        self.model_tree.select_ids(ids, self.selected_faces)
        if len(ids) == 1 and not self.selected_faces:
            self.property_panel.show_entity(ids[0])
        elif not ids and len(self.selected_faces) == 1:
            self.property_panel.show_face(self.selected_faces[0])
        else:
            self.property_panel.show_entity(None)
        self.delete_action.setEnabled(bool(ids))
        self._update_points_preview()
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

    @contextmanager
    def batch_updates(self):
        """Coalesce the refreshes of the commands run inside into one (the actions that run
        several commands: connectors on curves, imports, drags)."""
        self._batch_depth += 1
        try:
            yield
        finally:
            self._batch_depth -= 1
            if self._batch_depth == 0 and self._batched is not None:
                ids, whole = self._batched
                self._batched = None
                if whole:
                    self._on_model_changed(None)
                else:
                    self._on_model_changed(_Batched(sorted(ids)))

    def _on_model_changed(self, change) -> None:
        if self._batch_depth > 0:
            ids, whole = self._batched or (set(), False)
            if change is None:
                whole = True
            else:
                ids.update(change.ids())
            self._batched = (ids, whole)
            return
        had_geometry = self.scene.has_geometry()
        self._pick_cycle = None  # what lies under the last click may have changed
        if self.scene.patch_faces is not None:
            if self.session.is_3d():
                self.scene.patch_faces = self._patch_face_colors(previous=self.scene.patch_faces)
            else:
                self.scene.patch_faces = None
                self.scene.show_patch_faces(None)
        if change is None:
            self.scene.rebuild()
        else:
            self.scene.update_entities(change.ids())
            self._report_connector_warnings(change.ids())
            self._report_domain_warnings(change.ids())
        self.model_tree.populate(self.session)
        if self.scene.patch_colors is not None:
            self.scene.patch_colors = self._patch_colors()
            self.scene.update_entities([])
        self.set_selection(self.selected_ids, faces=self.selected_faces)
        self._update_edit_actions()
        self._auto_navigation()
        if not had_geometry and len(self.session.store) > 0:
            self.scene.fit()
        if self.high_order_action.isChecked():
            self.refresh_high_order()
        self._boundary_timer.start()

    def _refresh_boundary(self) -> None:
        """The Boundary panel after a change: recomputed (a 3D report reuses the meshes drawn
        in the view), unless the last 3D report was slow: then it waits for Refresh. Not while
        a long operation runs in the worker thread (it would mesh the same volumes at once)."""
        if self._long_running:
            self._boundary_timer.start()
            return
        if self.session.is_3d() and self._report_seconds > SLOW_REPORT_SECONDS:
            self.boundary_dock.mark_stale()
            return
        try:
            self.validate_boundary()
        except (ValidationError, ValueError, RuntimeError) as error:
            log.warning("boundary check: %s", error)

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

    def _report_domain_warnings(self, ids) -> None:
        """What the mesher of a changed unstructured domain reported, and inflation layers that
        did not all grow (and why)."""
        from opengrid.model.unstructured import layer_shortfall

        for entity_id in ids:
            if entity_id not in self.session.store:
                continue
            entity = self.session.store.get(entity_id)
            if entity.kind != "unstructured_domain":
                continue
            try:
                mesh = self.session.unstructured_mesh(entity.name)
            except (ValidationError, ValueError, RuntimeError):
                continue  # (its error is shown where it is drawn)
            for warning in mesh["warnings"]:
                log.warning("%s: %s", entity.name, warning)
            text = layer_shortfall(self.session, entity, mesh)
            if text is not None:
                log.warning("%s", text)

    def _refresh_all(self) -> None:
        self.scene.rebuild()
        self.model_tree.populate(self.session)
        self.set_selection([])
        self._update_edit_actions()
        self._update_tool_feedback(None)
        self._auto_navigation()
        action = getattr(self, "high_order_action", None)
        if action is not None and action.isChecked():
            self.refresh_high_order()

    # ------------------------------------------------------ high-order overlay
    def _toggle_high_order(self, on: bool) -> None:
        if not on:
            self.scene.show_high_order(None)
            return
        if self.sender() is self.high_order_action and not self._high_order_prompted:
            from PySide6.QtWidgets import QInputDialog

            order, ok = QInputDialog.getInt(
                self, "High-order elements", "Order:", self._high_order_settings[0], 1, 6
            )
            if ok:
                choices = ["curved boundary elements", "straight-sided elements"]
                choice, ok = QInputDialog.getItem(
                    self, "High-order elements", "Elements:", choices, 0, False
                )
                if ok:
                    self._high_order_settings = (order, choice == choices[0])
            if not ok:
                self.high_order_action.blockSignals(True)
                self.high_order_action.setChecked(False)
                self.high_order_action.blockSignals(False)
                return
        self.refresh_high_order()

    _high_order_prompted = False  # True while show_high_order sets the settings itself
    _last_high_order = None

    def show_high_order(self, order: int = 4, curved: bool = True):
        """Overlay the order-p elements (curved onto the geometry, or straight-sided); the
        overlay follows every model change until hidden. Returns the mesh, or None when it
        cannot be built (reported)."""
        self._high_order_settings = (order, curved)
        self._high_order_prompted = True
        try:
            if self.high_order_action.isChecked():
                return self.refresh_high_order()
            self.high_order_action.setChecked(True)  # refreshes through _toggle_high_order
        finally:
            self._high_order_prompted = False
        return self._last_high_order

    def hide_high_order(self) -> None:
        self.high_order_action.setChecked(False)

    def refresh_high_order(self):
        import time

        order, curved = self._high_order_settings
        started = time.perf_counter()
        mesh = self._run_long(
            "Building the high-order elements…", lambda: self.session.high_order(order, curved)
        )
        self._last_high_order = mesh
        if mesh is None:
            self.scene.show_high_order(None)
            return None
        self.scene.show_high_order(mesh)
        report = mesh.curving
        message = f"order {order}: {mesh.elements} elements, {len(mesh.points)} nodes"
        if curved:
            message += f", {report['moved']} boundary nodes on the geometry"
            if report.get("min_scaled_jacobian") is not None:
                message += (
                    f", scaled Jacobian of the curved elements "
                    f"{report['min_scaled_jacobian']:.3f} "
                    f"(straight-sided {report['min_scaled_jacobian_straight']:.3f})"
                )
        else:
            message += ", straight-sided"
        message += f" ({time.perf_counter() - started:.2f} s)"
        self.statusBar().showMessage(message)
        log.info("high-order overlay: %s", message)
        return mesh

    # ------------------------------------------------------------- navigation
    def _model_is_3d(self) -> bool:
        if self.session.is_3d():
            return True
        return any(entity.kind in SURFACE_KINDS for entity in self.session.store)

    def _auto_navigation(self) -> None:
        """Follow the model with the navigation mode until the user chose one."""
        if self._navigation_chosen:
            return
        wanted = "2d"
        if self._model_is_3d():
            wanted = "3d"
        if wanted != self.navigation_mode:
            self.set_navigation_mode(wanted)
            if wanted == "3d" and self.scene.has_geometry():
                self.scene.fit()

    def _on_navigation_toggled(self, checked: bool) -> None:
        mode = "2d"
        if checked:
            mode = "3d"
        if mode != self.navigation_mode:
            self.set_navigation_mode(mode, chosen=True)
        else:
            self._navigation_chosen = True

    def set_navigation_mode(self, mode: str, chosen: bool = False) -> None:
        """Mouse bindings "2d" or "3d" (NAVIGATION); ``chosen``: the user's pick, kept
        whatever the model becomes."""
        if mode not in NAVIGATION:
            raise ValueError(f"unknown navigation mode {mode!r} (2d or 3d)")
        self.navigation_mode = mode
        if chosen:
            self._navigation_chosen = True
        iren = self.viewport.iren
        if iren is not None:
            set_navigation_style(iren, mode)
        self.scene.default_view = {"2d": "xy", "3d": "iso"}[mode]
        self.zoom_window_action.setEnabled(mode == "2d")
        self.navigation_3d_action.blockSignals(True)
        self.navigation_3d_action.setChecked(mode == "3d")
        self.navigation_3d_action.blockSignals(False)
        self.statusBar().showMessage(NAVIGATION_HINTS[mode], 8000)

    def actions_with_single_keys(self) -> list[QAction]:
        """The actions bound to one unmodified key (tools, views, camera)."""
        out = []
        for action in (
            *self.tool_actions.values(),
            self.fit_action,
            self.zoom_selection_action,
            *self.view_actions.values(),
            self.perspective_action,
            self.hide_action,
        ):
            sequence = action.shortcut()
            if (
                sequence.count() == 1
                and sequence[0].keyboardModifiers() == Qt.KeyboardModifier.NoModifier
            ):
                out.append(action)
        return out

    def view_from(self, key: str) -> None:
        """A view preset: "xy", "xz", "yz", "iso", or "-xy", "-xz", "-yz" for the opposite
        side; a plain preset pressed again turns to the opposite side."""
        negative = None
        name = key
        if key.startswith("-"):
            negative = True
            name = key[1:]
        self.scene.view_preset(name, negative)
        shown = self.scene.view_name()
        if shown is None:
            shown = "iso"
        message = f"View: {VIEW_NAMES.get(shown, shown)}"
        if self.navigation_mode == "3d" and shown != "iso":
            message += (
                " (3D navigation stays on: left-drag rotates; the orbit button gives the 2D "
                "bindings)"
            )
        self.statusBar().showMessage(message, 6000)

    def zoom_to_selection(self) -> None:
        if not self.selected_ids or not self.scene.zoom_to(self.selected_ids):
            self.scene.reset_camera()

    def set_perspective(self, enabled: bool) -> None:
        self.scene.set_perspective(enabled)

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
        self.setWindowTitle(f"{name}{dirty} — OpenGrid {version_label()}")

    def is_modified(self) -> bool:
        return len(self.session.commands.journal.lines) != self._saved_journal_length

    def _on_panel_message(self, text: str, level: str) -> None:
        """A message of the property panel: "status" (the status bar only), "info", "warning"
        or "error" (also logged; warnings and errors are shown above the view)."""
        if level == "error":
            log.error("%s", text)
        elif level == "warning":
            log.warning("%s", text)
        elif level == "info":
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

    def mesh_selected_surfaces(self, spacing: float | None = None) -> list[str] | None:
        """Create > Mesh surfaces: domains on the selected surfaces (a box, a cylinder, a
        wing), the spacing asked (a thirtieth of the selection's size suggested)."""
        surfaces = [
            n for n in self.selected_names() if self.session.store.get(n).kind in SURFACE_KINDS
        ]
        if not surfaces:
            self.statusBar().showMessage(
                self._selection_hint("the surfaces to mesh (a box's faces, a cylinder, a wing)"),
                6000,
            )
            return None
        if spacing is None:
            points = []
            for name in surfaces:
                geometry = self.session.geometry(name)
                points.extend(geometry.eval(u, v) for u in (0.0, 0.5, 1.0) for v in (0.0, 0.5, 1.0))
            extent = float(np.max(np.ptp(np.asarray(points, dtype=float), axis=0)))
            suggested = max(extent / 30.0, 1e-6)
            spacing, ok = QInputDialog.getDouble(
                self,
                "Mesh surfaces",
                f"Spacing along the edges and in the mesh ({len(surfaces)} surface(s)):",
                suggested,
                1e-12,
                1e12,
                6,
            )
            if not ok:
                return None
        with self.batch_updates():
            names = self._run(lambda: self.session.mesh_surfaces(surfaces, spacing=spacing))
            if names:
                self._run_long("Meshing the surfaces…", lambda: self._warm_caches(names))
        if not names:
            return None
        self.set_selection([self.session.store.get(n).id for n in names])
        log.info("surface meshes: %s", ", ".join(names))
        self.statusBar().showMessage(
            f"{len(names)} surface mesh(es) at spacing {spacing:g}, selected: Volume from faces "
            "encloses them",
            8000,
        )
        return names

    def create_volume_on_selection(self) -> str | None:
        faces = self._selected_of_kind("domain") + self._selected_of_kind("unstructured_domain")
        if len(faces) < 2:
            self.statusBar().showMessage("Select the domains enclosing the volume", 5000)
            return None
        with self.batch_updates():
            name = self._run(lambda: self.session.create_volume(faces))
            if name:
                self._run_long(f"Meshing {name}…", lambda: self._warm_caches([name]))
        if name:
            log.info("created %s: %s", name, self.session.volume_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def extrude_selection(self) -> str | None:
        """Extrude the selected domain a quarter of its size along its normal: a structured
        domain becomes a block (hexahedra), an unstructured one a volume of prisms and
        hexahedra; both are edited afterwards in the property panel."""
        domains = self._selected_of_kind("domain") + self._selected_of_kind("unstructured_domain")
        if len(domains) != 1:
            self.statusBar().showMessage(
                self._selection_hint(
                    "one domain to extrude (structured: a block; unstructured: "
                    "a volume of prisms and hexahedra)"
                ),
                5000,
            )
            return None
        structured = self.session.store.get(domains[0]).kind == "domain"
        if structured:
            grid = self.session.domain_grid(domains[0])
            points = grid.reshape(-1, 3)
            normal = np.cross(grid[0, -1] - grid[0, 0], grid[-1, 0] - grid[0, 0])
        else:
            mesh = self.session.unstructured_mesh(domains[0])
            points = np.asarray(mesh["points"], dtype=float)
            normal = np.zeros(3)
            triangles = np.asarray(mesh["triangles"], dtype=np.int64).reshape(-1, 3)
            if len(triangles):
                p = points[triangles]
                normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]).sum(axis=0)
            quads = np.asarray(mesh["quads"], dtype=np.int64).reshape(-1, 4)
            if len(quads):
                p = points[quads]
                normal = normal + np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]).sum(axis=0)
        size = float(np.linalg.norm(np.ptp(points, axis=0)))
        if np.linalg.norm(normal) == 0:
            normal = np.array([0.0, 0.0, 1.0])
        vector = tuple(float(v) for v in 0.25 * size * normal / np.linalg.norm(normal))
        dimension = self.dimension_box.value()
        if structured:
            name = self._run_meshed(
                lambda: self.session.create_extruded_block(
                    domains[0], "translate", vector=vector, dimension=dimension
                )
            )
            if name:
                log.info("created %s: %s", name, self.session.block_status(name))
        else:
            name = self._run_meshed(
                lambda: self.session.create_extruded_volume(
                    domains[0], "translate", vector=vector, dimension=dimension
                )
            )
            if name:
                log.info("created %s: %s", name, self.session.volume_status(name))
        if name:
            self.set_selection([self.session.store.get(name).id])
        return name

    def _patch_members_from_selection(self) -> list:
        """Connectors (2D boundaries, side walls of extruded volumes), structured domains
        (faces of blocks and volumes), extruded volumes (both their ends) and the selected
        faces of blocks and volumes."""
        from opengrid.model.faces import face_identifier

        members = self._selected_of_kind("connector") + self._selected_of_kind("domain")
        for face in self.selected_faces:  # by identifier: "block-1.imin" in the journal
            members.append(face_identifier(self.session.store, face))
        for name in self._selected_of_kind("volume"):
            if self.session.store.get(name).params["faces"] is None:
                members.append({"volume": name, "face": "bottom"})
                members.append({"volume": name, "face": "top"})
        return members

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
        """Structured domain bounded by the selected curves; those without points get the
        'points' box's, in the same undo step."""
        connectors = self._selected_of_kind("connector")
        curves = self._selected_unmeshed_curves()
        if not connectors and not curves:
            self.statusBar().showMessage(self._selection_hint("the curves of the domain"), 5000)
            return None
        if curves:
            dimension = self.dimension_box.value()
            name = self._run(
                lambda: self.session.structured_domain_from_curves(
                    curves, dimension, connectors=connectors
                )
            )
        else:
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
        """Unstructured domain bounded by the selected curves; those without points get the
        'points' box's, all in one undo step."""
        connectors = self._selected_of_kind("connector")
        curves = self._selected_unmeshed_curves()
        if not connectors and not curves:
            self.statusBar().showMessage(
                self._selection_hint("the curves of the boundary loops"), 5000
            )
            return None
        if curves:
            dimension = self.dimension_box.value()
            name = self._run(
                lambda: self.session.unstructured_domain_from_curves(
                    curves, dimension, connectors=connectors
                )
            )
        else:
            name = self._run(lambda: self.session.create_unstructured_domain(connectors=connectors))
        if name:
            log.info("created %s: %s", name, self.session.domain_status(name))
            self.set_selection([self.session.store.get(name).id])
        return name

    def run_smoother(self, name: str) -> bool:
        self.smoother_dock_widget.raise_()
        started = self.smoother_dock.start(name)
        if started:
            self.statusBar().showMessage(f"smoothing {name}…", 5000)
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
            status = self.session.group_status(name)
            log.debug("%s: %s", name, status)
            if "NOT converged" in status:
                log.warning(
                    "%s: the smoothing did not converge in %d iterations (the grid shown is its "
                    "last iterate; see the Smoother panel)",
                    name,
                    result.iterations,
                )
            else:
                log.info("%s smoothed (%d iterations)", name, result.iterations)
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
                        if d not in self.session.store:
                            continue  # regenerated meanwhile (a sketch block changed kind)
                        for edge in self.session.store.get(d).params["edges"]:
                            ids += [item["connector"] for item in edge]
                ids = sorted({i for i in ids if i in self.session.store})
            self.scene.update_entities(ids)
        self.model_tree.populate(self.session)
        self.property_panel.refresh()
        self.refresh_quality()

    def set_quality_colouring(self, enabled: bool) -> None:
        if enabled:
            self.quality_dock_widget.raise_()  # (its metric and threshold)
            self.scene.quality_metric = self.quality_dock.metric_key()
            self.scene.quality_threshold = self.quality_dock.threshold_value()
        else:
            self.scene.quality_metric = None
        self.quality_dock.colour.blockSignals(True)
        self.quality_dock.colour.setChecked(enabled)
        self.quality_dock.colour.blockSignals(False)
        self.scene.rebuild()

    def set_quality_threshold(self, value) -> None:
        """Only the cells worse than ``value`` filled (None: every cell coloured); turns the
        colouring on."""
        self.scene.quality_threshold = value
        if value is not None and not self.quality_action.isChecked():
            self.quality_action.setChecked(True)  # (rebuilds)
        elif self.scene.quality_metric is not None:
            self.scene.rebuild()
        self.refresh_quality()  # (its statistics count the cells worse than the threshold)

    def _on_metric_changed(self, key: str) -> None:
        if self.scene.quality_metric is not None:
            self.scene.quality_metric = key
            self.scene.quality_threshold = self.quality_dock.threshold_value()
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
            path, _filter = self._save_dialog(self, "Export VTKHDF", "", "VTKHDF (*.vtkhdf *.hdf)")
            if not path:
                return None
        summary = self._run_long("Exporting VTKHDF…", lambda: self.session.export_vtkhdf(path))
        if summary:
            self.trace.note(f"og.export_vtkhdf({str(path)!r})")
            log.info("exported %s: %s", path, summary)
        return summary

    def export_high_order(
        self,
        path: str | None = None,
        order: int | None = None,
        nodes: str = "equispaced",
        curved: bool | None = None,
    ) -> dict | None:
        """Order-p elements of every mesh cell (boundary elements curved onto the geometry, or
        straight-sided): Gmsh, CGNS, VTKHDF Lagrange cells, or .h5 with the chosen node set."""
        if order is None or curved is None:
            from PySide6.QtWidgets import QInputDialog

        if order is None:
            order, ok = QInputDialog.getInt(
                self, "High-order elements", "Order (CGNS: up to 4):", 2, 1, 6
            )
            if not ok:
                return None
        if curved is None:
            choices = ["curved boundary elements", "straight-sided elements"]
            choice, ok = QInputDialog.getItem(
                self, "High-order elements", "Elements:", choices, 0, False
            )
            if not ok:
                return None
            curved = choice == choices[0]
        if path is None:
            path, chosen = self._save_dialog(
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
        summary = self._run_long(
            "Exporting the high-order elements…",
            lambda: self.session.export_high_order(path, order=order, nodes=nodes, curved=curved),
        )
        if summary:
            self.trace.note(
                f"og.export_high_order({str(path)!r}, order={order}, nodes={nodes!r}, "
                f"curved={curved})"
            )
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
            path, chosen = self._save_dialog(
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
        summary = self._run_long(
            "Exporting the curved high-order elements…",
            lambda: self.session.export_curved_high_order(path, order=order, nodes=nodes),
        )
        if summary:
            self.trace.note(
                f"og.export_curved_high_order({str(path)!r}, order={order}, nodes={nodes!r})"
            )
            log.info("exported %s: %s", path, summary)
        return summary

    def export_plot3d(self, path: str | None = None, binary: bool = True) -> bool:
        if not self.session.names("domain"):
            self.statusBar().showMessage("No domains to export", 5000)
            return False
        if path is None:
            path, _filter = self._save_dialog(self, "Export Plot3D", "", "Plot3D (*.xyz *.x *.p3d)")
            if not path:
                return False
        self._run_long("Exporting Plot3D…", lambda: self.session.export_plot3d(path, binary=binary))
        self.trace.note(f"og.export_plot3d({str(path)!r}, binary={binary})")
        log.info("exported %s", path)
        return True

    # ------------------------------------------------------- patches and zones
    def _selected_of_kind(self, kind: str) -> list[str]:
        return [n for n in self.selected_names() if self.session.store.get(n).kind == kind]

    def convert_selected_domain(self) -> str | None:
        domains = self._selected_of_kind("domain") + self._selected_of_kind("unstructured_domain")
        if len(domains) != 1:
            self.statusBar().showMessage(self._selection_hint("one domain to convert"), 5000)
            return None
        store = self.session.store
        groups = [d.name for d in store.dependents(domains[0]) if d.kind == "smoothing_group"]
        from opengrid.model.sketch import output_owner

        sketched = store.get(domains[0]).id in output_owner(store)  # its sketch decides
        if groups and not sketched and not self._confirm_leave_groups(domains[0], groups):
            return None
        leave = bool(groups)
        name = self._run(lambda: self.session.convert_domain(domains[0], leave_groups=leave))
        if name:
            kind = self.session.store.get(name).kind.replace("_", " ")
            log.info("%s is now a %s", name, kind)
            self.set_selection([self.session.store.get(name).id])
            self.statusBar().showMessage(f"{name} is now a {kind}", 5000)
        return name

    def _confirm_leave_groups(self, domain: str, groups: list[str]) -> bool:
        """Ask before a converted domain leaves its smoothing groups; tests replace this."""
        answer = QMessageBox.question(
            self,
            "Convert domain",
            f"{domain} is smoothed in {', '.join(groups)}. Take it out of the group and "
            "convert it?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )
        return answer == QMessageBox.StandardButton.Yes

    def pair_selection_as_periodic(self) -> list[str] | None:
        """Two selected sides (two connectors; 3D: two faces, or one extruded volume for its
        two ends) become a periodic pair; two selected periodic patches are paired."""
        patches = [
            n
            for n in self._selected_of_kind("patch")
            if self.session.store.get(n).params["patch_type"] == "periodic"
        ]
        if len(patches) == 2:
            name = self._run(lambda: self.session.pair_periodic(patches[0], patches[1]))
            if name:
                self._after_pairing(name)
                return [patches[0], patches[1]]
            return None
        members = self._patch_members_from_selection()
        if len(members) != 2:
            self.statusBar().showMessage(
                self._selection_hint(
                    "the two sides: two curves (3D: two faces, or one extruded volume)"
                ),
                6000,
            )
            return None
        names = self._run(lambda: self.session.make_periodic_pair([members[0]], [members[1]]))
        if names:
            self._after_pairing(names[0])
        return names

    def _after_pairing(self, source: str) -> None:
        entity = self.session.store.get(source)
        transform = entity.params["periodic"]
        t = transform["translation"]
        log.info("%s paired with %s", source, self.session.store.get(transform["partner"]).name)
        self._show_boundary_result()
        self.set_selection([entity.id])  # the panel shows the pairing (axis, angle, ...)
        self.statusBar().showMessage(
            f"{source}: periodic with translation ({t[0]:.4g}, {t[1]:.4g}, {t[2]:.4g}); the "
            "property panel edits the pairing (rotation, translation)",
            8000,
        )

    def extrude_model(self, thickness: float | None = None, ends: str | None = None):
        """Create > Extrude 2D model: one cell thick along z, the fronts and backs periodic
        (or symmetry planes)."""
        from opengrid.app.property_panel import mean_connector_spacing

        if self.session.is_3d():
            self.statusBar().showMessage("The model is 3D already", 5000)
            return None
        if not self.session.names("domain") and not self.session.names("unstructured_domain"):
            self.statusBar().showMessage("No domains to extrude", 5000)
            return None
        if thickness is None:
            spacing = mean_connector_spacing(
                self.session, [e.id for e in self.session.store.entities("connector")]
            )
            suggested = 1.0
            if spacing is not None:
                suggested = spacing
            thickness, ok = QInputDialog.getDouble(
                self,
                "Extrude 2D model",
                "Thickness along z (one cell; the mean boundary spacing is suggested):",
                suggested,
                1e-12,
                1e12,
                6,
            )
            if not ok:
                return None
        if ends is None:
            ends, ok = QInputDialog.getItem(
                self,
                "Extrude 2D model",
                "Front and back:",
                ["periodic", "symmetry"],
                0,
                False,
            )
            if not ok:
                return None
        with self.batch_updates():
            names = self._run(lambda: self.session.extrude_model(thickness, ends=ends))
            if names:
                self._run_long("Meshing the extrusion…", lambda: self._warm_caches(names))
        if names:
            log.info("extruded the 2D model: %s", ", ".join(names))
            self.statusBar().showMessage(
                f"Extruded {len(names) - 2} domain(s) by {thickness:g}; front / back {ends}", 8000
            )
        return names

    def shapes_menu(self, parent=None) -> QMenu:
        """Create > Shapes: standard shapes from a small dialog each."""
        from opengrid.app.shapes_dialog import SHAPES

        menu = QMenu("&Shapes", parent)
        for kind, spec in SHAPES.items():
            action = menu.addAction(f"{spec.title}…")
            action.setToolTip(spec.hint)
            action.triggered.connect(lambda _checked=False, k=kind: self.create_shape(k))
        return menu

    def create_shape(self, kind: str, values: dict | None = None):
        """A standard shape (Create > Shapes), from the dialog or from ``values`` (the
        Session creator's arguments); the new entities are selected. Returns their names."""
        from PySide6.QtWidgets import QDialog

        from opengrid.app.shapes_dialog import SHAPES, ShapeDialog

        spec = SHAPES[kind]
        if values is None:
            dialog = ShapeDialog(spec, self)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return None
            values = dialog.values()
        method = getattr(self.session, spec.method)
        names = self._run(lambda: method(**values))
        if not names:
            return None
        if isinstance(names, str):
            names = [names]
        self.set_selection([self.session.store.get(n).id for n in names])
        log.info("%s: %s", spec.title, ", ".join(names))
        self.statusBar().showMessage(f"{spec.title}: {', '.join(names)} (selected)", 6000)
        return names

    def boundary_type_menu(self, parent=None) -> QMenu:
        """A submenu with one entry per patch type: the selection joins a patch of it."""
        menu = QMenu("Boundary type", parent)
        for patch_type in PATCH_TYPES:
            action = menu.addAction(patch_type)
            action.setStatusTip(
                f"The selected boundaries join a {patch_type} patch (made when there is none)"
            )
            action.triggered.connect(
                lambda _checked=False, t=patch_type: self.set_boundary_type_on_selection(t)
            )
        return menu

    def set_boundary_type_on_selection(
        self, patch_type: str, name: str | None = None
    ) -> str | None:
        """The selected connectors (3D: domains, extruded volumes) join a patch of
        ``patch_type`` (the one called ``name`` when given); the patch colouring is switched
        on and the Boundary panel raised to show it."""
        members = self._patch_members_from_selection()
        if not members:
            self.statusBar().showMessage(
                self._selection_hint("the boundary curves (3D: domains, block faces or volumes)"),
                5000,
            )
            return None
        name = self._run(lambda: self.session.set_boundary_type(members, patch_type, name=name))
        if not name:
            return None
        count = len(members)
        log.info("%s: %d boundary member(s) now %s", name, count, patch_type)
        self._show_boundary_result()
        self.statusBar().showMessage(
            f"{count} boundary member(s) in patch {name} ({patch_type})", 6000
        )
        return name

    def _set_type_from_dock(self, patch_type: str, name: str) -> None:
        if not name:
            name = None
        self.set_boundary_type_on_selection(patch_type, name)

    def clear_boundary_type_on_selection(self) -> list[str] | None:
        """The selected boundaries leave their patches (the Boundary panel's Clear type)."""
        members = self._patch_members_from_selection()
        if not members:
            self.statusBar().showMessage(
                self._selection_hint("the boundaries to clear (curves; 3D: faces)"), 5000
            )
            return None
        changed = self._run(lambda: self.session.clear_boundary_type(members))
        if changed is None:
            return None
        self.statusBar().showMessage(f"{len(members)} boundary member(s) now have no type", 6000)
        return changed

    def _show_boundary_result(self) -> None:
        """Patch colouring on and the Boundary panel in front, after a boundary edit."""
        if not self.patch_color_action.isChecked():
            self.patch_color_action.setChecked(True)
        self.boundary_dock_widget.raise_()

    def create_patch_from_selection(self, name: str, patch_type: str = "wall") -> str | None:
        connectors = self._patch_members_from_selection()
        if not name:
            self.statusBar().showMessage("Enter a patch name", 5000)
            return None
        created = self._run(lambda: self.session.create_patch(name, patch_type, connectors))
        if created:
            log.info("created patch %s (%s) with %d curve(s)", created, patch_type, len(connectors))
            self.validate_boundary()
        return created

    def add_selection_to_patch(self, name: str) -> bool:
        connectors = self._patch_members_from_selection()
        if not connectors:
            self.statusBar().showMessage(
                self._selection_hint("curves, domains or extruded volumes"), 5000
            )
            return False
        done = self._run(lambda: self.session.add_to_patch(name, connectors))
        self.validate_boundary()
        return bool(done)

    def remove_selection_from_patch(self, name: str) -> bool:
        connectors = self._patch_members_from_selection()
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
        import time

        started = time.perf_counter()
        report = self.session.boundary_report()
        if self.session.is_3d():
            self._report_seconds = time.perf_counter() - started
        self.boundary_dock.show_report(report, self.session.store)
        if self.scene.patch_colors is not None:
            self.scene.patch_colors = self._patch_colors(report)
            if self.session.is_3d():
                self.scene.patch_faces = self._patch_face_colors()
            self.scene.rebuild()
        return report

    def _patch_colors(self, report=None) -> dict[int, str]:
        """Connector colours: interface, unassigned boundary, or the patch's colour (3D: the
        side-wall connectors of extruded volumes)."""
        colors = {}
        if self.session.is_3d():
            from opengrid.model.volumes import analyze_3d

            for entry in analyze_3d(self.session)["entries"]:
                for key in entry["keys"]:
                    if key[0] == "connector":
                        colors[key[1]] = UNASSIGNED_COLOR
        else:
            topology = self.session.topology()
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
        self.scene.patch_faces = None
        if enabled:
            self.scene.patch_colors = self._patch_colors()
            if self.session.is_3d():
                self.scene.patch_faces = self._patch_face_colors()
        else:
            self.scene.patch_colors = None
            self.scene.show_patch_faces(None)
        self.scene.rebuild()

    def _patch_face_colors(self, previous=None) -> list:
        """3D: [(face, colour)] of every boundary face: its patch's colour, or the unassigned
        colour. The boundary faces are found by meshing the model (analyze_3d), or taken from
        ``previous`` (after an edit of the patches: cheap)."""
        from opengrid.api import _face_member
        from opengrid.model.volumes import analyze_3d

        coloured = []
        for index, patch in enumerate(self.session.store.entities("patch")):
            color = PATCH_PALETTE[index % len(PATCH_PALETTE)]
            for member in patch.params["members"]:
                coloured.append((dict(member), color))
        assigned = [m for m, _c in coloured]
        if previous is None:
            faces = [_face_member(e) for e in analyze_3d(self.session)["entries"]]
        else:
            faces = [m for m, _c in previous]
        for member in faces:
            if member not in assigned:
                assigned.append(member)
                coloured.append((member, UNASSIGNED_COLOR))
        return coloured

    def export_cgns(self, path: str | None = None) -> dict | None:
        if path is None:
            path, _filter = self._save_dialog(self, "Export CGNS", "", "CGNS (*.cgns)")
            if not path:
                return None
        summary = self._run_long("Exporting CGNS…", lambda: self.session.export_cgns(path))
        if summary:
            self.trace.note(f"og.export_cgns({str(path)!r})")
            log.info("exported %s: %s", path, summary)
        self.validate_boundary()
        return summary

    def export_gmsh(self, path: str | None = None) -> dict | None:
        if path is None:
            path, _filter = self._save_dialog(self, "Export Gmsh", "", "Gmsh (*.msh)")
            if not path:
                return None
        summary = self._run_long("Exporting Gmsh…", lambda: self.session.export_gmsh(path))
        if summary:
            self.trace.note(f"og.export_gmsh({str(path)!r})")
            log.info("exported %s: %s", path, summary)
        self.validate_boundary()
        return summary

    # -------------------------------------------------------------- connectors
    def create_connectors_on_selection(self) -> list[str]:
        """Points ('points' box) on each selected curve that has none yet."""
        curves = self._selected_unmeshed_curves()
        if not curves:
            if self._selected_curves():
                self._has_points_message()
                return []
            self.statusBar().showMessage(self._selection_hint("one or more curves"), 5000)
            return []
        names = self._run(
            lambda: self.session.create_connectors_on_curves(curves, self.dimension_box.value())
        )
        if names:
            log.info("created %s", ", ".join(names))
            self.set_selection([self.session.store.get(n).id for n in names])
        return names or []

    def _update_points_preview(self) -> None:
        """Where the selected curves without points would get them (the 'points' box),
        before a domain is made on them: pink dots."""
        import contextlib

        arrays = []
        for curve in self._selected_unmeshed_curves():
            with contextlib.suppress(ValidationError, ValueError):
                arrays.append(self.session.preview_points(curve, self.dimension_box.value()))
        self.scene.show_points_preview(arrays)
        if arrays:
            self.statusBar().showMessage(
                f"Preview: {self.dimension_box.value()} points on each selected curve (the "
                "'points' box); they get them with Structured / Unstructured domain or Points "
                "on curves",
                6000,
            )

    def _show_edit_preview(self, arrays) -> None:
        """The points a curve's panel edits would give, before Apply (or the plain preview
        again when the edits are undone)."""
        if arrays:
            self.scene.show_points_preview(arrays)
            return
        self._update_points_preview()

    def inflation_on_selection(self, values: dict | None = None) -> list[str] | None:
        """Inflation layers from the selected curves (with points) into the unstructured domains
        they bound: the dialog's settings, or ``values`` (for Session.set_inflation; with
        "remove": True they stop growing layers). One undo step."""
        from opengrid.app.inflation_dialog import InflationDialog
        from opengrid.model.layer_spacing import suggested_first_height

        curves = self._selected_of_kind("connector")
        if self.session.is_3d():  # 3D: surfaces (the domains bounding volumes)
            curves += self._selected_of_kind("unstructured_domain") + self._selected_of_kind(
                "domain"
            )
        if not curves:
            self.statusBar().showMessage(
                self._selection_hint(
                    "the curves the layers grow from (curves with points; 3D: surfaces)"
                ),
                6000,
            )
            return None
        if values is None:
            store = self.session.store
            current = None
            ids = {store.get(c).id for c in curves}
            for domain in store.entities("unstructured_domain"):
                inflation = domain.params["generator"]["inflation"]
                if inflation is not None and ids & set(inflation["connectors"]):
                    current = inflation
                    break
            from opengrid.model.topology_edit import _domain_connectors

            refs = []
            for name in curves:
                entity = store.get(name)
                if entity.kind == "connector":
                    refs.append(entity.id)
                else:
                    refs += _domain_connectors(entity)
            spacing = min(
                self.session.connectors.nodes(r).length
                / max(store.get(r).params["dimension"] - 1, 1)
                for r in refs
            )
            dialog = InflationDialog(curves, current, suggested_first_height(spacing), self)
            if dialog.exec() != dialog.DialogCode.Accepted:
                return None
            values = {"remove": True}
            if not dialog.removed:
                values = dialog.values()
        names = self._run(lambda: self.session.set_inflation(curves, **values))
        if names:  # (the session logs its warnings)
            text = f"layers from {', '.join(curves)} in {', '.join(names)}"
            if values.get("remove"):
                text = f"no layers from {', '.join(curves)} in {', '.join(names)}"
            log.info("%s", text)
            self.statusBar().showMessage(text, 6000)
        return names

    def give_points(self, curve: str) -> list[str]:
        """Points ('points' box) on one curve without points (its panel's button)."""
        names = self._run(
            lambda: self.session.create_connectors_on_curves([curve], self.dimension_box.value())
        )
        if names:
            log.info("%s: %d points", names[0], self.dimension_box.value())
            self.set_selection([self.session.store.get(names[0]).id])
        return names or []

    def create_chain_connector_on_selection(self) -> str | None:
        """One curve of points along the selected curves (without points yet), end to end."""
        curves = self._selected_unmeshed_curves()
        if not curves:
            if self._selected_curves():
                self.statusBar().showMessage(
                    "The selected curves have points already: join them instead "
                    "(Topology > Join two curves)",
                    6000,
                )
                return None
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
        """Delete the selection. A sketch's domain removes its block from the sketch (the
        sketch itself when every block goes); anything else takes what depends on it along
        after a confirmation that lists it; patches, zones and groups just lose members."""
        from opengrid.model.sketch import output_owner

        names = self.selected_names()
        if not names:
            return
        store = self.session.store
        owners = output_owner(store)
        blocks: dict[str, set[int]] = {}
        plain = []
        for name in names:
            entity = store.get(name)
            owner = owners.get(entity.id)
            if owner is None or owner[0].name in names:
                plain.append(name)
                continue
            sketch, role = owner
            if role[:1] == "d" and role[1:].isdigit():
                blocks.setdefault(sketch.name, set()).add(int(role[1:]))
            else:
                self.statusBar().showMessage(
                    f"{name} is generated by sketch {sketch.name}: delete the domain of the "
                    "block it belongs to, or the sketch",
                    8000,
                )
                return
        for sketch_name, indices in list(blocks.items()):
            if len(indices) >= len(store.get(sketch_name).params["blocks"]):
                plain.append(sketch_name)
                del blocks[sketch_name]
        # A curve is its points and its shape: both go (the shape is never listed alone).
        shown = list(plain)
        index = CurveIndex(store)
        for name in shown:
            curve = index.geometry_of.get(store.get(name).id)
            if curve is not None and curve not in owners and store.get(curve).name not in plain:
                plain.append(store.get(curve).name)
        extra = []
        if plain:
            try:
                extra = self.session.deletion_dependents(plain)
            except ValidationError as error:
                log.warning("%s", error)
                self.statusBar().showMessage(str(error).splitlines()[0], 8000)
                return
            if extra and not self._confirm_deletion(shown, extra):
                return
        with self.batch_updates():
            for sketch_name, indices in blocks.items():
                removed = self._run(
                    lambda s=sketch_name, i=indices: self.session.remove_sketch_blocks(s, sorted(i))
                )
                if removed:
                    log.info(
                        "%s: block(s) %s removed", sketch_name, ", ".join(map(str, sorted(indices)))
                    )
            count = 0
            if plain:
                count = self._run(lambda: self.session.delete(plain, dependents=True))
        if plain and not count:
            return
        text = f"Deleted {', '.join(shown)}"
        if extra:
            text += f" and {len(extra)} dependent entities"
        if blocks:
            text = "Removed " + ", ".join(
                f"block(s) {', '.join(map(str, sorted(i)))} of {s}" for s, i in blocks.items()
            )
            if plain:
                text += f"; deleted {', '.join(shown)}"
        log.info("%s", text)
        self.statusBar().showMessage(text, 8000)

    def _confirm_deletion(self, names, extra) -> bool:
        """Ask before entities that were not selected go too (the dependents of the
        selection); tests replace this."""
        shown = ", ".join(extra[:12])
        if len(extra) > 12:
            shown += f" and {len(extra) - 12} more"
        answer = QMessageBox.question(
            self,
            "Delete",
            f"Deleting {', '.join(names)} also deletes {len(extra)} entities that depend on "
            f"them:\n{shown}\n\nDelete them all?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        return answer == QMessageBox.StandardButton.Yes

    def _selection_hint(self, needed: str) -> str:
        """Status text when the selection does not fit an action: what is selected (and what
        it is), then what is needed."""
        names = self.selected_names()
        if not names:
            return f"Select {needed} first"
        shown = []
        for name in names[:3]:
            kind = self.session.store.get(name).kind.replace("_", " ")
            shown.append(f"{name} is a {kind}")
        if len(names) > 3:
            shown.append(f"{len(names) - 3} more")
        return f"{', '.join(shown)}: select {needed}"

    def _run_long(self, label: str, function, cancel=None):
        """Run a long operation in a worker thread behind a busy dialog (app/worker.py), the
        window refreshed once afterwards; errors reported as in _run."""
        from opengrid.app.worker import run_blocking

        self._hover_timer.stop()
        self.scene.set_hover(None)
        self._long_running = True
        try:
            with self.batch_updates():
                return run_blocking(self, label, function, cancel)
        except (ValidationError, KeyError, ValueError, RuntimeError, OSError) as error:
            text = str(error)
            log.warning("%s", text)
            self.statusBar().showMessage(text.splitlines()[0], 8000)
            return None
        finally:
            self._long_running = False

    def _run_meshed(self, function):
        """A command that creates a meshed entity: the mesh is computed in a worker thread
        before the window draws it."""
        with self.batch_updates():
            name = self._run(function)
            if name:
                self._run_long(f"Meshing {name}…", lambda: self._warm_caches([name]))
        return name

    def _warm_caches(self, names=None) -> None:
        """Evaluate the meshes the scene will draw (worker thread): domains, blocks, volumes;
        problems are left to the refresh, which reports them."""
        if names is None:
            names = [e.name for e in self.session.store if e.kind in MESHED_KINDS]
        for name in names:
            if name not in self.session.store:
                continue
            kind = self.session.store.get(name).kind
            try:
                if kind == "domain":
                    self.session.domain_grid(name, compute=False)
                elif kind == "unstructured_domain":
                    self.session.unstructured_mesh(name)
                elif kind == "block":
                    self.session.block_grid(name)
                elif kind == "volume":
                    self.session.volume_mesh(name)
            except (ValidationError, ValueError, RuntimeError):
                continue

    def _run(self, function):
        """Run a model operation (its commands refresh the window once); report validation
        errors instead of raising."""
        try:
            with self.batch_updates():
                return function()
        except (ValidationError, KeyError, ValueError) as error:
            text = str(error)
            log.error("%s", text)  # (shown above the view)
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
            if self.session.store.get(n).kind in (*CURVE_KINDS, "connector", *SURFACE_KINDS)
        ]
        if not geometry:
            self.statusBar().showMessage("Select curves or surfaces with the sketch", 5000)
            return None
        name = self._run_long(
            "Assigning the geometry…", lambda: self.session.assign_sketch_geometry(sketch, geometry)
        )
        if name:
            log.info("%s: assigned %s", name, ", ".join(geometry))
            self._after_sketch(name)
        return name

    def _selected_curves(self) -> list[str]:
        """The shapes (database curves) of the selected curves, with or without points."""
        from opengrid.model.connectors import CURVE_KINDS

        store = self.session.store
        index = CurveIndex(store)
        names = []
        for entity_id in self.selected_ids:
            entity = store.get(entity_id)
            if entity.kind in CURVE_KINDS:
                names.append(entity.name)
            elif entity_id in index.geometry_of:
                names.append(store.get(index.geometry_of[entity_id]).name)
        return names

    def _selected_unmeshed_curves(self) -> list[str]:
        """The selected curves that have no points yet (their database curves)."""
        from opengrid.model.connectors import CURVE_KINDS

        return [n for n in self.selected_names() if self.session.store.get(n).kind in CURVE_KINDS]

    def _has_points_message(self) -> None:
        self.statusBar().showMessage(
            "The selected curves have points already: change them in the Properties panel", 6000
        )

    def template_from_selection(self, kind: str) -> str | None:
        """A 2D template from the selected curves (the body, or the boundary)."""
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select the curves of the body first", 5000)
            return None
        key = "body"
        if kind == "butterfly":
            key = "boundary"
        name = self._run_long(
            "Building the template…", lambda: self.session.create_template(kind, **{key: curves})
        )
        if name:
            log.info("template %s: %s", kind, name)
            self._after_sketch(name)
        return name

    def auto_topology_from_selection(self) -> str | None:
        curves = self._selected_curves()
        if not curves:
            self.statusBar().showMessage("Select the curves bounding the region first", 5000)
            return None
        name = self._run_long("Computing the topology…", lambda: self.session.auto_topology(curves))
        if name:
            log.info("automatic topology: %s", name)
            self.scene.rebuild()
            self._after_sketch(name)
        return name

    def regenerate_sketch(self) -> str | None:
        sketch = self._one("sketch", "sketch")
        if sketch is None:
            return None
        name = self._run_long(
            "Regenerating the sketch…", lambda: self.session.regenerate_sketch(sketch)
        )
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
        name = self._one("connector", "curve with points")
        if name is None:
            return None
        n = len(self.session.connector_nodes(name).points)
        if index is None:
            index = self._ask_int("Split curve", f"At point (1 … {n - 2}):", 1, n - 2)
            if index is None:
                return None
        names = self._run(lambda: self.session.split_connector(name, index))
        self._after_topology(names)
        return names

    def join_connectors(self) -> str | None:
        names = self._two("connector", "curves (with points) sharing an end")
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

    def _check_for_update(self) -> None:
        """Say once, in the status bar, when OpenGrid's code changed on disk since it started:
        this window still runs the code it started with."""
        if self._code_stamp is None or code_stamp() <= self._code_stamp:
            return
        self._code_stamp = None
        notice = QLabel("OpenGrid was updated: restart it to use the new version")
        notice.setStyleSheet("color: #ffd166;")
        self.statusBar().addPermanentWidget(notice)
        log.info("OpenGrid's code changed on disk: restart OpenGrid to use the new version")

    def open_document(self, path: str) -> bool:
        """Open a project handed over by the system (Finder double-click, command line):
        unsaved changes are confirmed first."""
        if not self._confirm_discard():
            return False
        opened = self.open_project(path)
        self.raise_()
        self.activateWindow()
        return opened

    def open_project(self, path: str | None = None) -> bool:
        if path is None:
            if not self._confirm_discard():
                return False
            path, _filter = self._open_dialog(self, "Open project", "", PROJECT_FILTER)
            if not path:
                return False
        try:
            self.session.load(path)
        except (ValidationError, OSError) as error:
            self._report_error("Cannot open project", error)
            return False
        self._saved_journal_length = len(self.session.commands.journal.lines)
        self._run_long("Loading the project's meshes…", self._warm_caches)
        self._refresh_all()
        self.scene.fit()
        self._add_recent(path)
        log.info("opened %s", path)
        return True

    def save_project(self) -> bool:
        if self.session.project_path is None:
            return self.save_project_as()
        return self._save_to(self.session.project_path)

    def save_project_as(self, path: str | None = None) -> bool:
        if path is None:
            path, _filter = self._save_dialog(self, "Save project", "", PROJECT_FILTER)
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
        self.trace.note(f"og.save({str(path)!r})")
        self._saved_journal_length = len(self.session.commands.journal.lines)
        self._update_title()
        self._add_recent(path)
        log.info("saved %s", path)
        return True

    def import_airfoil(self, path: str | None = None, mode: str | None = None) -> list[str]:
        if path is None:
            path, _filter = self._open_dialog(
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
            path, _filter = self._open_dialog(
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
            with self.batch_updates():
                from opengrid.app.worker import run_blocking

                self._hover_timer.stop()
                result = run_blocking(
                    self, "Importing CAD…", lambda: self.session.import_cad(path, domains=domains)
                )
        except (ValidationError, RuntimeError) as error:
            self._report_error("Cannot import CAD", error)
            return None
        for line in result.report().splitlines():
            log.info("%s", line)
        self.scene.fit()
        return result

    def import_csv(self, path: str | None = None, mode: str | None = None) -> list[str]:
        if path is None:
            path, _filter = self._open_dialog(
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
            path, _filter = self._save_dialog(
                self, "Save journal", "journal.py", "Python script (*.py)"
            )
            if not path:
                return
        self.session.save_journal(path)
        log.info("journal saved to %s", path)

    def run_script(self, code: str) -> tuple[str, str]:
        """Run code written in the Script tab in this window's model, statement by statement,
        with the console's names (``og``, ``window`` and what was defined there). Each
        statement that runs joins the trace as written; the first that fails stops: it and
        the rest are returned (the code not run) with a message."""
        import ast
        import contextlib
        import io

        try:
            tree = ast.parse(code)
        except SyntaxError as error:
            return code, f"Not run: syntax error at line {error.lineno}: {error.msg}"
        lines = code.splitlines(keepends=True)
        namespace = self.console.user_namespace()
        start = 0  # where the next statement's text starts (the comments above go with it)
        ran = 0
        with self.batch_updates():
            for statement in tree.body:
                end = statement.end_lineno
                source = "".join(lines[start:end])
                printed = io.StringIO()
                try:
                    with (
                        self.trace.typed(source),
                        contextlib.redirect_stdout(printed),
                        contextlib.redirect_stderr(printed),
                    ):
                        exec(compile(ast.Module([statement], []), "<script>", "exec"), namespace)
                except Exception as error:
                    self._log_printed(printed.getvalue())
                    message = f"Stopped at line {statement.lineno}: {type(error).__name__}"
                    if str(error):
                        message += f": {str(error).splitlines()[0]}"
                    log.error("script: %s", message)
                    return "".join(lines[start:]), message
                self._log_printed(printed.getvalue())
                start = end
                ran += 1
        rest = "".join(lines[start:])
        if rest.strip():  # comments after the last statement
            self.trace.note(rest.strip("\n"))
        message = f"Ran {ran} statement(s)"
        if namespace.get("og") is not self.session:
            message += " (og no longer names this window's model: og = window.session)"
        log.info("script: %s", message)
        return "", message

    def _log_printed(self, text: str) -> None:
        for line in text.splitlines():
            log.info("script: %s", line)

    def save_trace(self, path: str | None = None) -> str | None:
        """Save the Python trace of this window (everything done in it, File > New and Open
        included) as a script."""
        if path is None:
            path, _filter = self._save_dialog(
                self, "Save the script", "opengrid_session.py", "Python script (*.py)"
            )
            if not path:
                return None
        Path(path).write_text(self.trace.text(), encoding="utf-8")
        log.info("script saved to %s", path)
        return str(path)

    def _report_error(self, title: str, error: Exception) -> None:
        dialog = self.isVisible() and not self._off_screen
        log.error("%s: %s", title, error, extra={"on_screen": not dialog})
        if dialog:
            QMessageBox.critical(self, title, str(error))

    # ----------------------------------------------------------------- settings
    def _restore_settings(self) -> None:
        s = self.settings
        geometry = s.value("window/geometry")
        if geometry is not None:
            self.restoreGeometry(geometry)
        state = s.value(WINDOW_STATE_KEY)
        if state is not None:
            self.restoreState(state)
        self.zoom_window_action.setChecked(_setting_bool(s.value("view/zoom_window"), True))
        self.nodes_action.setChecked(_setting_bool(s.value("view/nodes"), True))
        self.color_action.setChecked(_setting_bool(s.value("view/stretching"), True))
        self.grid_action.setChecked(_setting_bool(s.value("view/grid"), False))
        dimension = s.value("tools/dimension")
        if dimension is not None:
            self.dimension_box.setValue(int(dimension))
        self._last_dir = str(s.value("files/last_dir", ""))
        recent = s.value("files/recent", [])
        if isinstance(recent, str):
            recent = [recent]
        self._recent = [str(p) for p in (recent or []) if p]
        self._update_recent_menu()

    def _save_settings(self) -> None:
        s = self.settings
        s.setValue("window/geometry", self.saveGeometry())
        s.setValue(WINDOW_STATE_KEY, self.saveState())
        s.setValue("view/zoom_window", self.zoom_window_action.isChecked())
        s.setValue("view/nodes", self.nodes_action.isChecked())
        s.setValue("view/stretching", self.color_action.isChecked())
        s.setValue("view/grid", self.grid_action.isChecked())
        s.setValue("tools/dimension", self.dimension_box.value())
        s.setValue("files/last_dir", self._last_dir)
        s.setValue("files/recent", self._recent)
        s.sync()

    def _add_recent(self, path) -> None:
        path = str(Path(path).resolve())
        self._recent = [path, *[p for p in self._recent if p != path]][:RECENT_FILES]
        self._update_recent_menu()

    def _update_recent_menu(self) -> None:
        menu = getattr(self, "recent_menu", None)
        if menu is None:
            return
        menu.clear()
        for path in self._recent:
            action = menu.addAction(Path(path).name)
            action.setStatusTip(path)
            action.setToolTip(path)
            action.triggered.connect(lambda _checked=False, p=path: self.open_document(p))
        menu.setEnabled(bool(self._recent))

    def _open_dialog(
        self, _parent, title: str, _directory: str = "", filter: str = "", *args, **kwargs
    ):
        """QFileDialog.getOpenFileName from the last directory used, which it remembers."""
        path, chosen = QFileDialog.getOpenFileName(
            self, title, self._last_dir, filter, *args, **kwargs
        )
        if path:
            self._last_dir = str(Path(path).parent)
        return path, chosen

    def _save_dialog(
        self, _parent, title: str, suggested: str = "", filter: str = "", *args, **kwargs
    ):
        """QFileDialog.getSaveFileName in the last directory used (``suggested``: a file name
        there), which it remembers."""
        start = self._last_dir
        if suggested:
            start = str(Path(self._last_dir) / Path(suggested).name)
        path, chosen = QFileDialog.getSaveFileName(self, title, start, filter, *args, **kwargs)
        if path:
            self._last_dir = str(Path(path).parent)
        return path, chosen

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt override)
        self._save_settings()
        logging.getLogger("opengrid").removeHandler(self.log_panel.handler)
        logging.getLogger("opengrid").removeHandler(self.message_bar.handler)
        if self._on_model_changed in self.session.commands._listeners:
            self.session.commands.remove_listener(self._on_model_changed)
        if self.trace._on_change in self.session.commands._listeners:
            self.session.commands.remove_listener(self.trace._on_change)
        self.console.shutdown()
        if self._trace_dir is not None:
            import shutil

            shutil.rmtree(self._trace_dir, ignore_errors=True)
        self.viewport.close()
        super().closeEvent(event)
