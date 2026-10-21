"""Boundary panel: every boundary with its patch and type, the periodic pairs and whether
they match, the interfaces (shared, no type needed) and the empty patches; one status line
says whether the boundary is ready to export. The selection gets a type (a named patch when a
name is given), becomes a periodic pair, or loses its type. The panel follows a 2D model at
once; a 3D model's report meshes every block and volume, so it waits for Refresh.
"""

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QGridLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from opengrid.model.patches import PATCH_TYPES, partner_of

# Categorical colours for patches (cycled); unassigned boundary connectors and interfaces have
# their own colours.
PATCH_PALETTE = (
    "#4c8bf5",
    "#2ee59d",
    "#ffb86b",
    "#c678dd",
    "#56b6c2",
    "#e5c07b",
    "#98c379",
    "#f78c6c",
)
UNASSIGNED_COLOR = "#ff3b30"
INTERFACE_COLOR = "#7f8c8d"
OK_COLOR = "#2ee59d"
ENTITY_ROLE = Qt.ItemDataRole.UserRole + 1
COLUMNS = ("boundary", "patch", "type", "status")


def transform_text(pairing: dict) -> str:
    """'translation (2, 0, 0)' or 'rotation 60° about (0, 1, 0)'."""
    if pairing["rotation_angle"] != 0.0:
        axis = ", ".join(f"{c:.4g}" for c in pairing["rotation_axis"])
        return f"rotation {math.degrees(pairing['rotation_angle']):.4g}\u00b0 about ({axis})"
    translation = ", ".join(f"{c:.4g}" for c in pairing["translation"])
    return f"translation ({translation})"


def match_text(matched) -> str:
    """What a consistent pair matched: connectors node by node (2D) or nodes and faces (3D)."""
    nodes = getattr(matched, "nodes", None)
    if nodes is not None:
        return f"match: {nodes} nodes, {matched.faces} faces"
    count = len(matched)
    unit = "curves"
    if count == 1:
        unit = "curve"
    return f"match: {count} {unit} node by node"


class BoundaryDock(QWidget):
    type_requested = Signal(str, str)  # patch type, patch name ("" = the patch of that type)
    pair_requested = Signal()
    clear_requested = Signal()
    zone_requested = Signal(str)
    refresh_requested = Signal()
    entity_chosen = Signal(int)  # a row clicked: the entity to select

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.status = QLabel("No boundaries yet")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        # The controls in three short rows (the panel must not need to be wide):
        # type + name, the three actions on the selection, zone + refresh.
        grid = QGridLayout()
        self.patch_type = QComboBox()
        self.patch_type.addItems(list(PATCH_TYPES))
        self.patch_type.setToolTip(
            "The boundary type the selected boundaries get (wall, inlet, ...); periodic sides "
            "are paired with Pair as periodic instead"
        )
        self.name = QLineEdit()
        self.name.setPlaceholderText("patch name (optional)")
        self.name.setToolTip(
            "Leave empty to use the patch named after the type; give a name to keep several "
            "patches of one type apart (wall-top, wall-bottom)"
        )
        self.set_type = QPushButton("Set type")
        self.set_type.setToolTip("The selected boundaries join a patch of this type")
        self.pair = QPushButton("Pair as periodic")
        self.pair.setToolTip(
            "The two selected sides (two curves; 3D: two faces, or one extruded volume for "
            "its ends) become a periodic pair"
        )
        self.clear = QPushButton("Clear type")
        self.clear.setToolTip("The selected boundaries leave their patches")
        self.new_zone = QPushButton("New zone")
        self.new_zone.setToolTip("A named group of the selected domains (the name field names it)")
        self.refresh = QPushButton("Refresh")
        self.refresh.setToolTip(
            "Check the boundary again (a 2D model is checked at every change; a 3D model's "
            "check meshes every block and volume, so it runs on request)"
        )
        grid.addWidget(QLabel("Type"), 0, 0)
        grid.addWidget(self.patch_type, 0, 1)
        grid.addWidget(self.name, 0, 2)
        grid.addWidget(self.set_type, 1, 0)
        grid.addWidget(self.pair, 1, 1)
        grid.addWidget(self.clear, 1, 2, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(self.new_zone, 2, 0)
        grid.addWidget(self.refresh, 2, 2, Qt.AlignmentFlag.AlignRight)
        grid.setColumnStretch(2, 1)  # the name field takes the spare width, not the buttons
        layout.addLayout(grid)
        self.set_type.clicked.connect(
            lambda: self.type_requested.emit(
                self.patch_type.currentText(), self.name.text().strip()
            )
        )
        self.pair.clicked.connect(self.pair_requested.emit)
        self.clear.clicked.connect(self.clear_requested.emit)
        self.new_zone.clicked.connect(lambda: self.zone_requested.emit(self.name.text().strip()))
        self.refresh.clicked.connect(self.refresh_requested.emit)
        self.tree = QTreeWidget()
        self.tree.setColumnCount(len(COLUMNS))
        self.tree.setHeaderLabels(list(COLUMNS))
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.tree.setWordWrap(True)
        self.tree.setToolTip("Click a row to select that boundary, pair or patch in the viewport")
        self.tree.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(self.tree, 1)
        self.detail = QLabel("")  # the clicked row's full status (long reasons get cut in a column)
        self.detail.setWordWrap(True)
        self.detail.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.detail)
        self.report = None
        self.stale = False

    # ------------------------------------------------------------------ report
    def show_report(self, report, store) -> None:
        """Fill the rows from a BoundaryReport and the store (the patches' pairings)."""
        self.report = report
        self.stale = False
        self.tree.clear()
        by_where: dict[str, list[str]] = {}
        for issue in report.issues:
            by_where.setdefault(issue.where, []).append(issue.message)
        patches = {p.name: p for p in store.entities("patch")}
        boundaries = QTreeWidgetItem(self.tree, [f"Boundaries ({len(report.boundary)})"])
        for label in report.boundary:
            patch = report.assignment.get(label)
            kind = ""
            patch_text = "none"
            status = "ok"
            color = OK_COLOR
            if patch is None:
                status = "no type yet: select it and press Set type"
                color = UNASSIGNED_COLOR
            else:
                patch_text = patch
                kind = patches[patch].params["patch_type"]
                if by_where.get(label):
                    status = "; ".join(by_where.pop(label))
                    color = UNASSIGNED_COLOR
            item = QTreeWidgetItem(boundaries, [label, patch_text, kind, status])
            item.setData(0, ENTITY_ROLE, report.ids.get(label))
            item.setForeground(3, QBrush(QColor(color)))
        boundaries.setExpanded(True)
        pairs = []
        unpaired = []
        for patch in store.entities("patch"):
            if patch.params["patch_type"] != "periodic":
                continue
            pairing = patch.params.get("periodic")
            if pairing is not None:
                pairs.append((patch, store.get(pairing["partner"])))
            elif partner_of(patch, store) is None:
                unpaired.append(patch)
        matched = {(s.name, t.name): m for s, t, m in report.periodic}
        broken = 0
        if pairs or unpaired:
            group = QTreeWidgetItem(self.tree, [f"Periodic pairs ({len(pairs)})"])
            for source, target in pairs:
                problems = by_where.pop(source.name, []) + by_where.pop(target.name, [])
                color = OK_COLOR
                if (source.name, target.name) in matched and not problems:
                    status = match_text(matched[(source.name, target.name)])
                else:
                    status = "; ".join(problems)
                    if not status:
                        status = "not checked"
                    color = UNASSIGNED_COLOR
                    broken += 1
                item = QTreeWidgetItem(
                    group,
                    [
                        f"{source.name} \u2194 {target.name}",
                        transform_text(source.params["periodic"]),
                        "periodic",
                        status,
                    ],
                )
                item.setData(0, ENTITY_ROLE, source.id)
                item.setForeground(3, QBrush(QColor(color)))
            for patch in unpaired:
                by_where.pop(patch.name, None)
                item = QTreeWidgetItem(
                    group,
                    [
                        patch.name,
                        "",
                        "periodic",
                        "no partner: select both sides and press Pair as periodic",
                    ],
                )
                item.setData(0, ENTITY_ROLE, patch.id)
                item.setForeground(3, QBrush(QColor(UNASSIGNED_COLOR)))
                broken += 1
            group.setExpanded(True)
        if report.interfaces:
            group = QTreeWidgetItem(self.tree, [f"Interfaces ({len(report.interfaces)})"])
            for label in report.interfaces:
                item = QTreeWidgetItem(
                    group, [label, "", "interface", "shared by two domains: no type needed"]
                )
                item.setForeground(3, QBrush(QColor(INTERFACE_COLOR)))
            group.setExpanded(True)
        empty = [p for p in store.entities("patch") if not p.params["members"]]
        if empty:
            group = QTreeWidgetItem(self.tree, [f"Empty patches ({len(empty)})"])
            for patch in empty:
                item = QTreeWidgetItem(
                    group,
                    [
                        patch.name,
                        "",
                        patch.params["patch_type"],
                        "holds nothing: delete it or set it on a selection",
                    ],
                )
                item.setData(0, ENTITY_ROLE, patch.id)
                item.setForeground(3, QBrush(QColor(INTERFACE_COLOR)))
            group.setExpanded(True)
        others = [(where, text) for where, texts in by_where.items() for text in texts]
        if others:
            group = QTreeWidgetItem(self.tree, [f"Other problems ({len(others)})"])
            for where, text in others:
                item = QTreeWidgetItem(group, [where, "", "", text])
                item.setForeground(3, QBrush(QColor(UNASSIGNED_COLOR)))
            group.setExpanded(True)
        for column in range(len(COLUMNS) - 1):
            self.tree.resizeColumnToContents(column)
        for k in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(k)
            for j in range(top.childCount()):
                child = top.child(j)
                for column in range(len(COLUMNS)):
                    child.setToolTip(column, child.text(3))
        self.detail.setText("")
        self._set_status(report, broken, len(others))

    def _set_status(self, report, broken: int, others: int) -> None:
        total = len(report.boundary)
        unit = "boundaries"
        if total == 1:
            unit = "boundary"
        if report.ok:
            text = f"All {total} {unit} have a type"
            if report.periodic:
                text += f", {len(report.periodic)} periodic pair(s) match"
            if total == 0:
                text = "No boundaries yet (make a domain first)"
            else:
                text += ": ready to export"
            color = OK_COLOR
        else:
            parts = []
            if report.unassigned:
                parts.append(f"{len(report.unassigned)} of {total} {unit} have no type")
            if broken:
                parts.append(f"{broken} periodic pair(s) do not match")
            if others:
                parts.append(f"{others} other problem(s)")
            if not parts:
                parts.append(f"{len(report.issues)} problem(s)")
            text = "; ".join(parts)
            color = UNASSIGNED_COLOR
        self.status.setText(text)
        self.status.setStyleSheet(f"color: {color}; font-weight: bold")

    def mark_stale(self) -> None:
        """A 3D model changed: the report is not recomputed by itself."""
        self.stale = True
        self.status.setText("The model changed: press Refresh (it meshes the blocks and volumes)")
        self.status.setStyleSheet("")

    def rows(self, group: str) -> list[list[str]]:
        """The rows of a group ("Boundaries", "Periodic pairs", ...), for tests."""
        for k in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(k)
            if top.text(0).split(" (")[0] == group:
                return [
                    [top.child(j).text(c) for c in range(len(COLUMNS))]
                    for j in range(top.childCount())
                ]
        return []

    def _on_item_clicked(self, item, _column) -> None:
        self.detail.setText(item.text(3))
        entity_id = item.data(0, ENTITY_ROLE)
        if entity_id is not None:
            self.entity_chosen.emit(int(entity_id))
