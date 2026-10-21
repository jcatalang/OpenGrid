"""Model tree: entities grouped as curves, surfaces, points, domains, blocks, patches, zones.

A curve and the connector of its points are one row (``model/curves.py``)."""

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from opengrid.model.curves import CurveIndex
from opengrid.model.entities import SURFACE_KINDS
from opengrid.validation import ValidationError

# Top-level categories, in display order. Points and Surfaces show only when there are some.
CATEGORIES = ("Curves", "Surfaces", "Points", "Domains", "Blocks", "Patches", "Zones")
HIDDEN_WHEN_EMPTY = ("Surfaces", "Points")
ENTITY_ID_ROLE = Qt.ItemDataRole.UserRole


def connector_info(session, entity, index: CurveIndex | None = None) -> str:
    """'line, 21 pts' (the curve's kind when it is one with its points), with a warning sign
    if the points have warnings or exceed the stretching threshold."""
    if index is None:
        index = CurveIndex(session.store)
    text = index.describe(entity)
    nodes = session.connector_nodes(entity.name)
    if nodes.warnings or nodes.max_stretching_ratio > session.stretching_threshold:
        text += "  ⚠"
    return text


def domain_info(session, entity) -> str:
    """'ni x nj', with a warning sign while smoothing is pending, not converged, or cells are
    inverted."""
    ni, nj = session.domain_dimensions(entity.name)
    text = f"{ni} x {nj}"
    status = session.domain_status(entity.name)
    if "pending" in status or "NOT" in status:
        text += "  ⚠"
    return text


class ModelTree(QTreeWidget):
    """Tree of model entities grouped by category; emits the selected entity ids."""

    selection_changed = Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Single keys the window binds (tools, views): the tree's type-ahead search would take
        # them first; it leaves these to the window's shortcuts.
        self.passthrough_keys: set[int] = set()
        self.setHeaderLabels(["Entity", "Type"])
        self.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
        self._category_items = {}
        self._entity_items: dict[int, QTreeWidgetItem] = {}
        for name in CATEGORIES:
            item = QTreeWidgetItem(self, [name, ""])
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            item.setExpanded(True)
            self._category_items[name] = item
        self.itemSelectionChanged.connect(self._emit_selection)

    def event(self, event) -> bool:
        if (
            event.type() == QEvent.Type.ShortcutOverride
            and int(event.key()) in self.passthrough_keys
            and not (event.modifiers() & ~Qt.KeyboardModifier.KeypadModifier)
        ):
            event.ignore()
            return False
        return super().event(event)

    def category_names(self) -> list[str]:
        """Names of the top-level category nodes, in display order."""
        names = []
        for index in range(self.topLevelItemCount()):
            names.append(self.topLevelItem(index).text(0))
        return names

    def category_item(self, name: str) -> QTreeWidgetItem:
        return self._category_items[name]

    def populate(self, session) -> None:
        """Rebuild the categories (entity ids kept as item data; a curve's row stands for
        its database curve and its connector)."""
        self.blockSignals(True)
        self._entity_items.clear()
        curves = self._category_items["Curves"]
        surfaces = self._category_items["Surfaces"]
        points = self._category_items["Points"]
        domains = self._category_items["Domains"]
        blocks = self._category_items["Blocks"]
        patches = self._category_items["Patches"]
        zones = self._category_items["Zones"]
        categories = (curves, surfaces, points, domains, blocks, patches, zones)
        for category in categories:
            category.takeChildren()
        index = CurveIndex(session.store)
        for entity in session.store:
            if entity.kind == "connector":
                info = connector_info(session, entity, index)
                item = QTreeWidgetItem(curves, [entity.name, info])
            elif entity.kind in SURFACE_KINDS:
                item = QTreeWidgetItem(surfaces, [entity.name, entity.kind.replace("_", " ")])
            elif entity.kind == "point":
                item = QTreeWidgetItem(points, [entity.name, "point"])
            elif entity.kind == "domain":
                item = QTreeWidgetItem(domains, [entity.name, domain_info(session, entity)])
            elif entity.kind == "unstructured_domain":
                mesh = session.unstructured_mesh(entity.name)
                parts = []
                if len(mesh["triangles"]) or not len(mesh["quads"]):
                    parts.append(f"{len(mesh['triangles'])} tris")
                if len(mesh["quads"]):
                    parts.append(f"{len(mesh['quads'])} quads")
                item = QTreeWidgetItem(domains, [entity.name, ", ".join(parts)])
            elif entity.kind == "smoothing_group":
                info = f"group of {len(entity.params['domains']) + len(entity.params['blocks'])}"
                try:
                    status = session.group_status(entity.name)
                    if "pending" in status or "NOT" in status:
                        info += "  ⚠"
                except ValidationError:
                    info += "  ⚠"
                item = QTreeWidgetItem(domains, [entity.name, info])
            elif entity.kind == "block":
                info = "invalid"
                try:
                    nk, nj, ni, _ = session.block_grid(entity.name).shape
                    info = f"{ni} x {nj} x {nk}"
                except ValidationError:
                    pass
                item = QTreeWidgetItem(blocks, [entity.name, info])
            elif entity.kind == "volume":
                info = "invalid"
                try:
                    mesh = session.volume_mesh(entity.name)
                    info = f"{len(mesh['tets'])} tets"
                    for key, short in (("pyramids", "pyr"), ("prisms", "pri"), ("hexes", "hex")):
                        if len(mesh[key]):
                            info += f", {len(mesh[key])} {short}"
                except ValidationError:
                    pass
                item = QTreeWidgetItem(blocks, [entity.name, info])
            elif entity.kind == "sketch":
                info = f"sketch, {len(entity.params['blocks'])} blocks"
                item = QTreeWidgetItem(blocks, [entity.name, info])
            elif entity.kind == "patch":
                item = QTreeWidgetItem(patches, [entity.name, patch_summary(entity, session.store)])
            elif entity.kind == "zone":
                count = len(entity.params["members"])
                item = QTreeWidgetItem(zones, [entity.name, f"{count} domain(s)"])
            else:  # a database curve: listed when no connector stands for it
                if not index.is_listed(entity):
                    continue
                item = QTreeWidgetItem(curves, [entity.name, index.describe(entity)])
            item.setData(0, ENTITY_ID_ROLE, entity.id)
            self._entity_items[entity.id] = item
        for curve, connector in index.connector_of.items():
            self._entity_items[curve] = self._entity_items[connector]  # selecting either
        for category in categories:
            category.setText(1, f"{category.childCount()}")
            name = category.text(0)
            category.setHidden(name in HIDDEN_WHEN_EMPTY and category.childCount() == 0)
        self.resizeColumnToContents(0)
        self.blockSignals(False)

    def entity_names(self, category: str = "Curves") -> list[str]:
        parent = self._category_items[category]
        return [parent.child(i).text(0) for i in range(parent.childCount())]

    def select_ids(self, ids) -> None:
        self.blockSignals(True)
        self.clearSelection()
        for entity_id in ids:
            item = self._entity_items.get(entity_id)
            if item is not None:
                item.setSelected(True)
                self.scrollToItem(item)
        self.blockSignals(False)

    def _emit_selection(self) -> None:
        ids = []
        for item in self.selectedItems():
            entity_id = item.data(0, ENTITY_ID_ROLE)
            if entity_id is not None:
                ids.append(int(entity_id))
        self.selection_changed.emit(sorted(ids))


def patch_summary(entity, store) -> str:
    """The tree's text for a patch: its type, the partner of a periodic one (or NOT PAIRED),
    and what it holds ("wall, 2 curves"; 3D: faces)."""
    from opengrid.model.patches import partner_of

    members = entity.params["members"]
    unit = "curve"
    if any("connector" not in m for m in members):
        unit = "face"
    if len(members) != 1:
        unit += "s"
    kind = entity.params["patch_type"]
    if kind == "periodic":
        partner = partner_of(entity, store)
        if partner is None:
            kind = "periodic, NOT PAIRED"
        else:
            kind = f"periodic \u2194 {partner}"
    return f"{kind}, {len(members)} {unit}"
