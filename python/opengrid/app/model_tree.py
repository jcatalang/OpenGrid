"""Model tree: entities grouped as database, connectors, domains, blocks, patches, zones."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

# Top-level categories, in display order (hierarchy, see docs/SPEC.md).
from opengrid.validation import ValidationError

CATEGORIES = ("Database", "Connectors", "Domains", "Blocks", "Patches", "Zones")
ENTITY_ID_ROLE = Qt.ItemDataRole.UserRole


def connector_info(session, entity) -> str:
    """'<n> pts', with a warning sign if the connector has warnings or exceeds the threshold."""
    text = f"{entity.params['dimension']} pts"
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

    def category_names(self) -> list[str]:
        """Names of the top-level category nodes, in display order."""
        names = []
        for index in range(self.topLevelItemCount()):
            names.append(self.topLevelItem(index).text(0))
        return names

    def category_item(self, name: str) -> QTreeWidgetItem:
        return self._category_items[name]

    def populate(self, session) -> None:
        """Rebuild the Database and Connectors categories (entity ids kept as item data)."""
        self.blockSignals(True)
        self._entity_items.clear()
        database = self._category_items["Database"]
        connectors = self._category_items["Connectors"]
        domains = self._category_items["Domains"]
        blocks = self._category_items["Blocks"]
        patches = self._category_items["Patches"]
        zones = self._category_items["Zones"]
        for category in (database, connectors, domains, blocks, patches, zones):
            category.takeChildren()
        for entity in session.store:
            if entity.kind == "connector":
                item = QTreeWidgetItem(connectors, [entity.name, connector_info(session, entity)])
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
            elif entity.kind == "block":
                info = "invalid"
                try:
                    nk, nj, ni, _ = session.block_grid(entity.name).shape
                    info = f"{ni} x {nj} x {nk}"
                except ValidationError:
                    pass
                item = QTreeWidgetItem(blocks, [entity.name, info])
            elif entity.kind == "patch":
                count = len(entity.params["members"])
                item = QTreeWidgetItem(
                    patches, [entity.name, f"{entity.params['patch_type']}, {count}"]
                )
            elif entity.kind == "zone":
                count = len(entity.params["members"])
                item = QTreeWidgetItem(zones, [entity.name, f"{count} domain(s)"])
            else:
                item = QTreeWidgetItem(database, [entity.name, entity.kind])
            item.setData(0, ENTITY_ID_ROLE, entity.id)
            self._entity_items[entity.id] = item
        for category in (database, connectors, domains, blocks, patches, zones):
            category.setText(1, f"{category.childCount()}")
        self.resizeColumnToContents(0)
        self.blockSignals(False)

    def entity_names(self, category: str = "Database") -> list[str]:
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
