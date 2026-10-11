"""Model tree: entities grouped as database, connectors, domains, blocks, patches, zones."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

# Top-level categories, in display order (hierarchy, see docs/SPEC.md).
CATEGORIES = ("Database", "Connectors", "Domains", "Blocks", "Patches", "Zones")
ENTITY_ID_ROLE = Qt.ItemDataRole.UserRole


def connector_info(session, entity) -> str:
    """'<n> pts', with a warning sign if the connector has warnings or exceeds the threshold."""
    text = f"{entity.params['dimension']} pts"
    nodes = session.connector_nodes(entity.name)
    if nodes.warnings or nodes.max_stretching_ratio > session.stretching_threshold:
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
        database.takeChildren()
        connectors.takeChildren()
        for entity in session.store:
            if entity.kind == "connector":
                item = QTreeWidgetItem(connectors, [entity.name, connector_info(session, entity)])
            else:
                item = QTreeWidgetItem(database, [entity.name, entity.kind])
            item.setData(0, ENTITY_ID_ROLE, entity.id)
            self._entity_items[entity.id] = item
        database.setText(1, f"{database.childCount()}")
        connectors.setText(1, f"{connectors.childCount()}")
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
