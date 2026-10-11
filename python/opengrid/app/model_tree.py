"""Model tree: entities grouped as database, connectors, domains, blocks, patches, zones."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

# Top-level categories, in display order (hierarchy, see docs/SPEC.md).
CATEGORIES = ("Database", "Connectors", "Domains", "Blocks", "Patches", "Zones")
ENTITY_ID_ROLE = Qt.ItemDataRole.UserRole


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

    def populate(self, store) -> None:
        """Rebuild the Database category from the store (ids kept as item data)."""
        self.blockSignals(True)
        database = self._category_items["Database"]
        database.takeChildren()
        self._entity_items.clear()
        for entity in store:
            item = QTreeWidgetItem(database, [entity.name, entity.kind])
            item.setData(0, ENTITY_ID_ROLE, entity.id)
            self._entity_items[entity.id] = item
        database.setText(1, f"{len(store)}")
        self.resizeColumnToContents(0)
        self.blockSignals(False)

    def entity_names(self) -> list[str]:
        database = self._category_items["Database"]
        return [database.child(i).text(0) for i in range(database.childCount())]

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
