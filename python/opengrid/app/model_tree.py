"""Model tree: entities grouped as database, connectors, domains, blocks, patches, zones."""

from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

# Top-level categories, in display order (hierarchy, see docs/SPEC.md).
CATEGORIES = ("Database", "Connectors", "Domains", "Blocks", "Patches", "Zones")


class ModelTree(QTreeWidget):
    """Tree of model entities grouped by category. Populated from the model in Phase 1+."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setHeaderLabels(["Entity", "Info"])
        self.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
        self._category_items = {}
        for name in CATEGORIES:
            item = QTreeWidgetItem(self, [name, ""])
            item.setExpanded(True)
            self._category_items[name] = item

    def category_names(self) -> list[str]:
        """Names of the top-level category nodes, in display order."""
        names = []
        for index in range(self.topLevelItemCount()):
            names.append(self.topLevelItem(index).text(0))
        return names

    def category_item(self, name: str) -> QTreeWidgetItem:
        return self._category_items[name]
