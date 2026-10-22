"""Model tree: entities grouped as curves, surfaces, points, domains, blocks and zones.

A curve and the connector of its points are one row (``model/curves.py``). Boundary types are
shown where they are given: on the curves (2D) and on the faces of blocks and volumes, listed
under them (3D, ``model/faces.py``); there is no category of patches."""

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from opengrid.model.curves import CurveIndex
from opengrid.model.entities import SURFACE_KINDS
from opengrid.validation import ValidationError

# Top-level categories, in display order. Points and Surfaces show only when there are some.
CATEGORIES = ("Curves", "Surfaces", "Points", "Domains", "Blocks", "Zones")
HIDDEN_WHEN_EMPTY = ("Surfaces", "Points")
ENTITY_ID_ROLE = Qt.ItemDataRole.UserRole
FACE_ROLE = Qt.ItemDataRole.UserRole + 1  # a face of a block or volume (a patch member)


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


def boundary_types(session) -> dict:
    """{member key: text} of every boundary given a type: 'wall', 'periodic ↔ c-right'."""
    from opengrid.model.patches import partner_of

    texts = {}
    store = session.store
    for patch in store.entities("patch"):
        kind = patch.params["patch_type"]
        if kind == "periodic":
            partner = partner_of(patch, store)
            if partner is None:
                kind = "periodic, NOT PAIRED"
            else:
                labels = [member_label(m, store) for m in store.get(partner).params["members"]]
                kind = f"periodic \u2194 {', '.join(labels)}"
        for member in patch.params["members"]:
            texts[member_key(member)] = kind
    return texts


def member_key(member: dict) -> tuple:
    return tuple(sorted(member.items()))


def member_label(member: dict, store) -> str:
    from opengrid.model.faces import face_identifier

    return face_identifier(store, member)


class ModelTree(QTreeWidget):
    """Tree of model entities grouped by category; emits the selected entity ids and faces."""

    selection_changed = Signal(list, list)  # entity ids, faces (patch members)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Single keys the window binds (tools, views): the tree's type-ahead search would take
        # them first; it leaves these to the window's shortcuts.
        self.passthrough_keys: set[int] = set()
        self.setHeaderLabels(["Entity", "Type"])
        self.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
        self._category_items = {}
        self._entity_items: dict[int, QTreeWidgetItem] = {}
        self._face_items: dict[tuple, QTreeWidgetItem] = {}
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
        from opengrid.model.faces import faces_of

        self.blockSignals(True)
        expanded = {
            self._entity_items[i].text(0)
            for i in self._entity_items
            if self._entity_items[i].isExpanded()
        }
        self._entity_items.clear()
        self._face_items.clear()
        curves = self._category_items["Curves"]
        surfaces = self._category_items["Surfaces"]
        points = self._category_items["Points"]
        domains = self._category_items["Domains"]
        blocks = self._category_items["Blocks"]
        zones = self._category_items["Zones"]
        categories = (curves, surfaces, points, domains, blocks, zones)
        for category in categories:
            category.takeChildren()
        index = CurveIndex(session.store)
        types = boundary_types(session)
        for entity in session.store:
            if entity.kind == "connector":
                info = connector_info(session, entity, index)
                kind = types.get(member_key({"connector": entity.id}))
                if kind is not None:
                    info += f"  \u00b7 {kind}"
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
                continue  # boundary types are shown on the curves and faces they are given to
            elif entity.kind == "zone":
                count = len(entity.params["members"])
                item = QTreeWidgetItem(zones, [entity.name, f"{count} domain(s)"])
            else:  # a database curve: listed when no connector stands for it
                if not index.is_listed(entity):
                    continue
                item = QTreeWidgetItem(curves, [entity.name, index.describe(entity)])
            item.setData(0, ENTITY_ID_ROLE, entity.id)
            self._entity_items[entity.id] = item
            for face in faces_of(session.store, entity):  # a block's or volume's faces
                kind = types.get(member_key(face.member), "no type")
                child = QTreeWidgetItem(item, [face.short, f"face \u00b7 {kind}"])
                child.setData(0, FACE_ROLE, face.member)
                child.setToolTip(0, face.identifier)
                self._face_items[member_key(face.member)] = child
            if entity.name in expanded:
                item.setExpanded(True)
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

    def select_ids(self, ids, faces=()) -> None:
        self.blockSignals(True)
        self.clearSelection()
        for entity_id in ids:
            item = self._entity_items.get(entity_id)
            if item is not None:
                item.setSelected(True)
                self.scrollToItem(item)
        for face in faces:
            item = self._face_items.get(member_key(face))
            if item is not None:
                item.parent().setExpanded(True)
                item.setSelected(True)
                self.scrollToItem(item)
        self.blockSignals(False)

    def face_rows(self, entity_id: int) -> list[list[str]]:
        """The face rows under a block or volume, for tests."""
        item = self._entity_items[entity_id]
        return [[item.child(k).text(0), item.child(k).text(1)] for k in range(item.childCount())]

    def _emit_selection(self) -> None:
        ids = []
        faces = []
        for item in self.selectedItems():
            entity_id = item.data(0, ENTITY_ID_ROLE)
            face = item.data(0, FACE_ROLE)
            if entity_id is not None:
                ids.append(int(entity_id))
            elif face is not None:
                faces.append(dict(face))
        self.selection_changed.emit(sorted(ids), faces)


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
