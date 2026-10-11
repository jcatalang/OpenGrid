"""Property panel: view and edit the parameters of the selected entity.

Edits are applied through ``Session.modify`` / ``Session.rename``, i.e. as undoable commands.
Input that cannot be parsed or that the model rejects is reported in the panel; nothing is
changed in that case.
"""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from opengrid.api import Session
from opengrid.model.connectors import DISTRIBUTION_TYPES, SPACING_MODES
from opengrid.model.entities import SPLINE_METHODS
from opengrid.validation import ValidationError

# Field editors per kind: (parameter, editor type). Types: vec3, opt_vec3, number, integer,
# points, numbers, method.
FIELDS = {
    "point": (("xyz", "vec3"),),
    "line": (("p0", "vec3"), ("p1", "vec3")),
    "arc": (
        ("center", "vec3"),
        ("radius", "number"),
        ("start_angle", "number"),
        ("end_angle", "number"),
        ("normal", "vec3"),
        ("x_axis", "vec3"),
    ),
    "polyline": (("points", "points"),),
    "spline": (
        ("method", "method"),
        ("points", "points"),
        ("start_tangent", "opt_vec3"),
        ("end_tangent", "opt_vec3"),
    ),
    "nurbs": (
        ("degree", "integer"),
        ("control_points", "points"),
        ("weights", "numbers"),
        ("knots", "numbers"),
    ),
}


CONNECTOR_FIELDS = (
    ("dimension", "dimension"),
    ("mode", "mode"),
    ("distribution", "distribution"),
    ("breakpoints", "breakpoints"),
    ("propagate_spacing", "flag"),
)


class ParseError(ValueError):
    pass


def format_optional(value) -> str:
    if value is None:
        return ""
    return format_number(value)


def parse_optional(text: str, field: str) -> float | None:
    if not text.strip() or text.strip() == "-":
        return None
    values = parse_numbers(text, field)
    if len(values) != 1:
        raise ParseError(f"{field}: expected one number")
    return values[0]


def format_breakpoints(breakpoints) -> str:
    """One line per breakpoint: position index type begin end ratio ('-' = none)."""
    lines = []
    for item in breakpoints:
        dist = item["distribution"]
        values = []
        for key in ("begin_spacing", "end_spacing", "growth_ratio"):
            if dist[key] is None:
                values.append("-")
            else:
                values.append(format_number(dist[key]))
        lines.append(
            f"{format_number(item['position'])} {item['index']} {dist['type']} {' '.join(values)}"
        )
    return "\n".join(lines)


def parse_breakpoints(text: str) -> tuple:
    result = []
    for number, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        where = f"breakpoints line {number}"
        if len(fields) not in (3, 6):
            raise ParseError(f"{where}: expected 'position index type [begin end ratio]'")
        position = parse_numbers(fields[0], where)[0]
        try:
            index = int(fields[1])
        except ValueError:
            raise ParseError(f"{where}: index {fields[1]!r} is not an integer") from None
        dist = {"type": fields[2], "begin_spacing": None, "end_spacing": None, "growth_ratio": None}
        if len(fields) == 6:
            dist["begin_spacing"] = parse_optional(fields[3], where)
            dist["end_spacing"] = parse_optional(fields[4], where)
            dist["growth_ratio"] = parse_optional(fields[5], where)
        result.append({"position": position, "index": index, "distribution": dist})
    return tuple(result)


class DistributionEditor(QWidget):
    """Type + begin/end spacing + growth ratio (empty fields mean 'not specified')."""

    def __init__(self, value: dict, parent=None):
        super().__init__(parent)
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.type = QComboBox()
        self.type.addItems(list(DISTRIBUTION_TYPES))
        self.type.setCurrentText(value["type"])
        self.begin = QLineEdit(format_optional(value["begin_spacing"]))
        self.end = QLineEdit(format_optional(value["end_spacing"]))
        self.ratio = QLineEdit(format_optional(value["growth_ratio"]))
        for edit in (self.begin, self.end, self.ratio):
            edit.setPlaceholderText("none")
        layout.addRow("type", self.type)
        layout.addRow("begin spacing", self.begin)
        layout.addRow("end spacing", self.end)
        layout.addRow("growth ratio", self.ratio)

    def value(self) -> dict:
        return {
            "type": self.type.currentText(),
            "begin_spacing": parse_optional(self.begin.text(), "begin spacing"),
            "end_spacing": parse_optional(self.end.text(), "end spacing"),
            "growth_ratio": parse_optional(self.ratio.text(), "growth ratio"),
        }


def format_number(value: float) -> str:
    return repr(float(value))


def format_vec3(value) -> str:
    if value is None:
        return ""
    return ", ".join(format_number(c) for c in value)


def parse_numbers(text: str, field: str) -> list[float]:
    values = []
    for part in text.replace(",", " ").split():
        try:
            values.append(float(part))
        except ValueError:
            raise ParseError(f"{field}: {part!r} is not a number") from None
    return values


def parse_vec3(text: str, field: str) -> tuple:
    values = parse_numbers(text, field)
    if len(values) not in (2, 3):
        raise ParseError(f"{field}: expected 2 or 3 numbers, got {len(values)}")
    return tuple(values)


def parse_points(text: str, field: str) -> list[tuple]:
    points = []
    for number, line in enumerate(text.splitlines(), start=1):
        if line.strip():
            points.append(parse_vec3(line, f"{field} line {number}"))
    return points


class PropertyPanel(QWidget):
    """Editor for one entity. Emits ``message(text, is_error)`` after apply attempts."""

    message = Signal(str, bool)

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self.entity_id: int | None = None
        self._editors: dict[str, tuple[str, QWidget]] = {}

        self._layout = QVBoxLayout(self)
        self._title = QLabel("No selection")
        self._layout.addWidget(self._title)
        self._form_host = QWidget()
        self._form = QFormLayout(self._form_host)
        self._layout.addWidget(self._form_host)
        self._error = QLabel("")
        self._error.setWordWrap(True)
        self._error.setStyleSheet("color: #ff6b6b;")
        self._layout.addWidget(self._error)
        buttons = QHBoxLayout()
        self.apply_button = QPushButton("Apply")
        self.revert_button = QPushButton("Revert")
        self.apply_button.clicked.connect(self.apply)
        self.revert_button.clicked.connect(self.refresh)
        buttons.addWidget(self.apply_button)
        buttons.addWidget(self.revert_button)
        self._layout.addLayout(buttons)
        self._layout.addStretch(1)
        self.show_entity(None)

    def show_entity(self, entity_id: int | None) -> None:
        self.entity_id = entity_id
        self.refresh()

    def refresh(self) -> None:
        while self._form.rowCount():
            self._form.removeRow(0)
        self._editors.clear()
        self._error.setText("")
        if self.entity_id is None or self.entity_id not in self.session.store:
            self.entity_id = None
            self._title.setText("No selection")
            self.apply_button.setEnabled(False)
            self.revert_button.setEnabled(False)
            return
        entity = self.session.store.get(self.entity_id)
        self._title.setText(f"<b>{entity.kind}</b>  (id {entity.id})")
        name_edit = QLineEdit(entity.name)
        self._form.addRow("name", name_edit)
        self._editors["name"] = ("name", name_edit)
        if entity.kind == "connector":
            self._connector_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        for field, kind in FIELDS[entity.kind]:
            editor = self._make_editor(kind, entity.params[field])
            self._form.addRow(field, editor)
            self._editors[field] = (kind, editor)
        if entity.kind != "point":
            length = self.session.length(entity.name)
            self._form.addRow("length", QLabel(f"{length:.10g}"))
        self.apply_button.setEnabled(True)
        self.revert_button.setEnabled(True)

    def _connector_rows(self, entity) -> None:
        params = entity.params
        pieces = []
        for segment in params["segments"]:
            curve = self.session.store.get(segment["curve"]).name
            pieces.append(f"{curve} [{segment['t0']:g} → {segment['t1']:g}]")
        segments = QLabel("\n".join(pieces))
        segments.setWordWrap(True)
        self._form.addRow("segments", segments)
        for field, kind in CONNECTOR_FIELDS:
            editor = self._make_editor(kind, params[field])
            self._form.addRow(field.replace("_", " "), editor)
            self._editors[field] = (kind, editor)
        nodes = self.session.connector_nodes(entity.name)
        self._form.addRow("length", QLabel(f"{nodes.length:.10g}"))
        self._form.addRow(
            "begin / end spacing", QLabel(f"{nodes.begin_spacing:.6g} / {nodes.end_spacing:.6g}")
        )
        ratio = QLabel(f"{nodes.max_stretching_ratio:.4g}")
        if nodes.max_stretching_ratio > self.session.stretching_threshold:
            ratio.setText(
                f"{nodes.max_stretching_ratio:.4g}  (> {self.session.stretching_threshold:g})"
            )
            ratio.setStyleSheet("color: #ff6b6b;")
        self._form.addRow("max stretching", ratio)
        for end, text in nodes.propagated.items():
            self._form.addRow(f"{end} spacing from", QLabel(text))
        if nodes.warnings:
            warnings = QLabel("\n".join(nodes.warnings))
            warnings.setWordWrap(True)
            warnings.setStyleSheet("color: #ff6b6b;")
            self._form.addRow("warnings", warnings)

    def _make_editor(self, kind: str, value) -> QWidget:
        if kind == "dimension":
            editor = QSpinBox()
            editor.setRange(2, 10_000_000)
            editor.setValue(int(value))
            return editor
        if kind == "mode":
            editor = QComboBox()
            editor.addItems(list(SPACING_MODES))
            editor.setCurrentText(value)
            return editor
        if kind == "distribution":
            return DistributionEditor(value)
        if kind == "flag":
            editor = QCheckBox()
            editor.setChecked(bool(value))
            return editor
        if kind == "breakpoints":
            editor = QPlainTextEdit(format_breakpoints(value))
            editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
            editor.setPlaceholderText("position index type [begin end ratio], '-' = none")
            editor.setMinimumHeight(60)
            return editor
        if kind == "vec3" or kind == "opt_vec3":
            editor = QLineEdit(format_vec3(value))
            if kind == "opt_vec3":
                editor.setPlaceholderText("none")
            return editor
        if kind == "number":
            return QLineEdit(format_number(value))
        if kind == "integer":
            editor = QSpinBox()
            editor.setRange(1, 15)
            editor.setValue(int(value))
            return editor
        if kind == "method":
            editor = QComboBox()
            editor.addItems(list(SPLINE_METHODS))
            editor.setCurrentText(value)
            return editor
        editor = QPlainTextEdit()
        editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        if kind == "points":
            editor.setPlainText("\n".join(format_vec3(p) for p in value))
        else:
            editor.setPlainText("\n".join(format_number(v) for v in value))
        editor.setMinimumHeight(90)
        return editor

    def read_changes(self) -> tuple[str | None, dict]:
        """Parse the editors; returns (new name or None, changed parameters)."""
        entity = self.session.store.get(self.entity_id)
        new_name = None
        changes = {}
        for field, (kind, editor) in self._editors.items():
            if field == "name":
                text = editor.text().strip()
                if text != entity.name:
                    new_name = text
                continue
            value = self._read_editor(field, kind, editor)
            if value != entity.params[field]:
                changes[field] = value
        return new_name, changes

    def _read_editor(self, field: str, kind: str, editor: QWidget):
        if kind == "dimension":
            return editor.value()
        if kind == "mode":
            return editor.currentText()
        if kind == "distribution":
            return editor.value()
        if kind == "flag":
            return editor.isChecked()
        if kind == "breakpoints":
            return parse_breakpoints(editor.toPlainText())
        if kind == "vec3":
            return _canonical_vec(parse_vec3(editor.text(), field))
        if kind == "opt_vec3":
            if not editor.text().strip():
                return None
            return _canonical_vec(parse_vec3(editor.text(), field))
        if kind == "number":
            values = parse_numbers(editor.text(), field)
            if len(values) != 1:
                raise ParseError(f"{field}: expected one number")
            return values[0]
        if kind == "integer":
            return editor.value()
        if kind == "method":
            return editor.currentText()
        if kind == "points":
            return tuple(_canonical_vec(p) for p in parse_points(editor.toPlainText(), field))
        return tuple(parse_numbers(editor.toPlainText(), field))

    def apply(self) -> bool:
        """Apply edits as commands; returns True if the model changed."""
        if self.entity_id is None:
            return False
        try:
            new_name, changes = self.read_changes()
        except ParseError as error:
            self._show_error(str(error))
            return False
        name = self.session.store.get(self.entity_id).name
        try:
            if changes:
                self.session.modify(name, **changes)
            if new_name is not None:
                self.session.rename(name, new_name)
        except ValidationError as error:
            self._show_error(str(error))
            return False
        if not changes and new_name is None:
            self.message.emit("No changes", False)
            return False
        self.message.emit(f"Updated {self.session.store.get(self.entity_id).name}", False)
        return True

    def _show_error(self, text: str) -> None:
        self._error.setText(text)
        self.message.emit(text, True)


def _canonical_vec(values) -> tuple:
    coords = [float(v) for v in values]
    if len(coords) == 2:
        coords.append(0.0)
    return tuple(coords)
