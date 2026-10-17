"""Create > Shapes: one small dialog per standard shape (rectangle, circle, ellipse, NACA
airfoil, box, cylinder, sphere, wing) whose fields map onto a Session creator."""

from dataclasses import dataclass

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QVBoxLayout,
)


class ShapeParseError(ValueError):
    pass


@dataclass(frozen=True)
class Field:
    key: str
    label: str
    default: str
    kind: str  # "number", "int", "vec2", "vec3", "text", "bool"


@dataclass(frozen=True)
class ShapeSpec:
    title: str
    method: str  # the Session method
    fields: tuple
    hint: str


NAME = Field("name", "name (optional)", "", "text")
SHAPES = {
    "rectangle": ShapeSpec(
        "Rectangle",
        "create_rectangle",
        (
            Field("center", "centre x y", "0 0", "vec2"),
            Field("width", "width", "1", "number"),
            Field("height", "height", "1", "number"),
            NAME,
        ),
        "Four lines: bottom, right, top, left (one connector each, or a sketch block on them)",
    ),
    "circle": ShapeSpec(
        "Circle",
        "create_circle",
        (
            Field("center", "centre x y", "0 0", "vec2"),
            Field("radius", "radius", "0.5", "number"),
            NAME,
        ),
        "A full circle starting on the +x side",
    ),
    "ellipse": ShapeSpec(
        "Ellipse",
        "create_ellipse",
        (
            Field("center", "centre x y", "0 0", "vec2"),
            Field("a", "semi-axis along x", "1", "number"),
            Field("b", "semi-axis along y", "0.5", "number"),
            NAME,
        ),
        "An exact ellipse (rational NURBS) starting on the +x side",
    ),
    "airfoil": ShapeSpec(
        "NACA airfoil",
        "create_airfoil",
        (
            Field("code", "NACA 4-digit code", "0012", "text"),
            Field("chord", "chord", "1", "number"),
            Field("leading_edge", "leading edge x y", "0 0", "vec2"),
            Field("angle", "angle of attack [deg]", "0", "number"),
            Field("points", "points per side", "101", "int"),
            NAME,
        ),
        "Two splines, upper and lower, from the leading edge to the sharp trailing edge; "
        "select both (and an outer boundary) for an unstructured domain",
    ),
    "box": ShapeSpec(
        "Box",
        "create_box",
        (
            Field("center", "centre x y z", "0 0 0", "vec3"),
            Field("size", "size x y z", "1 1 1", "vec3"),
            NAME,
        ),
        "Six plane patches, normals outward",
    ),
    "cylinder": ShapeSpec(
        "Cylinder",
        "create_cylinder",
        (
            Field("center", "bottom centre x y z", "0 0 0", "vec3"),
            Field("radius", "radius", "0.5", "number"),
            Field("height", "height", "1", "number"),
            Field("axis", "axis x y z", "0 0 1", "vec3"),
            Field("caps", "end caps", "yes", "bool"),
            NAME,
        ),
        "A surface of revolution (the caps are discs, degenerate at their centre)",
    ),
    "sphere": ShapeSpec(
        "Sphere",
        "create_sphere",
        (
            Field("center", "centre x y z", "0 0 0", "vec3"),
            Field("radius", "radius", "0.5", "number"),
            NAME,
        ),
        "A meridian half circle revolved about the z axis through the centre",
    ),
    "wing": ShapeSpec(
        "Wing",
        "create_wing",
        (
            Field("code", "NACA 4-digit code", "0012", "text"),
            Field("root_chord", "root chord", "1", "number"),
            Field("tip_chord", "tip chord", "0.5", "number"),
            Field("span", "span (along y)", "3", "number"),
            Field("sweep", "sweep [deg]", "0", "number"),
            Field("dihedral", "dihedral [deg]", "0", "number"),
            Field("points", "points per side", "101", "int"),
            NAME,
        ),
        "Root at y = 0 and tip sections in the xz plane, ruled upper and lower surfaces",
    ),
}


def parse_field(field: Field, text: str):
    text = text.strip()
    try:
        if field.kind == "text":
            if field.key == "name" and text == "":
                return None
            return text
        if field.kind == "bool":
            return text.lower() in ("yes", "y", "true", "1", "on")
        if field.kind == "number":
            return float(text)
        if field.kind == "int":
            return int(text)
        parts = [float(p) for p in text.replace(",", " ").split()]
        if field.kind == "vec2" and len(parts) == 2:
            return (parts[0], parts[1])
        if field.kind == "vec3" and len(parts) == 3:
            return (parts[0], parts[1], parts[2])
        raise ValueError(f"{len(parts)} numbers given")
    except ValueError as error:
        raise ShapeParseError(f"{field.label}: {text!r} is not a valid value ({error})") from None


def parse_values(spec: ShapeSpec, texts: dict) -> dict:
    """{key: value} from {key: text}; fields not given take their defaults."""
    values = {}
    for field in spec.fields:
        text = texts.get(field.key, field.default)
        values[field.key] = parse_field(field, text)
    return values


class ShapeDialog(QDialog):
    def __init__(self, spec: ShapeSpec, parent=None):
        super().__init__(parent)
        self.spec = spec
        self.setWindowTitle(spec.title)
        layout = QVBoxLayout(self)
        hint = QLabel(spec.hint)
        hint.setWordWrap(True)
        layout.addWidget(hint)
        form = QFormLayout()
        self.edits = {}
        for field in spec.fields:
            if field.kind == "bool":
                box = QCheckBox()
                box.setChecked(parse_field(field, field.default))
                form.addRow(field.label, box)
                self.edits[field.key] = box
                continue
            edit = QLineEdit(field.default)
            if field.key == "name":
                edit.setText("")
                edit.setPlaceholderText(f"{spec.method.removeprefix('create_')}-<n>")
            form.addRow(field.label, edit)
            self.edits[field.key] = edit
        layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def texts(self) -> dict:
        out = {}
        for key, widget in self.edits.items():
            if isinstance(widget, QCheckBox):
                if widget.isChecked():
                    out[key] = "yes"
                else:
                    out[key] = "no"
            else:
                out[key] = widget.text()
        return out

    def values(self) -> dict:
        return parse_values(self.spec, self.texts())

    def _accept(self) -> None:
        try:
            self.values()
        except ShapeParseError as error:
            QMessageBox.warning(self, self.spec.title, str(error))
            return
        self.accept()
