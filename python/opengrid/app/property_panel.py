"""Property panel: view and edit the parameters of the selected entity.

Edits are applied through ``Session.modify`` / ``Session.rename``, i.e. as undoable commands.
Input that cannot be parsed or that the model rejects is reported in the panel; nothing is
changed in that case.
"""

import math

from PySide6.QtCore import Qt, Signal
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
from opengrid.model.blocks import EXTRUSION_MODES
from opengrid.model.connectors import CURVE_KINDS, DISTRIBUTION_TYPES, EXPLICIT, SPACING_MODES
from opengrid.model.domains import (
    CONTROLS,
    DEFAULT_ELLIPTIC,
    DEFAULT_WALL,
    EDGE_NAMES,
    MAX_ACCELERATION,
    METHODS,
    SOLVERS,
)
from opengrid.model.entities import SPLINE_METHODS
from opengrid.model.patches import PATCH_TYPES, member_label
from opengrid.model.periodic import DEFAULT_AXIS
from opengrid.model.unstructured import (
    ALGORITHMS,
    DEFAULT_GENERATOR,
    ELEMENT_TYPES,
    SMOOTHING_METHODS,
)
from opengrid.model.volumes import BACKENDS, DEFAULT_VOLUME_INFLATION, PYRAMID_MODES
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
    "plane": (("origin", "vec3"), ("u_vector", "vec3"), ("v_vector", "vec3")),
    "ruled_surface": (),
    "extruded_surface": (("direction", "vec3"),),
    "revolved_surface": (
        ("axis_point", "vec3"),
        ("axis", "vec3"),
        ("start_angle", "number"),
        ("end_angle", "number"),
    ),
    "nurbs_surface": (("knots_u", "numbers"), ("knots_v", "numbers")),
}


class ExtrusionEditor(QWidget):
    """Extrusion of a block: mode, vector / axis, angle (degrees in the editor, radians in the
    model), distance, number of points and distribution."""

    def __init__(self, value: dict, parent=None):
        super().__init__(parent)
        self.domain = value["domain"]
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.mode = QComboBox()
        self.mode.addItems(list(EXTRUSION_MODES))
        self.mode.setCurrentText(value["mode"])
        self.vector = QLineEdit(format_vec3(value["vector"]))
        self.axis_point = QLineEdit(format_vec3(value["axis_point"]))
        self.angle = QLineEdit(format_number(math.degrees(value["angle"])))
        self.distance = QLineEdit(format_number(value["distance"]))
        self.dimension = QSpinBox()
        self.dimension.setRange(2, 100000)
        self.dimension.setValue(value["dimension"])
        self.distribution = DistributionEditor(
            value["distribution"], length=1.0, dimension=value["dimension"], unit=" of the length"
        )
        self.edges = QLineEdit(", ".join(value["edges"]))
        self.edges.setToolTip(
            "hyperbolic: conditions of the base edges imin, imax, jmin, jmax "
            "(free, constant_x, constant_y, constant_z, periodic)"
        )
        self.dissipation = QLineEdit(format_number(value["dissipation"]))
        layout.addRow("mode", self.mode)
        layout.addRow("vector / axis", self.vector)
        layout.addRow("axis point", self.axis_point)
        layout.addRow("angle [deg]", self.angle)
        layout.addRow("distance", self.distance)
        layout.addRow("points", self.dimension)
        layout.addRow("distribution", self.distribution)
        layout.addRow("edges", self.edges)
        layout.addRow("dissipation", self.dissipation)

    def value(self) -> dict:
        def vec(edit, label):
            values = parse_numbers(edit.text(), label)
            if len(values) not in (2, 3):
                raise ParseError(f"{label}: expected 2 or 3 numbers")
            if len(values) == 2:
                values.append(0.0)
            return tuple(values)

        def one(edit, label):
            values = parse_numbers(edit.text(), label)
            if len(values) != 1:
                raise ParseError(f"{label}: expected one number")
            return values[0]

        distribution = self.distribution.value()
        distribution.pop("mode", None)
        return {
            "domain": self.domain,
            "mode": self.mode.currentText(),
            "vector": vec(self.vector, "vector"),
            "axis_point": vec(self.axis_point, "axis point"),
            "angle": math.radians(one(self.angle, "angle")),
            "distance": one(self.distance, "distance"),
            "dimension": self.dimension.value(),
            "distribution": distribution,
            "edges": tuple(e.strip() for e in self.edges.text().split(",") if e.strip()),
            "dissipation": one(self.dissipation, "dissipation"),
        }


class PeriodicEditor(QWidget):
    """Periodic pairing of a patch: the partner (another periodic patch, or none), the rotation
    axis, centre and angle (degrees here, radians in the model) and the translation (empty: the
    translation that maps the centroid of this patch's nodes onto the partner's)."""

    def __init__(self, value: dict | None, partners, paired_from=(), parent=None):
        super().__init__(parent)
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.partner = QComboBox()
        self.partner.addItem("(none)", None)
        for name, entity_id in partners:
            self.partner.addItem(name, entity_id)
        axis = DEFAULT_AXIS
        center = (0.0, 0.0, 0.0)
        angle = 0.0
        translation = ""
        if value is not None:
            index = self.partner.findData(value["partner"])
            if index >= 0:
                self.partner.setCurrentIndex(index)
            axis = value["rotation_axis"]
            center = value["rotation_center"]
            angle = value["rotation_angle"]
            translation = format_vec3(value["translation"])
        self.axis = QLineEdit(format_vec3(axis))
        self.center = QLineEdit(format_vec3(center))
        self.angle = QLineEdit(format_number(math.degrees(angle)))
        self.translation = QLineEdit(translation)
        self.translation.setPlaceholderText("automatic (centroid to centroid)")
        layout.addRow("partner", self.partner)
        if paired_from:
            note = QLabel(f"paired from {', '.join(paired_from)} (pair from one side only)")
            note.setWordWrap(True)
            layout.addRow(note)
        layout.addRow("rotation axis", self.axis)
        layout.addRow("rotation centre", self.center)
        layout.addRow("angle [deg]", self.angle)
        layout.addRow("translation", self.translation)

    def value(self) -> dict | None:
        """None without a partner; otherwise the pairing, ``translation`` None when the
        automatic one is wanted."""
        partner = self.partner.currentData()
        if partner is None:
            return None
        translation = None
        if self.translation.text().strip():
            translation = parse_vec3(self.translation.text(), "translation")
        angle = parse_numbers(self.angle.text(), "angle")
        if len(angle) != 1:
            raise ParseError("angle: expected one number (degrees)")
        return {
            "partner": partner,
            "rotation_center": parse_vec3(self.center.text(), "rotation centre"),
            "rotation_axis": parse_vec3(self.axis.text(), "rotation axis"),
            "rotation_angle": math.radians(angle[0]),
            "translation": translation,
        }


FIRST_HEIGHT_FRACTION = 0.05  # suggested first layer: a twentieth of the boundary spacing


def mean_connector_spacing(session, connector_ids) -> float | None:
    """Mean node interval of the given connectors (None without any)."""
    import numpy as np

    lengths = []
    for connector_id in connector_ids:
        if connector_id not in session.store:
            continue
        points = np.asarray(session.connector_nodes(session.store.get(connector_id).name).points)
        if len(points) >= 2:
            lengths.extend(np.linalg.norm(np.diff(points, axis=0), axis=1).tolist())
    if not lengths:
        return None
    return float(sum(lengths) / len(lengths))


def _first_height_placeholder(boundary_spacing) -> str:
    if boundary_spacing is None:
        return "required with walls"
    return f"suggested {FIRST_HEIGHT_FRACTION * boundary_spacing:.3g}"


def _inflation_hint(boundary_spacing) -> QLabel:
    if boundary_spacing is None:
        text = "first height: a fraction of the boundary spacing (1/20 is a usual start)"
    else:
        first = FIRST_HEIGHT_FRACTION * boundary_spacing
        text = (
            f"boundary spacing is {boundary_spacing:.3g}: first height {first:.3g} suggested "
            f"(used when empty); the layers stop at the core cell size unless a total "
            "thickness is given"
        )
    label = QLabel(text)
    label.setWordWrap(True)
    return label


def _suggested_first_height(editor) -> float | None:
    if editor._boundary_spacing is None:
        return None
    return FIRST_HEIGHT_FRACTION * editor._boundary_spacing


class VolumeEditor(QWidget):
    """Volume generator: backend, max size (empty: none), grading, pyramids, quality
    optimisation and floor, and boundary layers (wall faces; empty: none)."""

    def __init__(self, value: dict, store=None, parent=None, boundary_spacing=None):
        super().__init__(parent)
        self._boundary_spacing = boundary_spacing  # mean boundary node spacing: suggestions
        self.store = store
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.backend = QComboBox()
        self.backend.addItems(list(BACKENDS))
        self.backend.setCurrentText(value["backend"])
        self.max_size = QLineEdit(format_optional(value["max_size"]))
        self.max_size.setPlaceholderText("none")
        self.grading = QLineEdit(format_number(value["grading"]))
        self.pyramids = QComboBox()
        self.pyramids.addItems(list(PYRAMID_MODES))
        self.pyramids.setCurrentText(value["pyramids"])
        layout.addRow("backend", self.backend)
        layout.addRow("max size", self.max_size)
        layout.addRow("grading", self.grading)
        layout.addRow("pyramids", self.pyramids)
        self.optimize = QCheckBox("optimise quality")
        self.optimize.setChecked(value["optimize"])
        layout.addRow(self.optimize)
        floor = value["quality_floor"]
        self.floor_jacobian = QLineEdit(format_number(floor["min_jacobian"]))
        self.floor_skewness = QLineEdit(format_number(floor["max_skewness"]))
        layout.addRow("min scaled Jacobian", self.floor_jacobian)
        layout.addRow("max skewness", self.floor_skewness)
        inflation = value["inflation"]
        if inflation is None:
            inflation = dict(DEFAULT_VOLUME_INFLATION)
        names = []
        for ref in inflation["faces"]:
            if store is not None and ref in store:
                names.append(store.get(ref).name)
        self.walls = QLineEdit(", ".join(names))
        self.walls.setPlaceholderText("wall faces (none: no layers)")
        self.first_height = QLineEdit(format_optional(inflation["first_height"]))
        self.layer_growth = QLineEdit(format_number(inflation["growth_rate"]))
        self.layers = QSpinBox()
        self.layers.setRange(1, 1000)
        self.layers.setValue(inflation["layers"])
        self.thickness = QLineEdit(format_optional(inflation["total_thickness"]))
        self.thickness.setPlaceholderText("none")
        self.transition = QCheckBox("stop at core size")
        self.transition.setChecked(inflation["transition"])
        self.slide = QCheckBox("slide on planar faces")
        self.slide.setChecked(inflation["slide"])
        layout.addRow("inflation walls", self.walls)
        self.first_height.setPlaceholderText(_first_height_placeholder(boundary_spacing))
        layout.addRow("first height", self.first_height)
        layout.addRow("", _inflation_hint(boundary_spacing))
        layout.addRow("layer growth", self.layer_growth)
        layout.addRow("max layers", self.layers)
        layout.addRow("total thickness", self.thickness)
        layout.addRow(self.transition)
        layout.addRow(self.slide)

    def value(self) -> dict:
        grading = parse_numbers(self.grading.text(), "grading")
        if len(grading) != 1:
            raise ParseError("grading: expected one number")
        return {
            "backend": self.backend.currentText(),
            "max_size": parse_optional(self.max_size.text(), "max size"),
            "grading": grading[0],
            "pyramids": self.pyramids.currentText(),
            "inflation": self._inflation(),
            "optimize": self.optimize.isChecked(),
            "quality_floor": {
                "min_jacobian": self._one(self.floor_jacobian, "min scaled Jacobian"),
                "max_skewness": self._one(self.floor_skewness, "max skewness"),
            },
        }

    @staticmethod
    def _one(field, label):
        values = parse_numbers(field.text(), label)
        if len(values) != 1:
            raise ParseError(f"{label}: expected one number")
        return values[0]

    def suggested_first_height(self) -> float | None:
        return _suggested_first_height(self)

    def _inflation(self):
        names = [n.strip() for n in self.walls.text().replace(";", ",").split(",") if n.strip()]
        if not names:
            return None
        ids = []
        for name in names:
            if self.store is None or name not in self.store:
                raise ParseError(f"inflation walls: no domain {name!r}")
            ids.append(self.store.get(name).id)
        first = parse_optional(self.first_height.text(), "first height")
        if first is None:
            first = self.suggested_first_height()  # the hint's value when the field is empty
        if first is None:
            raise ParseError("first height: required when inflation walls are given")
        growth = parse_numbers(self.layer_growth.text(), "layer growth")
        if len(growth) != 1:
            raise ParseError("layer growth: expected one number")
        return {
            "faces": tuple(ids),
            "first_height": first,
            "growth_rate": growth[0],
            "layers": self.layers.value(),
            "total_thickness": parse_optional(self.thickness.text(), "total thickness"),
            "transition": self.transition.isChecked(),
            "slide": self.slide.isChecked(),
        }


class BlockSmoothingEditor(QWidget):
    """3D elliptic smoothing of a block: control (none = off), wall faces with a first spacing,
    maximum iterations."""

    def __init__(self, value: dict | None, parent=None):
        super().__init__(parent)
        self.original = value
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.control = QComboBox()
        self.control.addItems(["none", "laplace", "thomas_middlecoff"])
        walls = ""
        spacing = ""
        iterations = 500
        if value is not None:
            self.control.setCurrentText(value["control"])
            walls = ", ".join(w["face"] for w in value["walls"])
            if value["walls"] and value["walls"][0]["spacing"] is not None:
                spacing = format_number(value["walls"][0]["spacing"])
            iterations = value["max_iterations"]
        self.walls = QLineEdit(walls)
        self.walls.setPlaceholderText("wall faces, e.g. kmin")
        self.spacing = QLineEdit(spacing)
        self.spacing.setPlaceholderText("from the block edges")
        self.iterations = QSpinBox()
        self.iterations.setRange(1, 1_000_000)
        self.iterations.setValue(iterations)
        layout.addRow("control", self.control)
        layout.addRow("wall faces", self.walls)
        layout.addRow("wall spacing", self.spacing)
        layout.addRow("iterations", self.iterations)

    def value(self) -> dict | None:
        if self.control.currentText() == "none":
            return None
        spacing = parse_optional(self.spacing.text(), "wall spacing")
        faces = [f.strip() for f in self.walls.text().split(",") if f.strip()]
        previous = {}
        if self.original is not None:
            previous = {w["face"]: w for w in self.original["walls"]}
        walls = []
        for face in faces:
            wall = dict(
                previous.get(
                    face, {"face": face, "decay": 0.1, "relaxation": 0.1, "max_source": 100.0}
                )
            )
            wall["spacing"] = spacing
            walls.append(wall)
        result = {
            "control": self.control.currentText(),
            "walls": tuple(walls),
            "max_iterations": self.iterations.value(),
        }
        if self.original is not None:
            result = {**self.original, **result}
        else:
            result.update({"omega": 1.0, "tolerance": 1e-6})
        return result


CONNECTOR_FIELDS = (
    ("dimension", "dimension"),
    ("mode", "mode"),
    ("distribution", "distribution"),
    ("breakpoints", "breakpoints"),
    ("propagate_spacing", "flag"),
)


class ParseError(ValueError):
    pass


def format_group_walls(walls, store) -> str:
    lines = []
    for wall in walls:
        key = "connector"
        if "domain" in wall:
            key = "domain"
        line = store.get(wall[key]).name
        layers = wall.get("layers", 0)
        if wall["spacing"] is not None:
            line += f" {format_number(wall['spacing'])}"
        elif layers > 0:
            line += " -"
        if layers > 0:
            line += f" {layers} {format_number(wall.get('growth', 1.2))}"
        lines.append(line)
    return "\n".join(lines)


def parse_group_walls(text: str, key: str = "connector") -> tuple:
    """'name [spacing [layers [growth]]]' per line -> wall dicts (decay etc. at their
    defaults; spacing "-" = the initial height); ``key`` is "connector" (2D groups) or
    "domain" (block faces, 3D groups)."""
    walls = []
    for raw in text.splitlines():
        parts = raw.split()
        if not parts:
            continue
        if len(parts) > 4:
            raise ParseError(
                f"walls: expected '{key} [spacing [layers [growth]]]', got {raw.strip()!r}"
            )
        wall = {key: parts[0], "spacing": None}
        if len(parts) >= 2:
            wall["spacing"] = parse_optional(parts[1], f"walls: {parts[0]}")
        if len(parts) >= 3:
            try:
                wall["layers"] = int(parts[2])
            except ValueError as error:
                raise ParseError(f"walls: {parts[0]}: layers must be a whole number") from error
        if len(parts) == 4:
            wall["growth"] = parse_numbers(parts[3], f"walls: {parts[0]}")[0]
        walls.append(wall)
    return tuple(walls)


def format_slide(slide, store) -> str:
    if slide == "auto":
        return "auto"
    if not slide:
        return "none"
    return ", ".join(store.get(c).name for c in slide)


def parse_slide(text: str):
    text = text.strip()
    if text in ("", "auto"):
        return "auto"
    if text == "none":
        return ()
    return tuple(part.strip() for part in text.split(",") if part.strip())


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


def format_attractions(attractions) -> str:
    """One per line: 'xi_line xi amplitude decay', 'eta_line eta amplitude decay' or
    'point xi eta amplitude decay'."""
    lines = []
    for item in attractions:
        if item["kind"] == "xi_line":
            values = [item["xi"], item["amplitude"], item["decay"]]
        elif item["kind"] == "eta_line":
            values = [item["eta"], item["amplitude"], item["decay"]]
        else:
            values = [item["xi"], item["eta"], item["amplitude"], item["decay"]]
        lines.append(" ".join([item["kind"], *(format_number(v) for v in values)]))
    return "\n".join(lines)


def parse_attractions(text: str) -> tuple:
    result = []
    for number, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        where = f"attractions line {number}"
        kind = fields[0]
        values = parse_numbers(" ".join(fields[1:]), where)
        if kind in ("xi_line", "eta_line") and len(values) == 3:
            item = {
                "kind": kind,
                "xi": None,
                "eta": None,
                "amplitude": values[1],
                "decay": values[2],
            }
            if kind == "xi_line":
                item["xi"] = values[0]
            else:
                item["eta"] = values[0]
        elif kind == "point" and len(values) == 4:
            item = {
                "kind": kind,
                "xi": values[0],
                "eta": values[1],
                "amplitude": values[2],
                "decay": values[3],
            }
        else:
            raise ParseError(
                f"{where}: expected 'xi_line|eta_line pos amp decay' or 'point xi eta amp decay'"
            )
        result.append(item)
    return tuple(result)


class GeneratorEditor(QWidget):
    """Structured-domain generator: TFI type, Hermite walls and elliptic smoothing."""

    def __init__(self, value: dict, parent=None):
        super().__init__(parent)
        self.original = value
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.method = QComboBox()
        self.method.addItems(list(METHODS))
        self.method.setCurrentText(value["method"])
        self.hermite_bottom = QCheckBox("j = 0")
        self.hermite_top = QCheckBox("j = jmax")
        self.hermite_bottom.setChecked(0 in value["hermite_walls"])
        self.hermite_top.setChecked(2 in value["hermite_walls"])
        walls_row = QWidget()
        walls_layout = QHBoxLayout(walls_row)
        walls_layout.setContentsMargins(0, 0, 0, 0)
        walls_layout.addWidget(self.hermite_bottom)
        walls_layout.addWidget(self.hermite_top)
        self.strength = QLineEdit(format_number(value["hermite_strength"]))
        layout.addRow("method", self.method)
        layout.addRow("hermite walls", walls_row)
        layout.addRow("hermite strength", self.strength)

        elliptic = value["elliptic"]
        self.elliptic_enabled = QCheckBox("elliptic smoothing")
        self.elliptic_enabled.setChecked(elliptic is not None)
        if elliptic is None:
            elliptic = DEFAULT_ELLIPTIC
        self.elliptic_original = elliptic
        self.control = QComboBox()
        self.control.addItems(list(CONTROLS))
        self.control.setCurrentText(elliptic["control"])
        self.solver = QComboBox()
        self.solver.addItems(list(SOLVERS))
        self.solver.setCurrentText(elliptic["solver"])
        self.omega = QLineEdit(format_number(elliptic["omega"]))
        self.iterations = QSpinBox()
        self.iterations.setRange(1, 10_000_000)
        self.iterations.setValue(elliptic["max_iterations"])
        self.tolerance = QLineEdit(format_number(elliptic["tolerance"]))
        self.acceleration = QSpinBox()
        self.acceleration.setRange(0, MAX_ACCELERATION)
        self.acceleration.setValue(elliptic["acceleration"])
        self.acceleration.setToolTip("Anderson acceleration depth (0: plain relaxation)")
        walls = {wall["edge"]: wall for wall in elliptic["walls"]}
        self.wall_boxes = {}
        self.wall_spacing = {}
        for edge, label in ((0, "wall control j = 0"), (2, "wall control j = jmax")):
            box = QCheckBox(label)
            box.setChecked(edge in walls)
            spacing = QLineEdit()
            spacing.setPlaceholderText("from side boundaries")
            if edge in walls:
                spacing.setText(format_optional(walls[edge]["spacing"]))
            self.wall_boxes[edge] = box
            self.wall_spacing[edge] = spacing
            layout.addRow(box, spacing)
        self.walls_original = walls
        self.attractions = QPlainTextEdit(format_attractions(elliptic["attractions"]))
        self.attractions.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.attractions.setPlaceholderText(
            "ttm control: 'xi_line 20 1.0 0.5' / 'point 10 5 1.0 0.5'"
        )
        self.attractions.setMaximumHeight(60)
        layout.insertRow(3, self.elliptic_enabled)
        layout.insertRow(4, "control", self.control)
        layout.insertRow(5, "solver", self.solver)
        layout.insertRow(6, "omega", self.omega)
        layout.insertRow(7, "max iterations", self.iterations)
        layout.insertRow(8, "tolerance", self.tolerance)
        layout.insertRow(9, "acceleration", self.acceleration)
        layout.addRow("attractions", self.attractions)

    def value(self) -> dict:
        walls = []
        if self.hermite_bottom.isChecked():
            walls.append(0)
        if self.hermite_top.isChecked():
            walls.append(2)
        strength = parse_numbers(self.strength.text(), "hermite strength")
        if len(strength) != 1:
            raise ParseError("hermite strength: expected one number")
        result = {
            "method": self.method.currentText(),
            "hermite_walls": tuple(walls),
            "hermite_strength": strength[0],
            "elliptic": None,
        }
        if not self.elliptic_enabled.isChecked():
            return result
        omega = parse_numbers(self.omega.text(), "omega")
        tolerance = parse_numbers(self.tolerance.text(), "tolerance")
        if len(omega) != 1 or len(tolerance) != 1:
            raise ParseError("omega and tolerance: expected one number each")
        wall_items = []
        for edge in (0, 2):
            if not self.wall_boxes[edge].isChecked():
                continue
            wall = dict(self.walls_original.get(edge, DEFAULT_WALL))
            wall["edge"] = edge
            wall["spacing"] = parse_optional(self.wall_spacing[edge].text(), "wall spacing")
            wall_items.append(wall)
        result["elliptic"] = {
            "control": self.control.currentText(),
            "attractions": parse_attractions(self.attractions.toPlainText()),
            "walls": tuple(wall_items),
            "solver": self.solver.currentText(),
            "omega": omega[0],
            "max_iterations": self.iterations.value(),
            "tolerance": tolerance[0],
            "acceleration": self.acceleration.value(),
        }
        return result


def format_sources(sources) -> str:
    """One per line: 'x0 y0 size radius' (point) or 'x0 y0 x1 y1 size radius' (line)."""
    lines = []
    for item in sources:
        values = list(item["p0"][:2])
        if item["p1"] != item["p0"]:
            values += list(item["p1"][:2])
        values += [item["size"], item["radius"]]
        lines.append(" ".join(format_number(v) for v in values))
    return "\n".join(lines)


def parse_sources(text: str) -> tuple:
    result = []
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        values = parse_numbers(line, f"sources line {number}")
        if len(values) == 4:
            p0 = (values[0], values[1], 0.0)
            result.append({"p0": p0, "p1": p0, "size": values[2], "radius": values[3]})
        elif len(values) == 6:
            result.append(
                {
                    "p0": (values[0], values[1], 0.0),
                    "p1": (values[2], values[3], 0.0),
                    "size": values[4],
                    "radius": values[5],
                }
            )
        else:
            raise ParseError(
                f"sources line {number}: expected 'x0 y0 size radius' or 'x0 y0 x1 y1 size radius'"
            )
    return tuple(result)


class UnstructuredEditor(QWidget):
    """Triangulation settings of an unstructured domain."""

    def __init__(self, value: dict, store=None, parent=None, boundary_spacing=None):
        super().__init__(parent)
        self._boundary_spacing = boundary_spacing  # mean boundary node spacing: suggestions
        self.store = store
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.algorithm = QComboBox()
        self.algorithm.addItems(list(ALGORITHMS))
        self.algorithm.setCurrentText(value["algorithm"])
        self.refine = QCheckBox("refine")
        self.refine.setChecked(value["refine"])
        self.min_angle = QLineEdit(format_number(value["min_angle"]))
        self.growth = QLineEdit(format_number(value["growth_rate"]))
        self.max_size = QLineEdit(format_optional(value["max_size"]))
        self.max_size.setPlaceholderText("none")
        self.max_points = QSpinBox()
        self.max_points.setRange(3, 2_000_000_000)
        self.max_points.setValue(value["max_points"])
        self.sources = QPlainTextEdit(format_sources(value["sources"]))
        self.sources.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.sources.setPlaceholderText("x0 y0 size radius  |  x0 y0 x1 y1 size radius")
        self.sources.setMaximumHeight(60)
        self.element_type = QComboBox()
        self.element_type.addItems(list(ELEMENT_TYPES))
        self.element_type.setCurrentText(value["element_type"])
        self.smoothing = QComboBox()
        self.smoothing.addItems(list(SMOOTHING_METHODS))
        self.smoothing.setCurrentText(value["smoothing"])
        self.smoothing_iterations = QSpinBox()
        self.smoothing_iterations.setRange(0, 1000)
        self.smoothing_iterations.setValue(value["smoothing_iterations"])
        self.edge_swap = QCheckBox("edge swapping")
        self.edge_swap.setChecked(value["edge_swap"])
        layout.addRow("element type", self.element_type)
        layout.addRow("algorithm", self.algorithm)
        layout.addRow(self.refine)
        layout.addRow("min angle [deg]", self.min_angle)
        layout.addRow("growth rate", self.growth)
        layout.addRow("max size", self.max_size)
        layout.addRow("max points", self.max_points)
        layout.addRow("sources", self.sources)
        layout.addRow("smoothing", self.smoothing)
        layout.addRow("iterations", self.smoothing_iterations)
        layout.addRow(self.edge_swap)
        # Boundary layers: wall connectors (empty: no inflation) and layer settings.
        inflation = value["inflation"]
        if inflation is None:
            inflation = {
                "connectors": (),
                "first_height": None,
                "growth_rate": 1.2,
                "layers": 30,
                "total_thickness": None,
                "transition": True,
            }
        names = []
        for ref in inflation["connectors"]:
            if store is not None and ref in store:
                names.append(store.get(ref).name)
        self.walls = QLineEdit(", ".join(names))
        self.walls.setPlaceholderText("wall connectors (none: no layers)")
        self.first_height = QLineEdit(format_optional(inflation["first_height"]))
        self.layer_growth = QLineEdit(format_number(inflation["growth_rate"]))
        self.layers = QSpinBox()
        self.layers.setRange(1, 1000)
        self.layers.setValue(inflation["layers"])
        self.thickness = QLineEdit(format_optional(inflation["total_thickness"]))
        self.thickness.setPlaceholderText("none")
        self.transition = QCheckBox("stop at core size")
        self.transition.setChecked(inflation["transition"])
        layout.addRow("inflation walls", self.walls)
        self.first_height.setPlaceholderText(_first_height_placeholder(boundary_spacing))
        layout.addRow("first height", self.first_height)
        layout.addRow("", _inflation_hint(boundary_spacing))
        layout.addRow("layer growth", self.layer_growth)
        layout.addRow("max layers", self.layers)
        layout.addRow("total thickness", self.thickness)
        layout.addRow(self.transition)
        # refine and min angle apply to the Delaunay mesher only.
        self.algorithm.currentTextChanged.connect(self._update_enabled)
        self._update_enabled(self.algorithm.currentText())

    def _update_enabled(self, algorithm: str) -> None:
        delaunay = algorithm == "delaunay"
        self.refine.setEnabled(delaunay)
        self.min_angle.setEnabled(delaunay)

    def value(self) -> dict:
        def one(edit, label):
            values = parse_numbers(edit.text(), label)
            if len(values) != 1:
                raise ParseError(f"{label}: expected one number")
            return values[0]

        return {
            "algorithm": self.algorithm.currentText(),
            "refine": self.refine.isChecked(),
            "min_angle": one(self.min_angle, "min angle"),
            "growth_rate": one(self.growth, "growth rate"),
            "max_size": parse_optional(self.max_size.text(), "max size"),
            "sources": parse_sources(self.sources.toPlainText()),
            "max_points": self.max_points.value(),
            "element_type": self.element_type.currentText(),
            "smoothing": self.smoothing.currentText(),
            "smoothing_iterations": self.smoothing_iterations.value(),
            "edge_swap": self.edge_swap.isChecked(),
            "inflation": self._inflation(),
        }

    def suggested_first_height(self) -> float | None:
        return _suggested_first_height(self)

    def _inflation(self):
        names = [n.strip() for n in self.walls.text().replace(";", ",").split(",") if n.strip()]
        if not names:
            return None
        ids = []
        for name in names:
            if self.store is None or name not in self.store:
                raise ParseError(f"inflation walls: no connector {name!r}")
            ids.append(self.store.get(name).id)
        first = parse_optional(self.first_height.text(), "first height")
        if first is None:
            first = self.suggested_first_height()  # the hint's value when the field is empty
        if first is None:
            raise ParseError("first height: required when inflation walls are given")
        growth = parse_numbers(self.layer_growth.text(), "layer growth")
        if len(growth) != 1:
            raise ParseError("layer growth: expected one number")
        return {
            "connectors": tuple(ids),
            "first_height": first,
            "growth_rate": growth[0],
            "layers": self.layers.value(),
            "total_thickness": parse_optional(self.thickness.text(), "total thickness"),
            "transition": self.transition.isChecked(),
        }


class DistributionEditor(QWidget):
    """Type + begin/end spacing + growth ratio (empty fields mean 'not specified'). Explicit
    positions (from split / join) are kept while the type stays 'explicit'. With a length and
    a point count it says the uniform spacing and, when the type changes, fills the fields
    with a clustering to start from: tanh and sine a quarter of the uniform spacing at the
    ends they cluster, geometric a growth ratio of 1.2."""

    CLUSTER_FRACTION = 0.25
    GROWTH_RATIO = 1.2

    def __init__(self, value: dict, parent=None, length=None, dimension=None, unit: str = ""):
        super().__init__(parent)
        layout = QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.type = QComboBox()
        self.type.addItems(list(DISTRIBUTION_TYPES))
        self._positions = value.get("positions")
        if value["type"] == EXPLICIT:
            self.type.addItem(EXPLICIT)
        self.type.setCurrentText(value["type"])
        self.begin = QLineEdit(format_optional(value["begin_spacing"]))
        self.end = QLineEdit(format_optional(value["end_spacing"]))
        self.ratio = QLineEdit(format_optional(value["growth_ratio"]))
        for edit in (self.begin, self.end, self.ratio):
            edit.setPlaceholderText("none")
        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        layout.addRow("type", self.type)
        layout.addRow("", self.hint)
        layout.addRow("begin spacing", self.begin)
        layout.addRow("end spacing", self.end)
        layout.addRow("growth ratio", self.ratio)
        self._length = None
        self._dimension = None
        self._unit = unit
        self.set_context(length, dimension)
        self.type.currentTextChanged.connect(self._on_type_changed)

    def set_context(self, length, dimension) -> None:
        """The curve length and point count the suggestions are made from."""
        self._length = length
        self._dimension = dimension
        self._update_hint()

    def uniform_spacing(self) -> float | None:
        if self._length is None or self._dimension is None or self._dimension < 2:
            return None
        return float(self._length) / (int(self._dimension) - 1)

    def suggestions(self) -> dict:
        """Suggested begin / end spacing and growth ratio for the current type."""
        kind = self.type.currentText()
        h = self.uniform_spacing()
        out = {"begin": None, "end": None, "ratio": None}
        if kind == "geometric":
            out["ratio"] = self.GROWTH_RATIO
        elif kind in ("tanh", "sine") and h is not None:
            out["begin"] = self.CLUSTER_FRACTION * h
            if kind == "tanh":
                out["end"] = self.CLUSTER_FRACTION * h
        return out

    def _update_hint(self) -> None:
        h = self.uniform_spacing()
        kind = self.type.currentText()
        text = ""
        if h is not None:
            text = f"uniform spacing {h:.4g}{self._unit}"
        if kind == "tanh":
            text += "; clusters at the ends that have a spacing (leave one empty for one side)"
        elif kind == "sine":
            text += "; clusters at the beginning (begin spacing)"
        elif kind == "geometric":
            text += "; grows from the beginning by the ratio (or give one spacing)"
        elif kind == "uniform":
            text += "; no clustering"
        self.hint.setText(text.lstrip("; "))
        for edit, key in ((self.begin, "begin"), (self.end, "end"), (self.ratio, "ratio")):
            suggested = self.suggestions()[key]
            placeholder = "none"
            if suggested is not None:
                placeholder = f"suggested {suggested:.4g}"
            edit.setPlaceholderText(placeholder)

    def _on_type_changed(self, _kind: str) -> None:
        """A new type starts from its suggested clustering (the old fields rarely apply)."""
        for edit, key in ((self.begin, "begin"), (self.end, "end"), (self.ratio, "ratio")):
            suggested = self.suggestions()[key]
            edit.setText("")
            if suggested is not None:
                edit.setText(format_number(suggested))
        self._update_hint()

    def value(self) -> dict:
        result = {
            "type": self.type.currentText(),
            "begin_spacing": parse_optional(self.begin.text(), "begin spacing"),
            "end_spacing": parse_optional(self.end.text(), "end spacing"),
            "growth_ratio": parse_optional(self.ratio.text(), "growth ratio"),
        }
        if result["type"] == EXPLICIT:
            result["positions"] = self._positions
        return result


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
    smooth_requested = Signal(str)

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
        self._build()
        for edit in self.findChildren(QLineEdit):
            edit.returnPressed.connect(self.apply)  # Enter applies (the widgets are new)

    def _build(self) -> None:
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
        if entity.kind == "domain":
            self._domain_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind == "unstructured_domain":
            self._unstructured_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind == "volume":
            self._volume_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind == "block":
            self._block_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind == "smoothing_group":
            self._smoothing_group_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind in ("patch", "zone"):
            self._group_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        if entity.kind == "sketch":
            self._sketch_rows(entity)
            self.apply_button.setEnabled(True)
            self.revert_button.setEnabled(True)
            return
        for field, kind in FIELDS[entity.kind]:
            editor = self._make_editor(kind, entity.params[field])
            self._form.addRow(field, editor)
            self._editors[field] = (kind, editor)
        if entity.kind in CURVE_KINDS:
            length = self.session.length(entity.name)
            self._form.addRow("length", QLabel(f"{length:.10g}"))
        self.apply_button.setEnabled(True)
        self.revert_button.setEnabled(True)

    def _unstructured_rows(self, entity) -> None:
        lines = []
        for k, loop in enumerate(entity.params["loops"]):
            names = []
            for item in loop:
                name = self.session.store.get(item["connector"]).name
                if item["reversed"]:
                    name += " (rev)"
                names.append(name)
            kind = "outer"
            if k > 0:
                kind = f"hole {k}"
            lines.append(f"{kind}: {', '.join(names)}")
        loops = QLabel("\n".join(lines))
        loops.setWordWrap(True)
        self._form.addRow("loops", loops)
        status = QLabel(self.session.domain_status(entity.name))
        status.setWordWrap(True)
        self._form.addRow("status", status)
        loop_connectors = [item["connector"] for loop in entity.params["loops"] for item in loop]
        editor = UnstructuredEditor(
            entity.params["generator"],
            self.session.store,
            boundary_spacing=mean_connector_spacing(self.session, loop_connectors),
        )
        self._form.addRow("generator", editor)
        self._editors["generator"] = ("generator", editor)
        from opengrid.model.sketch import output_owner

        owner = output_owner(self.session.store).get(entity.id)
        if owner is not None:
            hint = QLabel(
                f"a block of sketch {owner[0].name}: these settings go to the sketch (the "
                "triangulation settings are shared by its unstructured blocks; the inflation "
                "walls become sketch walls, stopping at the core cell size)"
            )
            hint.setWordWrap(True)
            self._form.addRow("", hint)

    def _volume_rows(self, entity) -> None:
        try:
            text = self.session.volume_status(entity.name)
        except ValidationError as error:
            text = str(error)
        status = QLabel(text)
        status.setWordWrap(True)
        self._form.addRow("status", status)
        if entity.params["faces"] is None:
            source = QLabel(self.session.store.get(entity.params["extrusion"]["domain"]).name)
            self._form.addRow("base domain", source)
            editor = ExtrusionEditor(entity.params["extrusion"])
            self._form.addRow("extrusion", editor)
            self._editors["extrusion"] = ("extrusion", editor)
            return
        faces = QLabel(", ".join(self.session.store.get(f).name for f in entity.params["faces"]))
        faces.setWordWrap(True)
        self._form.addRow("faces", faces)
        face_connectors = []
        for face in entity.params["faces"] or ():
            face_params = self.session.store.get(face).params
            for edge in face_params.get("edges", ()):
                face_connectors.extend(item["connector"] for item in edge)
            for loop in face_params.get("loops", ()):
                face_connectors.extend(item["connector"] for item in loop)
        editor = VolumeEditor(
            entity.params["generator"],
            self.session.store,
            boundary_spacing=mean_connector_spacing(self.session, face_connectors),
        )
        self._form.addRow("generator", editor)
        self._editors["generator"] = ("generator", editor)

    def _block_rows(self, entity) -> None:
        params = entity.params
        try:
            text = self.session.block_status(entity.name)
        except ValidationError as error:
            text = str(error)
        status = QLabel(text)
        status.setWordWrap(True)
        self._form.addRow("status", status)
        smoothing = BlockSmoothingEditor(params["smoothing"])
        self._form.addRow("smoothing", smoothing)
        self._editors["smoothing"] = ("smoothing", smoothing)
        if params["faces"] is not None:
            faces = QLabel(", ".join(self.session.store.get(f).name for f in params["faces"]))
            faces.setWordWrap(True)
            self._form.addRow("faces", faces)
            return
        source = QLabel(self.session.store.get(params["extrusion"]["domain"]).name)
        self._form.addRow("base domain", source)
        editor = ExtrusionEditor(params["extrusion"])
        self._form.addRow("extrusion", editor)
        self._editors["extrusion"] = ("extrusion", editor)

    def _smoothing_group_rows(self, entity) -> None:
        params = entity.params
        store = self.session.store
        try:
            text = self.session.group_status(entity.name)
        except ValidationError as error:
            text = str(error)
        status = QLabel(text)
        status.setWordWrap(True)
        self._form.addRow("status", status)
        members = "domains"
        wall_kind = "connector"
        if params["blocks"]:
            members = "blocks"
            wall_kind = "face (domain)"
        names = QPlainTextEdit("\n".join(store.get(d).name for d in params[members]))
        names.setPlaceholderText(f"one {members[:-1]} per line")
        names.setMinimumHeight(60)
        self._form.addRow(members, names)
        self._editors[members] = ("names", names)
        walls = QPlainTextEdit(format_group_walls(params["walls"], store))
        walls.setPlaceholderText(f"{wall_kind} [spacing] per line (no spacing: keep the initial)")
        walls.setMinimumHeight(50)
        self._form.addRow("walls", walls)
        self._editors["walls"] = ("group_walls", walls)
        slide = QLineEdit(format_slide(params["slide"], store))
        slide.setToolTip("auto (every non-wall boundary connector), none, or connector names")
        self._form.addRow("slide", slide)
        self._editors["slide"] = ("slide", slide)
        for field in ("omega", "tolerance"):
            editor = QLineEdit(format_number(params[field]))
            self._form.addRow(field, editor)
            self._editors[field] = ("number", editor)
        iterations = QSpinBox()
        iterations.setRange(1, 1_000_000)
        iterations.setValue(params["max_iterations"])
        self._form.addRow("max iterations", iterations)
        self._editors["max_iterations"] = ("spin", iterations)
        acceleration = QSpinBox()
        acceleration.setRange(0, MAX_ACCELERATION)
        acceleration.setValue(params["acceleration"])
        acceleration.setToolTip("Anderson acceleration depth (0: plain relaxation)")
        self._form.addRow("acceleration", acceleration)
        self._editors["acceleration"] = ("spin", acceleration)
        self.smooth_button = QPushButton("Run smoother")
        self.smooth_button.clicked.connect(lambda: self.smooth_requested.emit(entity.name))
        self._form.addRow("", self.smooth_button)

    def _sketch_rows(self, entity) -> None:
        params = entity.params
        report = QLabel(self.session.sketch_report(entity.name))
        report.setWordWrap(True)
        report.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._form.addRow("sketch", report)
        dimensions = QPlainTextEdit(
            "\n".join(f"{d['edge'][0]} {d['edge'][1]} {d['points']}" for d in params["dimensions"])
        )
        dimensions.setPlaceholderText("a b points (one edge of a family per line)")
        dimensions.setMinimumHeight(60)
        self._form.addRow("dimensions", dimensions)
        self._editors["dimensions"] = ("sketch_dimensions", dimensions)
        lines = []
        for wall in params["walls"]:
            ref = wall.get("edge", wall.get("face"))
            text = " ".join(str(v) for v in ref)
            layers = wall.get("layers", 0)
            if wall["spacing"] is not None:
                text += f" {format_number(wall['spacing'])}"
            elif layers > 0:
                text += " -"
            if layers > 0:
                text += f" {layers} {format_number(wall.get('growth', 1.2))}"
            lines.append(text)
        walls = QPlainTextEdit("\n".join(lines))
        walls.setPlaceholderText(
            "edge a b (3D: face a b c d) [first-cell height [layers [growth]]] per line"
        )
        walls.setMinimumHeight(50)
        self._form.addRow("walls", walls)
        self._editors["walls"] = ("sketch_walls", walls)
        own = [
            entity_id
            for _role, entity_id in params["outputs"]
            if entity_id in self.session.store
            and self.session.store.get(entity_id).kind == "connector"
        ]
        spacing = mean_connector_spacing(self.session, own)
        if spacing is not None:
            hint = QLabel(
                f"edge spacing is {spacing:.3g}: a first-cell height of "
                f"{FIRST_HEIGHT_FRACTION * spacing:.3g} is a usual start for a wall"
            )
            hint.setWordWrap(True)
            self._form.addRow("", hint)
        relax = QCheckBox("vertices on curves slide (2D)")
        relax.setChecked(params["relax"])
        self._form.addRow("relax", relax)
        self._editors["relax"] = ("flag", relax)
        unstructured = QLineEdit(", ".join(str(k) for k in params["unstructured"]))
        unstructured.setPlaceholderText("block indices meshed unstructured, e.g. 1, 3")
        self._form.addRow("unstructured blocks", unstructured)
        self._editors["unstructured"] = ("sketch_unstructured", unstructured)
        if params["unstructured"]:
            generator = params["generator"]
            if generator is None:
                generator = dict(DEFAULT_GENERATOR)
            sketch_connectors = [
                entity_id
                for _role, entity_id in params["outputs"]
                if entity_id in self.session.store
                and self.session.store.get(entity_id).kind == "connector"
            ]
            editor = UnstructuredEditor(
                generator,
                self.session.store,
                boundary_spacing=mean_connector_spacing(self.session, sketch_connectors),
            )
            editor.walls.setEnabled(False)
            editor.walls.setPlaceholderText("from the sketch's walls")
            self._form.addRow("unstructured generator", editor)
            self._editors["generator"] = ("sketch_generator", editor)

    def _group_rows(self, entity) -> None:
        key = "domain"
        if entity.kind == "patch":
            key = "connector"
            combo = QComboBox()
            combo.addItems(list(PATCH_TYPES))
            combo.setCurrentText(entity.params["patch_type"])
            self._form.addRow("patch type", combo)
            self._editors["patch_type"] = ("patch_type", combo)
        names = [member_label(item, self.session.store) for item in entity.params["members"]]
        members = QPlainTextEdit("\n".join(names))
        placeholder = f"one {key} name per line"
        if entity.kind == "patch":
            placeholder = "one connector, domain, block.face or volume.top / .bottom per line"
        members.setPlaceholderText(placeholder)
        self._form.addRow("members", members)
        self._editors["members"] = (f"members:{key}", members)
        if entity.kind == "patch" and entity.params["patch_type"] == "periodic":
            partners = []
            paired_from = []
            for patch in self.session.store.entities("patch"):
                if patch.id == entity.id or patch.params["patch_type"] != "periodic":
                    continue
                partners.append((patch.name, patch.id))
                pairing = patch.params["periodic"]
                if pairing is not None and pairing["partner"] == entity.id:
                    paired_from.append(patch.name)
            editor = PeriodicEditor(entity.params["periodic"], partners, paired_from)
            self._form.addRow("periodic", editor)
            self._editors["periodic"] = ("periodic", editor)

    def _domain_rows(self, entity) -> None:
        lines = []
        for k, edge in enumerate(entity.params["edges"]):
            names = []
            for item in edge:
                name = self.session.store.get(item["connector"]).name
                if item["reversed"]:
                    name += " (rev)"
                names.append(name)
            lines.append(f"{EDGE_NAMES[k]}: {', '.join(names)}")
        edges = QLabel("\n".join(lines))
        edges.setWordWrap(True)
        self._form.addRow("edges", edges)
        ni, nj = self.session.domain_dimensions(entity.name)
        self._form.addRow("dimensions", QLabel(f"{ni} x {nj}"))
        status = QLabel(self.session.domain_status(entity.name))
        status.setWordWrap(True)
        self._form.addRow("status", status)
        editor = GeneratorEditor(entity.params["generator"])
        self._form.addRow("generator", editor)
        self._editors["generator"] = ("generator", editor)
        buttons = QWidget()
        row = QHBoxLayout(buttons)
        row.setContentsMargins(0, 0, 0, 0)
        self.rotate_button = QPushButton("Rotate edges")
        self.rotate_button.clicked.connect(self._rotate)
        self.smooth_button = QPushButton("Run smoother")
        self.smooth_button.clicked.connect(lambda: self.smooth_requested.emit(entity.name))
        self.smooth_button.setEnabled(entity.params["generator"]["elliptic"] is not None)
        row.addWidget(self.rotate_button)
        row.addWidget(self.smooth_button)
        self._form.addRow("", buttons)

    def _rotate(self) -> None:
        name = self.session.store.get(self.entity_id).name
        try:
            self.session.rotate_domain(name)
        except ValidationError as error:
            self._show_error(str(error))

    def _connector_rows(self, entity) -> None:
        params = entity.params
        pieces = []
        for segment in params["segments"]:
            curve = self.session.store.get(segment["curve"]).name
            pieces.append(f"{curve} [{segment['t0']:g} → {segment['t1']:g}]")
        segments = QLabel("\n".join(pieces))
        segments.setWordWrap(True)
        self._form.addRow("segments", segments)
        nodes = self.session.connector_nodes(entity.name)
        for field, kind in CONNECTOR_FIELDS:
            context = None
            if kind == "distribution":
                context = (nodes.length, params["dimension"], "")
            editor = self._make_editor(kind, params[field], context)
            self._form.addRow(field.replace("_", " "), editor)
            self._editors[field] = (kind, editor)
        dimension_editor = self._editors["dimension"][1]
        distribution_editor = self._editors["distribution"][1]
        dimension_editor.valueChanged.connect(
            lambda n: distribution_editor.set_context(nodes.length, n)
        )
        self._boundary_rows(entity, "connector")
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

    def _boundary_rows(self, entity, key: str) -> None:
        """The patch the entity is in, and a combo to give it a boundary type (it then joins
        a patch of that type: Session.set_boundary_type)."""
        current = None
        for patch in self.session.store.entities("patch"):
            if {key: entity.id} in patch.params["members"]:
                current = patch
                break
        combo = QComboBox()
        combo.addItem("(none)")
        combo.addItems(list(PATCH_TYPES))
        label = QLabel("in no patch")
        if current is not None:
            combo.setCurrentText(current.params["patch_type"])
            label.setText(f"patch {current.name}")
        box = QWidget()
        row = QHBoxLayout(box)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(combo)
        row.addWidget(label, 1)
        self._form.addRow("boundary", box)
        self._editors["boundary"] = ("boundary", combo)
        self._boundary_before = combo.currentText()

    def _make_editor(self, kind: str, value, context=None) -> QWidget:
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
            if context is not None:
                length, dimension, unit = context
                return DistributionEditor(value, length=length, dimension=dimension, unit=unit)
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
            if kind == "boundary":
                continue  # applied as its own command (apply)
            value = self._read_editor(field, kind, editor)
            if value != entity.params[field]:
                changes[field] = value
        if entity.kind == "sketch" and changes:
            from opengrid.model.sketch import normalize_sketch, sketch_script_params

            merged = dict(sketch_script_params(entity.params, self.session.store))
            merged.update(changes)
            try:
                canonical = normalize_sketch(merged, self.session.store, entity.name)
            except ValidationError:
                return new_name, changes  # reported by apply
            changes = {k: v for k, v in changes.items() if canonical[k] != entity.params[k]}
        if entity.kind == "smoothing_group" and changes:
            # Names in the editors vs ids in the model: compare canonical parameters.
            from opengrid.model.groups import normalize_group

            merged = dict(entity.params)
            merged.update(changes)
            try:
                canonical = normalize_group(merged, self.session.store, entity.name)
            except ValidationError:
                return new_name, changes  # reported by apply
            changes = {k: v for k, v in changes.items() if canonical[k] != entity.params[k]}
        return new_name, changes

    def _read_editor(self, field: str, kind: str, editor: QWidget):
        if kind == "names":
            return tuple(line.strip() for line in editor.toPlainText().splitlines() if line.strip())
        if kind == "group_walls":
            key = "connector"
            if self.session.store.get(self.entity_id).params["blocks"]:
                key = "domain"
            return parse_group_walls(editor.toPlainText(), key)
        if kind == "slide":
            return parse_slide(editor.text())
        if kind in ("sketch_dimensions", "sketch_walls"):
            return parse_sketch_lines(
                editor.toPlainText(), kind, self.session.store.get(self.entity_id)
            )
        if kind == "sketch_unstructured":
            parts = [p for p in editor.text().replace(",", " ").split() if p]
            try:
                return tuple(sorted({int(p) for p in parts}))
            except ValueError:
                raise ParseError(
                    f"unstructured blocks: {editor.text()!r} is not a list of block indices"
                ) from None
        if kind == "sketch_generator":
            value = editor.value()
            value["inflation"] = None  # the sketch's walls give the layers
            return value
        if kind == "spin":
            return editor.value()
        if kind == "patch_type":
            return editor.currentText()
        if kind.startswith("members:"):
            key = kind.split(":", 1)[1]
            names = [line.strip() for line in editor.toPlainText().splitlines() if line.strip()]
            store = self.session.store
            members = []
            for name in names:
                if name in store:
                    entity = store.get(name)
                    member_key = key
                    if key == "connector" and entity.kind in ("domain", "unstructured_domain"):
                        member_key = "domain"
                    members.append({member_key: entity.id})
                    continue
                owner, _dot, face = name.rpartition(".")
                if key == "connector" and owner in store:
                    if store.get(owner).kind == "volume":
                        members.append({"volume": store.get(owner).id, "face": face})
                    else:
                        members.append({"block": store.get(owner).id, "face": face})
                    continue
                raise ParseError(f"members: no entity {name!r}")
            return tuple(members)
        if kind == "periodic":
            return self._periodic_value(editor)
        if kind in ("generator", "extrusion", "smoothing"):
            return editor.value()
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

    def _periodic_value(self, editor):
        """The pairing from the editor: the automatic translation filled in, and the stored
        pairing itself when nothing changed beyond the rounding of the editor's text."""
        from opengrid.model.periodic import automatic_translation

        value = editor.value()
        if value is None:
            return None
        entity = self.session.store.get(self.entity_id)
        if value["translation"] is None:
            if value["rotation_angle"] != 0.0:
                raise ParseError("periodic: give the translation together with a rotation")
            try:
                value["translation"] = automatic_translation(
                    self.session, entity, self.session.store.get(value["partner"])
                )
            except (ValueError, ValidationError) as error:
                raise ParseError(f"periodic: {error}") from None
        current = entity.params["periodic"]
        if current is None or current["partner"] != value["partner"]:
            return value
        same = math.isclose(current["rotation_angle"], value["rotation_angle"], abs_tol=1e-12)
        for key in ("rotation_axis", "rotation_center", "translation"):
            for a, b in zip(current[key], value[key], strict=True):
                same = same and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
        if same:
            return current
        return value

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
        boundary = None
        if "boundary" in self._editors:
            chosen = self._editors["boundary"][1].currentText()
            if chosen != self._boundary_before:
                boundary = chosen
        dimension = None
        generator = None
        entity = self.session.store.get(self.entity_id)
        if entity.kind == "unstructured_domain" and "generator" in changes:
            # Through set_unstructured_generator: a sketch's block sends it to the sketch.
            generator = changes.pop("generator")
        if entity.kind == "connector":
            from opengrid.model.sketch import output_owner

            owned = entity.id in output_owner(self.session.store)
            if "dimension" in changes or (owned and changes):
                # The points propagate through structured domains; a sketch's connector takes
                # its points and distribution from the sketch, so its edits go there.
                dimension = changes.pop("dimension", entity.params["dimension"])
        try:
            if generator is not None:
                self.session.set_unstructured_generator(name, generator)
                for warning in self.session.last_warnings:
                    self.message.emit(warning, False)
            if dimension is not None:
                # The connector's other edits go in the same step.
                followed = self.session.set_dimension(name, dimension, **changes)
                changes = {}
                for warning in self.session.last_warnings:
                    self.message.emit(warning, False)
                if len(followed) > 1:
                    self.message.emit(
                        f"{len(followed) - 1} opposite connector(s) followed: "
                        + ", ".join(n for n in followed if n != name),
                        False,
                    )
            if changes or new_name is not None:
                self.session.modify(name, new_name=new_name, **changes)  # one undo step
                if new_name is not None:
                    name = new_name
            if boundary is not None:
                self._apply_boundary(name, boundary)
        except ValidationError as error:
            self._show_error(str(error))
            return False
        if (
            not changes
            and new_name is None
            and boundary is None
            and dimension is None
            and generator is None
        ):
            self.message.emit("No changes", False)
            return False
        self.message.emit(f"Updated {self.session.store.get(self.entity_id).name}", False)
        return True

    def _apply_boundary(self, name: str, chosen: str) -> None:
        if chosen == "(none)":
            entity_id = self.session.store.get(name).id
            for patch in self.session.store.entities("patch"):
                if any(entity_id in item.values() for item in patch.params["members"]):
                    self.session.remove_from_patch(patch.name, [name])
            return
        self.session.set_boundary_type([name], chosen)

    def _show_error(self, text: str) -> None:
        self._error.setText(text)
        self.message.emit(text, True)


def parse_sketch_lines(text: str, kind: str, entity) -> tuple:
    """Sketch dimensions ("a b points") or walls ("a b [spacing [layers [growth]]]", 3D
    "a b c d [...]"; spacing "-": the initial height)."""
    size = 2
    key = "edge"
    if len(entity.params["blocks"][0]) == 8:
        size = 4
        key = "face"
    out = []
    for number, line in enumerate(text.splitlines(), start=1):
        parts = line.split()
        if not parts:
            continue
        try:
            if kind == "sketch_dimensions":
                if len(parts) != 3:
                    raise ValueError
                out.append({"edge": (int(parts[0]), int(parts[1])), "points": int(parts[2])})
                continue
            if not size <= len(parts) <= size + 3:
                raise ValueError
            ref = tuple(int(v) for v in parts[:size])
            wall = {key: ref, "spacing": None}
            if len(parts) > size and parts[size] != "-":
                wall["spacing"] = float(parts[size])
            if len(parts) > size + 1:
                wall["layers"] = int(parts[size + 1])
            if len(parts) > size + 2:
                wall["growth"] = float(parts[size + 2])
            out.append(wall)
        except ValueError:
            raise ParseError(f"line {number}: cannot read {line!r}") from None
    return tuple(out)


def _canonical_vec(values) -> tuple:
    coords = [float(v) for v in values]
    if len(coords) == 2:
        coords.append(0.0)
    return tuple(coords)
