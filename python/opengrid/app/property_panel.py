"""Property panel: view and edit the parameters of the selected entity.

Edits are applied through ``Session.modify`` / ``Session.rename``, i.e. as undoable commands.
Input that cannot be parsed or that the model rejects is reported in the panel; nothing is
changed in that case.
"""

import math

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
from opengrid.model.blocks import EXTRUSION_MODES
from opengrid.model.connectors import CURVE_KINDS, DISTRIBUTION_TYPES, SPACING_MODES
from opengrid.model.domains import (
    CONTROLS,
    DEFAULT_ELLIPTIC,
    DEFAULT_WALL,
    EDGE_NAMES,
    METHODS,
    SOLVERS,
)
from opengrid.model.entities import SPLINE_METHODS
from opengrid.model.patches import PATCH_TYPES, member_label
from opengrid.model.unstructured import ALGORITHMS, ELEMENT_TYPES, SMOOTHING_METHODS
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
        self.distribution = DistributionEditor(value["distribution"])
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


class VolumeEditor(QWidget):
    """Volume generator: backend, max size (empty: none), grading, pyramids, quality
    optimisation, and boundary layers (wall faces; empty: none)."""

    def __init__(self, value: dict, store=None, parent=None):
        super().__init__(parent)
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
        layout.addRow("first height", self.first_height)
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
        }

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

    def __init__(self, value: dict, store=None, parent=None):
        super().__init__(parent)
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
        layout.addRow("first height", self.first_height)
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
        if entity.kind in ("patch", "zone"):
            self._group_rows(entity)
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
        editor = UnstructuredEditor(entity.params["generator"], self.session.store)
        self._form.addRow("generator", editor)
        self._editors["generator"] = ("generator", editor)

    def _volume_rows(self, entity) -> None:
        try:
            text = self.session.volume_status(entity.name)
        except ValidationError as error:
            text = str(error)
        status = QLabel(text)
        status.setWordWrap(True)
        self._form.addRow("status", status)
        faces = QLabel(", ".join(self.session.store.get(f).name for f in entity.params["faces"]))
        faces.setWordWrap(True)
        self._form.addRow("faces", faces)
        editor = VolumeEditor(entity.params["generator"], self.session.store)
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
            placeholder = "one connector, domain or block.face per line"
        members.setPlaceholderText(placeholder)
        self._form.addRow("members", members)
        self._editors["members"] = (f"members:{key}", members)

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
                    if key == "connector" and entity.kind == "domain":
                        member_key = "domain"
                    members.append({member_key: entity.id})
                    continue
                block, _dot, face = name.rpartition(".")
                if key == "connector" and block in store:
                    members.append({"block": store.get(block).id, "face": face})
                    continue
                raise ParseError(f"members: no entity {name!r}")
            return tuple(members)
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
