"""Database (geometry) entity kinds: parameter schemas, validation and geometry construction.

Parameters are kept in a canonical form: floats, 3-component tuples for vectors, tuples of
those for point lists. Canonicalisation is idempotent, which makes project save/load/save
byte-identical. Values are stored as given (e.g. an arc normal is not normalised); the C++
constructors normalise internally.
"""

import math
from collections.abc import Callable
from typing import Any

from opengrid import _core
from opengrid.validation import Issue, ValidationError

DATABASE_KINDS = ("point", "line", "arc", "polyline", "spline", "nurbs")
KINDS = (*DATABASE_KINDS, "connector", "domain", "patch", "zone")
SPLINE_METHODS = ("natural", "clamped", "catmull_rom")

_SPLINE_ENUM = {
    "natural": _core.SplineMethod.NATURAL,
    "clamped": _core.SplineMethod.CLAMPED,
    "catmull_rom": _core.SplineMethod.CATMULL_ROM,
}


class _Checker:
    """Collects issues while converting raw parameter values."""

    def __init__(self, where: str):
        self.where = where
        self.issues: list[Issue] = []

    def fail(self, field: str, message: str) -> None:
        self.issues.append(Issue(f"{self.where}: {field}", message))

    def number(self, value: Any, field: str) -> float | None:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            self.fail(field, f"expected a number, got {value!r}")
            return None
        result = float(value)
        if not math.isfinite(result):
            self.fail(field, f"must be finite, got {result!r}")
            return None
        return result

    def integer(self, value: Any, field: str) -> int | None:
        if isinstance(value, bool) or not isinstance(value, int):
            self.fail(field, f"expected an integer, got {value!r}")
            return None
        return value

    def vec3(self, value: Any, field: str) -> tuple[float, float, float] | None:
        """A 3-vector; 2-vectors are accepted and get z = 0."""
        if isinstance(value, (str, bytes)) or not hasattr(value, "__len__"):
            self.fail(field, f"expected 2 or 3 coordinates, got {value!r}")
            return None
        if len(value) not in (2, 3):
            self.fail(field, f"expected 2 or 3 coordinates, got {len(value)}")
            return None
        coords = []
        for index, component in enumerate(value):
            number = self.number(_as_python_number(component), f"{field}[{index}]")
            if number is None:
                return None
            coords.append(number)
        if len(coords) == 2:
            coords.append(0.0)
        return (coords[0], coords[1], coords[2])

    def points(self, value: Any, field: str, minimum: int) -> tuple | None:
        if isinstance(value, (str, bytes)) or not hasattr(value, "__len__"):
            self.fail(field, f"expected a list of points, got {value!r}")
            return None
        result = []
        for index, item in enumerate(value):
            point = self.vec3(item, f"{field}[{index}]")
            if point is None:
                return None
            result.append(point)
        if len(result) < minimum:
            self.fail(field, f"at least {minimum} points are required, got {len(result)}")
            return None
        for index in range(1, len(result)):
            if result[index] == result[index - 1]:
                self.fail(f"{field}[{index}]", f"coincides with {field}[{index - 1}]")
                return None
        return tuple(result)

    def numbers(self, value: Any, field: str) -> tuple | None:
        if isinstance(value, (str, bytes)) or not hasattr(value, "__len__"):
            self.fail(field, f"expected a list of numbers, got {value!r}")
            return None
        result = []
        for index, item in enumerate(value):
            number = self.number(_as_python_number(item), f"{field}[{index}]")
            if number is None:
                return None
            result.append(number)
        return tuple(result)


def _as_python_number(value: Any) -> Any:
    """Convert NumPy scalars to Python numbers so type checks treat them uniformly."""
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        return value.item()
    return value


def _point(c: _Checker, raw: dict) -> dict:
    return {"xyz": c.vec3(raw["xyz"], "xyz")}


def _line(c: _Checker, raw: dict) -> dict:
    return {"p0": c.vec3(raw["p0"], "p0"), "p1": c.vec3(raw["p1"], "p1")}


def _arc(c: _Checker, raw: dict) -> dict:
    return {
        "center": c.vec3(raw["center"], "center"),
        "radius": c.number(_as_python_number(raw["radius"]), "radius"),
        "normal": c.vec3(raw.get("normal", (0.0, 0.0, 1.0)), "normal"),
        "x_axis": c.vec3(raw.get("x_axis", (1.0, 0.0, 0.0)), "x_axis"),
        "start_angle": c.number(_as_python_number(raw["start_angle"]), "start_angle"),
        "end_angle": c.number(_as_python_number(raw["end_angle"]), "end_angle"),
    }


def _polyline(c: _Checker, raw: dict) -> dict:
    return {"points": c.points(raw["points"], "points", 2)}


def _spline(c: _Checker, raw: dict) -> dict:
    method = raw.get("method", "natural")
    if method not in SPLINE_METHODS:
        c.fail("method", f"must be one of {', '.join(SPLINE_METHODS)}, got {method!r}")
    start = raw.get("start_tangent")
    end = raw.get("end_tangent")
    if start is not None:
        start = c.vec3(start, "start_tangent")
    if end is not None:
        end = c.vec3(end, "end_tangent")
    return {
        "points": c.points(raw["points"], "points", 2),
        "method": method,
        "start_tangent": start,
        "end_tangent": end,
    }


def _nurbs(c: _Checker, raw: dict) -> dict:
    control_points = c.points(raw["control_points"], "control_points", 2)
    weights = raw.get("weights")
    if weights is None and control_points is not None:
        weights = [1.0] * len(control_points)
    return {
        "degree": c.integer(raw["degree"], "degree"),
        "knots": c.numbers(raw["knots"], "knots"),
        "control_points": control_points,
        "weights": c.numbers(weights, "weights"),
    }


_SCHEMAS: dict[str, tuple[tuple[str, ...], tuple[str, ...], Callable[[_Checker, dict], dict]]] = {
    # kind: (required keys, optional keys, converter)
    "point": (("xyz",), (), _point),
    "line": (("p0", "p1"), (), _line),
    "arc": (("center", "radius", "start_angle", "end_angle"), ("normal", "x_axis"), _arc),
    "polyline": (("points",), (), _polyline),
    "spline": (("points",), ("method", "start_tangent", "end_tangent"), _spline),
    "nurbs": (("degree", "knots", "control_points"), ("weights",), _nurbs),
}


_CONNECTOR_KEYS = (
    ("segments", "dimension"),
    ("mode", "distribution", "breakpoints", "propagate_spacing"),
)


def normalize_params(kind: str, raw: dict, where: str, store=None) -> dict:
    """Validate raw parameters of an entity and return them in canonical form.

    Raises ValidationError listing every problem (unknown kind, missing or unexpected keys,
    malformed values, and geometric errors reported by the C++ constructors). Connectors
    reference other entities and need the ``store``.
    """
    if kind in ("patch", "zone"):
        from opengrid.model.patches import normalize_patch, normalize_zone

        if store is None:
            raise ValueError(f"{kind} parameters need the store")
        if kind == "patch":
            return normalize_patch(raw, store, where)
        return normalize_zone(raw, store, where)
    if kind == "domain":
        from opengrid.model.domains import normalize_domain

        if store is None:
            raise ValueError("domain parameters need the store to resolve connectors")
        return normalize_domain(raw, store, where)
    if kind == "connector":
        from opengrid.model.connectors import normalize_connector

        checker = _Checker(where)
        required, optional = _CONNECTOR_KEYS
        for key in required:
            if key not in raw:
                checker.fail(key, "missing")
        for key in raw:
            if key not in required and key not in optional:
                checker.fail(key, "unexpected parameter for a connector")
        if checker.issues:
            raise ValidationError(checker.issues)
        if store is None:
            raise ValueError("connector parameters need the store to resolve curves")
        return normalize_connector(raw, store, where)
    if kind not in _SCHEMAS:
        raise ValidationError([Issue(where, f"unknown entity type {kind!r}")])
    required, optional, converter = _SCHEMAS[kind]
    checker = _Checker(where)
    for key in required:
        if key not in raw:
            checker.fail(key, "missing")
    for key in raw:
        if key not in required and key not in optional:
            checker.fail(key, f"unexpected parameter for a {kind}")
    if checker.issues:
        raise ValidationError(checker.issues)
    params = converter(checker, raw)
    if checker.issues:
        raise ValidationError(checker.issues)
    if kind != "point":
        try:
            build_geometry(kind, params)
        except ValueError as error:
            raise ValidationError([Issue(where, str(error))]) from None
    return params


def build_geometry(kind: str, params: dict):
    """Construct the C++ curve for canonical parameters (None for a point)."""
    if kind == "point":
        return None
    if kind in ("connector", "domain", "patch", "zone"):
        raise ValueError(f"{kind}s have no single curve")
    if kind == "line":
        return _core.Line(params["p0"], params["p1"])
    if kind == "arc":
        return _core.Arc(
            params["center"],
            params["radius"],
            params["normal"],
            params["x_axis"],
            params["start_angle"],
            params["end_angle"],
        )
    if kind == "polyline":
        return _core.Polyline(list(params["points"]))
    if kind == "spline":
        return _core.CubicSpline(
            list(params["points"]),
            _SPLINE_ENUM[params["method"]],
            params["start_tangent"],
            params["end_tangent"],
        )
    if kind == "nurbs":
        return _core.Nurbs(
            params["degree"],
            list(params["knots"]),
            list(params["control_points"]),
            list(params["weights"]),
        )
    raise ValueError(f"unknown entity type {kind!r}")
