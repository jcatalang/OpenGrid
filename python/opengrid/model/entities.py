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

SURFACE_KINDS = ("plane", "ruled_surface", "extruded_surface", "revolved_surface", "nurbs_surface")
DATABASE_KINDS = ("point", "line", "arc", "polyline", "spline", "nurbs", *SURFACE_KINDS)
DOMAIN_KINDS = ("domain", "unstructured_domain")
KINDS = (*DATABASE_KINDS, "connector", *DOMAIN_KINDS, "block", "volume", "patch", "zone")
# Surface parameters that reference curves (entity ids once canonical).
SURFACE_CURVES = {
    "ruled_surface": ("curve0", "curve1"),
    "extruded_surface": ("curve",),
    "revolved_surface": ("curve",),
}
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


def _plane(c: _Checker, raw: dict) -> dict:
    return {
        "origin": c.vec3(raw["origin"], "origin"),
        "u_vector": c.vec3(raw["u_vector"], "u_vector"),
        "v_vector": c.vec3(raw["v_vector"], "v_vector"),
    }


def _nurbs_surface(c: _Checker, raw: dict) -> dict:
    rows = raw["control_points"]
    control = None
    if not isinstance(rows, (list, tuple)) or not rows:
        c.fail("control_points", "expected a list of rows (one per u index) of points")
    else:
        control = []
        for k, row in enumerate(rows):
            converted = c.points(row, f"control_points[{k}]", 2)
            control.append(converted)
        if any(r is not None and len(r) != len(control[0] or ()) for r in control):
            c.fail("control_points", "all rows must have the same number of points")
        control = tuple(control)
    weights = raw.get("weights")
    if weights is None and control is not None:
        weights = [[1.0] * len(row or ()) for row in control]
    weight_rows = None
    if weights is not None:
        weight_rows = tuple(c.numbers(row, "weights") for row in weights)
    return {
        "degree_u": c.integer(raw["degree_u"], "degree_u"),
        "degree_v": c.integer(raw["degree_v"], "degree_v"),
        "knots_u": c.numbers(raw["knots_u"], "knots_u"),
        "knots_v": c.numbers(raw["knots_v"], "knots_v"),
        "control_points": control,
        "weights": weight_rows,
    }


def _curve_ref(c: _Checker, raw, field: str, store) -> int | None:
    from opengrid.model.connectors import CURVE_KINDS

    if isinstance(raw, bool) or not isinstance(raw, (str, int)) or raw not in store:
        c.fail(field, f"no curve {raw!r}")
        return None
    entity = store.get(raw)
    if entity.kind not in CURVE_KINDS:
        c.fail(field, f"{entity.name!r} is a {entity.kind}, not a curve")
        return None
    return entity.id


def _ruled(c: _Checker, raw: dict, store) -> dict:
    return {
        "curve0": _curve_ref(c, raw["curve0"], "curve0", store),
        "curve1": _curve_ref(c, raw["curve1"], "curve1", store),
    }


def _extruded(c: _Checker, raw: dict, store) -> dict:
    return {
        "curve": _curve_ref(c, raw["curve"], "curve", store),
        "direction": c.vec3(raw["direction"], "direction"),
    }


def _revolved(c: _Checker, raw: dict, store) -> dict:
    return {
        "curve": _curve_ref(c, raw["curve"], "curve", store),
        "axis_point": c.vec3(raw.get("axis_point", (0.0, 0.0, 0.0)), "axis_point"),
        "axis": c.vec3(raw.get("axis", (0.0, 0.0, 1.0)), "axis"),
        "start_angle": c.number(_as_python_number(raw.get("start_angle", 0.0)), "start_angle"),
        "end_angle": c.number(_as_python_number(raw["end_angle"]), "end_angle"),
    }


# kind: (required, optional, converter taking the store)
_SURFACE_SCHEMAS = {
    "ruled_surface": (("curve0", "curve1"), (), _ruled),
    "extruded_surface": (("curve", "direction"), (), _extruded),
    "revolved_surface": (("curve", "end_angle"), ("axis_point", "axis", "start_angle"), _revolved),
}


_SCHEMAS: dict[str, tuple[tuple[str, ...], tuple[str, ...], Callable[[_Checker, dict], dict]]] = {
    # kind: (required keys, optional keys, converter)
    "point": (("xyz",), (), _point),
    "line": (("p0", "p1"), (), _line),
    "arc": (("center", "radius", "start_angle", "end_angle"), ("normal", "x_axis"), _arc),
    "polyline": (("points",), (), _polyline),
    "spline": (("points",), ("method", "start_tangent", "end_tangent"), _spline),
    "nurbs": (("degree", "knots", "control_points"), ("weights",), _nurbs),
    "plane": (("origin", "u_vector", "v_vector"), (), _plane),
    "nurbs_surface": (
        ("degree_u", "degree_v", "knots_u", "knots_v", "control_points"),
        ("weights",),
        _nurbs_surface,
    ),
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
    if kind == "unstructured_domain":
        from opengrid.model.unstructured import normalize_unstructured

        if store is None:
            raise ValueError("unstructured domain parameters need the store")
        return normalize_unstructured(raw, store, where)
    if kind == "volume":
        from opengrid.model.volumes import normalize_volume

        if store is None:
            raise ValueError("volume parameters need the store")
        return normalize_volume(raw, store, where)
    if kind == "block":
        from opengrid.model.blocks import normalize_block

        if store is None:
            raise ValueError("block parameters need the store")
        return normalize_block(raw, store, where)
    if kind in _SURFACE_SCHEMAS:
        if store is None:
            raise ValueError(f"{kind} parameters need the store to resolve curves")
        required, optional, converter = _SURFACE_SCHEMAS[kind]
        checker = _Checker(where)
        for key in required:
            if key not in raw:
                checker.fail(key, "missing")
        for key in raw:
            if key not in required and key not in optional:
                checker.fail(key, f"unexpected parameter for a {kind}")
        if checker.issues:
            raise ValidationError(checker.issues)
        params = converter(checker, raw, store)
        if checker.issues:
            raise ValidationError(checker.issues)
        try:
            build_surface(kind, params, store.geometry)
        except ValueError as error:
            raise ValidationError([Issue(where, str(error))]) from None
        return params
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
    if kind in SURFACE_KINDS:
        try:
            build_surface(kind, params, None)
        except (ValueError, TypeError, IndexError) as error:
            raise ValidationError([Issue(where, str(error))]) from None
        return params
    if kind != "point":
        try:
            build_geometry(kind, params)
        except ValueError as error:
            raise ValidationError([Issue(where, str(error))]) from None
    return params


def build_surface(kind: str, params: dict, curve_of):
    """C++ surface for canonical parameters; ``curve_of(id)`` returns a referenced curve."""
    if kind == "plane":
        return _core.PlaneSurface(params["origin"], params["u_vector"], params["v_vector"])
    if kind == "ruled_surface":
        return _core.RuledSurface(curve_of(params["curve0"]), curve_of(params["curve1"]))
    if kind == "extruded_surface":
        return _core.ExtrudedSurface(curve_of(params["curve"]), params["direction"])
    if kind == "revolved_surface":
        return _core.RevolvedSurface(
            curve_of(params["curve"]),
            params["axis_point"],
            params["axis"],
            params["start_angle"],
            params["end_angle"],
        )
    if kind == "nurbs_surface":
        rows = params["control_points"]
        nu = len(rows)
        nv = len(rows[0])
        control = [p for row in rows for p in row]
        weights = [w for row in params["weights"] for w in row]
        return _core.NurbsSurface(
            params["degree_u"],
            params["degree_v"],
            list(params["knots_u"]),
            list(params["knots_v"]),
            nu,
            nv,
            control,
            weights,
        )
    raise ValueError(f"{kind} is not a surface")


def build_geometry(kind: str, params: dict):
    """Construct the C++ curve for canonical parameters (None for a point)."""
    if kind == "point":
        return None
    if kind in ("connector", "domain", "unstructured_domain", "block", "volume", "patch", "zone"):
        raise ValueError(f"{kind}s have no single curve")
    if kind in SURFACE_KINDS:
        raise ValueError(f"{kind}: use build_surface (surfaces reference curves)")
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
