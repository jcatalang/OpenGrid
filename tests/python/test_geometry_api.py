import math
from itertools import pairwise

import numpy as np
import pytest

import opengrid
from opengrid import _core


def test_circle_length_and_projection():
    og = opengrid.Session()
    name = og.create_circle((0.5, 0.0, 0.0), 20.0)
    assert abs(og.length(name) - 40.0 * math.pi) / (40.0 * math.pi) < 1e-12
    t, point, distance = og.project_point(name, (0.5 + 25.0 * math.cos(1.0), 25.0 * math.sin(1.0)))
    assert abs(t - 1.0 / math.tau) < 1e-12
    assert abs(distance - 5.0) < 1e-11
    assert abs(math.hypot(point[0] - 0.5, point[1]) - 20.0) < 1e-12


def test_arc_through_three_points():
    og = opengrid.Session()
    name = og.create_arc_3pt((2.0, 1.0), (1.0, 2.0), (0.0, 1.0))
    params = og.entity(name).params
    assert params["center"] == pytest.approx((1.0, 1.0, 0.0), abs=1e-15)
    assert params["radius"] == pytest.approx(1.0, abs=1e-15)
    assert params["end_angle"] == pytest.approx(math.pi, abs=1e-15)
    curve = og.geometry(name)
    assert curve.eval(0.5) == pytest.approx([1.0, 2.0, 0.0], abs=1e-15)
    # Clockwise order gives the opposite normal.
    cw = og.create_arc_3pt((0.0, 1.0), (1.0, 2.0), (2.0, 1.0))
    assert og.entity(cw).params["normal"][2] == pytest.approx(-1.0)


def test_arc_3pt_collinear_rejected():
    og = opengrid.Session()
    with pytest.raises(opengrid.ValidationError, match="collinear"):
        og.create_arc_3pt((0, 0), (1, 1), (2, 2))


def test_spline_arc_length_vs_reference():
    og = opengrid.Session()
    points = [(0.0, 0.0), (1.0, 0.8), (2.5, -0.3), (3.0, 1.2), (4.0, 0.0)]
    name = og.create_spline(points, method="natural")
    curve = og.geometry(name)
    # Independent reference: Richardson-extrapolated chord sums per segment (uses eval only).
    knots = curve.knots
    total = 0.0
    for a, b in pairwise(knots):
        sums = []
        for n in (4096, 8192):
            xyz = curve.eval_many(np.linspace(a, b, n + 1))
            sums.append(np.sum(np.linalg.norm(np.diff(xyz, axis=0), axis=1)))
        total += (4.0 * sums[1] - sums[0]) / 3.0
    assert abs(curve.length() - total) / total < 1e-10


def test_arc_length_reparametrization_uniform_spacing():
    og = opengrid.Session()
    name = og.create_spline([(0, 0), (1, 1), (3, 0), (4, 2)], method="catmull_rom")
    curve = og.geometry(name)
    total = curve.length()
    for k in range(11):
        s = total * k / 10
        t = curve.t_at_length(s)
        assert abs(curve.length_between(0.0, t) - s) / total < 1e-10


def test_projection_points_on_and_off_spline():
    og = opengrid.Session()
    curve = og.geometry(og.create_spline([(0, 0), (1, 1), (2, 0)]))
    for t in np.linspace(0.05, 0.95, 19):
        on = curve.eval(float(t))
        t_found, _point, distance = curve.project(on)
        assert distance < 1e-12 * curve.length()
        assert abs(t_found - t) < 1e-9


@pytest.mark.parametrize(
    ("kind", "kwargs", "message"),
    [
        ("line", {"p0": (0, 0), "p1": (0, 0)}, "coincide"),
        ("polyline", {"points": [(0, 0), (1, 0), (1, 0)]}, r"points\[2\]"),
        ("spline", {"points": [(0, 0), (1, 1)], "method": "akima"}, "method"),
        ("spline", {"points": [(0, 0), (1, 1)], "method": "clamped"}, "tangent"),
        ("arc", {"center": (0, 0), "radius": -1, "start_angle": 0, "end_angle": 1}, "radius"),
        (
            "nurbs",
            {"degree": 2, "knots": [0, 1], "control_points": [(0, 0), (1, 1), (2, 0)]},
            "knots",
        ),
        ("point", {"xyz": (0, "a")}, r"xyz\[1\]"),
    ],
)
def test_invalid_parameters_are_reported(kind, kwargs, message):
    og = opengrid.Session()
    create = getattr(og, f"create_{kind}")
    with pytest.raises(opengrid.ValidationError, match=message):
        create(**kwargs)
    assert len(og.store) == 0
    assert not og.commands.can_undo


def test_core_classes_exposed():
    assert issubclass(_core.Line, _core.Curve)
    assert issubclass(_core.Nurbs, _core.Curve)
