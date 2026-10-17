"""Closed curves stay closed when their seam is dragged, a loop that does not quite close
says so, and the region between an outer and an inner loop meshes as one domain with a hole
(the user could not mesh a freehand loop after moving its control points, and asked how to
mesh between two curves)."""

import math

import numpy as np
import pytest

import opengrid
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def wavy_loop(og, close_tolerance=0.05, turns=1.0):
    t = np.linspace(0, 2 * math.pi * turns, 121)
    r = 1.0 + 0.15 * np.cos(3 * t)
    pts = [(float(r[k] * math.cos(t[k])), float(r[k] * math.sin(t[k])), 0.0) for k in range(121)]
    return og.create_freehand(pts, 0.005, 60.0, close_tolerance=close_tolerance)


def test_seam_of_a_closed_curve_moves_as_one_point():
    og = opengrid.Session()
    (name,) = wavy_loop(og)
    control = og.entity(name).params["control_points"]
    assert control[0] == control[-1]
    last = len(control) - 1
    og.move_handle(name, "control_points", 0, (control[0][0] + 0.1, control[0][1] + 0.1, 0.0))
    moved = og.entity(name).params["control_points"]
    assert moved[0] == moved[-1] == (control[0][0] + 0.1, control[0][1] + 0.1, 0.0)
    og.move_handle(name, "control_points", last, (control[0][0], control[0][1] - 0.2, 0.0))
    moved = og.entity(name).params["control_points"]
    assert moved[0] == moved[-1] == (control[0][0], control[0][1] - 0.2, 0.0)
    assert moved[1] == control[1]
    domain = og.unstructured_domain_from_curves([name], 41)  # still a closed loop
    assert len(og.unstructured_mesh(domain)["triangles"]) > 50
    # An open curve: only the dragged end moves.
    og.create_polyline([(0, 0), (1, 0), (2, 0), (3, 0)], name="open")
    og.move_handle("open", "points", 0, (0, 1, 0))
    assert og.entity("open").params["points"][-1] == (3.0, 0.0, 0.0)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_a_loop_that_does_not_quite_close_says_so():
    og = opengrid.Session()
    (name,) = wavy_loop(og, close_tolerance=0.0, turns=0.98)  # released short of the start
    with pytest.raises(ValidationError, match=r"its two ends are .* apart: close the loop"):
        og.unstructured_domain_from_curves([name], 41)


def test_region_between_an_outer_and_an_inner_loop():
    og = opengrid.Session()
    (inner,) = wavy_loop(og)
    og.create_line((-3, -3), (3, -3), name="b")
    og.create_line((3, -3), (3, 3), name="r")
    og.create_line((3, 3), (-3, 3), name="t")
    og.create_line((-3, 3), (-3, -3), name="l")
    domain = og.unstructured_domain_from_curves([inner, "b", "r", "t", "l"], 31)
    loops = og.entity(domain).params["loops"]
    assert len(loops) == 2  # the square outside, the loop a hole
    outer = {og.store.get(i["connector"]).name for i in loops[0]}
    assert len(outer) == 4
    mesh = og.unstructured_mesh(domain)
    points = np.asarray(mesh["points"])
    inside = np.hypot(points[:, 0], points[:, 1]) < 0.7  # the hole's heart is empty
    assert not inside.any()
    assert len(mesh["triangles"]) > 500
