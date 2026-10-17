"""Phase 10c: live mouse editing (drag handles), freehand curves, domains from curve networks."""

import math
import time

import numpy as np
import pytest

import opengrid
from opengrid.model.fitting import find_corners, fit_stroke
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def replays(og) -> bool:
    return to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def distance_to(pieces, points) -> float:
    curves = [
        opengrid._core.Nurbs(
            p["degree"], list(p["knots"]), list(p["control_points"]), list(p["weights"])
        )
        for p in pieces
    ]
    return max(min(c.project(tuple(q))[2] for c in curves) for q in points)


def stroke_length(points) -> float:
    return float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))


@pytest.mark.parametrize("case", ["quarter_circle", "naca", "noisy_s"])
def test_freehand_fit_within_tolerance(case):
    rng = np.random.default_rng(0)
    if case == "quarter_circle":
        t = np.linspace(0, math.pi / 2, 200)
        points = np.stack([np.cos(t), np.sin(t), 0 * t], axis=1)
    elif case == "naca":
        x = np.linspace(0, 1, 300)
        y = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
        points = np.stack([x, y, 0 * x], axis=1)
    else:
        t = np.linspace(0, 1, 400)
        points = np.stack([t, 0.3 * np.sin(2 * np.pi * t), 0 * t], axis=1)
        points[:, :2] += 0.0005 * rng.standard_normal((400, 2))
    pieces = fit_stroke(points, tolerance=0.005)
    assert len(pieces) == 1
    assert distance_to(pieces, points) <= 0.005 * stroke_length(points)
    if case == "quarter_circle":
        assert len(pieces[0]["control_points"]) <= 12
    assert pieces[0]["control_points"][0] == tuple(points[0])  # end points interpolated


def test_freehand_corners_split_the_stroke():
    a = np.stack([np.linspace(0, 1, 60), np.zeros(60), np.zeros(60)], axis=1)
    b = np.stack([np.ones(60), np.linspace(0, 1, 60), np.zeros(60)], axis=1)[1:]
    c = np.stack([np.linspace(1, 0, 60), np.ones(60), np.zeros(60)], axis=1)[1:]
    gentle = np.stack(
        [-0.5 * np.linspace(0, 1, 60), 1 + 0.5 * np.linspace(0, 1, 60), np.zeros(60)], axis=1
    )[1:]
    stroke = np.concatenate([a, b, c, gentle])
    assert len(find_corners(stroke)) == 2  # two 90 degree turns; the 45 degree turn is no corner
    pieces = fit_stroke(stroke)
    assert len(pieces) == 3
    assert distance_to(pieces, stroke) <= 0.005 * stroke_length(stroke)
    assert len(fit_stroke(stroke, corner_angle=None)) == 1


def test_create_freehand_is_one_undo_step_and_replays():
    og = opengrid.Session()
    t = np.linspace(0, math.pi, 100)
    names = og.create_freehand(np.stack([np.cos(t), np.sin(t), 0 * t], axis=1), name="arcish")
    assert names == ["arcish"]
    assert og.entity("arcish").kind == "nurbs"
    assert replays(og)
    og.undo()
    assert og.names() == []


def lattice(og, columns=(0.0, 1.0, 2.0, 3.0), rows=(0.0, 1.0, 2.0)):
    names = []
    for k, x in enumerate(columns):
        y = np.linspace(-0.2, rows[-1] + 0.2, 11)
        pts = np.stack([x + 0.05 * np.sin(2 * y + k), y, 0 * y], axis=1)
        names.append(og.create_spline([tuple(p) for p in pts], name=f"v{k}"))
    for k, y0 in enumerate(rows):
        x = np.linspace(-0.2, columns[-1] + 0.2, 13)
        pts = np.stack([x, y0 + 0.05 * np.cos(2 * x + k), 0 * x], axis=1)
        names.append(og.create_spline([tuple(p) for p in pts], name=f"h{k}"))
    return names


def test_curve_lattice_gives_conforming_structured_domains(tmp_path):
    og = opengrid.Session()
    curves = lattice(og)
    domains = og.domains_from_curves(curves, spacing=0.05)
    assert len(domains) == 6 and og.last_warnings == []
    assert len(og.topology().interfaces()) == 7
    for d in domains:
        assert og.domain_quality(d)["scaled_jacobian"].min() > 0.9
    path = tmp_path / "lattice.cgns"
    og.create_patch("all", "wall", og.unassigned_boundary_connectors())
    og.export_cgns(path)
    assert path.stat().st_size > 0
    assert replays(og)
    og.undo()  # the patch
    og.undo()  # the network: one step
    assert og.names("domain") == []


def test_curve_network_reports_faces_without_four_corners():
    og = opengrid.Session()
    for name, a, b in (
        ("s0", (0, 0), (2, 0)),
        ("s1", (2, 0), (2, 2)),
        ("s2", (2, 2), (0, 2)),
        ("s3", (0, 2), (0, 0)),
        ("cut", (1.5, -0.2), (2.2, 0.5)),  # cuts the corner (2, 0): a triangle and a pentagon
    ):
        og.create_line(a, b, name=name)
    curves = ["s0", "s1", "s2", "s3", "cut"]
    with pytest.raises(ValidationError, match="5 corner"):
        og.domains_from_curves(curves, spacing=0.1)
    created = og.domains_from_curves(curves, spacing=0.1, unstructured=True)
    assert len(created) == 2
    assert all(og.entity(n).kind == "unstructured_domain" for n in created)


def test_single_closed_loop_of_four_curves():
    og = opengrid.Session()
    og.create_spline([(0, 0), (1, -0.1), (2, 0)], name="b")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_spline([(2, 1), (1, 1.2), (0, 1)], name="t")
    og.create_line((0, 1), (0, 0), name="l")
    (domain,) = og.domains_from_curves(["b", "r", "t", "l"], spacing=0.05)
    ni, nj = og.domain_dimensions(domain)
    assert ni > 30 and nj > 15


def test_gui_drag_freehand_nurbs_and_network(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    og.create_spline([(0, 0), (1, 0.1), (2, 0)], name="b")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_line((2, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    og.create_connector("b", 200, name="cb")
    og.create_connector("t", 200, name="ct")
    og.create_connector("r", 100, name="cr")
    og.create_connector("l", 100, name="cl")
    og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="d")
    journal = len(og.journal_text().splitlines())
    # Drag the middle spline point: live preview, one command on release.
    window.set_selection([og.entity("b").id])
    assert len(window.visible_handles()) == 3
    assert window.handle_press((1, 0.1, 0))
    before = og.domain_grid("d").copy()
    start = time.perf_counter()
    for k in range(10):
        window.handle_drag((1, 0.1 + 0.02 * k, 0))
    per_frame = (time.perf_counter() - start) / 10
    assert per_frame < 0.1  # >= 10 frames per second
    assert np.array_equal(og.domain_grid("d"), before)  # previews are not committed
    assert window.handle_release((1, 0.3, 0)) == ["b"]
    assert og.entity("b").params["points"][1] == (1.0, 0.3, 0.0)
    lines = og.journal_text().splitlines()
    assert len(lines) == journal + 1 and lines[-1].startswith("og.move_handle('b'")
    og.undo()
    assert og.entity("b").params["points"][1] == (1.0, 0.1, 0.0)
    # Freehand stroke.
    window.set_tool("freehand")
    assert window.handle_press((3, 0, 0))
    for k in range(1, 40):
        window.handle_drag((3 + 0.025 * k, 0.3 * math.sin(k / 6), 0))
    created = window.handle_release((4, 0.3 * math.sin(40 / 6), 0))
    assert len(created) == 1 and og.entity(created[0]).kind == "nurbs"
    # NURBS by control points.
    window.set_tool("nurbs")
    before_names = set(og.names("nurbs"))
    for point in ((50, 0), (60, 10), (70, 0), (80, 10)):
        window.handle_click((*point, 0))
    window.finish_tool()
    (made,) = set(og.names("nurbs")) - before_names
    assert og.entity(made).params["degree"] == 3
    # Domains from the selected curves.
    window.set_tool("select")
    names = lattice(og, columns=(10.0, 11.0), rows=(0.0, 1.0))
    window.set_selection([og.entity(n).id for n in names])
    created = window.domains_from_selected_curves()
    assert len(created) == 1
    window.close()


def test_move_handle_keeps_a_curve_network_connected():
    og = opengrid.Session()
    curves = lattice(og)
    domains = og.domains_from_curves(curves, spacing=0.05)
    # Move the middle point of the middle horizontal stroke well up: its crossings with the
    # vertical strokes move and the connectors follow; every domain stays closed.
    point = og.entity("h1").params["points"][6]
    og.move_handle("h1", "points", 6, (point[0], point[1] + 0.3, 0.0))
    assert og.entity("h1").params["points"][6][1] == pytest.approx(point[1] + 0.3)
    for d in domains:
        assert og.domain_quality(d)["scaled_jacobian"].min() > 0.5
    assert len(og.topology().interfaces()) == 7
    assert og.journal_text().splitlines()[-1].startswith("og.move_handle('h1', 'points', 6")
    assert replays(og)
    og.undo()
    assert og.entity("h1").params["points"][6] == point
    with pytest.raises(ValidationError, match="no movable field"):
        og.move_handle("h1", "control_points", 0, (0, 0, 0))


def test_freehand_stroke_closes_within_the_close_tolerance():
    og = opengrid.Session()
    t = np.linspace(0.0, 2 * np.pi, 120)[:-2]  # the hand stops short of the start
    pts = np.stack([np.cos(t), np.sin(t), 0 * t], axis=1)
    gap = float(np.linalg.norm(pts[-1] - pts[0]))
    (open_,) = og.create_freehand(pts, close_tolerance=0.5 * gap)
    (closed,) = og.create_freehand(pts, close_tolerance=2.0 * gap)
    c_open = og.geometry(open_)
    c_closed = og.geometry(closed)
    assert np.linalg.norm(np.subtract(c_open.eval(1.0), c_open.eval(0.0))) > 0.5 * gap
    assert np.linalg.norm(np.subtract(c_closed.eval(1.0), c_closed.eval(0.0))) == 0.0
    assert distance_to([og.entity(closed).params], pts) < 0.005 * 2 * np.pi  # the fit tolerance
    # A scribble shorter than four tolerances is never closed into a loop.
    (scribble,) = og.create_freehand(pts[:3] * 0.01, close_tolerance=1.0)
    c = og.geometry(scribble)
    assert np.linalg.norm(np.subtract(c.eval(1.0), c.eval(0.0))) > 0.0


def spiral_with_tail():
    """One stroke: almost a full turn, then a tail that crosses its own beginning."""
    t = np.linspace(0.0, 2 * np.pi - 0.3, 200)
    radius = 1.0 - 0.1 * t / (2 * np.pi)
    spiral = np.stack([radius * np.cos(t), radius * np.sin(t), 0 * t], axis=1)
    steps = np.linspace(0.0, 1.0, 25)[1:, None]
    tail = spiral[-1] + steps * (np.array([1.0, 0.8, 0.0]) - spiral[-1])
    return spiral, np.vstack([spiral, tail])


def test_self_crossing_stroke_bounds_a_region_and_an_open_one_does_not():
    og = opengrid.Session()
    spiral, stroke = spiral_with_tail()
    (curve,) = og.create_freehand(stroke)
    (domain,) = og.domains_from_curves([curve], spacing=0.1, unstructured=True)
    assert og.last_warnings == []  # the tails beyond the crossing are dropped silently
    assert "triangles" in og.domain_status(domain)
    assert replays(og)
    og2 = opengrid.Session()
    (curve,) = og2.create_freehand(spiral)
    with pytest.raises(ValidationError, match="enclose no region"):
        og2.domains_from_curves([curve], spacing=0.1, unstructured=True)
    assert og2.names("unstructured_domain") == []


def test_loose_curves_and_bridges_are_reported_not_meshed():
    og = opengrid.Session()
    square = (("s0", (0, 0), (3, 0)), ("s1", (3, 0), (3, 3)), ("s2", (3, 3), (0, 3)))
    for name, a, b in (*square, ("s3", (0, 3), (0, 0)), ("loose", (1, 1), (1.5, 1.5))):
        og.create_line(a, b, name=name)
    (domain,) = og.domains_from_curves(["s0", "s1", "s2", "s3", "loose"], spacing=0.25)
    assert og.entity(domain).kind == "domain"
    assert og.last_warnings == [
        "loose: bounds no region (an open curve, or one whose ends meet nothing): ignored"
    ]
    # An inner square joined to the outer one by a line: the ring between them has that
    # line on both sides of its boundary; it is reported, the inner square is meshed.
    inner = (("i0", (1, 1), (2, 1)), ("i1", (2, 1), (2, 2)), ("i2", (2, 2), (1, 2)))
    for name, a, b in (*inner, ("i3", (1, 2), (1, 1)), ("bridge", (2, 1.5), (3, 1.5))):
        og.create_line(a, b, name=name)
    og.delete(["loose", domain])
    curves = ["s0", "s1", "s2", "s3", "i0", "i1", "i2", "i3", "bridge"]
    made = og.domains_from_curves(curves, spacing=0.25, unstructured=True)
    assert len(made) == 1 and og.entity(made[0]).kind == "domain"
    assert any("runs along a curve twice" in warning for warning in og.last_warnings)


def test_insert_and_remove_points_of_a_polyline_and_a_spline():
    og = opengrid.Session()
    og.create_polyline([(0, 0), (1, 0), (1, 1)], name="pl")
    og.insert_handle("pl", (0.5, 0.1))
    assert og.entity("pl").params["points"] == (
        (0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
    )
    og.remove_handle("pl", 1)
    og.remove_handle("pl", 1)
    assert og.entity("pl").params["points"] == ((0.0, 0.0, 0.0), (1.0, 1.0, 0.0))
    with pytest.raises(ValidationError, match="at least two"):
        og.remove_handle("pl", 0)
    og.create_spline([(0, 0), (1, 0.5), (2, 0)], name="sp")
    before = og.geometry("sp")
    og.insert_handle("sp", (0.5, 0.3))
    points = og.entity("sp").params["points"]
    assert len(points) == 4
    assert float(before.project(points[1])[2]) < 1e-9  # the new point lies on the old curve
    after = og.geometry("sp")
    along = [float(after.project(p)[0]) for p in points]
    assert along == sorted(along)  # in order along the curve
    assert replays(og)
    og.undo()
    assert len(og.entity("sp").params["points"]) == 3  # one step


def test_nurbs_control_point_by_knot_insertion_keeps_the_curve():
    from opengrid.model.knots import removable_control_points

    og = opengrid.Session()
    t = np.linspace(0.0, 2 * np.pi, 40)
    (name,) = og.create_freehand(np.stack([np.cos(t), np.sin(t), 0 * t], axis=1), name="nb")
    s = np.linspace(0.0, 1.0, 101)
    reference = og.geometry(name).eval_many(s)
    first = og.entity(name).params
    og.insert_handle(name, (0.0, 1.0))
    p = og.entity(name).params
    assert len(p["control_points"]) == len(first["control_points"]) + 1
    assert len(p["knots"]) == len(p["control_points"]) + p["degree"] + 1
    assert np.max(np.linalg.norm(og.geometry(name).eval_many(s) - reference, axis=1)) < 1e-12
    # The inserted knot's control point comes out again exactly ...
    (new_knot,) = set(p["knots"]) - set(first["knots"])
    removable = removable_control_points(p["degree"], p["knots"])
    (index,) = [i for i, r in removable.items() if p["knots"][r] == new_knot]
    og.remove_handle(name, index)
    assert og.last_warnings == []
    assert np.max(np.linalg.norm(og.geometry(name).eval_many(s) - reference, axis=1)) < 1e-12
    # ... a genuine control point does not: the curve changes a little, and says so.
    og.remove_handle(name, 3)
    assert len(og.entity(name).params["control_points"]) == len(first["control_points"]) - 1
    assert any("not exact" in w for w in og.last_warnings)
    with pytest.raises(ValidationError, match="cannot be removed"):
        og.remove_handle(name, 0)
    with pytest.raises(ValidationError, match="curve's end"):
        og.insert_handle(name, (1.0, 0.0))  # the stroke's start
    assert replays(og)


def test_inserted_spline_point_keeps_a_curve_network_connected():
    og = opengrid.Session()
    curves = lattice(og)
    domains = og.domains_from_curves(curves, spacing=0.05)
    interfaces = len(og.topology().interfaces())
    grids = {d: og.domain_grid(d).copy() for d in domains}
    og.insert_handle(curves[0], (0.1, 0.6))  # the first vertical spline, between crossings
    assert len(og.topology().interfaces()) == interfaces
    for d in domains:
        assert og.domain_quality(d)["scaled_jacobian"].min() > 0.9
        assert np.max(np.abs(og.domain_grid(d) - grids[d])) < 0.05  # the curve barely moved
    assert replays(og)
