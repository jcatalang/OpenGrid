"""Phase 10e: topology templates and automatic 2D topology (cross fields)."""

import math

import numpy as np
import pytest

import opengrid
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def naca(x, t=0.12):
    return 5 * t * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)


def airfoil(og, name="foil", x0=0.0, y0=0.0, chord=1.0, aoa=0.0):
    x = 0.5 * (1 - np.cos(np.linspace(0, math.pi, 61)))
    c, s = math.cos(aoa), math.sin(aoa)

    def at(px, py):
        return (x0 + chord * (c * px + s * py), y0 + chord * (-s * px + c * py))

    og.create_spline([at(a, naca(a)) for a in x], name=f"{name}-u")
    og.create_spline([at(a, -naca(a)) for a in x], name=f"{name}-l")
    return [f"{name}-u", f"{name}-l"]


def channel(og, bump=False, width=3.0):
    if bump:
        x = np.linspace(0, width, 41)
        og.create_spline(
            [(a, 0.1 * np.exp(-(((a - width / 2) / 0.3) ** 2))) for a in x], name="bottom"
        )
    else:
        og.create_line((0, 0), (width, 0), name="bottom")
    og.create_line((width, 0), (width, 1), name="right")
    og.create_line((width, 1), (0, 1), name="top")
    og.create_line((0, 1), (0, 0), name="left")
    return {"bottom": "bottom", "right": "right", "top": "top", "left": "left"}


def worst(og) -> float:
    if og.names("block"):
        return min(float(og.block_quality(b)["scaled_jacobian"].min()) for b in og.names("block"))
    return min(float(og.domain_quality(d)["scaled_jacobian"].min()) for d in og.names("domain"))


def smoothed(og) -> float:
    og.smooth_group(og.names("smoothing_group")[0])
    return worst(og)


def ellipse(og):
    t = np.linspace(0, 2 * math.pi, 73)[:-1]
    points = [(0.5 * math.cos(a), 0.2 * math.sin(a)) for a in t]
    og.create_spline([*points, points[0]], method="catmull_rom", name="ell")
    return ["ell"]


@pytest.mark.parametrize(
    "kind",
    ["ogrid", "butterfly", "cgrid", "hgrid", "body_in_channel"],
)
def test_2d_templates_mesh_their_golden_geometry(kind):
    og = opengrid.Session()
    if kind == "ogrid":
        name = og.create_template("ogrid", body=ellipse(og), wall=0.002)
    elif kind == "butterfly":
        og.create_circle((0, 0), 1.0, name="c")
        name = og.create_template("butterfly", boundary=["c"])
    elif kind == "cgrid":
        name = og.create_template("cgrid", body=airfoil(og), wall=0.001)
    elif kind == "hgrid":
        name = og.create_template("hgrid", **channel(og, bump=True), blocks=(3, 1), wall=0.005)
    else:
        og.create_circle((1.0, 0.5), 0.2, name="cyl")
        name = og.create_template("body_in_channel", body=["cyl"], **channel(og), wall=0.005)
    assert og.entity(name).kind == "sketch"
    assert smoothed(og) >= 0.4


@pytest.mark.parametrize("kind", ["pipe", "cylinder_in_channel", "bump_channel"])
def test_3d_templates_mesh_their_golden_geometry(kind):
    og = opengrid.Session()
    if kind == "pipe":
        og.create_circle((0, 0), 0.5, name="c")
        og.create_template("pipe", boundary=["c"], length=2.0, span_points=11)
    elif kind == "cylinder_in_channel":
        og.create_circle((1.0, 0.5), 0.2, name="cyl")
        og.create_template(
            "cylinder_in_channel",
            body=["cyl"],
            **channel(og),
            wall=0.005,
            length=1.0,
            span_points=7,
        )
    else:
        og.create_template(
            "bump_channel", **channel(og, bump=True), blocks=(3, 1), wall=0.005, length=1.0
        )
    assert og.names("block")
    assert smoothed(og) >= 0.4


def test_sphere_in_box_template():
    og = opengrid.Session()
    og.create_template("sphere_in_box", radius=0.5, wall=0.01)
    assert len(og.names("block")) == 32
    assert smoothed(og) >= 0.4


def test_template_is_one_undo_step_and_replays():
    og = opengrid.Session()
    body = airfoil(og)
    before = len(og.names())
    og.create_template("cgrid", body=body, wall=0.001, name="c")
    assert og.entity("c-group").kind == "smoothing_group"
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.undo()
    assert len(og.names()) == before


def test_ogrid_refuses_a_cusp_and_bad_options_are_reported():
    og = opengrid.Session()
    body = airfoil(og)
    with pytest.raises(ValidationError, match="cgrid template"):
        og.create_template("ogrid", body=body)
    with pytest.raises(ValidationError, match="unexpected options"):
        og.create_template("cgrid", body=body, colour="red")
    with pytest.raises(ValidationError, match="unknown template"):
        og.create_template("star", body=body)


# --------------------------------------------------------------------- automatic topology
def box(og, x0, y0, x1, y1):
    og.create_line((x0, y0), (x1, y0), name="b0")
    og.create_line((x1, y0), (x1, y1), name="b1")
    og.create_line((x1, y1), (x0, y1), name="b2")
    og.create_line((x0, y1), (x0, y0), name="b3")
    return ["b0", "b1", "b2", "b3"]


def airfoil_in_box(og):
    return box(og, -2, -2, 3, 2) + airfoil(og)


def two_element(og):
    main = airfoil(og, "main")
    flap = airfoil(og, "flap", x0=1.05, y0=-0.08, chord=0.35, aoa=math.radians(25))
    return box(og, -2, -2, 3.5, 2) + main + flap


def plate(og):
    curves = box(og, 0, 0, 4, 2)
    for k, (x, y) in enumerate(((1, 1), (2, 0.8), (3, 1.2))):
        og.create_circle((x, y), 0.3, name=f"hole{k}")
        curves.append(f"hole{k}")
    return curves


def bracket(og):
    w, h, r = 1.2, 0.8, 0.12
    og.create_line((r, 0), (w - r, 0), name="s0")
    og.create_arc((w - r, r), r, -math.pi / 2, 0, name="f0")
    og.create_line((w, r), (w, h - r), name="s1")
    og.create_arc((w - r, h - r), r, 0, math.pi / 2, name="f1")
    og.create_line((w - r, h), (r, h), name="s2")
    og.create_arc((r, h - r), r, math.pi / 2, math.pi, name="f2")
    og.create_line((0, h - r), (0, r), name="s3")
    og.create_arc((r, r), r, math.pi, 1.5 * math.pi, name="f3")
    og.create_circle((0.6, 0.4), 0.15, name="hole")
    return ["s0", "f0", "s1", "f1", "s2", "f2", "s3", "f3", "hole"]


CASES = {
    "airfoil_in_box": airfoil_in_box,
    "two_element": two_element,
    "plate_three_holes": plate,
    "bracket_section": bracket,
}


@pytest.mark.parametrize("case", list(CASES))
def test_cross_field_index_theorem(case):
    og = opengrid.Session()
    layout = og.cross_field(CASES[case](og))
    field = layout.field
    assert field.index_sum() == pytest.approx(field.euler())
    assert all(abs(abs(s.index) - 0.25) < 1e-12 for s in field.singularities)


@pytest.mark.parametrize(
    "case",
    [
        "airfoil_in_box",
        pytest.param(
            "two_element",
            marks=pytest.mark.xfail(
                reason="a face in the slot between main element and flap gets 5 corners",
                strict=True,
            ),
        ),
        "plate_three_holes",
        "bracket_section",
    ],
)
def test_automatic_topology_gives_four_sided_conforming_blocks(case):
    og = opengrid.Session()
    curves = CASES[case](og)
    sketch = og.auto_topology(curves)
    blocks = og.entity(sketch).params["blocks"]
    assert blocks and all(len(b) == 4 for b in blocks)
    # Conforming: every generated domain is a valid structured domain and the group smooths.
    assert len(og.names("domain")) == len(blocks)
    assert smoothed(og) > 0.0


def test_automatic_topology_undo_and_replay():
    og = opengrid.Session()
    curves = airfoil_in_box(og)
    before = len(og.names())
    og.auto_topology(curves, name="auto")
    assert og.names("polyline")  # the separatrices
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.undo()
    assert len(og.names()) == before


def test_gui_templates_and_automatic_topology(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    body = ellipse(og)
    window.set_selection([og.entity(n).id for n in body])
    name = window.template_from_selection("ogrid")
    assert og.entity(name).kind == "sketch"
    window.smoother_dock.wait()
    curves = airfoil_in_box(og)
    window.set_selection([og.entity(n).id for n in curves])
    auto = window.auto_topology_from_selection()
    assert og.entity(auto).kind == "sketch"
    window.smoother_dock.wait()
    window.close()
