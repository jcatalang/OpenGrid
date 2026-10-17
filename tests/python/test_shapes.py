"""Standard shapes: NACA 4-digit sections, exact ellipses, rectangles, boxes, cylinders,
spheres and wings, from the Session and from Create > Shapes."""

import math

import numpy as np
import pytest

import opengrid
from opengrid.model.shapes import box_planes, ellipse_nurbs, frame, naca4, wing_sections
from opengrid.project.ogp import to_json


def test_naca4_sections():
    upper, lower = naca4("0012", 101)
    assert upper[0] == lower[0] == (0.0, 0.0) and upper[-1] == lower[-1] == (1.0, 0.0)
    xs = [p[0] for p in upper]
    assert xs == sorted(xs) and abs(xs[50] - 0.5) < 1e-12  # cosine spaced
    thickest = max(p[1] for p in upper)
    assert abs(thickest - 0.06) < 1e-3  # 12 % thick: 0.06 each side
    assert all(abs(u[1] + lo[1]) < 1e-12 for u, lo in zip(upper, lower, strict=True))  # symmetric
    upper, lower = naca4("2412", 51)
    camber = [(u[1] + lo[1]) / 2 for u, lo in zip(upper, lower, strict=True)]
    assert abs(max(camber) - 0.02) < 2e-3  # 2 % camber at 40 % chord
    assert upper[-1] == lower[-1] == (1.0, 0.0)
    _upper, lower_open = naca4("0012", 11, closed_trailing_edge=False)
    assert lower_open[-1][1] < -1e-3
    with pytest.raises(ValueError, match="four digits"):
        naca4("12", 11)
    with pytest.raises(ValueError, match="at least 3"):
        naca4("0012", 2)


def test_ellipse_frame_and_box():
    degree, knots, control, weights = ellipse_nurbs((1, 2), 3.0, 1.0)
    assert degree == 2 and len(knots) == 12 and len(control) == 9 == len(weights)
    assert control[0] == (4.0, 2.0, 0.0) and control[2] == (1.0, 3.0, 0.0)
    u, v, w = frame((0, 0, 2))
    assert w == (0.0, 0.0, 1.0) and u == (1.0, 0.0, 0.0) and v == (0.0, 1.0, 0.0)
    u, v, w = frame((1, 0, 0))
    assert w == (1.0, 0.0, 0.0) and abs(np.linalg.norm(np.cross(u, v) - np.asarray(w))) < 1e-12
    with pytest.raises(ValueError, match="zero"):
        frame((0, 0, 0))
    faces = box_planes((0, 0, 0), (2, 4, 6))
    assert [f[0] for f in faces] == ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    for suffix, _origin, uvec, vvec in faces:
        normal = np.cross(uvec, vvec)
        normal /= np.linalg.norm(normal)
        axis = "xyz".index(suffix[0])
        sign = {"min": -1.0, "max": 1.0}[suffix[1:]]
        assert normal[axis] == pytest.approx(sign)  # outward
    sections = wing_sections("0012", 1.0, 0.5, 2.0, 30.0, 0.0, 21)
    tip_upper = sections["tip"][0]
    assert tip_upper[0] == pytest.approx((2.0 * math.tan(math.radians(30.0)), 2.0, 0.0))
    assert tip_upper[-1][0] - tip_upper[0][0] == pytest.approx(0.5)


def test_session_shapes_and_replay():
    og = opengrid.Session()
    assert og.create_rectangle((0, 0), 2, 1) == [
        "rectangle-1-bottom",
        "rectangle-1-right",
        "rectangle-1-top",
        "rectangle-1-left",
    ]
    assert og.create_rectangle()[0] == "rectangle-2-bottom"
    assert og.entity("rectangle-1-bottom").params["p0"] == (-1.0, -0.5, 0.0)
    og.create_ellipse((0, 0), 2.0, 1.0, name="e")
    curve = og.geometry("e")
    for k in range(41):
        x, y, _z = curve.eval(k / 40)
        assert abs((x / 2.0) ** 2 + y**2 - 1.0) < 1e-12
    assert og.create_airfoil("2412", chord=2.0, leading_edge=(1, 1), angle=5.0, name="af") == [
        "af-upper",
        "af-lower",
    ]
    upper = og.geometry("af-upper")
    lower = og.geometry("af-lower")
    assert tuple(upper.eval(0.0)) == tuple(lower.eval(0.0)) == (1.0, 1.0, 0.0)
    trailing = upper.eval(1.0)
    assert tuple(trailing) == tuple(lower.eval(1.0))
    assert trailing[0] == pytest.approx(1.0 + 2.0 * math.cos(math.radians(5.0)))
    assert trailing[1] == pytest.approx(1.0 - 2.0 * math.sin(math.radians(5.0)))  # nose up
    assert og.create_box(name="bx") == [
        f"bx-{s}" for s in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
    ]
    names = og.create_cylinder((0, 0, 1), 0.5, 2.0, name="cy")
    assert names == [
        "cy-generatrix",
        "cy-side",
        "cy-bottom-radius",
        "cy-bottom",
        "cy-top-radius",
        "cy-top",
    ]
    assert og.entity("cy-generatrix").params["p1"] == (0.5, 0.0, 3.0)
    assert og.create_cylinder(caps=False, axis=(1, 0, 0)) == [
        "cylinder-1-generatrix",
        "cylinder-1-side",
    ]
    assert og.create_sphere((1, 2, 3), 2.0, name="sp") == ["sp-meridian", "sp"]
    meridian = og.geometry("sp-meridian")
    assert tuple(meridian.eval(0.0)) == pytest.approx((1.0, 2.0, 5.0))  # the +z pole
    assert tuple(meridian.eval(1.0)) == pytest.approx((1.0, 2.0, 1.0))
    wing = og.create_wing("0012", 1.0, 0.5, 3.0, sweep=20.0, dihedral=5.0, name="wg")
    assert wing[-4:] == ["wg-upper", "wg-lower", "wg-root", "wg-tip"] and len(wing) == 8
    assert og.entity("wg-upper").kind == "ruled_surface"
    tip = og.geometry("wg-tip")  # the flat cap: inside the tip section, in its plane y = 3
    for u, v in ((0.3, 0.5), (0.7, 0.2), (0.5, 0.9)):
        point = tip.eval(u, v)
        assert abs(point[1] - 3.0) < 1e-12
    assert og.create_wing(caps=False, name="open")[-2:] == ["open-upper", "open-lower"]
    # An airfoil in a far-field rectangle meshes as one domain with the airfoil as a hole.
    og.create_airfoil(name="naca")
    og.create_rectangle((0.5, 0), 20, 20, name="far")
    domain = og.unstructured_domain_from_curves(
        ["naca-upper", "naca-lower", "far-bottom", "far-right", "far-top", "far-left"], 61
    )
    assert len(og.entity(domain).params["loops"]) == 2
    assert len(og.unstructured_mesh(domain)["triangles"]) > 1000
    og.undo()  # the domain
    og.undo()  # the rectangle
    assert "far-top" not in og.store and "naca-upper" in og.store
    og.undo()
    assert "naca-upper" not in og.store  # the airfoil went as one step
    for _ in range(3):
        og.redo()
    assert domain in og.store
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_shapes_from_the_gui(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow
    from opengrid.app.shapes_dialog import SHAPES, ShapeDialog, ShapeParseError, parse_values

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    names = window.create_shape(
        "airfoil",
        {
            "code": "0012",
            "chord": 1.0,
            "leading_edge": (0, 0),
            "angle": 0.0,
            "points": 51,
            "name": None,
        },
    )
    assert names == ["airfoil-1-upper", "airfoil-1-lower"]
    assert window.selected_names() == names
    (create_menu,) = [a.menu() for a in window.menuBar().actions() if "Create" in a.text()]
    assert "Shapes" in [a.text().replace("&", "") for a in create_menu.actions()]
    # The dialog parses its fields; defaults give a valid shape.
    dialog = ShapeDialog(SHAPES["cylinder"], window)
    values = dialog.values()
    assert values["center"] == (0.0, 0.0, 0.0) and values["caps"] is True and values["name"] is None
    dialog.edits["radius"].setText("abc")
    with pytest.raises(ShapeParseError, match="radius"):
        dialog.values()
    dialog.edits["radius"].setText("2")
    dialog.edits["axis"].setText("1, 0, 0")
    dialog.edits["name"].setText("drum")
    values = dialog.values()
    assert values["axis"] == (1.0, 0.0, 0.0) and values["name"] == "drum"
    assert window.create_shape("cylinder", values)[1] == "drum-side"
    assert parse_values(SHAPES["wing"], {})["tip_chord"] == 0.5
    window.close()
