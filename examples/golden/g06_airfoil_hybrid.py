"""Golden case 6: NACA 0012 hybrid mesh: quad boundary layers (first height 1e-5, growth 1.2,
up to 40 layers) from the airfoil and triangles (constrained Delaunay) to a circular farfield.

Checks: quads and triangles present, first layer height 1e-5 to 1e-9, area equal to the
discrete farfield polygon minus the discrete airfoil polygon, all cells valid.
"""

import math

import numpy as np

import opengrid

FIRST = 1e-5


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def build() -> opengrid.Session:
    og = opengrid.Session()
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    for sign, side in ((1, "upper"), (-1, "lower")):
        og.create_spline(
            [(x, sign * naca0012(x)) for x in xs], "clamped", start_tangent=(0, sign), name=side
        )
        og.create_connector(
            side, 81, {"type": "tanh", "begin_spacing": 3e-3, "end_spacing": 2e-3}, name=f"c-{side}"
        )
    og.create_circle((0.5, 0.0), 20.0, name="far")
    og.create_connector("far", 101, name="c-far")
    walls = ["c-upper", "c-lower"]
    og.create_unstructured_domain(
        connectors=[*walls, "c-far"],
        generator={
            "growth_rate": 1.25,
            "inflation": {
                "connectors": walls,
                "first_height": FIRST,
                "growth_rate": 1.2,
                "layers": 40,
            },
        },
        name="flow",
    )
    og.create_patch("airfoil", "wall", walls)
    og.create_patch("farfield", "farfield", ["c-far"])
    return og


def polygon_area(points) -> float:
    x, y = points[:, 0], points[:, 1]
    return 0.5 * abs(float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def check(og) -> dict:
    mesh = og.unstructured_mesh("flow")
    q = og.domain_quality("flow")
    assert len(mesh["quads"]) > 0 and len(mesh["triangles"]) > 0
    assert q["scaled_jacobian"].min() > 0
    upper = og.connector_nodes("c-upper").points
    lower = og.connector_nodes("c-lower").points
    airfoil = np.concatenate([upper, lower[::-1][1:-1]])
    far = og.connector_nodes("c-far").points[:-1]
    expected = polygon_area(far) - polygon_area(airfoil)
    area = float(q["area"].sum())
    assert abs(area - expected) < 1e-9 * expected
    # First layer: the wall node closest to mid-chord and the layer node above it.
    points = mesh["points"]
    wall = upper[len(upper) // 2]
    distances = np.linalg.norm(points[:, :2] - wall[:2], axis=1)
    order = np.argsort(distances)
    first = float(distances[order[1]])
    assert abs(first - FIRST) < 1e-9, first
    assert og.boundary_report().ok
    return {
        "quads": len(mesh["quads"]),
        "triangles": len(mesh["triangles"]),
        "area": area,
        "first_height": first,
    }


if __name__ == "__main__":
    print(check(build()))
