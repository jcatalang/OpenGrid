"""Golden case 3: NACA 0012 O-grid and C-grid (structured, elliptic smoothing with wall control
on the airfoil: first-cell height 1e-4 chord, orthogonality).

Checks: all cells valid (scaled Jacobian > 0), smoothing converged, mean first-cell height on
the airfoil within 5 % of 1e-4.
"""

import math

import numpy as np

import opengrid

WALL = 1e-4
N = 61


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def airfoil(og, prefix):
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 80)) for i in range(81)]
    for sign, side in ((1, "upper"), (-1, "lower")):
        og.create_spline(
            [(x, sign * naca0012(x)) for x in xs],
            "clamped",
            start_tangent=(0, sign),
            name=f"{prefix}{side}",
        )
        og.create_connector(
            f"{prefix}{side}",
            N,
            {"type": "tanh", "begin_spacing": 4e-3, "end_spacing": 2e-3},
            name=f"{prefix}c-{side}",
        )


def generator(method):
    return {
        "method": method,
        "elliptic": {
            "control": "thomas_middlecoff",
            "walls": [{"edge": 0, "spacing": WALL}],
            "max_iterations": 20000,
        },
    }


def build() -> opengrid.Session:
    og = opengrid.Session()
    # O-grid: airfoil inside a circle of radius 15, cut along y = 0 behind the trailing edge.
    airfoil(og, "o-")
    og.create_circle((0.5, 0.0), 15.0, name="o-far")
    og.create_line((1.0, 0.0), (15.5, 0.0), name="o-cut")
    og.create_connector("o-far", 2 * N - 1, name="o-c-far")
    normal = {"type": "tanh", "begin_spacing": WALL, "end_spacing": 1.5}
    og.create_connector("o-cut", 61, normal, name="o-c-cut")
    og.create_domain(
        [[("o-c-lower", True), "o-c-upper"], ["o-c-cut"], ["o-c-far"], [("o-c-cut", True)]],
        generator=generator("hermite_tfi"),
        name="o-grid",
    )
    # C-grid: wake cut from the trailing edge to x = 16, C-shaped farfield of radius 15.
    airfoil(og, "c-")
    og.create_line((1.0, 0.0), (16.0, 0.0), name="c-wake")
    og.create_line((16.0, 0.0), (16.0, 15.0), name="c-out-top")
    og.create_line((16.0, -15.0), (16.0, 0.0), name="c-out-bottom")
    og.create_line((1.0, 15.0), (16.0, 15.0), name="c-top")
    og.create_line((1.0, -15.0), (16.0, -15.0), name="c-bottom")
    og.create_arc((1.0, 0.0), 15.0, math.pi / 2, 3 * math.pi / 2, name="c-front")
    wake = {"type": "tanh", "begin_spacing": 2e-3, "end_spacing": 1.5}
    og.create_connector("c-wake", 41, wake, name="c-c-wake")
    og.create_connector("c-out-top", 61, normal, name="c-c-out-top")
    reverse_normal = {"type": "tanh", "begin_spacing": 1.5, "end_spacing": WALL}
    og.create_connector("c-out-bottom", 61, reverse_normal, name="c-c-out-bottom")
    og.create_connector("c-top", 41, wake, name="c-c-top")
    og.create_connector("c-bottom", 41, wake, name="c-c-bottom")
    og.create_connector("c-front", 2 * N - 1, name="c-c-front")
    og.create_domain(
        [
            [("c-c-wake", True), ("c-c-lower", True), "c-c-upper", "c-c-wake"],
            ["c-c-out-top"],
            [("c-c-top", True), "c-c-front", "c-c-bottom"],
            ["c-c-out-bottom"],
        ],
        generator=generator("tfi"),
        name="c-grid",
    )
    return og


def check(og) -> dict:
    result = {}
    for name, airfoil_nodes in (("o-grid", slice(1, -1)), ("c-grid", slice(41, -41))):
        smoothing = og.smooth_domain(name)
        assert smoothing.converged, f"{name}: smoothing did not converge"
        grid = og.domain_grid(name)
        q = og.domain_quality(name)
        assert q["scaled_jacobian"].min() > 0
        height = np.linalg.norm(grid[1, airfoil_nodes] - grid[0, airfoil_nodes], axis=1)
        mean = float(height.mean())
        assert abs(mean - WALL) < 0.05 * WALL, f"{name}: mean first height {mean}"
        result[f"{name}_mean_first_height"] = mean
        result[f"{name}_min_sj"] = float(q["scaled_jacobian"].min())
    return result


if __name__ == "__main__":
    print(check(build()))
