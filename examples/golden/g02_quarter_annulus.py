"""Golden case 2: quarter annulus 1 <= r <= 2, 0 <= theta <= 90 deg, linear TFI.

Linear TFI between two concentric arcs and two radial lines is exactly the polar grid
r = 1 + xi, theta = eta * 90 deg; checked to 1e-13.
"""

import math

import numpy as np

import opengrid


def build() -> opengrid.Session:
    og = opengrid.Session()
    og.create_line((1, 0), (2, 0), name="bottom")
    og.create_arc((0, 0), 2.0, 0.0, math.pi / 2, name="outer")
    og.create_line((0, 2), (0, 1), name="top")
    og.create_arc((0, 0), 1.0, 0.0, math.pi / 2, name="inner")
    og.create_connector("bottom", 11, name="c-bottom")
    og.create_connector("outer", 31, name="c-outer")
    og.create_connector("top", 11, name="c-top")
    og.create_connector("inner", 31, name="c-inner")
    og.create_domain([["c-bottom"], ["c-outer"], ["c-top"], [("c-inner", True)]], name="annulus")
    return og


def check(og) -> dict:
    grid = og.domain_grid("annulus")
    nj, ni = grid.shape[:2]
    xi = np.linspace(0.0, 1.0, ni)
    eta = np.linspace(0.0, 1.0, nj)
    radius = np.hypot(grid[..., 0], grid[..., 1])
    angle = np.arctan2(grid[..., 1], grid[..., 0])
    r_error = np.abs(radius - (1.0 + xi)[None, :]).max()
    t_error = np.abs(angle - (eta * math.pi / 2)[:, None]).max()
    assert r_error < 1e-13 and t_error < 1e-13
    return {"radius_error": float(r_error), "angle_error": float(t_error)}


if __name__ == "__main__":
    print(check(build()))
