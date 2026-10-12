"""Golden case 8: annular sector duct 1 <= r <= 2, 0 <= theta <= 60 deg, 0 <= z <= 1: the
sector cross-section (linear TFI, the exact polar grid) extruded along z.

Checks: volume equal to the polygonal sector (straight cell edges) to 1e-12, all hexahedra
valid, the cross-section grid on the exact polar grid.
"""

import math

import numpy as np

import opengrid

THETA = math.pi / 3
N_THETA = 31


def build() -> opengrid.Session:
    og = opengrid.Session()
    og.create_line((1, 0, 0), (2, 0, 0), name="bottom")
    og.create_arc((0, 0, 0), 2.0, 0.0, THETA, name="outer")
    og.create_line(
        (2 * math.cos(THETA), 2 * math.sin(THETA), 0),
        (math.cos(THETA), math.sin(THETA), 0),
        name="top",
    )
    og.create_arc((0, 0, 0), 1.0, 0.0, THETA, name="inner")
    radial = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.01}
    og.create_connector("bottom", 21, radial, name="c-bottom")
    og.create_connector("outer", N_THETA, name="c-outer")
    og.create_connector(
        "top", 21, {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.01}, name="c-top"
    )
    og.create_connector("inner", N_THETA, name="c-inner")
    og.create_domain([["c-bottom"], ["c-outer"], ["c-top"], [("c-inner", True)]], name="section")
    og.create_extruded_block(
        "section",
        "translate",
        vector=(0, 0, 1),
        dimension=21,
        distribution={"type": "tanh", "begin_spacing": 0.02, "end_spacing": 0.02},
        name="sector",
    )
    return og


def check(og) -> dict:
    q = og.block_quality("sector")
    assert q["scaled_jacobian"].min() > 0
    volume = float(q["volume"].sum())
    n = N_THETA - 1
    polygon = 0.5 * n * math.sin(THETA / n) * (2.0**2 - 1.0**2)
    assert abs(volume - polygon) < 1e-12
    grid = og.domain_grid("section")
    angle = np.arctan2(grid[..., 1], grid[..., 0])
    expected = np.linspace(0.0, THETA, grid.shape[0])[:, None]
    assert np.abs(angle - expected).max() < 1e-13
    return {"volume": volume, "polygon_volume": polygon}


if __name__ == "__main__":
    print(check(build()))
