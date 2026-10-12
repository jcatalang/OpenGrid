"""Golden case 4: O-grid around a cylinder (radius 0.5, farfield radius 10) with a wall spacing
of 1e-4 (tanh radial distribution on the cut).

Between concentric circles linear TFI is the polar grid, so the first cell height equals the
requested spacing exactly and the grid lines meet the wall at 90 deg; checked to 1e-10.
"""

import numpy as np

import opengrid

WALL = 1e-4


def build() -> opengrid.Session:
    og = opengrid.Session()
    og.create_circle((0, 0), 0.5, name="cylinder")
    og.create_circle((0, 0), 10.0, name="far")
    og.create_line((0.5, 0), (10.0, 0), name="cut")
    og.create_connector("cylinder", 129, name="c-cylinder")
    og.create_connector("far", 129, name="c-far")
    og.create_connector(
        "cut", 81, {"type": "tanh", "begin_spacing": WALL, "end_spacing": 1.0}, name="c-cut"
    )
    og.create_domain(
        [[("c-cylinder", True)], ["c-cut"], ["c-far"], [("c-cut", True)]], name="ogrid"
    )
    og.create_patch("wall", "wall", ["c-cylinder"])
    og.create_patch("farfield", "farfield", ["c-far"])
    return og


def check(og) -> dict:
    grid = og.domain_grid("ogrid")  # (nj, ni, 3): j = 0 on the cylinder
    height = np.linalg.norm(grid[1] - grid[0], axis=1)
    radial = (grid[1] - grid[0]) / height[:, None]
    outward = grid[0] / np.linalg.norm(grid[0], axis=1)[:, None]
    angle_error = float(
        np.degrees(np.arccos(np.clip(np.sum(radial * outward, axis=1), -1, 1))).max()
    )
    assert np.abs(height - WALL).max() < 1e-10
    assert angle_error < 1e-5
    assert og.boundary_report().ok
    return {"first_height_error": float(np.abs(height - WALL).max()), "angle_error": angle_error}


if __name__ == "__main__":
    print(check(build()))
