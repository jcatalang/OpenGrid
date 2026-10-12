"""Golden case 1: unit square, uniform and stretched (structured TFI).

Checks: 20 x 20 cells each, total area 1, uniform cells perfect (scaled Jacobian 1), the
stretched grid's first spacing equals the request (1e-3) to 1e-10.
"""

import numpy as np

import opengrid


def square(og, x0, tag, distribution):
    corners = [(x0, 0), (x0 + 1, 0), (x0 + 1, 1), (x0, 1)]
    names = []
    for k in range(4):
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"{tag}{k}")
        dist = "uniform"
        if k in (0, 3):  # cluster towards the corner (x0, 0) on edges 0 and 3
            dist = distribution
        names.append(og.create_connector(f"{tag}{k}", 21, dist, name=f"c-{tag}{k}"))
    edges = [[names[0]], [names[1]], [(names[2], False)], [(names[3], False)]]
    return og.create_domain(edges, name=tag)


def build() -> opengrid.Session:
    og = opengrid.Session()
    square(og, 0.0, "uniform", "uniform")
    stretched = {"type": "tanh", "begin_spacing": 1e-3}
    square(og, 2.0, "stretched", stretched)
    return og


def check(og) -> dict:
    result = {}
    for name in ("uniform", "stretched"):
        q = og.domain_quality(name)
        assert q["area"].shape == (20, 20)
        assert abs(q["area"].sum() - 1.0) < 1e-13
        assert q["scaled_jacobian"].min() > 0
        result[f"{name}_min_sj"] = float(q["scaled_jacobian"].min())
    assert result["uniform_min_sj"] > 1 - 1e-13
    spacing = og.connector_nodes("c-stretched0").begin_spacing
    assert abs(spacing - 1e-3) < 1e-10
    result["stretched_first_spacing"] = spacing
    grid = og.domain_grid("stretched")
    assert np.abs(grid[0, 1] - grid[0, 0]).max() - 1e-3 < 1e-10
    return result


if __name__ == "__main__":
    print(check(build()))
