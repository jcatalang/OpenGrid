import math

import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.project import to_json


def _rectangle(og, w=2.0, h=1.0, ni=11, nj=6, x0=0.0, prefix=""):
    p = [(x0, 0), (x0 + w, 0), (x0 + w, h), (x0, h)]
    names = []
    for k, (a, b) in enumerate(zip(p, p[1:] + p[:1], strict=True)):
        line = og.create_line(a, b, name=f"{prefix}l{k}")
        dimension = nj
        if k % 2 == 0:
            dimension = ni
        names.append(og.create_connector(line, dimension, name=f"{prefix}c{k}"))
    return names


def _annulus(og, ni=41, nj=21, r0=1.0, r1=5.0, wall=None):
    og.create_circle((0, 0), r0, name="inner")
    og.create_circle((0, 0), r1, name="outer")
    og.create_line((r0, 0), (r1, 0), name="cut")
    og.create_connector("inner", ni, name="c-inner")
    og.create_connector("outer", ni, name="c-outer")
    dist = "uniform"
    if wall is not None:
        dist = {"type": "tanh", "begin_spacing": wall}
    og.create_connector("cut", nj, dist, name="c-cut")
    return ["c-inner", "c-outer", "c-cut"]


def test_rectangle_auto_assembly_and_tfi():
    og = opengrid.Session()
    names = _rectangle(og)
    d = og.create_structured_domain(connectors=list(reversed(names)))
    assert og.domain_dimensions(d) == (11, 6)
    grid = og.domain_grid(d)
    assert grid.shape == (6, 11, 3)
    assert np.allclose(grid[0, :, 1], 0.0)
    assert np.allclose(grid[:, -1, 0], 2.0)
    q = og.domain_quality(d)
    assert np.allclose(q["scaled_jacobian"], 1.0)
    assert np.allclose(q["area"].sum(), 2.0)


def test_assembly_chooses_corners_with_matching_counts():
    og = opengrid.Session()
    # Bottom split into two connectors, top a single one: 6 + 6 - 1 = 11 points.
    og.create_line((0, 0), (1, 0), name="b1")
    og.create_line((1, 0), (2, 0), name="b2")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_line((2, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    for name, n in [("b1", 6), ("b2", 6), ("r", 5), ("t", 11), ("l", 5)]:
        og.create_connector(name, n, name=f"c-{name}")
    d = og.create_structured_domain(connectors=["c-b1", "c-b2", "c-r", "c-t", "c-l"])
    assert og.domain_dimensions(d) == (11, 5)
    edges = og.entity(d).params["edges"]
    assert len(edges[0]) == 2  # the straight junction between b1 and b2 is not a corner


def test_o_type_assembly_puts_inner_loop_at_j0():
    og = opengrid.Session()
    names = _annulus(og)
    d = og.create_structured_domain(connectors=names)
    assert og.domain_dimensions(d) == (41, 21)
    grid = og.domain_grid(d)
    radius = np.hypot(grid[..., 0], grid[..., 1])
    assert np.allclose(radius[0], 1.0)
    assert np.allclose(radius[-1], 5.0)
    assert og.domain_quality(d)["scaled_jacobian"].min() > 0.0


@pytest.mark.parametrize(
    ("edges", "message"),
    [
        ([["c0"], ["c1"], ["c2"]], "exactly 4 edges"),
        ([["c0"], ["c1"], ["c2"], ["missing"]], "no entity 'missing'"),
        ([["c0"], ["c1"], ["c2"], ["l3"]], "is a line, not a connector"),
        ([["c0"], ["c1"], ["c2"], [("c3", True)]], "does not end where"),
        ([["c0"], ["c3"], ["c2"], ["c1"]], "does not end where"),
        ([[("c2", True)], [("c1", True)], [("c0", True)], [("c3", True)]], "clockwise"),
    ],
)
def test_explicit_edge_validation(edges, message):
    og = opengrid.Session()
    _rectangle(og)
    with pytest.raises(opengrid.ValidationError, match=message):
        og.create_domain(edges)


def test_clockwise_and_count_mismatch_reported():
    og = opengrid.Session()
    _rectangle(og)
    # Clockwise loop: c3 reversed (up the left side), c2 reversed, c1 reversed, c0 reversed.
    with pytest.raises(opengrid.ValidationError, match="clockwise"):
        og.create_domain([[("c3", True)], [("c2", True)], [("c1", True)], [("c0", True)]])
    og.set_dimension("c2", 9)
    with pytest.raises(
        opengrid.ValidationError,
        match=(
            r"edge 0 \(j = 0\) \[c0\] has 11 points but the opposite "
            r"edge 2 \(j = jmax\) \[c2\] has 9"
        ),
    ):
        og.create_domain([["c0"], ["c1"], ["c2"], ["c3"]])


def test_dependencies_and_rotation():
    og = opengrid.Session()
    names = _rectangle(og)
    d = og.create_structured_domain(connectors=names)
    with pytest.raises(opengrid.ValidationError, match="used by domain"):
        og.delete("c0")
    with pytest.raises(opengrid.ValidationError, match="would become invalid"):
        og.set_dimension("c0", 21)  # breaks the opposite counts
    with pytest.raises(opengrid.ValidationError, match="would become invalid"):
        og.modify("l1", p1=(2.5, 1.0))  # breaks the loop through connector c1
    og.modify("l1", p0=(2.0, 0.0))  # unchanged geometry is fine
    before = og.domain_grid(d)
    og.rotate_domain(d)
    assert og.domain_dimensions(d) == (6, 11)
    after = og.domain_grid(d)
    assert np.allclose(
        np.sort(before.reshape(-1, 3), axis=0), np.sort(after.reshape(-1, 3), axis=0)
    )
    assert og.domain_quality(d)["scaled_jacobian"].min() > 0
    og.undo()
    assert og.domain_dimensions(d) == (11, 6)


def test_shared_connector_points_match_exactly():
    og = opengrid.Session()
    left = _rectangle(og, prefix="a")
    og.create_line((2, 0), (4, 0), name="bl0")
    og.create_line((4, 0), (4, 1), name="bl1")
    og.create_line((4, 1), (2, 1), name="bl2")
    og.create_connector("bl0", 9, name="bc0")
    og.create_connector("bl1", 6, name="bc1")
    og.create_connector("bl2", 9, name="bc2")
    a = og.create_structured_domain(connectors=left)
    b = og.create_domain([["bc0"], ["bc1"], ["bc2"], [("ac1", True)]])
    ga = og.domain_grid(a)
    gb = og.domain_grid(b)
    assert np.array_equal(ga[:, -1, :], gb[:, 0, :])  # bitwise identical shared column


def test_hermite_and_elliptic_wall_control():
    og = opengrid.Session()
    names = _annulus(og, ni=61, nj=31, wall=1e-3)
    # Distort the outer boundary so the grid is not trivially orthogonal: ellipse-like.
    d = og.create_structured_domain(connectors=names)
    og.set_domain_generator(d, method="hermite_tfi")
    assert og.domain_quality(d)["scaled_jacobian"].min() > 0
    og.set_domain_generator(
        d,
        elliptic={
            "control": "thomas_middlecoff",
            "walls": [{"edge": 0, "spacing": 5e-4}],
            "max_iterations": 4000,
        },
    )
    assert og.domain_status(d).startswith("smoothing pending")
    result = og.smooth_domain(d)
    assert result.converged
    assert result.residuals[-1] <= 1e-6 * result.residuals[0]
    grid = og.domain_grid(d)
    spacing = np.linalg.norm(grid[1] - grid[0], axis=1)
    assert np.allclose(spacing[1:-1], 5e-4, rtol=0.05)
    assert og.domain_quality(d)["scaled_jacobian"].min() > 0
    assert "converged" in og.domain_status(d)


def test_smoothing_cache_and_determinism():
    og = opengrid.Session()
    names = _annulus(og)
    d = og.create_structured_domain(
        connectors=names, generator={"elliptic": {"control": "laplace", "max_iterations": 300}}
    )
    first = og.smooth_domain(d)
    assert og.smooth_domain(d) is first  # cached
    og.create_point((9, 9))  # unrelated edit keeps the cache
    assert og.smooth_domain(d) is first
    og.set_dimension("c-cut", 25)
    second = og.smooth_domain(d)
    assert second is not first
    og.undo()
    assert og.smooth_domain(d) is first  # back to the cached inputs
    other = opengrid.replay(og.journal_text())
    assert np.array_equal(other.domain_grid(d), og.domain_grid(d))


def test_smoothing_progress_and_cancel():
    og = opengrid.Session()
    names = _annulus(og, ni=81, nj=41)
    d = og.create_structured_domain(
        connectors=names,
        generator={
            "elliptic": {"control": "laplace", "max_iterations": 10**7, "tolerance": 1e-300}
        },
    )
    progress = _core.SmoothingProgress()
    progress.cancel()
    result = og.smooth_domain(d, progress)
    assert result.cancelled
    assert og.domain_status(d).startswith("smoothing pending")  # cancelled runs are not cached


def test_domains_round_trip_and_replay(tmp_path):
    og = opengrid.Session()
    names = _annulus(og)
    d = og.create_structured_domain(connectors=names, generator={"method": "hermite_tfi"})
    og.rename("c-cut", "cut-connector")
    path = tmp_path / "d.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert to_json(other.store) == to_json(og.store)
    assert np.array_equal(other.domain_grid(d), og.domain_grid(d))
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_generator_validation():
    og = opengrid.Session()
    names = _rectangle(og)
    d = og.create_structured_domain(connectors=names)
    for changes, message in [
        ({"method": "spline"}, "generator.method"),
        ({"hermite_walls": (1,)}, "edges must be 0 or 2"),
        ({"elliptic": {"control": "magic"}}, "control"),
        ({"elliptic": {"omega": 2.5}}, "below 2"),
        ({"elliptic": {"walls": [{"edge": 3}]}}, "edge must be 0 or 2"),
        ({"elliptic": {"attractions": [{"kind": "point", "xi": 1, "eta": 1}]}}, "amplitude"),
    ]:
        with pytest.raises(opengrid.ValidationError, match=message):
            og.set_domain_generator(d, **changes)


def test_quarter_annulus_tfi_exact_through_api():
    og = opengrid.Session()
    og.create_line((1, 0), (3, 0), name="ray0")
    og.create_arc((0, 0), 3.0, 0.0, math.pi / 2, name="outer")
    og.create_line((0, 1), (0, 3), name="ray1")  # same direction as ray0: same distribution
    og.create_arc((0, 0), 1.0, 0.0, math.pi / 2, name="inner")  # used reversed in edge 3
    radial = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.1}
    angular = {"type": "tanh", "begin_spacing": 0.02, "end_spacing": 0.05}
    og.create_connector("ray0", 21, radial, name="c-ray0")
    og.create_connector(
        "outer", 31, {"type": "tanh", "begin_spacing": 0.06, "end_spacing": 0.15}, name="c-outer"
    )
    og.create_connector("ray1", 21, radial, name="c-ray1")
    og.create_connector("inner", 31, angular, name="c-inner")
    # Opposite edges carry the same normalised distributions (spacings scale with the radius).
    d = og.create_domain([["c-ray0"], ["c-outer"], [("c-ray1", True)], [("c-inner", True)]])
    grid = og.domain_grid(d)
    radius = np.hypot(grid[..., 0], grid[..., 1])
    angle = np.arctan2(grid[..., 1], grid[..., 0])
    assert np.allclose(radius, radius[0][None, :], atol=1e-12)
    assert np.allclose(angle, angle[:, 0][:, None], atol=1e-12)
