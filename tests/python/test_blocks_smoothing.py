"""Phase 7b: 3D elliptic smoothing of blocks, smoothing on surfaces, hyperbolic extrusion."""

import math
import shutil
import subprocess

import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.project import to_json


def cylinder_ogrid(og=None, n_theta=81, n_radial=41, **extra) -> opengrid.Session:
    """Wall r = 0.5, z in [0, 1]; hyperbolic extrusion to r = 10.5."""
    if og is None:
        og = opengrid.Session()
    og.create_circle((0, 0, 0), 0.5, name="circle0")
    og.create_circle((0, 0, 1), 0.5, name="circle1")
    og.create_line((0.5, 0, 0), (0.5, 0, 1), name="seam")
    og.create_connector("circle0", n_theta, name="c-circle0")
    og.create_connector("circle1", n_theta, name="c-circle1")
    og.create_connector("seam", 11, name="c-seam")
    og.create_domain(
        [["c-circle0"], ["c-seam"], [("c-circle1", True)], [("c-seam", True)]], name="wall"
    )
    og.create_extruded_block(
        "wall",
        "hyperbolic",
        distance=10.0,
        dimension=n_radial,
        distribution={"type": "tanh", "begin_spacing": 1e-4},
        edges=("periodic", "periodic", "constant_z", "constant_z"),
        name="ogrid",
        **extra,
    )
    og.create_patch("p-wall", "wall", ["wall"])
    og.create_patch("p-far", "farfield", [{"block": "ogrid", "face": "kmax"}])
    og.create_patch(
        "p-ends",
        "symmetry",
        [{"block": "ogrid", "face": "jmin"}, {"block": "ogrid", "face": "jmax"}],
    )
    return og


def test_hyperbolic_cylinder_ogrid(tmp_path):
    og = cylinder_ogrid()
    grid = og.block_grid("ogrid")
    nk, nj, ni, _ = grid.shape
    assert (ni, nj, nk) == (81, 11, 41)
    fractions, _w = _core.distribute(41, "tanh", 1e-4, None, None)
    r = np.hypot(grid[..., 0], grid[..., 1])
    for k in range(nk):
        assert np.allclose(r[k], 0.5 + 10.0 * fractions[k], rtol=1e-9)
    assert np.allclose(grid[:, :, :, 2], grid[0, :, :, 2][None])  # z constant along k
    first = r[1] - r[0]
    assert np.allclose(first, 1e-3, rtol=1e-6)  # 1e-4 of the 10 distance
    assert og.block_quality("ogrid")["scaled_jacobian"].min() > 0.99
    report = og.boundary_report()
    assert report.ok, report.issues
    assert report.interfaces == ["ogrid.imin | ogrid.imax"]
    assert "hyperbolic extrusion" in og.block_status("ogrid")
    if _core.cgns_available():
        path = tmp_path / "o.cgns"
        assert og.export_cgns(path)["connections"] == 2
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
            assert "error" not in out.lower() and "warning" not in out.lower(), out
    replayed = opengrid.replay(og.journal_text())
    assert np.array_equal(replayed.block_grid("ogrid"), grid)


def bump_channel() -> opengrid.Session:
    """A domain in the plane y = 0 over a bump (x in [0, 3], z from the bump to 1), extruded
    along y; the bump is the block's jmin face."""
    og = opengrid.Session()
    xs = np.linspace(0.0, 3.0, 31)
    bump = [(x, 0.0, 0.15 * math.exp(-8.0 * (x - 1.5) ** 2)) for x in xs]
    og.create_spline(bump, name="bump")
    og.create_line((3, 0, bump[-1][2]), (3, 0, 1), name="outlet")
    og.create_line((3, 0, 1), (0, 0, 1), name="top")
    og.create_line((0, 0, 1), (0, 0, bump[0][2]), name="inlet")
    og.create_connector("bump", 41, name="c-bump")
    og.create_connector("outlet", 21, name="c-outlet")
    og.create_connector("top", 41, name="c-top")
    og.create_connector("inlet", 21, name="c-inlet")
    og.create_domain([["c-bump"], ["c-outlet"], ["c-top"], ["c-inlet"]], name="side")
    og.create_extruded_block("side", "translate", vector=(0, 1, 0), dimension=9, name="channel")
    return og


def wall_metrics(grid, spacing):
    """Largest relative first-spacing error and angle (degrees) off the wall normal on the
    interior of the jmin face."""
    worst_spacing = 0.0
    worst_angle = 0.0
    nk, _nj, ni, _ = grid.shape
    for k in range(1, nk - 1):
        for i in range(1, ni - 1):
            line = grid[k, 1, i] - grid[k, 0, i]
            ti = grid[k, 0, i + 1] - grid[k, 0, i - 1]
            tk = grid[k + 1, 0, i] - grid[k - 1, 0, i]
            n = np.cross(ti, tk)
            worst_spacing = max(worst_spacing, abs(np.linalg.norm(line) / spacing - 1.0))
            c = abs(line @ n) / (np.linalg.norm(line) * np.linalg.norm(n))
            worst_angle = max(worst_angle, math.degrees(math.acos(min(1.0, c))))
    return worst_spacing, worst_angle


def test_block_smoothing_with_wall_control():
    og = bump_channel()
    plain = og.block_grid("channel")
    _s, angle_before = wall_metrics(plain, 0.01)
    og.set_block_smoothing(
        "channel",
        {
            "control": "thomas_middlecoff",
            "walls": [{"face": "jmin", "spacing": 0.01, "decay": 0.3}],
            "max_iterations": 4000,
            "tolerance": 1e-9,
        },
    )
    grid = og.block_grid("channel")
    spacing_error, angle = wall_metrics(grid, 0.01)
    assert spacing_error < 0.02
    assert angle < 1.5 < angle_before
    assert og.block_quality("channel")["scaled_jacobian"].min() > 0
    assert "elliptic smoothing converged" in og.block_status("channel")
    # Faces are unchanged.
    assert np.array_equal(grid[:, 0], plain[:, 0])
    assert np.array_equal(grid[0], plain[0])
    og.undo()
    assert og.entity("channel").params["smoothing"] is None
    with pytest.raises(opengrid.ValidationError, match="face"):
        og.set_block_smoothing("channel", {"walls": [{"face": "bottom"}]})
    with pytest.raises(opengrid.ValidationError, match="control"):
        og.set_block_smoothing("channel", {"control": "ttm"})


def test_surface_domain_smoothing():
    og = opengrid.Session()
    t = math.pi / 2
    og.create_arc((0, 0, 0), 2.0, 0.0, t, name="a0")
    og.create_arc((0, 0, 3), 2.0, 0.0, t, name="a1")
    og.create_line((2, 0, 0), (2, 0, 3), name="l0")
    og.create_line((0, 2, 0), (0, 2, 3), name="l1")
    og.create_connector("a0", 31, {"type": "tanh", "begin_spacing": 0.005}, name="c-a0")
    og.create_connector("a1", 31, name="c-a1")
    og.create_connector("l0", 21, name="c-l0")
    og.create_connector("l1", 21, {"type": "tanh", "begin_spacing": 0.01}, name="c-l1")
    og.create_extruded_surface("a0", (0, 0, 3), name="cyl")
    og.create_domain(
        [["c-a0"], ["c-l1"], [("c-a1", True)], [("c-l0", True)]],
        name="patch",
        surface="cyl",
        generator={"method": "tfi", "elliptic": {"control": "thomas_middlecoff"}},
    )
    result = og.smooth_domain("patch")
    assert result.converged
    grid = og.domain_grid("patch")
    assert np.allclose(np.hypot(grid[..., 0], grid[..., 1]), 2.0, atol=1e-10)
    with pytest.raises(opengrid.ValidationError, match="planar"):
        og.set_domain_generator(
            "patch", elliptic={"control": "thomas_middlecoff", "walls": [{"edge": 0}]}
        )


def test_validation_and_round_trip(tmp_path):
    og = cylinder_ogrid(n_theta=41, n_radial=11)
    with pytest.raises(opengrid.ValidationError, match="edges"):
        og.create_extruded_block("wall", "hyperbolic", distance=1.0, edges=("free", "x"))
    with pytest.raises(opengrid.ValidationError, match="must not be zero"):
        og.create_extruded_block("wall", "hyperbolic")
    og.set_block_smoothing("ogrid", {"control": "laplace", "max_iterations": 50})
    path = tmp_path / "o.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert to_json(other.store) == to_json(og.store)
    assert np.array_equal(other.block_grid("ogrid"), og.block_grid("ogrid"))
