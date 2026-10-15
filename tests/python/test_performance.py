"""Phase 9c: parallel hot loops (deterministic for any thread count), accelerated elliptic
smoothing, the C++ volume-quality metrics and the Netgen backend's bulk transfer."""

import hashlib
import math

import numpy as np
import pytest
from quality_reference import volume_quality_reference
from test_inflation3d import body_in_box
from test_volumes import box, box_with_body

import opengrid
from opengrid.validation import ValidationError


@pytest.fixture
def threads():
    """Restores the thread count after a test changes it."""
    before = opengrid.threads()
    yield
    opengrid.set_threads(before)


def digest(*arrays) -> str:
    h = hashlib.sha1()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def cylinder(og, ni=121, nj=61, acceleration=5):
    og.create_circle((0, 0), 0.5, name="cyl")
    og.create_circle((0, 0), 20.0, name="far")
    og.create_line((0.5, 0), (20, 0), name="cut")
    og.create_connector("cyl", ni, name="c-cyl")
    og.create_connector("far", ni, name="c-far")
    og.create_connector("cut", nj, {"type": "tanh", "begin_spacing": 1e-3}, name="c-cut")
    og.create_domain([[("c-cyl", True)], ["c-cut"], ["c-far"], [("c-cut", True)]], name="d")
    og.set_domain_generator(
        "d",
        elliptic={
            "walls": [{"edge": 0, "spacing": 1e-3}],
            "acceleration": acceleration,
            "max_iterations": 20000,
        },
    )


def test_thread_control():
    assert opengrid.threads() >= 1
    with pytest.raises(ValueError):
        opengrid.set_threads(0)
    if opengrid._core.openmp_enabled():
        assert opengrid.threads() >= 1


@pytest.mark.parametrize("solver", ["line_sor", "point_sor"])
def test_elliptic_2d_is_independent_of_the_thread_count(threads, solver):
    results = []
    for n in (1, 3, 8):
        opengrid.set_threads(n)
        og = opengrid.Session()
        cylinder(og, ni=81, nj=41)
        og.set_domain_generator("d", elliptic={"solver": solver, "max_iterations": 300})
        r = og.smooth_domain("d")
        results.append((r.iterations, digest(r.grid, r.residuals)))
    assert results[0] == results[1] == results[2]


def test_anderson_acceleration_reaches_the_same_grid():
    plain = opengrid.Session()
    cylinder(plain, acceleration=0)
    plain.set_domain_generator("d", elliptic={"omega": 1.5, "tolerance": 1e-10})
    fast = opengrid.Session()
    cylinder(fast)
    a = plain.smooth_domain("d")
    b = fast.smooth_domain("d")
    assert a.converged and b.converged
    assert b.iterations * 5 < a.iterations
    g = a.grid
    radial = np.linalg.norm(np.diff(g, axis=0), axis=2)  # (nj - 1, ni)
    gap = np.linalg.norm(b.grid - g, axis=2)[1:]
    assert (gap / radial).max() < 1e-3
    spacing = np.linalg.norm(b.grid[1] - b.grid[0], axis=1)
    assert np.allclose(spacing, 1e-3, rtol=1e-5)


def test_acceleration_is_validated():
    og = opengrid.Session()
    cylinder(og)
    for bad in (-1, 51, 2.5, True):
        with pytest.raises(ValidationError):
            og.set_domain_generator("d", elliptic={"acceleration": bad})


def test_elliptic_3d_is_independent_of_the_thread_count(threads):
    results = []
    for n in (1, 4):
        opengrid.set_threads(n)
        og = opengrid.Session()
        faces = box(og, (0, 0, 0), (2, 1, 1), 9, "o", "structured")
        # Bulge one face so the smoother has something to do.
        og.create_block(
            faces=faces,
            name="b",
            smoothing={"control": "laplace", "walls": [{"face": "kmin", "spacing": 0.02}]},
        )
        grid = og.block_grid("b")
        status = og.block_status("b")
        results.append((status, digest(grid)))
    assert results[0] == results[1]


def test_volume_quality_matches_the_reference():
    mixed = 0
    for og, name in ((box_with_body(), "fluid"), (body_in_box("structured", 7, 4)[0], "fluid")):
        mesh = og.volume_mesh(name)
        mixed += sum(len(mesh[k]) for k in ("pyramids", "prisms", "hexes"))
        fast = og.volume_quality(name)
        reference = volume_quality_reference(mesh)
        assert set(fast) == set(reference)
        for key in reference:
            assert fast[key].shape == reference[key].shape
            # Angles near 0 come from arccos near 1, where a rounding of 1e-16 is ~1e-6 degrees.
            atol = 1e-9
            if key in ("min_angle", "max_angle", "non_orthogonality"):
                atol = 1e-5
            np.testing.assert_allclose(fast[key], reference[key], rtol=1e-9, atol=atol, err_msg=key)
    assert mixed > 0


def test_volume_mesh_is_independent_of_the_thread_count(threads):
    results = []
    for n in (1, 6):
        opengrid.set_threads(n)
        og = box_with_body()
        mesh = og.volume_mesh("fluid")
        q = og.volume_quality("fluid")
        results.append(digest(mesh["points"], mesh["tets"], mesh["pyramids"], *q.values()))
    assert results[0] == results[1]


def test_netgen_bulk_transfer_keeps_the_boundary():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (1, 1, 1), 9, "o", "unstructured")
    og.create_volume(faces, generator={"max_size": 0.15}, name="v")
    mesh = og.volume_mesh("v")
    assert len(mesh["tets"]) > 500
    assert abs(og.volume_quality("v")["volume"].sum() - 1.0) < 1e-12
    assert math.isclose(mesh["enclosed_volume"], 1.0, rel_tol=1e-12)
