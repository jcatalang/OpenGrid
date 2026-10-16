"""Phase 12c: rigorous validity (Bernstein-Bezier bounds of the Jacobian) and the repair of
curved elements."""

import numpy as np
import pytest
from test_curved_elements import ring, sphere_in_box

import opengrid
from opengrid.io import lagrange
from opengrid.model import repair, shape


@pytest.mark.parametrize("kind", ["triangle", "quad", "tetra", "hexa", "prism"])
def test_bezier_coefficients_of_affine_elements_are_the_constant_jacobian(kind):
    rng = np.random.default_rng(1)
    d = lagrange.DIM[kind]
    for p in (1, 2, 3, 4):
        nodes = lagrange.as_array(lagrange.nodes("gmsh", kind, p))
        matrix = np.eye(3)
        matrix[:d, :d] += 0.3 * rng.standard_normal((d, d))
        xyz = np.zeros((3, len(nodes), 3))
        xyz[:, :, :d] = nodes
        xyz = xyz @ matrix.T + rng.standard_normal(3)
        degrees = repair.jacobian_degree(kind, p)
        points, _ = repair._lattice(kind, degrees)
        det = shape.jacobians(kind, p, xyz, d, points)
        coefficients = repair.bezier_coefficients(kind, degrees, det)
        assert np.ptp(coefficients) < 1e-6 * np.abs(det).max()  # degree 11 at hexa p = 4
        status, bound = repair.validity(kind, p, xyz, d)
        assert (status == 1).all() and np.allclose(bound, 1.0, atol=1e-6)
        assert repair.rigorous(kind, p)


def test_bounds_reject_pushed_and_folded_quads_and_accept_curved_ones():
    grid = np.array([[0.5 * i, 0.5 * j, 0.0] for j in range(3) for i in range(3)])
    quad = grid[[0, 2, 8, 6, 1, 5, 7, 3, 4]]  # Gmsh quad9: corners, edges, centre
    for shift, expected in ((0.0, 1), (0.9, -1), (1.2, -1)):
        x = quad.copy()
        x[8] = [0.5 + 0.5 * shift, 0.5 + 0.5 * shift, 0.0]
        status, bound = repair.validity("quad", 2, x[None], 2)
        assert status[0] == expected
        if expected == 1:
            assert bound[0] == pytest.approx(1.0)
    # A genuinely curved mesh: every element provably valid, the bound below the samples.
    mesh = ring(16, first=1e-5, radial=60).high_order(4)
    conn = mesh.cells["quad"]
    status, bound = repair.validity("quad", 4, mesh.points[conn], 2)
    lattice, _ = repair._lattice("quad", repair.jacobian_degree("quad", 4))
    on_lattice = shape.jacobians("quad", 4, mesh.points[conn], 2, lattice)
    sampled = shape.jacobians("quad", 4, mesh.points[conn], 2)
    assert (status == 1).all()
    # The bound (relative to the lattice maximum) is below every sampled value of J.
    assert (bound * on_lattice.max(axis=1) <= sampled.min(axis=1) + 1e-9).all()
    assert bound.min() > 0.8


def test_pyramids_and_high_order_hexahedra_are_sampled():
    assert not repair.rigorous("pyramid", 2)
    assert repair.rigorous("hexa", 4) and not repair.rigorous("hexa", 6)
    assert repair.rigorous("tetra", 6)
    nodes = lagrange.as_array(lagrange.nodes("gmsh", "pyramid", 2))
    status, bound = repair.validity("pyramid", 2, nodes[None], 3)
    assert status[0] == 1 and bound[0] == pytest.approx(1.0)


def test_repair_recovers_quality_on_curved_layers():
    pytest.importorskip("netgen.meshing")
    og, _radius = sphere_in_box("unstructured")
    plain = og.high_order(3, repair=False)
    repaired = og.high_order(3)
    before = plain.curving["min_scaled_jacobian"]
    after = repaired.curving["min_scaled_jacobian"]
    report = repaired.curving["repair"]
    assert report["kept"] and report["elements"] >= 1 and report["free_nodes"] > 0
    assert after > before + 0.03
    assert report["min_scaled_jacobian_after"] > report["min_scaled_jacobian_before"]
    # The linear vertices and the geometry nodes did not move.
    assert np.array_equal(plain.points[: plain.vertices], repaired.points[: repaired.vertices])
    assert repaired.curving["validity"]["undetermined"] == 0
    for kind, conn in repaired.cells.items():
        status, _bound = repair.validity(kind, 3, repaired.points[conn], 3)
        assert (status == 1).all(), kind


def test_nothing_to_repair_costs_nothing_and_hopeless_folds_are_still_reported():
    og = ring(16, first=1e-5, radial=60)
    mesh = og.high_order(4)
    assert "repair" not in mesh.curving  # no element flagged
    import test_groups3d as g3

    block = g3.sphere_block(7)
    block.create_patch("wall", "wall", list(block.names("domain")))
    with pytest.raises(opengrid.ValidationError, match="folds after curving"):
        block.high_order(4)  # the corner cells cannot be untangled: reported, with the bound
