"""Phase 11e: curved high-order elements from the structured grids."""

import math

import h5py
import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.validation import ValidationError


def annulus(elements: int, p: int) -> opengrid.Session:
    og = opengrid.Session()
    og.create_arc((0, 0), 1.0, 0.0, math.pi / 2, name="inner")
    og.create_arc((0, 0), 2.0, 0.0, math.pi / 2, name="outer")
    og.create_line((1, 0), (2, 0), name="bottom")
    og.create_line((0, 1), (0, 2), name="left")
    n = elements * p + 1
    connectors = [
        og.create_connector(c, dimension=n, name=f"c-{c}")
        for c in ("inner", "outer", "bottom", "left")
    ]
    og.create_structured_domain(connectors=connectors, name="d")
    return og


def lagrange(nodes, t):
    out = np.ones((len(t), len(nodes)))
    for j in range(len(nodes)):
        for k in range(len(nodes)):
            if k != j:
                out[:, j] *= (t - nodes[k]) / (nodes[j] - nodes[k])
    return out


def arc_error(mesh) -> float:
    """Largest distance from the unit circle of the element edges on the inner arc."""
    p = mesh.order
    basis = lagrange(_core.nodes_1d(mesh.nodes, p + 1), np.linspace(-1, 1, 81))
    worst = 0.0
    for e in mesh.xyz:
        e = e.reshape(p + 1, p + 1, 3)
        for side in (e[0, :], e[-1, :], e[:, 0], e[:, -1]):
            if np.all(np.abs(np.linalg.norm(side[:, :2], axis=1) - 1.0) < 1e-6):
                points = basis @ side
                worst = max(worst, float(np.abs(np.linalg.norm(points[:, :2], axis=1) - 1).max()))
    return worst


def test_node_sets_and_quadrature():
    assert np.allclose(_core.nodes_1d("gll", 4), [-1, -1 / math.sqrt(5), 1 / math.sqrt(5), 1])
    assert np.allclose(_core.weights_1d("gll", 4), [1 / 6, 5 / 6, 5 / 6, 1 / 6])
    assert np.allclose(_core.nodes_1d("gauss", 3), [-math.sqrt(0.6), 0, math.sqrt(0.6)])
    for n in (2, 5, 9):
        x = _core.nodes_1d("gauss", n)
        w = _core.weights_1d("gauss", n)
        assert abs(np.sum(w * x ** (2 * n - 2)) - 2 / (2 * n - 1)) < 1e-13  # exact to 2n - 1
        x = _core.nodes_1d("gll", n)
        w = _core.weights_1d("gll", n)
        assert abs(np.sum(w * x ** (2 * n - 4)) - 2 / (2 * n - 3)) < 1e-13  # exact to 2n - 3


def test_curved_boundary_converges_with_order_and_gll_nodes_lie_on_the_arc():
    errors = {}
    for p in (1, 2, 3, 4):
        coarse = arc_error(annulus(4, p).curved_high_order("d", order=p))
        fine = arc_error(annulus(8, p).curved_high_order("d", order=p))
        errors[p] = coarse
        assert coarse / fine > 2 ** (p + 1) * 0.9  # at least the order h^(p + 1)
    assert errors[1] > 100 * errors[3] > 100 * 100 * errors[4] / 10
    mesh = annulus(4, 3).curved_high_order("d", order=3, nodes="gll")
    assert mesh.snapped > 0
    inner = []
    for e in mesh.xyz:
        side = e.reshape(4, 4, 3)
        for edge in (side[0, :], side[-1, :], side[:, 0], side[:, -1]):
            r = np.linalg.norm(edge[:, :2], axis=1)
            if np.all(np.abs(r - 1) < 1e-3):
                inner.append(r)
    assert inner and np.abs(np.concatenate(inner) - 1).max() < 1e-12
    assert mesh.scaled_jacobian.min() > 0.5


def test_order_must_divide_the_intervals():
    og = annulus(4, 3)  # 12 intervals
    with pytest.raises(ValidationError, match="not a multiple of the order 5"):
        og.curved_high_order("d", order=5)
    with pytest.raises(ValidationError, match="nodes must be one of"):
        og.curved_high_order("d", order=3, nodes="chebyshev")


def test_vtkhdf_lagrange_cells_match_the_element_mapping(tmp_path):
    vtk = pytest.importorskip("vtk")
    og = annulus(4, 3)
    mesh = og.curved_high_order("d", order=3)
    path = tmp_path / "ho.vtkhdf"
    info = og.export_curved_high_order(path, "d", order=3)
    assert info["cells"] == 16
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfCells() == 16
    assert grid.GetCellType(0) == 70  # VTK_LAGRANGE_QUADRILATERAL
    assert grid.GetNumberOfPoints() == (4 * 3 + 1) ** 2
    # One cell evaluated by VTK at (r, s) equals our mapping at (2 r - 1, 2 s - 1).
    cell = grid.GetCell(5)
    r, s = 0.3, 0.65
    x = [0.0, 0.0, 0.0]
    weights = [0.0] * cell.GetNumberOfPoints()
    cell.EvaluateLocation(vtk.reference(0), [r, s, 0.0], x, weights)
    nodes = _core.nodes_1d("equispaced", 4)
    lu = lagrange(nodes, np.array([2 * r - 1]))[0]
    lv = lagrange(nodes, np.array([2 * s - 1]))[0]
    e = mesh.xyz[5].reshape(4, 4, 3)
    ours = np.einsum("b,a,bac->c", lv, lu, e)
    assert np.allclose(x, ours, atol=1e-12)


def test_h5_export_with_gauss_points_and_3d_blocks(tmp_path):
    og = annulus(4, 2)  # 8 intervals
    og.create_extruded_block("d", "translate", vector=(0, 0, 1), dimension=5, name="b")
    mesh = og.curved_high_order("b", order=2, nodes="gauss")
    assert mesh.dim == 3 and mesh.xyz.shape == (4 * 4 * 2, 27, 3)
    assert mesh.scaled_jacobian.min() > 0.5
    # Gauss points lie inside their elements: none on the inner arc or the end planes.
    radius = np.linalg.norm(mesh.xyz[..., :2], axis=-1)
    assert radius.min() > 1.0 and 0.0 < mesh.xyz[..., 2].min() < mesh.xyz[..., 2].max() < 1.0
    path = tmp_path / "ho.h5"
    og.export_curved_high_order(path, "b", order=2, nodes="gauss")
    with h5py.File(path) as f:
        assert f.attrs["node_set"] == "gauss"
        assert f["nodes"].shape == mesh.xyz.shape
        assert list(f["block_names"].asstr()) == ["b"]
    with pytest.raises(ValueError, match="equispaced"):
        og.export_curved_high_order(tmp_path / "ho.vtkhdf", "b", order=2, nodes="gauss")
    # Equispaced hexahedra read back by VTK as Lagrange hexahedra.
    vtk = pytest.importorskip("vtk")
    og.export_curved_high_order(tmp_path / "hex.vtkhdf", "b", order=2)
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(tmp_path / "hex.vtkhdf"))
    reader.Update()
    assert reader.GetOutput().GetCellType(0) == 72  # VTK_LAGRANGE_HEXAHEDRON
    assert reader.GetOutput().GetNumberOfPoints() == 9 * 9 * 5


def test_gui_exports_high_order_elements(qtbot, tmp_path):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True, session=annulus(4, 2))
    qtbot.addWidget(window)
    summary = window.export_curved_high_order(str(tmp_path / "ho.vtkhdf"), order=2)
    assert summary["cells"] == 16
    assert window.export_curved_high_order(str(tmp_path / "ho.h5"), order=3) is None  # 8 intervals
    window.close()
