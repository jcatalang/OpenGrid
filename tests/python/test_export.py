import h5py
import numpy as np
import pytest

import opengrid
from opengrid.io import read_plot3d, write_plot3d

vtk = pytest.importorskip("vtk")
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402


def _two_blocks() -> opengrid.Session:
    """Two rectangles sharing the connector at x = 2."""
    og = opengrid.Session()
    corners = {"a": [(0, 0), (2, 0), (2, 1), (0, 1)], "b": [(2, 0), (5, 0), (5, 1), (2, 1)]}
    og.create_line((2, 0), (2, 1), name="shared")
    og.create_connector("shared", 6, {"type": "tanh", "begin_spacing": 0.05}, name="c-shared")
    for block, (dim) in (("a", 11), ("b", 16)):
        p = corners[block]
        og.create_line(p[0], p[1], name=f"{block}-bottom")
        og.create_line(p[2], p[3], name=f"{block}-top")
        og.create_connector(f"{block}-bottom", dim, name=f"c-{block}-bottom")
        og.create_connector(f"{block}-top", dim, name=f"c-{block}-top")
    og.create_line((0, 1), (0, 0), name="a-left")
    og.create_line((5, 0), (5, 1), name="b-right")
    og.create_connector("a-left", 6, {"type": "tanh", "end_spacing": 0.05}, name="c-a-left")
    og.create_connector("b-right", 6, {"type": "tanh", "begin_spacing": 0.05}, name="c-b-right")
    og.create_domain([["c-a-bottom"], ["c-shared"], ["c-a-top"], ["c-a-left"]], name="block-a")
    og.create_domain(
        [["c-b-bottom"], ["c-b-right"], ["c-b-top"], [("c-shared", True)]], name="block-b"
    )
    return og


def _read_vtkhdf(path):
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def test_vtkhdf_read_back_by_vtk(tmp_path):
    og = _two_blocks()
    path = tmp_path / "blocks.vtkhdf"
    summary = og.export_vtkhdf(path)
    assert summary["points_before_merge"] == 66 + 96
    assert summary["points"] == 66 + 96 - 6  # the shared column is merged
    assert summary["cells"] == 50 + 75
    out = _read_vtkhdf(path)
    assert out.GetNumberOfPoints() == summary["points"]
    assert out.GetNumberOfCells() == 125
    assert out.GetCellType(0) == vtk.VTK_QUAD
    block = vtk_to_numpy(out.GetCellData().GetArray("block_id"))
    assert np.bincount(block).tolist() == [50, 75]
    jac = vtk_to_numpy(out.GetCellData().GetArray("scaled_jacobian"))
    assert np.allclose(jac, 1.0)
    area = vtk_to_numpy(out.GetCellData().GetArray("area"))
    assert area.sum() == pytest.approx(5.0, rel=1e-12)
    names = out.GetFieldData().GetAbstractArray("block_names")
    assert [names.GetValue(i) for i in range(names.GetNumberOfValues())] == ["block-a", "block-b"]
    # Conforming: every interior edge is shared by exactly two cells.
    edges = {}
    for c in range(out.GetNumberOfCells()):
        ids = [out.GetCell(c).GetPointId(k) for k in range(4)]
        for k in range(4):
            edge = tuple(sorted((ids[k], ids[(k + 1) % 4])))
            edges[edge] = edges.get(edge, 0) + 1
    assert max(edges.values()) == 2
    boundary_edges = sum(1 for count in edges.values() if count == 1)
    assert boundary_edges == 2 * (10 + 15) + 2 * 5


def test_vtkhdf_layout(tmp_path):
    og = _two_blocks()
    path = tmp_path / "layout.vtkhdf"
    og.export_vtkhdf(path, domains=["block-a"])
    with h5py.File(path) as f:
        root = f["VTKHDF"]
        assert root.attrs["Type"] == b"UnstructuredGrid"
        assert list(root.attrs["Version"]) == [1, 0]
        for name in (
            "NumberOfPoints",
            "NumberOfCells",
            "NumberOfConnectivityIds",
            "Points",
            "Offsets",
            "Connectivity",
            "Types",
        ):
            assert name in root
        assert root["Points"].dtype == np.float64
        assert set(root["CellData"]) >= {
            "block_id",
            "cell_i",
            "cell_j",
            "skewness",
            "non_orthogonality",
        }


def test_vtkhdf_is_deterministic(tmp_path):
    og = _two_blocks()
    a = tmp_path / "a.vtkhdf"
    b = tmp_path / "b.vtkhdf"
    og.export_vtkhdf(a)
    og.export_vtkhdf(b)
    with h5py.File(a) as fa, h5py.File(b) as fb:
        for name in ("Points", "Connectivity", "CellData/skewness"):
            assert np.array_equal(fa["VTKHDF"][name][()], fb["VTKHDF"][name][()])


@pytest.mark.parametrize("binary", [True, False])
@pytest.mark.parametrize("three_d", [False, True])
def test_plot3d_round_trip(tmp_path, binary, three_d):
    og = _two_blocks()
    path = tmp_path / "grid.xyz"
    og.export_plot3d(path, binary=binary, three_d=three_d)
    grids = read_plot3d(path, three_d=three_d)
    assert len(grids) == 2
    for grid, name in zip(grids, ["block-a", "block-b"], strict=True):
        original = og.domain_grid(name)
        assert grid.shape == original.shape
        if binary:
            assert np.array_equal(grid, original)
        else:
            assert np.allclose(grid, original, rtol=0, atol=0)  # repr() round-trips exactly


def test_plot3d_single_grid_and_errors(tmp_path):
    grid = np.zeros((3, 4, 3))
    grid[..., 0] = np.arange(4)[None, :]
    grid[..., 1] = np.arange(3)[:, None]
    path = tmp_path / "one.x"
    write_plot3d(path, [grid], multi_grid=False)
    assert np.array_equal(read_plot3d(path, multi_grid=False)[0], grid)
    with pytest.raises(ValueError, match="exactly one block"):
        write_plot3d(path, [grid, grid], multi_grid=False)
    with pytest.raises(ValueError, match="nothing"):
        write_plot3d(path, [])
