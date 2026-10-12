"""Phase 8c: mixed 3D models (structured blocks and unstructured volumes) in one CGNS base."""

import shutil
import subprocess

import numpy as np
import pytest
from test_volumes import box

import opengrid
from opengrid import _core

pytest.importorskip("netgen.meshing")


def block_and_volume():
    """A structured block [0, 1]^3 and a tetrahedral volume [1, 2] x [0, 1]^2 sharing the
    structured face x = 1 (the volume puts pyramids on its quads)."""
    og = opengrid.Session()
    left = box(og, (0, 0, 0), (1, 1, 1), 6, "b", "structured")
    right = box(og, (1, 0, 0), (2, 1, 1), 6, "v", "unstructured", shared={(0, 0): left[1]})
    og.create_block(faces=left, name="block")
    og.create_volume(right, generator={"max_size": 0.25}, name="tets")
    og.create_patch("walls", "wall", [f for f in left + right if f != left[1]])
    return og, left[1]


def test_mixed_model_cgns(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    h5py = pytest.importorskip("h5py")
    og, shared = block_and_volume()
    report = og.boundary_report()
    assert report.ok, report.issues
    mesh = og.volume_mesh("tets")
    assert len(mesh["pyramids"]) == 25  # one per quad of the shared 6 x 6 face
    path = tmp_path / "mixed.cgns"
    summary = og.export_cgns(path)
    # 2 zones (both in the default zone family, plus the patch family); 5 block BCs + 5 volume
    # boundary sections; one vertex connection each way.
    assert summary == {"zones": 2, "families": 2, "boundaries": 10, "connections": 2}
    checker = shutil.which("cgnscheck")
    if checker:
        out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
        assert "error" not in out.lower() and "warning" not in out.lower(), out
    with h5py.File(path, "r") as f:
        base = f["OpenGrid"]
        block = base["block"]
        tets = base["tets"]
        assert block.attrs["label"].tobytes().rstrip(b"\0") == b"Zone_t"
        assert base["block/ZoneType/ data"][()].tobytes() == b"Structured"
        assert base["tets/ZoneType/ data"][()].tobytes() == b"Unstructured"

        def coords(zone):
            return np.stack(
                [zone[f"GridCoordinates/Coordinate{c}/ data"][()].reshape(-1) for c in "XYZ"],
                axis=1,
            )

        # Structured coordinates are stored (k, j, i) in HDF5, i fastest.
        grid = coords(block)
        ni = nj = 6
        conn = block[f"ZoneGridConnectivity/{shared}"]
        assert conn["GridConnectivityType/ data"][()].tobytes() == b"Abutting1to1"
        ijk = conn["PointList/ data"][()].reshape(-1, 3) - 1  # (i, j, k) triples
        donor = conn["PointListDonor/ data"][()].reshape(-1) - 1
        assert len(ijk) == ni * nj
        here = grid[ijk[:, 0] + ni * (ijk[:, 1] + nj * ijk[:, 2])]
        there = coords(tets)[donor]
        np.testing.assert_array_equal(here, there)
        np.testing.assert_allclose(here[:, 0], 1.0)
        back = tets[f"ZoneGridConnectivity/{shared}"]
        assert back["PointList/ data"][()].size == ni * nj
        assert back["PointListDonor/ data"][()].size == 3 * ni * nj


def test_mixed_model_vtkhdf_and_gmsh(tmp_path):
    meshio = pytest.importorskip("meshio")
    og, _shared = block_and_volume()
    volume = og.volume_mesh("tets")
    n_cells = 5**3 + len(volume["tets"]) + len(volume["pyramids"])
    # Points: block nodes + volume nodes, the 36 nodes of the shared face counted once.
    n_points = 6**3 + len(volume["points"]) - 36
    og.export_gmsh(tmp_path / "m.msh")
    mesh = meshio.read(tmp_path / "m.msh")
    cells = {
        b.type: len(b.data) for b in mesh.cells if b.type in ("hexahedron", "tetra", "pyramid")
    }
    assert sum(cells.values()) == n_cells and cells["hexahedron"] == 125
    assert len(np.unique(mesh.points, axis=0)) == n_points
    vtk = pytest.importorskip("vtk")
    og.export_vtkhdf(tmp_path / "m.vtkhdf")
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(tmp_path / "m.vtkhdf"))
    reader.Update()
    grid = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    assert grid.GetNumberOfCells() == n_cells
    assert grid.GetNumberOfPoints() == n_points
