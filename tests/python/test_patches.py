import importlib.util
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

import opengrid
from opengrid import _core
from opengrid.project import to_json

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def step_session() -> opengrid.Session:
    og = opengrid.Session()
    _load("03_backward_step").build(og)
    og.create_patch(
        "walls", "wall", ["c-in-bottom", "c-in-top", "c-step", "c-lo-bottom", "c-up-top"]
    )
    og.create_patch("inflow", "inlet", ["c-inlet"])
    og.create_patch("outflow", "outlet", ["c-out-lo", "c-out-up"])
    og.create_zone("channel", ["inlet-block"])
    og.create_zone("expansion", ["lower-block", "upper-block"])
    return og


def airfoil_session() -> opengrid.Session:
    og = opengrid.Session()
    demo = _load("03_airfoil_grids")
    demo.o_grid(og)
    demo.c_grid(og)
    for name in ("o-grid", "c-grid"):
        og.set_domain_generator(name, elliptic=None)  # TFI keeps the tests fast
    og.create_patch("o-airfoil", "wall", ["o-c-upper", "o-c-lower"])
    og.create_patch("o-far", "farfield", ["o-c-farfield"])
    og.create_patch("c-airfoil", "wall", ["c-c-upper", "c-c-lower"])
    og.create_patch("c-far", "farfield", ["c-c-top", "c-c-front", "c-c-bottom"])
    og.create_patch("c-out", "outlet", ["c-c-outflow-top", "c-c-outflow-bottom"])
    return og


def test_boundary_report_and_interfaces():
    og = step_session()
    report = og.boundary_report()
    assert report.ok, report.issues
    assert sorted(report.interfaces) == [
        "c-interface (inlet-block | upper-block)",
        "c-mid (lower-block | upper-block)",
    ]
    af = airfoil_session()
    report = af.boundary_report()
    assert report.ok, report.issues
    assert sorted(report.interfaces) == ["c-c-wake (cut in c-grid)", "o-c-cut (cut in o-grid)"]


def test_validator_messages():
    og = step_session()
    og.remove_from_patch("walls", "c-step")
    og.add_to_patch("inflow", "c-mid")
    og.add_to_patch("outflow", "c-out-up")  # already there: no-op
    og.add_to_patch("walls", "c-out-lo")  # now in two patches
    messages = [str(i) for i in og.boundary_report().issues]
    assert "c-step: lies on the boundary of lower-block but belongs to no patch" in messages
    assert any(
        m.startswith("c-mid: is an interface (lower-block and upper-block)") for m in messages
    )
    assert "c-out-lo: belongs to several patches: walls, outflow" in messages
    assert og.unassigned_boundary_connectors() == ["c-step"]
    with pytest.raises(opengrid.ValidationError, match="belongs to no patch"):
        og.export_gmsh("unused.msh")
    og.create_zone("again", ["lower-block"])
    assert any("belongs to several zones" in str(i) for i in og.boundary_report().issues)


def test_non_manifold_connector_reported():
    og = opengrid.Session()
    p = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for k in range(4):
        og.create_line(p[k], p[(k + 1) % 4], name=f"l{k}")
        og.create_connector(f"l{k}", 5, name=f"c{k}")
    og.create_domain([["c0"], ["c1"], ["c2"], ["c3"]], name="a")
    og.create_domain([["c0"], ["c1"], ["c2"], ["c3"]], name="b")
    og.create_line((1, 0), (2, 0), name="m0")
    og.create_line((2, 0), (2, 1), name="m1")
    og.create_line((2, 1), (1, 1), name="m2")
    for k in range(3):
        og.create_connector(f"m{k}", 5, name=f"n{k}")
    og.create_domain([["n0"], ["n1"], ["n2"], [("c1", True)]], name="c")
    assert any("c1: is used 3 times" in str(i) for i in og.boundary_report().issues)


def test_patch_and_zone_entities():
    og = step_session()
    with pytest.raises(opengrid.ValidationError, match="used by patch 'walls'"):
        og.delete("c-step")
    with pytest.raises(opengrid.ValidationError, match="is a domain, not a connector"):
        og.create_patch("bad", "wall", ["lower-block"])
    with pytest.raises(opengrid.ValidationError, match="type"):
        og.create_patch("bad", "slip", ["c-step"])
    with pytest.raises(opengrid.ValidationError, match="does not contain"):
        og.remove_from_patch("walls", "c-mid")
    og.rename("c-step", "step-face")
    assert og.boundary_report().ok
    og.set_patch_type("walls", "symmetry")
    og.undo()
    assert og.entity("walls").params["patch_type"] == "wall"


def test_round_trip_and_replay(tmp_path):
    og = step_session()
    og.save(tmp_path / "s.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "s.ogp")
    assert to_json(other.store) == to_json(og.store)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    assert '"patches": [' in (tmp_path / "s.ogp").read_text()


def _patch_face_counts(og):
    counts = {}
    for patch in og.store.entities("patch"):
        counts[patch.name] = sum(
            og.store.get(m["connector"]).params["dimension"] - 1 for m in patch.params["members"]
        )
    return counts


def test_gmsh_physical_groups(tmp_path):
    meshio = pytest.importorskip("meshio")
    og = step_session()
    path = tmp_path / "step.msh"
    og.export_gmsh(path)
    mesh = meshio.read(path)
    groups = {name: (int(tag), int(dim)) for name, (tag, dim) in mesh.field_data.items()}
    assert groups == {
        "walls": (1, 1),
        "inflow": (2, 1),
        "outflow": (3, 1),
        "channel": (1, 2),
        "expansion": (2, 2),
    }
    physical = mesh.cell_data["gmsh:physical"]
    counted = {}
    for block, tags in zip(mesh.cells, physical, strict=True):
        for tag in np.unique(tags):
            dim = 1
            if block.type == "quad":
                dim = 2
            name = next(n for n, (t, d) in groups.items() if t == tag and d == dim)
            counted[name] = counted.get(name, 0) + int(np.count_nonzero(tags == tag))
    expected = _patch_face_counts(og)
    expected["channel"] = 40 * 50
    expected["expansion"] = 2 * 160 * 50
    assert counted == expected
    # Conforming: merged nodes.
    assert len(mesh.points) == 41 * 51 + 2 * 161 * 51 - 51 - 161


def test_vtkhdf_collection(tmp_path):
    vtk = pytest.importorskip("vtk")
    og = step_session()
    path = tmp_path / "step.vtkhdf"
    summary = og.export_vtkhdf(path)
    assert summary["blocks"] == ["channel", "expansion", "walls", "inflow", "outflow"]
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    pdc = reader.GetOutputDataObject(0)
    assert pdc.IsA("vtkPartitionedDataSetCollection")
    names = [
        pdc.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
        for i in range(pdc.GetNumberOfPartitionedDataSets())
    ]
    assert names == summary["blocks"]
    counts = {
        name: pdc.GetPartitionedDataSet(i).GetPartition(0).GetNumberOfCells()
        for i, name in enumerate(names)
    }
    expected = _patch_face_counts(og)
    expected["channel"] = 40 * 50
    expected["expansion"] = 2 * 160 * 50
    assert counts == expected
    assembly = pdc.GetDataAssembly()
    patches = assembly.FindFirstNodeWithName("Patches")
    assert sorted(assembly.GetDataSetIndices(patches)) == [2, 3, 4]
    with h5py.File(path) as f:
        assert f["VTKHDF"].attrs["Type"] == b"PartitionedDataSetCollection"
        zone_ids = f["VTKHDF/expansion/CellData/zone_id"][()]
        assert set(zone_ids.tolist()) == {1}


def _cgns_children(group, label):
    return {
        name: node
        for name, node in group.items()
        if isinstance(node, h5py.Group) and node.attrs.get("label") == label.encode()
    }


def test_cgns_structure_and_cgnscheck(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    for og, name in ((step_session(), "step"), (airfoil_session(), "airfoil")):
        path = tmp_path / f"{name}.cgns"
        summary = og.export_cgns(path)
        checker = shutil.which("cgnscheck")
        if checker is not None:
            out = subprocess.run(
                [checker, str(path)], capture_output=True, text=True, check=False
            ).stdout
            assert "ERROR" not in out.upper()
            assert "warning" not in out.lower(), out
        with h5py.File(path) as f:
            base = f["OpenGrid"]
            zones = _cgns_children(base, "Zone_t")
            assert len(zones) == summary["zones"]
            families = _cgns_children(base, "Family_t")
            expected = {p.name for p in og.store.entities("patch")} | {
                z.name for z in og.store.entities("zone")
            }
            if name == "airfoil":
                expected.add("default")  # unzoned domains
            assert set(families) == expected
            transforms = {}
            bc_count = 0
            for zone_name, zone in zones.items():
                for conn_name, conn in _cgns_children(
                    zone.get("ZoneGridConnectivity", {}), "GridConnectivity1to1_t"
                ).items():
                    transforms[(zone_name, conn_name)] = tuple(
                        conn["Transform"][" data"][()].tolist()
                    )
                bc_count += len(_cgns_children(zone.get("ZoneBC", {}), "BC_t"))
            assert bc_count == summary["boundaries"]
        if name == "step":
            assert set(transforms.values()) == {(1, 2)}
            assert len(transforms) == 4  # two interfaces, both directions
        else:
            assert transforms[("o-grid", "o-c-cut")] == (1, 2)
            assert transforms[("c-grid", "c-c-wake")] == (-1, -2)
            assert len(transforms) == 4


def test_cgns_point_ranges_match_grid(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    og = airfoil_session()
    path = tmp_path / "a.cgns"
    og.export_cgns(path)
    grid = og.domain_grid("c-grid")
    with h5py.File(path) as f:
        zone = f["OpenGrid/c-grid"]
        conns = _cgns_children(zone["ZoneGridConnectivity"], "GridConnectivity1to1_t")
        for conn in conns.values():
            r = conn["PointRange"][" data"][()].reshape(-1).tolist()
            d = conn["PointRangeDonor"][" data"][()].reshape(-1).tolist()
            # Begin and end points of the range coincide with their donor images.
            assert np.array_equal(grid[r[1] - 1, r[0] - 1], grid[d[1] - 1, d[0] - 1])
            assert np.array_equal(grid[r[3] - 1, r[2] - 1], grid[d[3] - 1, d[2] - 1])
        wall = _cgns_children(zone["ZoneBC"], "BC_t")
        assert {n.split("_")[0] for n in wall} == {"c-airfoil", "c-far", "c-out"}


def test_cgns_name_length_limit(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    og = step_session()
    og.rename("walls", "a-very-long-patch-name-that-exceeds-cgns")
    with pytest.raises(opengrid.ValidationError, match="longer than 32 characters"):
        og.export_cgns(tmp_path / "x.cgns")
