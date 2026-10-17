import math
import re

import h5py
import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read

import opengrid
from opengrid import _core
from opengrid.project import to_json


def channel() -> opengrid.Session:
    """Two blocks stacked in y on [0, 2] x [0, 1]; the right side is periodic with the left
    (translation (2, 0)); one right connector runs downwards (opposite direction)."""
    og = opengrid.Session()
    lines = {
        "b": ((0, 0), (2, 0)),
        "m": ((0, 0.5), (2, 0.5)),
        "t": ((0, 1), (2, 1)),
        "l0": ((0, 0), (0, 0.5)),
        "l1": ((0, 0.5), (0, 1)),
        "r0": ((2, 0), (2, 0.5)),
        "r1": ((2, 1), (2, 0.5)),
    }
    for name, (a, b) in lines.items():
        og.create_line(a, b, name=name)
    spacing = {"type": "tanh", "begin_spacing": 0.02}
    for name in ("b", "m", "t"):
        og.create_connector(name, 21, name=f"c-{name}")
    og.create_connector("l0", 11, spacing, name="c-l0")
    og.create_connector("r0", 11, spacing, name="c-r0")
    og.create_connector("l1", 11, {"type": "tanh", "end_spacing": 0.02}, name="c-l1")
    og.create_connector("r1", 11, {"type": "tanh", "begin_spacing": 0.02}, name="c-r1")  # downwards
    og.create_domain([["c-b"], ["c-r0"], [("c-m", True)], [("c-l0", True)]], name="lower")
    og.create_domain([["c-m"], [("c-r1", True)], [("c-t", True)], [("c-l1", True)]], name="upper")
    og.create_patch("walls", "wall", ["c-b", "c-t"])
    og.create_patch("left", "periodic", ["c-l0", "c-l1"])
    og.create_patch("right", "periodic", ["c-r1", "c-r0"])
    return og


def sector() -> opengrid.Session:
    """Annular sector 0..60 degrees; the ray at 60 degrees is the rotated image of the ray at 0."""
    og = opengrid.Session()
    angle = math.pi / 3
    og.create_line((1, 0), (3, 0), name="ray0")
    og.create_arc((0, 0), 3.0, 0.0, angle, name="outer")
    og.create_line(
        (math.cos(angle), math.sin(angle)), (3 * math.cos(angle), 3 * math.sin(angle)), name="ray1"
    )
    og.create_arc((0, 0), 1.0, 0.0, angle, name="inner")
    og.create_connector("ray0", 11, name="c-ray0")
    og.create_connector("outer", 16, name="c-outer")
    og.create_connector("ray1", 11, name="c-ray1")
    og.create_connector("inner", 16, name="c-inner")
    og.create_domain(
        [["c-ray0"], ["c-outer"], [("c-ray1", True)], [("c-inner", True)]], name="sector"
    )
    og.create_patch("walls", "wall", ["c-inner", "c-outer"])
    og.create_patch("side0", "periodic", ["c-ray0"])
    og.create_patch("side1", "periodic", ["c-ray1"])
    og.pair_periodic("side0", "side1", translation=(0, 0, 0), rotation_angle=angle)
    return og


def test_pairing_validation():
    og = channel()
    report = og.boundary_report()
    assert not report.ok
    messages = [str(i) for i in report.issues]
    assert "left: periodic patch has no partner (use pair_periodic)" in messages
    og.pair_periodic("left", "right")
    params = og.entity("left").params["periodic"]
    assert params["translation"] == pytest.approx((2.0, 0.0, 0.0))
    report = og.boundary_report()
    assert report.ok, report.issues
    ((source, target, pairs),) = report.periodic
    assert (source.name, target.name) == ("left", "right")
    directions = {
        og.store.get(p.source).name: (og.store.get(p.target).name, p.same_direction) for p in pairs
    }
    assert directions == {"c-l0": ("c-r0", True), "c-l1": ("c-r1", False)}


def test_pairing_errors():
    og = channel()
    og.pair_periodic("left", "right", translation=(2.1, 0, 0))
    assert any(
        "c-l0 has no periodic image in right: its ends land 0.1 away from the ends of c-r0: "
        "the two sides are not images of each other under the translation" in str(i)
        for i in og.boundary_report().issues
    )
    og.pair_periodic("left", "right")
    og.pair_periodic("right", "left")
    assert any("from one side only" in str(i) for i in og.boundary_report().issues)
    og.unpair_periodic("right")
    with pytest.raises(opengrid.ValidationError, match="used by patch 'left'"):
        og.delete("right")
    og.set_patch_type("right", "wall")
    assert any(
        "partner right is a wall patch, not periodic" in str(i) for i in og.boundary_report().issues
    )
    with pytest.raises(opengrid.ValidationError, match="only periodic patches can be paired"):
        og.pair_periodic("walls", "left")


def test_periodic_round_trip_and_replay(tmp_path):
    og = channel()
    og.pair_periodic("left", "right")
    og.save(tmp_path / "p.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "p.ogp")
    assert to_json(other.store) == to_json(og.store)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    other.rename("right", "east")
    assert other.boundary_report().ok


def _cgns_nodes(group, label):
    return {
        n: g
        for n, g in group.items()
        if isinstance(g, h5py.Group) and g.attrs.get("label") == label.encode()
    }


def test_cgns_periodic(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    og = channel()
    og.pair_periodic("left", "right")
    path = tmp_path / "c.cgns"
    summary = og.export_cgns(path)
    assert summary["boundaries"] == 2  # walls only: periodic faces are connectivity
    assert summary["connections"] == 2 + 4  # c-m interface both ways + 2 periodic pairs both ways
    import shutil
    import subprocess

    checker = shutil.which("cgnscheck")
    if checker:
        out = subprocess.run(
            [checker, str(path)], capture_output=True, text=True, check=False
        ).stdout
        assert "warning" not in out.lower() and "error" not in out.lower(), out
    grids = {name: og.domain_grid(name) for name in ("lower", "upper")}
    found = 0
    with h5py.File(path) as f:
        for zone_name, zone in _cgns_nodes(f["OpenGrid"], "Zone_t").items():
            for _name, conn in _cgns_nodes(
                zone["ZoneGridConnectivity"], "GridConnectivity1to1_t"
            ).items():
                if "Periodic" not in conn.get("GridConnectivityProperty", {}):
                    continue
                found += 1
                translation = conn["GridConnectivityProperty/Periodic/Translation"][" data"][()]
                donor = conn[" data"][()].tobytes().decode("ascii").strip("\x00 ")
                r = conn["PointRange"][" data"][()].reshape(-1)
                d = conn["PointRangeDonor"][" data"][()].reshape(-1)
                g = grids[zone_name]
                dg = grids[donor]
                for (i, j), (di, dj) in (
                    ((r[0], r[1]), (d[0], d[1])),
                    ((r[2], r[3]), (d[2], d[3])),
                ):
                    shifted = g[j - 1, i - 1, :2] + translation
                    assert np.allclose(shifted, dg[dj - 1, di - 1, :2], atol=1e-6)
    assert found == 4


def test_cgns_rejects_rotational_periodicity(tmp_path):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    og = sector()
    assert og.boundary_report().ok
    with pytest.raises(opengrid.ValidationError, match="rotational periodicity"):
        og.export_cgns(tmp_path / "s.cgns")


@pytest.mark.parametrize("case", ["channel", "sector"])
def test_gmsh_periodic(tmp_path, case):
    meshio = pytest.importorskip("meshio")
    og = sector()
    if case == "channel":
        og = channel()
        og.pair_periodic("left", "right")
    path = tmp_path / f"{case}.msh"
    summary = og.export_gmsh(path)
    assert summary["periodic_links"] == 1
    mesh = meshio.read(path)
    (link,) = mesh.gmsh_periodic
    dim, (slave, master), affine, nodes = link  # meshio's layout
    assert dim == 1
    assert (int(slave), int(master)) == (3, 2)  # target curve 3 (right/side1), master 2
    matrix = np.asarray(affine).reshape(4, 4)
    nodes = np.asarray(nodes)
    master_points = np.c_[mesh.points[nodes[:, 1]], np.ones(len(nodes))]
    mapped = (matrix @ master_points.T).T[:, :3]
    assert np.allclose(mapped, mesh.points[nodes[:, 0]], atol=1e-12)
    if find_gmsh():
        log = gmsh_read(path)
        assert "Error" not in log and "Warning" not in log, log
        assert f"{len(mesh.points)} nodes" in log


def test_vtkhdf_periodic_fields(tmp_path):
    og = channel()
    og.pair_periodic("left", "right")
    path = tmp_path / "c.vtkhdf"
    og.export_vtkhdf(path)
    with h5py.File(path) as f:
        assert f["VTKHDF/left/FieldData/periodic_partner"][0] == b"right"
        affine = f["VTKHDF/left/FieldData/periodic_affine"][()].reshape(4, 4)
        assert np.allclose(affine[:3, 3], [2, 0, 0])
        assert f["VTKHDF/right/FieldData/periodic_source"][0] == b"left"


def test_gmsh_executable_reads_exports(tmp_path):
    if find_gmsh() is None:
        pytest.skip("Gmsh not found (set OPENGRID_GMSH)")
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "tp", Path(__file__).resolve().parent / "test_patches.py"
    )
    tp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tp)
    for og, name in ((tp.step_session(), "step"), (tp.airfoil_session(), "airfoil")):
        path = tmp_path / f"{name}.msh"
        summary = og.export_gmsh(path)
        log = gmsh_read(path)
        assert "Error" not in log and "Warning" not in log, log
        assert re.search(rf"\b{summary['nodes']} nodes", log)
        assert re.search(rf"\b{summary['elements']} elements", log)
