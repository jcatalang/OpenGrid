"""3D periodic pairing: the ends of extruded volumes and block faces, matched node by node and
face by face, and the pairs in the exports (CGNS, Gmsh, VTKHDF, high order) and the GUI."""

import math

import h5py
import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read
from test_volumes import disc

import opengrid
from opengrid import _core
from opengrid.model.periodic import affine_matrix, apply, inverse
from opengrid.model.periodic_faces import FaceMatch, match_points
from opengrid.project import to_json


def slab(og=None, mode="translate", layers=4, **settings) -> tuple[opengrid.Session, str]:
    """A quad-dominant disc extruded (translation by 0.5 in z unless told otherwise); patches:
    back (the domain), front (the volume's top) and rim (the connector)."""
    if og is None:
        og = opengrid.Session()
    disc(og, n=31, generator={"element_type": "quad_dominant"})
    if mode == "translate" and not settings:
        settings = {"vector": (0, 0, 0.5)}
    name = og.create_extruded_volume("disc", mode, dimension=layers + 1, **settings)
    og.create_patch("back", "periodic", ["disc"])
    og.create_patch("front", "periodic", [{"volume": name, "face": "top"}])
    og.create_patch("rim", "wall", ["cc"])
    return og, name


def sector_block(angle=math.pi / 3) -> tuple[opengrid.Session, str]:
    """The rectangle x in [1, 3], y in [0, 1] (z = 0) turned about the y axis: the block's kmax
    face is its kmin face rotated by ``angle``."""
    og = opengrid.Session()
    corners = ((1, 0), (3, 0), (3, 1), (1, 1))
    for k in range(4):
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"l{k}")
    for k, n in enumerate((9, 6, 9, 6)):
        og.create_connector(f"l{k}", n, name=f"c{k}")
    og.create_domain([["c0"], ["c1"], ["c2"], ["c3"]], name="base")
    name = og.create_extruded_block(
        "base", "rotate", vector=(0, 1, 0), axis_point=(0, 0, 0), angle=angle, dimension=7
    )
    og.create_patch("side0", "periodic", [{"block": name, "face": "kmin"}])
    og.create_patch("side1", "periodic", [{"block": name, "face": "kmax"}])
    og.create_patch(
        "walls", "wall", [{"block": name, "face": f} for f in ("imin", "imax", "jmin", "jmax")]
    )
    return og, name


def _cgns_nodes(group, label):
    return {
        n: g
        for n, g in group.items()
        if isinstance(g, h5py.Group) and g.attrs.get("label") == label.encode()
    }


def _zone_points(zone) -> np.ndarray:
    return np.stack(
        [zone[f"GridCoordinates/Coordinate{c}/ data"][()] for c in "XYZ"], axis=-1
    ).reshape(-1, 3)


def _periodic_property(connection) -> dict:
    group = connection["GridConnectivityProperty/Periodic"]
    keys = ("RotationCenter", "RotationAngle", "Translation")
    return {key: group[f"{key}/ data"][()] for key in keys}


def test_match_points_finds_images_within_the_tolerance_only():
    rng = np.random.default_rng(3)
    targets = rng.random((500, 3)) * 10
    image = targets[::-1] + 1e-10
    found, distance = match_points(image, targets, 1e-8)
    assert np.array_equal(found, np.arange(500)[::-1])
    assert distance.max() < 2e-10
    found, _ = match_points(targets + 1e-6, targets, 1e-8)
    assert (found == -1).all()
    # Points on cell boundaries are found in a neighbouring cell.
    on_edge = np.array([[2e-8, 2e-8, 2e-8], [4e-8, 0.0, 0.0]])
    found, _ = match_points(on_edge - 5e-9, on_edge, 1e-8)
    assert found.tolist() == [0, 1]


def test_extruded_volume_ends_pair_by_translation(tmp_path):
    og, _name = slab()
    report = og.boundary_report()
    messages = [str(i) for i in report.issues]
    assert "back: periodic patch has no partner (use pair_periodic)" in messages
    og.pair_periodic("back", "front")
    params = og.entity("back").params["periodic"]
    assert params["translation"] == pytest.approx((0.0, 0.0, 0.5))
    assert params["rotation_axis"] == (0.0, 0.0, 1.0)
    report = og.boundary_report()
    assert report.ok, report.issues
    ((source, target, match),) = report.periodic
    mesh2d = og.unstructured_mesh("disc")
    assert (source.name, target.name) == ("back", "front")
    assert match == FaceMatch(
        len(mesh2d["points"]), len(mesh2d["triangles"]) + len(mesh2d["quads"])
    )
    # CGNS: the periodic faces are connectivity (both ways, point lists), not BCs.
    if _core.cgns_available():
        path = tmp_path / "slab.cgns"
        summary = og.export_cgns(path)
        assert summary["boundaries"] == 1 and summary["connections"] == 2
        with h5py.File(path) as f:
            (zone,) = _cgns_nodes(f["OpenGrid"], "Zone_t").values()
            points = _zone_points(zone)
            connections = _cgns_nodes(zone["ZoneGridConnectivity"], "GridConnectivity_t")
            assert len(connections) == 2
            shifts = []
            for connection in connections.values():
                prop = _periodic_property(connection)
                assert np.allclose(prop["RotationAngle"], 0) and np.allclose(
                    prop["RotationCenter"], 0
                )
                mine = connection["PointList/ data"][()].reshape(-1) - 1
                theirs = connection["PointListDonor/ data"][()].reshape(-1) - 1
                assert len(mine) == len(mesh2d["points"])
                assert np.allclose(points[mine] + prop["Translation"], points[theirs], atol=1e-6)
                shifts.append(float(prop["Translation"][2]))
            assert sorted(shifts) == pytest.approx([-0.5, 0.5])
    # Gmsh: a $Periodic block whose affine matrix maps the master nodes onto the slaves.
    meshio = pytest.importorskip("meshio")
    msh = tmp_path / "slab.msh"
    assert og.export_gmsh(msh)["periodic_links"] == 1
    mesh = meshio.read(msh)
    ((dim, _tags, affine, nodes),) = mesh.gmsh_periodic
    assert dim == 2
    nodes = np.asarray(nodes)
    assert len(nodes) == len(mesh2d["points"])
    master = np.c_[mesh.points[nodes[:, 1]], np.ones(len(nodes))]
    mapped = (np.asarray(affine).reshape(4, 4) @ master.T).T[:, :3]
    assert np.allclose(mapped, mesh.points[nodes[:, 0]], atol=1e-12)
    if find_gmsh():
        log = gmsh_read(msh)
        assert "Error" not in log and "Warning" not in log, log
    # VTKHDF: the pairing as field data on the patch blocks.
    vtkhdf = tmp_path / "slab.vtkhdf"
    og.export_vtkhdf(vtkhdf)
    with h5py.File(vtkhdf) as f:
        assert f["VTKHDF/back/FieldData/periodic_partner"][0] == b"front"
        matrix = f["VTKHDF/back/FieldData/periodic_affine"][()].reshape(4, 4)
        assert np.allclose(matrix[:3, 3], [0, 0, 0.5])
        assert f["VTKHDF/front/FieldData/periodic_source"][0] == b"back"
    # High order: the links extend to the edge and face nodes of the periodic faces.
    high = og.high_order(2)
    ((s, t, transform, targets, sources),) = high.periodic
    assert (high.patch_names[s], high.patch_names[t]) == ("back", "front")
    edges = set()
    for cells in (mesh2d["triangles"], mesh2d["quads"]):
        for cell in np.asarray(cells).tolist():
            for k in range(len(cell)):
                a, b = cell[k], cell[(k + 1) % len(cell)]
                edges.add((min(a, b), max(a, b)))
    assert len(targets) == len(mesh2d["points"]) + len(edges) + len(mesh2d["quads"])
    assert np.allclose(apply(high.points[sources], transform), high.points[targets], atol=1e-12)
    assert og.export_high_order(tmp_path / "slab2.msh", order=2)["periodic_links"] == 1
    if _core.cgns_available():
        assert og.export_high_order(tmp_path / "slab2.cgns", order=2)["connections"] == 2
    # Replay and project round trip keep the pairing (with its axis).
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.save(tmp_path / "slab.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "slab.ogp")
    assert to_json(other.store) == to_json(og.store)
    assert other.boundary_report().ok


def test_rotational_periodicity_of_a_block_sector(tmp_path):
    angle = math.pi / 3
    og, name = sector_block(angle)
    grid = og.block_grid(name)
    pairing = {
        "rotation_center": (0.0, 0.0, 0.0),
        "rotation_axis": (0.0, 1.0, 0.0),
        "rotation_angle": angle,
        "translation": (0.0, 0.0, 0.0),
    }
    assert np.allclose(apply(grid[0].reshape(-1, 3), pairing), grid[-1].reshape(-1, 3), atol=1e-12)
    og.pair_periodic(
        "side0", "side1", translation=(0, 0, 0), rotation_angle=angle, rotation_axis=(0, 1, 0)
    )
    report = og.boundary_report()
    assert report.ok, report.issues
    assert report.periodic[0][2] == FaceMatch(9 * 6, 8 * 5)
    if _core.cgns_available():
        path = tmp_path / "sector.cgns"
        summary = og.export_cgns(path)
        assert summary["boundaries"] == 4 and summary["connections"] == 2
        with h5py.File(path) as f:
            (zone,) = _cgns_nodes(f["OpenGrid"], "Zone_t").values()
            connections = _cgns_nodes(zone["ZoneGridConnectivity"], "GridConnectivity1to1_t")
            assert len(connections) == 2
            angles = []
            for connection in connections.values():
                assert tuple(connection["Transform/ data"][()].reshape(-1)) == (1, 2, 3)
                prop = _periodic_property(connection)
                assert np.allclose(prop["RotationCenter"], 0) and np.allclose(
                    prop["Translation"], 0
                )
                assert prop["RotationAngle"][0] == 0 and prop["RotationAngle"][2] == 0
                angles.append(float(prop["RotationAngle"][1]))
                r = connection["PointRange/ data"][()].reshape(-1)
                d = connection["PointRangeDonor/ data"][()].reshape(-1)
                transform = dict(pairing, rotation_angle=float(prop["RotationAngle"][1]))
                for (i, j, k), (di, dj, dk) in (
                    ((r[0], r[1], r[2]), (d[0], d[1], d[2])),
                    ((r[3], r[4], r[5]), (d[3], d[4], d[5])),
                ):
                    image = apply(grid[k - 1, j - 1, i - 1], transform)
                    assert np.allclose(image, grid[dk - 1, dj - 1, di - 1], atol=1e-6)
            assert sorted(angles) == pytest.approx([-angle, angle], abs=1e-6)
    meshio = pytest.importorskip("meshio")
    msh = tmp_path / "sector.msh"
    assert og.export_gmsh(msh)["periodic_links"] == 1
    mesh = meshio.read(msh)
    ((dim, _tags, affine, nodes),) = mesh.gmsh_periodic
    nodes = np.asarray(nodes)
    assert dim == 2 and len(nodes) == 9 * 6
    master = np.c_[mesh.points[nodes[:, 1]], np.ones(len(nodes))]
    mapped = (np.asarray(affine).reshape(4, 4) @ master.T).T[:, :3]
    assert np.allclose(mapped, mesh.points[nodes[:, 0]], atol=1e-12)
    assert np.allclose(np.asarray(affine).reshape(4, 4), affine_matrix(pairing))
    # The inverse transform maps the images back.
    back = inverse(pairing)
    assert np.allclose(apply(apply(grid[0], pairing), back), grid[0], atol=1e-12)
    # A wrong angle: the report says how many nodes have no image and how far the first is.
    og.pair_periodic(
        "side0",
        "side1",
        translation=(0, 0, 0),
        rotation_angle=1.01 * angle,
        rotation_axis=(0, 1, 0),
    )
    messages = [str(i) for i in og.boundary_report().issues]
    assert any("of its 54 nodes have no periodic image in side1" in m for m in messages), messages
    assert any("from the nearest node of side1" in m for m in messages)


def test_rotated_extrusion_ends_pair_by_rotation(tmp_path):
    angle = math.pi / 6
    og, _name = slab(mode="rotate", vector=(1, 0, 0), axis_point=(0, -2, 0), angle=angle)
    og.pair_periodic(
        "back",
        "front",
        translation=(0, 0, 0),
        rotation_angle=angle,
        rotation_axis=(1, 0, 0),
        rotation_center=(0, -2, 0),
    )
    report = og.boundary_report()
    assert report.ok, report.issues
    if _core.cgns_available():
        path = tmp_path / "turned.cgns"
        assert og.export_cgns(path)["connections"] == 2
        with h5py.File(path) as f:
            (zone,) = _cgns_nodes(f["OpenGrid"], "Zone_t").values()
            points = _zone_points(zone)
            for connection in _cgns_nodes(
                zone["ZoneGridConnectivity"], "GridConnectivity_t"
            ).values():
                prop = _periodic_property(connection)
                transform = {
                    "rotation_center": prop["RotationCenter"],
                    "rotation_axis": (1.0, 0.0, 0.0),
                    "rotation_angle": float(prop["RotationAngle"][0]),
                    "translation": prop["Translation"],
                }
                mine = connection["PointList/ data"][()].reshape(-1) - 1
                theirs = connection["PointListDonor/ data"][()].reshape(-1) - 1
                assert np.allclose(apply(points[mine], transform), points[theirs], atol=1e-5)
    # Pairing refusals and mismatches are reported, never repaired.
    with pytest.raises(ValueError, match="give the translation explicitly"):
        og.pair_periodic("back", "front", rotation_angle=angle)
    with pytest.raises(opengrid.ValidationError, match=r"rotation_axis.*must not be zero"):
        og.pair_periodic("back", "front", translation=(0, 0, 0), rotation_axis=(0, 0, 0))
    og.set_patch_type("rim", "periodic")
    og.pair_periodic(
        "back", "rim", translation=(0, 0, 0), rotation_angle=angle, rotation_axis=(1, 0, 0)
    )
    messages = [str(i) for i in og.boundary_report().issues]
    assert any("back: has" in m and "nodes, rim has" in m for m in messages), messages
    assert any("faces, rim has" in m for m in messages)


def test_gui_pairs_periodic_patches_in_the_property_panel(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og, _name = slab(window.session)
    window.set_selection([og.entity("back").id])
    _kind, editor = window.property_panel._editors["periodic"]
    assert editor.partner.currentData() is None
    editor.partner.setCurrentText("front")
    assert window.property_panel.apply()
    params = og.entity("back").params["periodic"]
    assert og.store.get(params["partner"]).name == "front"
    assert params["translation"] == pytest.approx((0.0, 0.0, 0.5))
    assert og.boundary_report().ok
    # Applying again without a change is no change; the target says where it is paired from.
    window.set_selection([og.entity("back").id])
    assert not window.property_panel.apply()
    window.set_selection([og.entity("front").id])
    _kind, editor = window.property_panel._editors["periodic"]
    assert editor.partner.currentData() is None
    # Unpairing from the source.
    window.set_selection([og.entity("back").id])
    _kind, editor = window.property_panel._editors["periodic"]
    editor.partner.setCurrentIndex(0)
    assert window.property_panel.apply()
    assert og.entity("back").params["periodic"] is None
    window.close()
