"""Phase 5c: quad-dominant / all-quad unstructured domains, smoothing, CGNS unstructured zones."""

import shutil
import subprocess
from collections import Counter

import h5py
import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read

import opengrid
from opengrid import _core


def annulus(n_outer=61, n_inner=41, **generator) -> opengrid.Session:
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="hole")
    og.create_circle((0, 0), 5.0, name="outer")
    og.create_connector("hole", n_inner, name="c-hole")
    og.create_connector("outer", n_outer, name="c-outer")
    og.create_unstructured_domain(
        connectors=["c-hole", "c-outer"], generator=generator, name="ring"
    )
    og.create_patch("wall", "wall", ["c-hole"])
    og.create_patch("far", "farfield", ["c-outer"])
    return og


def hybrid(**generator) -> opengrid.Session:
    """Structured O-ring (r in [1, 2]) and an unstructured region (r in [2, 6])."""
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="wall")
    og.create_circle((0, 0), 2.0, name="mid")
    og.create_circle((0, 0), 6.0, name="far")
    og.create_line((1, 0), (2, 0), name="cut")
    og.create_connector("wall", 81, name="c-wall")
    og.create_connector("mid", 81, name="c-mid")
    og.create_connector("far", 61, name="c-far")
    og.create_connector("cut", 21, {"type": "tanh", "begin_spacing": 0.01}, name="c-cut")
    og.create_structured_domain(connectors=["c-wall", "c-mid", "c-cut"], name="ring")
    og.create_unstructured_domain(connectors=["c-mid", "c-far"], generator=generator, name="out")
    og.create_patch("p-wall", "wall", ["c-wall"])
    og.create_patch("p-far", "farfield", ["c-far"])
    return og


def directed_edges(mesh):
    edges = Counter()
    for cells in (mesh["triangles"], mesh["quads"]):
        for cell in cells:
            n = len(cell)
            for k in range(n):
                edges[(int(cell[k]), int(cell[(k + 1) % n]))] += 1
    return edges


def test_quad_dominant_and_all_quad_annulus():
    reference = annulus().domain_quality("ring")["area"].sum()
    for element_type in ("quad_dominant", "quads"):
        og = annulus(element_type=element_type)
        mesh = og.unstructured_mesh("ring")
        q = og.domain_quality("ring")
        assert q["area"].sum() == pytest.approx(reference, rel=1e-12)
        assert q["scaled_jacobian"].min() > 0
        assert q["area"].shape == (1, len(mesh["triangles"]) + len(mesh["quads"]))
        edges = directed_edges(mesh)
        assert max(edges.values()) == 1
        # Boundary: exactly the connector nodes, each boundary edge once.
        outer = og.connector_nodes("c-outer").points
        assert np.array_equal(mesh["points"][:60], outer[:-1])
        for k in range(60):
            assert edges[(k, (k + 1) % 60)] == 1
        if element_type == "quads":
            assert len(mesh["triangles"]) == 0
        else:
            assert len(mesh["quads"]) > 4 * len(mesh["triangles"])
        assert mesh["warnings"] == []
    status = annulus(element_type="quads").domain_status("ring")
    assert "quads" in status and "triangles" not in status


def test_all_quads_needs_even_loops():
    og = annulus(n_outer=62, element_type="quads")
    with pytest.raises(opengrid.ValidationError, match="even number of segments"):
        og.unstructured_mesh("ring")


def test_generator_validation_and_round_trip(tmp_path):
    og = annulus()
    for key, value in (
        ("element_type", "hexes"),
        ("smoothing", "winslow"),
        ("smoothing_iterations", -1),
        ("edge_swap", "yes"),
    ):
        with pytest.raises(opengrid.ValidationError, match=key):
            og.set_unstructured_generator("ring", **{key: value})
    og.set_unstructured_generator(
        "ring", element_type="quads", smoothing="angle_based", smoothing_iterations=3
    )
    path = tmp_path / "q.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert np.array_equal(
        other.unstructured_mesh("ring")["quads"], og.unstructured_mesh("ring")["quads"]
    )
    replayed = opengrid.replay(og.journal_text())
    assert replayed.entity("ring").params["generator"] == og.entity("ring").params["generator"]


def test_smoothing_improves_mean_quality_and_keeps_worst():
    plain = annulus().domain_quality("ring")
    for method in ("laplacian", "length_weighted", "angle_based"):
        q = annulus(smoothing=method).domain_quality("ring")
        assert q["scaled_jacobian"].mean() > plain["scaled_jacobian"].mean()
        assert q["min_angle"].min() >= plain["min_angle"].min() - 1e-9
        assert q["area"].sum() == pytest.approx(plain["area"].sum(), rel=1e-12)


def _children(group, label):
    return {
        name: node
        for name, node in group.items()
        if isinstance(node, h5py.Group) and node.attrs.get("label") == label.encode()
    }


def _check_cgns(path):
    checker = shutil.which("cgnscheck")
    if checker is None:
        return
    out = subprocess.run([checker, str(path)], capture_output=True, text=True, check=False).stdout
    assert "error" not in out.lower() and "warning" not in out.lower(), out


@pytest.mark.parametrize("element_type", ["triangles", "quad_dominant", "quads"])
def test_cgns_hybrid(tmp_path, element_type):
    if not _core.cgns_available():
        pytest.skip("built without CGNS")
    og = hybrid(element_type=element_type)
    path = tmp_path / "h.cgns"
    summary = og.export_cgns(path)
    assert summary == {"zones": 2, "families": 3, "boundaries": 2, "connections": 4}
    _check_cgns(path)
    mesh = og.unstructured_mesh("out")
    ring = og.domain_grid("ring")
    with h5py.File(path) as f:
        base = f["OpenGrid"]
        zone = base["out"]
        assert zone["ZoneType"][" data"][()].tobytes().decode().strip("\x00") == "Unstructured"
        x = zone["GridCoordinates/CoordinateX"][" data"][()]
        y = zone["GridCoordinates/CoordinateY"][" data"][()]
        assert np.array_equal(x, mesh["points"][:, 0]) and np.array_equal(y, mesh["points"][:, 1])
        sections = _children(zone, "Elements_t")
        cells = 0
        for name, kind, nodes, count in (
            ("Triangles", 5, 3, len(mesh["triangles"])),
            ("Quadrilaterals", 7, 4, len(mesh["quads"])),
        ):
            if count == 0:
                assert name not in sections
                continue
            assert sections[name][" data"][()][0] == kind  # TRI_3 = 5, QUAD_4 = 7
            assert len(sections[name]["ElementConnectivity"][" data"][()]) == nodes * count
            cells += count
        # Farfield BC: a BAR_2 section (ElementType 3) whose range the BC_t points to.
        bc = base["out/ZoneBC/p-far"]
        r = bc["PointRange"][" data"][()].reshape(-1)
        bars = sections["p-far"]
        assert bars[" data"][()][0] == 3
        assert bars["ElementRange"][" data"][()].tolist() == r.tolist() == [cells + 1, cells + 60]
        location = bc["GridLocation"][" data"][()].tobytes().decode().strip("\x00")
        assert location == "EdgeCenter"
        far = set(map(tuple, np.round(og.connector_nodes("c-far").points[:, :2], 12)))
        ids = bars["ElementConnectivity"][" data"][()] - 1
        assert set(map(tuple, np.round(mesh["points"][ids, :2], 12))) == far
        # Interface: vertex point lists in both zones refer to the same coordinates.
        conn = _children(zone["ZoneGridConnectivity"], "GridConnectivity_t")
        ((_name, c),) = conn.items()
        points = c["PointList"][" data"][()].reshape(-1) - 1
        donor = c["PointListDonor"][" data"][()].reshape(-1, 2) - 1
        assert len(points) == 81
        assert np.array_equal(mesh["points"][points], ring[donor[:, 1], donor[:, 0]])
        back = _children(base["ring/ZoneGridConnectivity"], "GridConnectivity_t")
        assert len(back) == 1


def periodic_channel(element_type="triangles") -> opengrid.Session:
    """Unstructured channel [0, 2] x [0, 1], periodic in x (the right side runs downwards)."""
    og = opengrid.Session()
    corners = {"b": ((0, 0), (2, 0)), "r": ((2, 1), (2, 0)), "t": ((2, 1), (0, 1))}
    corners["l"] = ((0, 1), (0, 0))
    for name, (a, b) in corners.items():
        og.create_line(a, b, name=name)
    og.create_connector("b", 41, name="c-b")
    og.create_connector("t", 41, name="c-t")
    og.create_connector("l", 21, {"type": "tanh", "begin_spacing": 0.02}, name="c-l")
    og.create_connector("r", 21, {"type": "tanh", "begin_spacing": 0.02}, name="c-r")
    og.create_unstructured_domain(
        connectors=["c-b", "c-r", "c-t", "c-l"],
        generator={"element_type": element_type},
        name="channel",
    )
    og.create_patch("walls", "wall", ["c-b", "c-t"])
    og.create_patch("left", "periodic", ["c-l"])
    og.create_patch("right", "periodic", ["c-r"])
    og.pair_periodic("left", "right")
    return og


@pytest.mark.parametrize("element_type", ["triangles", "quads"])
def test_unstructured_periodic_cgns_and_gmsh(tmp_path, element_type):
    og = periodic_channel(element_type)
    assert og.boundary_report().ok
    mesh = og.unstructured_mesh("channel")
    if element_type == "quads":
        assert len(mesh["triangles"]) == 0
    if _core.cgns_available():
        path = tmp_path / "p.cgns"
        summary = og.export_cgns(path)
        assert summary["connections"] == 2
        _check_cgns(path)
        with h5py.File(path) as f:
            conns = _children(f["OpenGrid/channel/ZoneGridConnectivity"], "GridConnectivity_t")
            assert len(conns) == 2
            for conn in conns.values():
                t = conn["GridConnectivityProperty/Periodic/Translation"][" data"][()]
                points = conn["PointList"][" data"][()].reshape(-1) - 1
                donor = conn["PointListDonor"][" data"][()].reshape(-1) - 1
                shifted = mesh["points"][points, :2] + t
                assert np.allclose(shifted, mesh["points"][donor, :2], atol=1e-6)
    msh = tmp_path / "p.msh"
    summary = og.export_gmsh(msh)
    assert summary["periodic_links"] == 1
    if find_gmsh():
        log = gmsh_read(msh)
        assert "Error" not in log and "Warning" not in log, log
