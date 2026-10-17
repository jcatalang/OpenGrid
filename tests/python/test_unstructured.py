import math
from collections import Counter

import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read

import opengrid
from opengrid.project import to_json


def annulus(n_outer=61, n_inner=41, generator=None) -> opengrid.Session:
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="hole")
    og.create_circle((0, 0), 5.0, name="outer")
    og.create_connector("hole", n_inner, name="c-hole")
    og.create_connector("outer", n_outer, name="c-outer")
    og.create_unstructured_domain(
        connectors=["c-hole", "c-outer"], generator=generator, name="ring"
    )
    return og


def polygon_area(n, r):
    return 0.5 * (n - 1) * r * r * math.sin(2 * math.pi / (n - 1))


def edge_counts(triangles):
    counts = Counter()
    for t in triangles:
        for k in range(3):
            counts[(int(t[k]), int(t[(k + 1) % 3]))] += 1
    return counts


def test_annulus_area_quality_and_conformity():
    og = annulus()
    mesh = og.unstructured_mesh("ring")
    q = og.domain_quality("ring")
    assert q["area"].sum() == pytest.approx(polygon_area(61, 5) - polygon_area(41, 1), rel=1e-12)
    assert q["area"].min() > 0
    assert q["min_angle"].min() >= 20.7
    counts = edge_counts(mesh["triangles"])
    assert max(counts.values()) == 1
    # Boundary: outer loop counter-clockwise (indices 0..59), hole clockwise (60..99).
    params = og.entity("ring").params
    assert [item["reversed"] for loop in params["loops"] for item in loop] == [False, True]
    outer = list(range(60))
    hole = list(range(60, 100))
    for loop in (outer, hole):
        for a, b in zip(loop, loop[1:] + loop[:1], strict=True):
            assert counts[(a, b)] == 1
            assert counts[(b, a)] == 0
    # Boundary points are exactly the connector nodes.
    nodes = og.connector_nodes("c-outer").points
    assert np.array_equal(mesh["points"][:60], nodes[:-1])


def test_deterministic_and_cached():
    a = annulus().unstructured_mesh("ring")
    og = annulus()
    b = og.unstructured_mesh("ring")
    assert np.array_equal(a["points"], b["points"])
    assert np.array_equal(a["triangles"], b["triangles"])
    assert og.unstructured_mesh("ring") is b
    og.create_point((9, 9))
    assert og.unstructured_mesh("ring") is b


def test_size_controls():
    coarse = annulus(generator={"growth_rate": 1.5})
    fine = annulus(generator={"max_size": 0.2})
    sourced = annulus(generator={"sources": [{"p0": (3, 0), "size": 0.03, "radius": 0.1}]})
    n = {
        name: len(s.unstructured_mesh("ring")["triangles"])
        for name, s in (("coarse", coarse), ("fine", fine), ("sourced", sourced))
    }
    base = len(annulus().unstructured_mesh("ring")["triangles"])
    assert n["coarse"] < base < n["fine"]
    assert n["sourced"] > base
    mesh = sourced.unstructured_mesh("ring")
    centroids = mesh["points"][mesh["triangles"]].mean(axis=1)
    near = np.hypot(centroids[:, 0] - 3, centroids[:, 1]) < 0.1
    q = sourced.domain_quality("ring")["area"].reshape(-1)
    assert q[near].max() < 0.01  # equilateral area at size ~0.036 is ~0.0006


def test_loop_validation_and_assembly():
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="hole")
    og.create_circle((0, 0), 5.0, name="outer")
    og.create_circle((10, 0), 1.0, name="away")
    og.create_line((0, 0), (1, 0), name="stub")
    for name, n in (("hole", 21), ("outer", 41), ("away", 21), ("stub", 3)):
        og.create_connector(name, n, name=f"c-{name}")
    with pytest.raises(opengrid.ValidationError, match="outer loop must run counter-clockwise"):
        og.create_unstructured_domain(loops=[[("c-outer", True)]])
    with pytest.raises(opengrid.ValidationError, match="holes\\) must run clockwise"):
        og.create_unstructured_domain(loops=[["c-outer"], ["c-hole"]])
    with pytest.raises(opengrid.ValidationError, match="lies outside the outer loop"):
        og.create_unstructured_domain(loops=[["c-outer"], [("c-away", True)]])
    with pytest.raises(opengrid.ValidationError, match="is not closed"):
        og.create_unstructured_domain(loops=[["c-stub"]])
    with pytest.raises(opengrid.ValidationError, match="must form closed loops"):
        og.create_unstructured_domain(connectors=["c-outer", "c-stub"])
    with pytest.raises(opengrid.ValidationError, match="growth_rate"):
        og.create_unstructured_domain(loops=[["c-outer"]], generator={"growth_rate": 0.9})
    name = og.create_unstructured_domain(connectors=["c-outer"])
    assert og.domain_quality(name)["area"].sum() == pytest.approx(polygon_area(41, 5), rel=1e-12)


def test_dependencies_round_trip_and_replay(tmp_path):
    og = annulus()
    with pytest.raises(opengrid.ValidationError, match="used by unstructured_domain 'ring'"):
        og.delete("c-hole")
    with pytest.raises(opengrid.ValidationError, match="would become invalid"):
        og.modify("hole", radius=6.0)  # hole outside the outer loop
    og.set_unstructured_generator("ring", min_angle=30.0)
    assert og.domain_quality("ring")["min_angle"].min() >= 25.0
    og.save(tmp_path / "u.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "u.ogp")
    assert to_json(other.store) == to_json(og.store)
    assert np.array_equal(
        other.unstructured_mesh("ring")["triangles"], og.unstructured_mesh("ring")["triangles"]
    )
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def hybrid() -> opengrid.Session:
    """Structured O-ring (r in [1, 2]) and an unstructured region (r in [2, 6]) sharing the
    circle r = 2."""
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
    og.create_unstructured_domain(connectors=["c-mid", "c-far"], name="outside")
    og.create_patch("p-wall", "wall", ["c-wall"])
    og.create_patch("p-farfield", "farfield", ["c-far"])
    return og


def test_hybrid_conformity_and_exports(tmp_path):
    og = hybrid()
    report = og.boundary_report()
    assert report.ok, report.issues
    assert sorted(report.interfaces) == ["c-cut (cut in ring)", "c-mid (ring | outside)"]
    ring = og.domain_grid("ring")
    tri = og.unstructured_mesh("outside")
    # The unstructured boundary on r = 2 uses the very same nodes as the ring's top row.
    shared = og.connector_nodes("c-mid").points
    assert np.array_equal(ring[-1], shared[::-1]) or np.array_equal(ring[-1], shared)
    meshio = pytest.importorskip("meshio")
    path = tmp_path / "hybrid.msh"
    summary = og.export_gmsh(path)
    mesh = meshio.read(path)
    expected_points = 81 * 21 - 21 + len(tri["points"]) - 80  # cut and shared circle merged
    assert len(mesh.points) == expected_points == summary["nodes"]
    types = Counter()
    for block in mesh.cells:
        types[block.type] += len(block.data)
    assert types == {"quad": 80 * 20, "triangle": len(tri["triangles"]), "line": 80 + 60}
    # Conforming: every interior edge shared by exactly two cells.
    edges = Counter()
    for block in mesh.cells:
        if block.type == "line":
            continue
        for cell in block.data:
            for k in range(len(cell)):
                edges[tuple(sorted((int(cell[k]), int(cell[(k + 1) % len(cell)]))))] += 1
    assert max(edges.values()) == 2
    assert sum(1 for c in edges.values() if c == 1) == 80 + 60  # wall and farfield only
    if find_gmsh():
        log = gmsh_read(path)
        assert "Error" not in log and "Warning" not in log, log
    vtk = pytest.importorskip("vtk")
    vpath = tmp_path / "hybrid.vtkhdf"
    og.export_vtkhdf(vpath)
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(vpath))
    reader.Update()
    block = reader.GetOutputDataObject(0).GetPartitionedDataSet(0).GetPartition(0)
    assert block.GetNumberOfCells() == 80 * 20 + len(tri["triangles"])
    assert {block.GetCellType(0), block.GetCellType(block.GetNumberOfCells() - 1)} == {
        vtk.VTK_QUAD,
        vtk.VTK_TRIANGLE,
    }


def test_unsupported_exports_are_reported(tmp_path):
    og = hybrid()
    with pytest.raises(opengrid.ValidationError, match="Plot3D holds structured blocks only"):
        og.export_plot3d(tmp_path / "h.xyz")


def test_zone_with_unstructured_domain():
    og = hybrid()
    og.create_zone("fluid", ["ring", "outside"])
    assert og.boundary_report().ok


def test_advancing_front_conformity_area_and_quality():
    og = annulus(generator={"algorithm": "advancing_front", "growth_rate": 1.15})
    mesh = og.unstructured_mesh("ring")
    q = og.domain_quality("ring")
    assert q["area"].sum() == pytest.approx(polygon_area(61, 5) - polygon_area(41, 1), rel=1e-12)
    assert q["area"].min() > 0
    assert q["min_angle"].min() >= 25.0
    assert q["min_angle"].mean() >= 45.0
    counts = edge_counts(mesh["triangles"])
    assert max(counts.values()) == 1
    assert np.array_equal(mesh["points"][:60], og.connector_nodes("c-outer").points[:-1])
    assert mesh["warnings"] == []
    again = annulus(generator={"algorithm": "advancing_front", "growth_rate": 1.15})
    assert np.array_equal(again.unstructured_mesh("ring")["triangles"], mesh["triangles"])


def test_algorithm_switch_undo_and_validation(tmp_path):
    og = annulus()
    delaunay = og.unstructured_mesh("ring")["triangles"]
    og.set_unstructured_generator("ring", algorithm="advancing_front")
    front = og.unstructured_mesh("ring")["triangles"]
    assert len(front) != len(delaunay) or not np.array_equal(front, delaunay)
    og.undo()
    assert np.array_equal(og.unstructured_mesh("ring")["triangles"], delaunay)
    og.redo()
    path = tmp_path / "af.ogp"
    og.save(path)
    loaded = opengrid.Session()
    loaded.load(path)
    assert loaded.entity("ring").params["generator"]["algorithm"] == "advancing_front"
    assert np.array_equal(loaded.unstructured_mesh("ring")["triangles"], front)
    with pytest.raises(opengrid.ValidationError, match="algorithm"):
        og.set_unstructured_generator("ring", algorithm="paving")
    # Projects written before the option existed default to the Delaunay mesher.
    text = path.read_text().replace('"algorithm": "advancing_front", ', "")
    assert '"algorithm"' not in text
    old = tmp_path / "old.ogp"
    old.write_text(text)
    previous = opengrid.Session()
    previous.load(old)
    assert previous.entity("ring").params["generator"]["algorithm"] == "delaunay"


def test_advancing_front_sources():
    og = annulus(generator={"algorithm": "advancing_front"})
    coarse = len(og.unstructured_mesh("ring")["triangles"])
    og.set_unstructured_generator("ring", sources=[{"p0": (3, 0), "size": 0.05, "radius": 0.3}])
    mesh = og.unstructured_mesh("ring")
    assert len(mesh["triangles"]) > coarse
    q = og.domain_quality("ring")
    assert q["min_angle"].min() >= 20.7


def test_unstructured_domain_from_curves_is_one_step_and_names_open_ends():
    from opengrid.validation import ValidationError

    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="c")
    name = og.unstructured_domain_from_curves("c", 61)
    assert og.names("connector") == ["connector-1"]
    assert og.entity("connector-1").params["dimension"] == 61
    assert "triangles" in og.domain_status(name)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.undo()
    assert og.names("connector") == [] and og.names("unstructured_domain") == []
    # Curves and existing connectors mixed: a half disc.
    og = opengrid.Session()
    og.create_spline([(1, 0), (0, 1), (-1, 0)], name="top")
    og.create_line((-1, 0), (1, 0), name="bottom")
    og.create_connector("bottom", 21, name="c-bottom")
    name = og.unstructured_domain_from_curves(["top"], 31, connectors=["c-bottom"])
    assert "triangles" in og.domain_status(name)
    with pytest.raises(ValidationError) as error:
        og.unstructured_domain_from_curves(["top"], 31)
    text = str(error.value)
    assert "closed loops" in text and "(on top)" in text and "meets 0 other" in text
    assert og.names("unstructured_domain") == [name]
