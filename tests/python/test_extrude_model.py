"""Session.make_periodic_pair (two sides into a periodic pair in one step) and
Session.extrude_model (a 2D model as a one-cell-thick 3D model with periodic front and back)."""

import math

import pytest

import opengrid
from opengrid import _core
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def channel_with_dome(og):
    """A structured channel [0, 2] x [0, 1] with a half disc on its top sharing the top
    connector (the dome is unstructured); the bottom is a wall, the arc far field."""
    for name, (a, b) in (
        ("bottom", ((0, 0), (2, 0))),
        ("right", ((2, 0), (2, 1))),
        ("top", ((2, 1), (0, 1))),
        ("left", ((0, 1), (0, 0))),
    ):
        og.create_line(a, b, name=name)
    og.create_arc((1, 1), 1.0, 0.0, math.pi, name="arc")
    for name, n in (("bottom", 21), ("right", 11), ("top", 21), ("left", 11), ("arc", 31)):
        og.create_connector(name, n, name=f"c{name[0]}")
    og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="chan")
    og.create_unstructured_domain(connectors=["ca", "ct"], name="dome")
    og.set_boundary_type(["cb"], "wall")
    og.set_boundary_type(["ca"], "farfield")
    return og


def test_make_periodic_pair_from_two_sides():
    og = channel_with_dome(opengrid.Session())
    names = og.make_periodic_pair(["cl"], ["cr"])
    assert names == ["periodic-1a", "periodic-1b"]
    source = og.entity("periodic-1a").params
    assert source["patch_type"] == "periodic"
    assert og.entity("periodic-1b").params["patch_type"] == "periodic"
    assert og.store.get(source["periodic"]["partner"]).name == "periodic-1b"
    assert source["periodic"]["translation"] == pytest.approx((2.0, 0.0, 0.0))
    report = og.boundary_report()
    assert report.ok, report.issues
    assert len(report.periodic) == 1
    og.undo()
    assert "periodic-1a" not in og.store and "periodic-1b" not in og.store
    og.redo()
    assert og.boundary_report().ok
    with pytest.raises(ValueError, match="give the translation"):
        og.make_periodic_pair(["cl"], ["cr"], rotation_angle=0.5)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_extrude_model_one_cell_thick_with_periodic_ends(tmp_path):
    og = channel_with_dome(opengrid.Session())
    og.make_periodic_pair(["cl"], ["cr"], name="side")
    og.create_zone("fluid", ["chan", "dome"])
    assert og.boundary_report().ok
    names = og.extrude_model(0.1)
    assert names == ["chan-ext", "dome-ext", "front", "back"]
    assert og.entity("chan-ext").kind == "block" and og.entity("dome-ext").kind == "volume"
    assert og.is_3d()
    block = og.entity("chan-ext").id
    volume = og.entity("dome-ext").id
    # The 2D patches moved to the side faces: block faces by edge, the volume's side walls
    # stay connectors; the periodic pairing of the sides is kept.
    assert og.entity("wall").params["members"] == ({"block": block, "face": "jmin"},)
    assert og.entity("farfield").params["members"] == ({"connector": og.entity("ca").id},)
    side = og.entity("side-a").params
    assert side["members"] == ({"block": block, "face": "imin"},)
    assert og.entity("side-b").params["members"] == ({"block": block, "face": "imax"},)
    assert side["periodic"]["translation"] == pytest.approx((2.0, 0.0, 0.0))
    front = og.entity("front").params
    back = og.entity("back").params
    assert front["patch_type"] == "periodic" and back["patch_type"] == "periodic"
    assert front["members"] == ({"domain": og.entity("chan").id}, {"domain": og.entity("dome").id})
    assert back["members"] == ({"block": block, "face": "kmax"}, {"volume": volume, "face": "top"})
    assert og.store.get(front["periodic"]["partner"]).name == "back"
    assert front["periodic"]["translation"] == pytest.approx((0.0, 0.0, 0.1))
    members = [og.store.get(m["domain"]).name for m in og.entity("fluid").params["members"]]
    assert members == ["chan", "dome", "chan-ext", "dome-ext"]
    report = og.boundary_report()
    assert report.ok, report.issues
    assert len(report.periodic) == 2
    assert any("ct" in i for i in report.interfaces)  # the shared connector: block | volume
    grid = og.block_grid("chan-ext")
    assert grid.shape[0] == 2 and abs(grid[-1, 0, 0, 2] - 0.1) < 1e-12
    assert len(og.volume_mesh("dome-ext")["prisms"]) > 0
    if _core.cgns_available():
        summary = og.export_cgns(tmp_path / "flat.cgns")
        assert summary["connections"] >= 4  # the sides' pair and the ends' pair, both ways
    og.undo()
    assert not og.is_3d() and "front" not in og.store
    assert og.entity("wall").params["members"] == ({"connector": og.entity("cb").id},)
    og.redo()
    assert og.is_3d()
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.save(tmp_path / "flat.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "flat.ogp")
    assert to_json(other.store) == to_json(og.store)
    with pytest.raises(ValidationError, match="3D already"):
        og.extrude_model(0.1)


def test_extrude_model_symmetry_ends_and_refusals():
    og = channel_with_dome(opengrid.Session())
    og.set_boundary_type(["cl", "cr"], "inlet")
    names = og.extrude_model(0.2, layers=2, ends="symmetry")
    assert names[-2:] == ["front", "back"]
    assert og.entity("front").params["patch_type"] == "symmetry"
    assert og.entity("front").params["periodic"] is None
    assert og.block_grid("chan-ext").shape[0] == 3
    assert og.boundary_report().ok
    og.undo()
    with pytest.raises(ValidationError, match=r"periodic.*symmetry"):
        og.extrude_model(0.1, ends="walls")
    with pytest.raises(ValidationError, match="positive"):
        og.extrude_model(-1.0)
    # A structured domain's edge in two patches cannot become one block face.
    og2 = opengrid.Session()
    for name, (a, b) in (
        ("b1", ((0, 0), (1, 0))),
        ("b2", ((1, 0), (2, 0))),
        ("r", ((2, 0), (2, 1))),
        ("t", ((2, 1), (0, 1))),
        ("l", ((0, 1), (0, 0))),
    ):
        og2.create_line(a, b, name=name)
    for name, n in (("b1", 6), ("b2", 6), ("r", 6), ("t", 11), ("l", 6)):
        og2.create_connector(name, n, name=f"c-{name}")
    og2.create_domain([["c-b1", "c-b2"], ["c-r"], ["c-t"], ["c-l"]], name="d")
    og2.set_boundary_type(["c-b1"], "wall")
    og2.set_boundary_type(["c-b2"], "inlet")
    with pytest.raises(ValidationError, match="split the domain"):
        og2.extrude_model(0.1)
