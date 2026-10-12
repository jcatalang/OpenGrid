"""Topology editing (Phase 9b): split / join of connectors, domains and blocks, O-grids."""

import itertools
import math

import numpy as np
import pytest
from test_volumes import box

import opengrid
from opengrid.project import to_json
from opengrid.validation import ValidationError


def rectangle(og, n_bottom=21, n_side=11, distribution=None):
    og.create_line((0, 0), (2, 0), name="b")
    og.create_line((2, 0), (2, 1), name="r")
    og.create_line((2, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    if distribution is None:
        distribution = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.2}
    og.create_connector("b", n_bottom, distribution, name="cb")
    og.create_connector("r", n_side, name="cr")
    og.create_connector("t", n_bottom, name="ct")
    og.create_connector("l", n_side, name="cl")
    return og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="d")


def replays(og) -> bool:
    return to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def point_set(grid) -> np.ndarray:
    return np.array(sorted(map(tuple, np.round(np.asarray(grid).reshape(-1, 3), 12))))


def test_split_and_join_connector_keep_nodes():
    og = opengrid.Session()
    rectangle(og)
    before = og.connector_nodes("cb").points.copy()
    grid = og.domain_grid("d").copy()
    assert og.split_connector("cb", 7) == ["cb", "cb-b"]
    a = og.connector_nodes("cb").points
    b = og.connector_nodes("cb-b").points
    assert len(a) == 8 and len(b) == 14
    assert np.abs(np.concatenate([a, b[1:]]) - before).max() < 1e-13
    edge0 = og.entity("d").params["edges"][0]
    assert [og.store.get(i["connector"]).name for i in edge0] == ["cb", "cb-b"]
    assert np.abs(og.domain_grid("d") - grid).max() < 1e-13
    assert og.join_connectors("cb-b", "cb") == "cb-b"  # joined into the first
    assert "cb" not in og.store
    assert np.abs(og.connector_nodes("cb-b").points - before).max() < 1e-13
    assert replays(og)
    og.undo()
    og.undo()
    assert np.array_equal(og.connector_nodes("cb").points, before)
    with pytest.raises(ValidationError):
        og.split_connector("cb", 0)
    with pytest.raises(ValidationError):
        og.join_connectors("cb", "ct")  # no shared end


def test_split_connector_updates_loops_inflation_and_patches():
    og = opengrid.Session()
    og.create_circle((0, 0, 0), 0.5, name="body")
    og.create_circle((0, 0, 0), 5.0, name="far")
    og.create_connector("body", 81, name="wall")
    og.create_connector("far", 41, name="outer")
    og.create_unstructured_domain(
        connectors=["outer", "wall"],
        generator={"inflation": {"connectors": ["wall"], "first_height": 1e-3, "layers": 5}},
        name="u",
    )
    og.create_patch("w", "wall", ["wall"])
    cells = len(og.unstructured_mesh("u")["quads"])
    og.split_connector("wall", 30)
    params = og.entity("u").params
    ids = {og.entity("wall").id, og.entity("wall-b").id}
    assert set(params["generator"]["inflation"]["connectors"]) == ids
    assert {i["connector"] for loop in params["loops"] for i in loop} >= ids
    assert {m["connector"] for m in og.entity("w").params["members"]} == ids
    assert len(og.unstructured_mesh("u")["quads"]) == cells
    og.join_connectors("wall", "wall-b")
    assert {m["connector"] for m in og.entity("w").params["members"]} == {og.entity("wall").id}
    assert replays(og)


def test_split_domain_along_i_and_j():
    og = opengrid.Session()
    rectangle(og)
    grid = og.domain_grid("d").copy()
    names = og.split_domain("d", i=5)
    assert names == ["d", "d-b", "d-cut"]
    low = og.domain_grid("d")
    high = og.domain_grid("d-b")
    assert low.shape == (11, 6, 3) and high.shape == (11, 16, 3)
    # Linear TFI on a rectangle is reproduced exactly by the parts.
    assert np.abs(np.concatenate([low, high[:, 1:]], axis=1) - grid).max() < 1e-13
    og.split_domain("d-b", j=4)
    assert og.domain_dimensions("d-b") == (16, 5)
    assert og.domain_dimensions("d-b-b") == (16, 7)
    assert replays(og)
    with pytest.raises(ValidationError):
        og.split_domain("d", i=0)
    with pytest.raises(ValidationError):
        og.split_domain("d", i=2, j=2)


def test_split_domain_keeps_walls_on_their_parts():
    og = opengrid.Session()
    rectangle(og)
    og.set_domain_generator("d", method="hermite_tfi", hermite_walls=(0, 2))
    og.split_domain("d", j=5)
    assert og.entity("d").params["generator"]["hermite_walls"] == (0,)
    assert og.entity("d-b").params["generator"]["hermite_walls"] == (2,)


def test_join_domains_back_and_errors():
    og = opengrid.Session()
    rectangle(og)
    grid = og.domain_grid("d").copy()
    og.split_domain("d", i=8)
    og.split_domain("d", j=3)
    assert og.join_domains("d", "d-b-2") == "d"
    assert og.join_domains("d", "d-b") == "d"
    assert np.abs(og.domain_grid("d") - grid).max() < 1e-13
    # The cut connectors were deleted with the joins.
    assert not [n for n in og.names("connector") if "cut" in n]
    assert replays(og)
    og.create_line((5, 0), (6, 0), name="x0")
    og.create_line((6, 0), (6, 1), name="x1")
    og.create_line((6, 1), (5, 1), name="x2")
    og.create_line((5, 1), (5, 0), name="x3")
    for k in range(4):
        og.create_connector(f"x{k}", 5, name=f"cx{k}")
    og.create_domain([["cx0"], ["cx1"], ["cx2"], ["cx3"]], name="far")
    with pytest.raises(ValidationError):
        og.join_domains("d", "far")


def test_join_adjacent_box_faces_in_every_orientation():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (2, 1, 1), 5, "o", "structured")
    for a, b in itertools.combinations(faces, 2):
        session = opengrid.replay(og.journal_text())
        try:
            joined = session.join_domains(a, b)
        except ValidationError as error:
            assert "shares no whole edge" in str(error)
            continue
        expected = np.unique(
            np.concatenate([point_set(og.domain_grid(a)), point_set(og.domain_grid(b))]), axis=0
        )
        assert np.array_equal(point_set(session.domain_grid(joined)), expected)


def test_split_and_join_block():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (2, 1, 1), 9, "o", "structured")
    og.create_block(faces=faces, name="blk")
    og.create_patch("walls", "wall", faces)
    grid = og.block_grid("blk").copy()
    for axis in ("i", "j", "k"):
        session = opengrid.replay(og.journal_text())
        low, high, middle = session.split_block("blk", axis, 3)
        assert session.entity(middle).kind == "domain"
        both = np.concatenate(
            [point_set(session.block_grid(low)), point_set(session.block_grid(high))]
        )
        assert np.array_equal(np.unique(both, axis=0), point_set(grid))
        cells = sum(np.prod(np.array(session.block_grid(n).shape[:3]) - 1) for n in (low, high))
        assert cells == 8**3
        assert session.boundary_report().ok
        assert session.join_blocks(low, high) == low
        assert np.array_equal(point_set(session.block_grid(low)), point_set(grid))
        assert len(session.names("domain")) == 6
        assert session.boundary_report().ok
        assert replays(session)


def circle_domain(og, n=17):
    for k in range(4):
        a0 = -3 * math.pi / 4 + k * math.pi / 2
        og.create_arc((0, 0, 0), 1.0, a0, a0 + math.pi / 2, name=f"a{k}")
        og.create_connector(f"a{k}", n, name=f"c{k}")
    return og.create_domain([[f"c{k}"] for k in range(4)], name="pipe")


def test_ogrid_domain_improves_a_circle():
    og = opengrid.Session()
    circle_domain(og)
    single = og.domain_quality("pipe")["scaled_jacobian"].min()
    og.create_zone("z", ["pipe"])
    names = og.ogrid_domain("pipe", fraction=0.55, radial_points=7, radial_spacing=0.02)
    assert names == ["pipe-0", "pipe-1", "pipe-2", "pipe-3", "pipe-core"]
    worst = min(og.domain_quality(n)["scaled_jacobian"].min() for n in names)
    assert single < 0.15 and worst > 0.7
    assert {m["domain"] for m in og.entity("z").params["members"]} == {
        og.entity(n).id for n in names
    }
    area = sum(og.domain_quality(n)["area"].sum() for n in names)
    polygon = 0.5 * 64 * math.sin(2 * math.pi / 64)  # 64 boundary segments
    assert abs(area - polygon) < 1e-9
    spacing = og.connector_nodes("pipe-r0").points
    assert abs(np.linalg.norm(spacing[1] - spacing[0]) - 0.02) < 1e-8
    assert replays(og)
    og.undo()
    assert og.names("domain") == ["pipe"]


def test_ogrid_domain_without_core_leaves_a_hole():
    og = opengrid.Session()
    circle_domain(og)
    names = og.ogrid_domain("pipe", fraction=0.4, core=False)
    assert len(names) == 4
    with pytest.raises(ValidationError):
        og.ogrid_domain("pipe-0", fraction=1.5)


def test_ogrid_block_all_hexahedral_ball():
    og = opengrid.Session()
    faces = box(og, (0, 0, 0), (1, 1, 1), 9, "o", "structured")
    og.create_block(faces=faces, name="cube")
    names = og.ogrid_block("cube", fraction=0.4, radial_points=5)
    assert len(names) == 7
    volume = 0.0
    for name in names:
        q = og.block_quality(name)
        assert q["scaled_jacobian"].min() > 0.3
        volume += q["volume"].sum()
    assert abs(volume - 1.0) < 1e-12
    assert replays(og)
    og2 = opengrid.Session()
    faces = box(og2, (0, 0, 0), (1, 1, 1), 9, "o", "structured")
    og2.create_block(faces=faces, name="cube")
    shell = og2.ogrid_block("cube", fraction=0.4, core=False)
    assert len(shell) == 6 and "cube-core" not in og2.store


def test_gui_topology_actions(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    rectangle(og)
    window.set_selection([og.entity("cb").id])
    assert window.split_connector(index=5) == ["cb", "cb-b"]
    window.set_selection([og.entity("cb").id, og.entity("cb-b").id])
    assert window.join_connectors() == "cb"
    window.set_selection([og.entity("d").id])
    assert window.split_domain("i", 6)[:2] == ["d", "d-b"]
    window.set_selection([og.entity("d").id, og.entity("d-b").id])
    assert window.join_domains() == "d"
    window.set_selection([og.entity("d").id])
    assert len(window.ogrid(0.5)) == 5
    window.set_selection([])
    assert window.join_blocks() is None
    window.close()
