"""Phase 10a: smoothing groups (multi-block elliptic smoothing as one grid)."""

import hashlib
import math

import numpy as np
import pytest

import opengrid
from opengrid.model.groups import interface_kinks
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def bump(og, ni=41, nj=21):
    xs = np.linspace(0.0, 3.0, 41)
    og.create_spline([(x, 0.25 * math.sin(math.pi * x / 3) ** 4, 0.0) for x in xs], name="b")
    og.create_line((3, 0), (3, 1), name="r")
    og.create_line((3, 1), (0, 1), name="t")
    og.create_line((0, 1), (0, 0), name="l")
    og.create_connector("b", ni, name="cb")
    og.create_connector("t", ni, name="ct")
    up = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.1}
    down = {"type": "tanh", "begin_spacing": 0.1, "end_spacing": 0.01}
    og.create_connector("r", nj, up, name="cr")
    og.create_connector("l", nj, down, name="cl")
    return og.create_domain([["cb"], ["cr"], ["ct"], ["cl"]], name="d")


def circle_ogrid(og, n=17, fraction=0.55):
    for k in range(4):
        a0 = -3 * math.pi / 4 + k * math.pi / 2
        og.create_arc((0, 0, 0), 1.0, a0, a0 + math.pi / 2, name=f"a{k}")
        og.create_connector(f"a{k}", n, name=f"c{k}")
    og.create_domain([[f"c{k}"] for k in range(4)], name="pipe")
    return og.ogrid_domain("pipe", fraction=fraction, radial_points=9, radial_spacing=0.03)


def replays(og) -> bool:
    return to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_split_domain_smoothed_as_a_group_matches_the_single_domain():
    tolerance = 1e-11
    single = opengrid.Session()
    bump(single)
    single.set_domain_generator(
        "d",
        elliptic={
            "control": "laplace",
            "walls": [{"edge": 0, "spacing": 0.01}],
            "tolerance": tolerance,
            "max_iterations": 20000,
        },
    )
    reference = single.domain_grid("d")
    og = opengrid.Session()
    bump(og)
    low, high, _cut = og.split_domain("d", i=17)
    walls = [{"connector": c, "spacing": 0.01} for c in ("cb", "cb-b")]
    group = og.create_smoothing_group(
        [low, high], walls=walls, slide=(), tolerance=tolerance, max_iterations=20000
    )
    result = og.smooth_group(group)
    assert result.converged
    assert result.node_counts["irregular"] == 0
    assert result.roles[og.entity("d-cut").id] == "interface"
    whole = np.concatenate([og.domain_grid(low), og.domain_grid(high)[:, 1:]], axis=1)
    assert np.abs(whole - reference).max() < 1e-9
    # The interface leaves no trace: shared nodes are identical in both parts.
    assert np.array_equal(og.domain_grid(low)[:, -1], og.domain_grid(high)[:, 0])
    assert result.max_kink < 0.1  # what the single-domain grid itself has there
    assert "in smoothing group" in og.domain_status(low)


def test_ogrid_group_smooths_across_interfaces_and_slides_on_the_circle():
    og = opengrid.Session()
    names = circle_ogrid(og)
    before = {og.entity(n).id: og.domain_grid(n).copy() for n in names}
    worst_before = min(og.domain_quality(n)["scaled_jacobian"].min() for n in names)
    group = og.create_smoothing_group(names, tolerance=1e-10)
    result = og.smooth_group(group)
    assert result.converged
    assert result.node_counts["irregular"] == 4
    assert result.node_counts["sliding"] == 4 * 15
    params = og.entity(group).params
    kink_before = max(interface_kinks(og.store, params, before, result.singular).values())
    assert kink_before > 10.0
    assert result.max_kink < 0.5
    worst_after = min(og.domain_quality(n)["scaled_jacobian"].min() for n in names)
    assert worst_after >= worst_before
    points = np.concatenate([og.domain_grid(n).reshape(-1, 3) for n in names])
    radius = np.linalg.norm(points[:, :2], axis=1)
    assert radius.max() <= 1.0 + 1e-12
    on_circle = radius > 0.999
    assert on_circle.sum() > 0
    assert np.abs(radius[on_circle] - 1.0).max() < 1e-10
    # The result depends on the topology, not on the initial core size.
    other = opengrid.Session()
    other_names = circle_ogrid(other, fraction=0.4)
    other.smooth_group(other.create_smoothing_group(other_names, tolerance=1e-10))
    for a, b in zip(names, other_names, strict=True):
        assert np.abs(og.domain_grid(a) - other.domain_grid(b)).max() < 1e-8
    assert "4 irregular" in og.group_status(group)


def test_group_is_thread_independent():
    digests = []
    before = opengrid.threads()
    try:
        for threads in (1, 4):
            opengrid.set_threads(threads)
            og = opengrid.Session()
            names = circle_ogrid(og)
            og.smooth_group(og.create_smoothing_group(names))
            points = np.concatenate([og.domain_grid(n) for n in names])
            digests.append(hashlib.sha256(points.tobytes()).hexdigest())
    finally:
        opengrid.set_threads(before)
    assert digests[0] == digests[1]


def test_group_roles_undo_journal_and_project(tmp_path):
    og = opengrid.Session()
    bump(og)
    low, high, cut = og.split_domain("d", i=17)
    group = og.create_smoothing_group([low, high], walls=["cb", "cb-b"], name="g")
    roles = {og.store.get(c).name: role for c, role in og.smooth_group(group).roles.items()}
    assert roles == {
        "cb": "wall",
        "cb-b": "wall",
        "ct": "sliding",
        "ct-b": "sliding",
        "cr": "sliding",
        "cl": "sliding",
        cut: "interface",
    }
    og.modify("g", slide=("ct",))
    roles = {og.store.get(c).name: r for c, r in og.smooth_group("g").roles.items()}
    assert roles["ct"] == "sliding" and roles["ct-b"] == "fixed" and roles["cr"] == "fixed"
    assert replays(og)
    path = tmp_path / "group.ogp"
    og.save(path)
    loaded = opengrid.Session()
    loaded.load(path)
    assert to_json(loaded.store) == to_json(og.store)
    assert np.array_equal(loaded.domain_grid(low), og.domain_grid(low))
    og.undo()
    og.undo()
    assert "g" not in og.store
    assert og.domain_status(low).startswith("tfi")


def test_connectors_used_outside_the_group_stay_fixed():
    og = opengrid.Session()
    bump(og)
    low, high, cut = og.split_domain("d", i=17)
    group = og.create_smoothing_group([low], name="g")
    result = og.smooth_group(group)
    assert result.roles[og.entity(cut).id] == "outside"
    assert np.array_equal(og.domain_grid(low)[:, -1], og.domain_grid(high)[:, 0])


def test_group_validation():
    og = opengrid.Session()
    bump(og)
    low, high, _cut = og.split_domain("d", i=17)
    og.create_smoothing_group([low], name="g0")
    with pytest.raises(ValidationError, match="unexpected key"):
        og.modify("g0", colour="red")
    og.delete("g0")
    with pytest.raises(ValidationError, match="not a domain"):
        og.create_smoothing_group(["cb"])
    with pytest.raises(ValidationError, match="not used by a member domain"):
        og.create_smoothing_group([low], walls=["cb-b"])
    with pytest.raises(ValidationError, match="is a wall"):
        og.create_smoothing_group([low], walls=["cb"], slide=("cb",))
    with pytest.raises(ValidationError, match="at least one domain"):
        og.create_smoothing_group([])
    with pytest.raises(ValidationError, match="acceleration"):
        og.create_smoothing_group([low], acceleration=-1)
    og.create_smoothing_group([low], name="g1")
    og.create_smoothing_group([low, high], name="g2")
    with pytest.raises(ValidationError, match="several smoothing groups"):
        og.domain_grid(low)
    og.delete(high)  # a group member: the group loses it (it keeps low)
    assert og.entity("g2").params["domains"] == (og.entity(low).id,)


def test_o_type_domain_smooths_across_its_cut():
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="inner")
    og.create_circle((0, 0), 6.0, name="outer")
    og.create_line((1, 0), (6, 0), name="cut")
    og.create_connector("inner", 41, name="c-inner")
    og.create_connector("outer", 41, name="c-outer")
    og.create_connector("cut", 21, {"type": "tanh", "begin_spacing": 0.02}, name="c-cut")
    og.create_structured_domain(connectors=["c-inner", "c-outer", "c-cut"], name="ring")
    group = og.create_smoothing_group(["ring"], walls=["c-inner"], slide=())
    result = og.smooth_group(group)
    assert result.converged
    assert result.roles[og.entity("c-cut").id] == "interface"
    grid = og.domain_grid("ring")
    # The cut is a grid line like any other: its two copies stay identical and its nodes moved
    # off the straight connector (smooth periodic grid lines across it).
    if np.array_equal(grid[:, 0], grid[:, -1]):
        line = grid[:, 0]
    else:
        assert np.array_equal(grid[0], grid[-1])
        line = grid[0]
    cut = og.connector_nodes("c-cut").points
    offsets = np.min(np.linalg.norm(line[:, None, :] - cut[None, :, :], axis=2), axis=1)
    assert offsets.max() > 1e-3
    assert result.max_kink < 0.5


def test_gui_group_action_dock_and_panel(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    names = circle_ogrid(og)
    ids = [og.entity(n).id for n in names]
    window.set_selection(ids)
    name = window.create_group_on_selection()
    assert name == "smoothing_group-1"
    window.smoother_dock.wait()
    qtbot.waitUntil(lambda: not window.smoother_dock.running, timeout=30000)
    assert "converged" in window.smoother_dock.status.text()
    assert "converged" in og.group_status(name)
    assert name in window.model_tree.entity_names("Domains")
    # Selecting a member shows its group in the dock.
    window.set_selection([ids[0]])
    assert window.smoother_dock.domain_name == name
    window.set_selection([og.entity(name).id])
    panel = window.property_panel
    assert panel.read_changes() == (None, {})
    panel._editors["slide"][1].setText("none")
    assert panel.apply()
    assert og.entity(name).params["slide"] == ()
    window.close()


def test_groups_follow_topology_edits():
    og = opengrid.Session()
    bump(og)
    low, high, _cut = og.split_domain("d", i=17)
    og.create_smoothing_group([low, high], walls=[{"connector": "cb", "spacing": 0.01}], name="g")
    third = og.split_domain(high, i=10)[1]
    members = [og.store.get(d).name for d in og.entity("g").params["domains"]]
    assert members == [low, high, third]
    parts = og.split_connector("cb", 8)
    walls = [og.store.get(w["connector"]).name for w in og.entity("g").params["walls"]]
    assert sorted(walls) == sorted(parts)
    assert all(w["spacing"] == 0.01 for w in og.entity("g").params["walls"])
    og.join_connectors(*parts)
    assert len(og.entity("g").params["walls"]) == 1
    created = og.ogrid_domain(low, fraction=0.5)
    members = [og.store.get(d).name for d in og.entity("g").params["domains"]]
    assert sorted(members) == sorted([high, third, *created])
    assert og.smooth_group("g").converged
    assert replays(og)
