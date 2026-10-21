import math

import numpy as np
import pytest
from helpers import naca0012_selig

import opengrid
from opengrid import _core
from opengrid.project import to_json


def _airfoil_session() -> opengrid.Session:
    og = opengrid.Session()
    points = naca0012_selig(61)
    le = 60
    upper = list(reversed(points[: le + 1]))
    lower = points[le:]
    og.create_spline(upper, method="clamped", start_tangent=(0, 1), name="upper")
    og.create_spline(lower, method="clamped", start_tangent=(0, -1), name="lower")
    og.create_line((1.0, 0.0), (5.0, 0.0), name="wake")
    return og


def test_distribute_binding():
    s, warnings = _core.distribute(41, "tanh", 1e-4, 1e-2)
    assert s[0] == 0.0
    assert s[-1] == 1.0
    assert np.all(np.diff(s) > 0)
    assert abs((s[1] - s[0]) - 1e-4) / 1e-4 < 1e-8
    assert abs((s[-1] - s[-2]) - 1e-2) / 1e-2 < 1e-8
    assert warnings == []
    with pytest.raises(ValueError, match="unknown distribution type"):
        _core.distribute(5, "cosine")


@pytest.mark.parametrize("kind", ["tanh", "sine"])
def test_achieved_spacings_through_api(kind):
    og = _airfoil_session()
    name = og.create_connector(
        "upper", 121, {"type": kind, "begin_spacing": 1e-4, "end_spacing": 2e-3}
    )
    nodes = og.connector_nodes(name)
    assert len(nodes.points) == 121
    assert abs(nodes.begin_spacing - 1e-4) / 1e-4 < 1e-8
    assert abs(nodes.end_spacing - 2e-3) / 2e-3 < 1e-8
    assert np.all(np.diff(nodes.arc_positions) > 0)
    assert nodes.warnings == []
    assert nodes.points[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-15)
    assert nodes.points[-1] == pytest.approx([1.0, 0.0, 0.0], abs=1e-12)
    # Nodes lie on the curve.
    curve = og.geometry("upper")
    for point in nodes.points[::10]:
        assert curve.project(tuple(point))[2] < 1e-12


def test_geometric_and_uniform():
    og = _airfoil_session()
    g = og.create_connector("wake", 41, {"type": "geometric", "begin_spacing": 0.01})
    nodes = og.connector_nodes(g)
    spacing = np.diff(nodes.arc_positions)
    ratios = spacing[1:] / spacing[:-1]
    assert np.allclose(ratios, ratios[0], rtol=1e-9)
    assert abs(nodes.begin_spacing - 0.01) / 0.01 < 1e-10
    u = og.create_connector("wake", 5)
    assert np.allclose(np.diff(og.connector_nodes(u).arc_positions), 1.0)


def test_chain_connector_ordering():
    og = _airfoil_session()
    # Lower surface runs LE -> TE, upper LE -> TE: the chain must reverse one of them.
    name = og.create_chain_connector(["upper", "lower", "wake"], 201)
    segments = og.entity(name).params["segments"]
    ids = {og.entity(n).id: n for n in ("upper", "lower", "wake")}
    order = [(ids[s["curve"]], s["t0"], s["t1"]) for s in segments]
    assert order[0] == ("wake", 1.0, 0.0)  # starts at the free end of the wake
    assert [o[0] for o in order] == ["wake", "upper", "lower"] or [o[0] for o in order] == [
        "wake",
        "lower",
        "upper",
    ]
    nodes = og.connector_nodes(name)
    total = og.length("upper") + og.length("lower") + og.length("wake")
    assert nodes.length == pytest.approx(total, rel=1e-12)
    assert nodes.points[0] == pytest.approx([5.0, 0.0, 0.0])


def test_breakpoints():
    og = _airfoil_session()
    name = og.create_connector(
        "wake",
        51,
        distribution={"type": "tanh", "begin_spacing": 0.002, "end_spacing": 0.05},
        breakpoints=[
            {
                "position": 0.25,
                "index": 20,
                "distribution": {"type": "geometric", "begin_spacing": 0.05},
            }
        ],
    )
    nodes = og.connector_nodes(name)
    s = nodes.arc_positions
    assert s[20] == pytest.approx(1.0, abs=1e-14)
    assert abs((s[20] - s[19]) - 0.05) / 0.05 < 1e-8
    assert abs((s[21] - s[20]) - 0.05) / 0.05 < 1e-8
    assert abs((s[1] - s[0]) - 0.002) / 0.002 < 1e-8


def test_parametric_mode():
    og = _airfoil_session()
    a = og.connector_nodes(og.create_connector("wake", 11, "uniform"))
    p = og.connector_nodes(og.create_connector("wake", 11, "uniform", mode="parametric"))
    assert np.allclose(a.points, p.points, atol=1e-13)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"segments": "missing", "dimension": 5}, "no entity 'missing'"),
        ({"segments": "pt", "dimension": 5}, "is a point, not a curve"),
        ({"segments": "wake", "dimension": 1}, "dimension"),
        ({"segments": "wake", "dimension": 5, "mode": "chord"}, "mode"),
        ({"segments": [("wake", 0.5, 0.5)], "dimension": 5}, "t0, t1"),
        ({"segments": ["wake", "upper"], "dimension": 5}, r"segments\[1\] does not start"),
        ({"segments": "wake", "dimension": 5, "distribution": "cosine"}, "distribution.type"),
        (
            {
                "segments": "wake",
                "dimension": 5,
                "distribution": {"type": "uniform", "begin_spacing": 0.1},
            },
            "uniform takes no",
        ),
        (
            {"segments": "wake", "dimension": 5, "distribution": {"type": "geometric"}},
            "exactly one",
        ),
        (
            {
                "segments": "wake",
                "dimension": 5,
                "distribution": {"type": "tanh", "growth_ratio": 1.1},
            },
            "not used by tanh",
        ),
        (
            {
                "segments": "wake",
                "dimension": 5,
                "distribution": {"type": "tanh", "begin_spacing": -1},
            },
            "positive",
        ),
        (
            {"segments": "wake", "dimension": 5, "breakpoints": [{"position": 0.5, "index": 4}]},
            r"index 4 must increase",
        ),
        (
            {"segments": "wake", "dimension": 9, "breakpoints": [{"position": 1.5, "index": 4}]},
            "position 1.5",
        ),
        ({"segments": "wake", "dimension": 9, "breakpoints": [{"index": 4}]}, "needs 'position'"),
    ],
)
def test_connector_validation(kwargs, message):
    og = _airfoil_session()
    og.create_point((0, 0), name="pt")
    with pytest.raises(opengrid.ValidationError, match=message):
        og.create_connector(**kwargs)


def test_infeasible_spacing_warns_never_nan():
    og = _airfoil_session()
    name = og.create_connector(
        "wake", 11, {"type": "tanh", "begin_spacing": 3.0, "end_spacing": 2.0}
    )
    nodes = og.connector_nodes(name)
    assert nodes.warnings
    assert np.all(np.isfinite(nodes.points))
    assert np.all(np.diff(nodes.arc_positions) > 0)
    assert "WARNING" in og.connector_report(name)


def test_spacing_propagation_and_undo():
    og = _airfoil_session()
    og.create_connector(
        "upper", 81, {"type": "tanh", "begin_spacing": 1e-3, "end_spacing": 5e-4}, name="c-upper"
    )
    wake = og.create_connector("wake", 41, {"type": "tanh", "end_spacing": 0.5}, name="c-wake")
    nodes = og.connector_nodes(wake)
    assert abs(nodes.begin_spacing - 5e-4) / 5e-4 < 1e-8
    assert "c-upper end" in nodes.propagated["begin"]
    # Changing the source spacing updates the neighbour; undo restores it.
    og.set_distribution("c-upper", "tanh", begin_spacing=1e-3, end_spacing=2e-4)
    assert abs(og.connector_nodes(wake).begin_spacing - 2e-4) / 2e-4 < 1e-8
    og.undo()
    assert abs(og.connector_nodes(wake).begin_spacing - 5e-4) / 5e-4 < 1e-8
    # Propagation can be switched off per connector.
    og.modify(wake, propagate_spacing=False)
    assert og.connector_nodes(wake).propagated == {}
    # Several sources: the smallest spacing wins.
    og.modify(wake, propagate_spacing=True)
    og.create_connector(
        "lower", 81, {"type": "tanh", "begin_spacing": 1e-3, "end_spacing": 1e-4}, name="c-lower"
    )
    assert abs(og.connector_nodes(wake).begin_spacing - 1e-4) / 1e-4 < 1e-8


def test_dependencies_protect_connectors():
    og = _airfoil_session()
    og.create_connector("wake", 11, name="c-wake")
    with pytest.raises(opengrid.ValidationError, match="used by connector 'c-wake'"):
        og.delete("wake")
    og.delete(["c-wake", "wake"])  # together is fine
    og.undo()
    chain = og.create_chain_connector(["upper", "wake"], 31)
    with pytest.raises(opengrid.ValidationError, match="would become invalid"):
        og.modify("wake", p0=(1.5, 0.0))  # breaks the chain upper -> wake
    assert og.entity("wake").params["p0"] == (1.0, 0.0, 0.0)
    og.modify("wake", p1=(6.0, 0.0))  # still connected: connector follows
    assert og.connector_nodes(chain).points[-1] == pytest.approx([6.0, 0.0, 0.0])
    og.rename("wake", "wake-line")
    assert og.connector_nodes(chain).points[-1] == pytest.approx([6.0, 0.0, 0.0])


def test_connectors_round_trip_and_replay(tmp_path):
    og = _airfoil_session()
    og.create_connector("upper", 81, {"type": "sine", "begin_spacing": 1e-3, "end_spacing": 5e-4})
    og.create_connector(
        [("wake", 1.0, 0.0)],
        31,
        {"type": "geometric", "growth_ratio": 1.1},
        breakpoints=[{"position": 0.5, "index": 10}],
        mode="parametric",
        propagate_spacing=False,
    )
    og.create_connectors_on_curves(["lower"], 21, "uniform")
    og.rename("upper", "suction")  # connector references survive renames (by id)
    path = tmp_path / "c.ogp"
    og.save(path)
    other = opengrid.Session()
    other.load(path)
    assert to_json(other.store) == to_json(og.store)
    for name in og.names("connector"):
        assert np.array_equal(other.connector_nodes(name).points, og.connector_nodes(name).points)
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)
    assert '"connectors": [' in path.read_text()


def test_undo_redo_connector_commands():
    og = _airfoil_session()
    before = og.store.snapshot()
    names = og.create_connectors_on_curves(
        ["upper", "lower"], 31, {"type": "tanh", "begin_spacing": 1e-3}
    )
    assert names == ["curve-1", "curve-2"]  # automatic names of points: curves
    og.undo()
    assert og.store.snapshot() == before
    og.redo()
    og.set_dimension("curve-1", 51)
    assert len(og.connector_nodes("curve-1").points) == 51
    og.undo()
    assert len(og.connector_nodes("curve-1").points) == 31


def test_stretching_ratio_report():
    og = _airfoil_session()
    name = og.create_connector("wake", 11, {"type": "geometric", "growth_ratio": 1.5})
    nodes = og.connector_nodes(name)
    assert nodes.max_stretching_ratio == pytest.approx(1.5, rel=1e-9)
    assert np.all(nodes.interval_ratios() >= 1.0)
    assert "exceeds 1.3" in og.connector_report(name)
    og.stretching_threshold = 2.0
    assert "exceeds" not in og.connector_report(name)


def test_closed_circle_connector():
    og = opengrid.Session()
    og.create_circle((0.5, 0), 20.0, name="farfield")
    nodes = og.connector_nodes(og.create_connector("farfield", 101))
    spacing = np.diff(nodes.arc_positions)
    assert np.allclose(spacing, 40.0 * math.pi / 100, rtol=1e-12)
    assert np.allclose(nodes.points[0], nodes.points[-1], atol=1e-12)
