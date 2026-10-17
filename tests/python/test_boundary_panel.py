"""The Boundary panel and what feeds it: periodic mismatch reasons in plain words, pair
names that belong together, interfaces refused as patch members, Clear type, the tree's patch
labels, the live rows of the panel and grid snapping for exact sides."""

import math

import pytest

import opengrid
from opengrid.validation import ValidationError


def quad(og, right_height=1.0, left_points=11, right_points=11, left_distribution=None):
    """An unstructured quad [0, 2] x [0, 1] whose right side may be taller (a trapezoid) or
    have another point count or distribution than the left: bottom, right, top, left."""
    corners = ((0, 0), (2, 0), (2, right_height), (0, 1))
    for k, name in enumerate(("bottom", "right", "top", "left")):
        og.create_line(corners[k], corners[(k + 1) % 4], name=name)
    counts = {"bottom": 21, "right": right_points, "top": 21, "left": left_points}
    for name, n in counts.items():
        kwargs = {}
        if name == "left" and left_distribution is not None:
            kwargs["distribution"] = left_distribution
        og.create_connector(name, n, name=f"c-{name}", **kwargs)
    og.create_unstructured_domain(connectors=["c-bottom", "c-right", "c-top", "c-left"], name="d")
    return og


def pair_issue(og):
    report = og.boundary_report()
    assert not report.ok
    (issue,) = [i for i in report.issues if "no periodic image" in i.message]
    return issue.message


def test_mismatch_reasons_in_plain_words():
    og = quad(opengrid.Session(), left_points=11, right_points=13)
    og.make_periodic_pair(["c-left"], ["c-right"])
    assert "has 11 points and the other side's connectors have 13" in pair_issue(og)
    og = quad(opengrid.Session(), right_height=1.2)
    og.make_periodic_pair(["c-left"], ["c-right"])
    message = pair_issue(og)
    assert "not images of each other under the translation" in message
    assert "snap to the reference grid" in message
    og = quad(
        opengrid.Session(),
        left_distribution={"type": "tanh", "begin_spacing": 0.02, "end_spacing": 0.02},
    )
    og.make_periodic_pair(["c-left"], ["c-right"])
    message = pair_issue(og)
    assert "ends coincide" in message and "same distribution" in message
    og = quad(opengrid.Session())
    og.set_boundary_type(["c-bottom", "c-top"], "wall")
    og.make_periodic_pair(["c-left"], ["c-right"])
    assert og.boundary_report().ok


def test_pair_names_belong_together():
    og = quad(opengrid.Session())
    assert og.make_periodic_pair(["c-left"], ["c-right"]) == ["periodic-1a", "periodic-1b"]
    assert og.make_periodic_pair(["c-bottom"], ["c-top"]) == ["periodic-2a", "periodic-2b"]
    og.undo()
    assert og.make_periodic_pair(["c-bottom"], ["c-top"], name="span") == ["span-a", "span-b"]
    og.undo()
    og.create_patch("span-a", "wall", ["c-bottom"])
    assert og.make_periodic_pair(["c-bottom"], ["c-top"], name="span") == ["span-2a", "span-2b"]
    from opengrid.model.patches import partner_of

    assert partner_of(og.entity("span-2b"), og.store) == "span-2a"
    assert partner_of(og.entity("periodic-1a"), og.store) == "periodic-1b"
    og.create_patch("lonely", "periodic", [])
    assert partner_of(og.entity("lonely"), og.store) is None


def test_shared_connectors_are_refused_as_patch_members():
    og = quad(opengrid.Session())
    og.create_line((2, 1), (2, 2), name="r2")  # a second quad on top shares c-top
    og.create_line((2, 2), (0, 2), name="t2")
    og.create_line((0, 2), (0, 1), name="l2")
    for name in ("r2", "t2", "l2"):
        og.create_connector(name, 11, name=f"c-{name}")
    og.create_unstructured_domain(connectors=["c-top", "c-r2", "c-t2", "c-l2"], name="d2")
    for call in (
        lambda: og.set_boundary_type(["c-top"], "interface"),
        lambda: og.create_patch("x", "wall", ["c-top"]),
        lambda: og.make_periodic_pair(["c-top"], ["c-bottom"]),
    ):
        with pytest.raises(ValidationError, match=r"shared by d and d2.*takes no boundary type"):
            call()
    og.create_patch("walls", "wall", ["c-bottom"])
    with pytest.raises(ValidationError, match="interface between them"):
        og.add_to_patch("walls", ["c-top"])
    assert og.names("patch") == ["walls"]
    report = og.boundary_report()
    assert report.interfaces == ["c-top (d | d2)"]
    assert "c-top" not in report.boundary


def test_clear_boundary_type_and_tree_labels():
    og = quad(opengrid.Session())
    og.set_boundary_type(["c-bottom", "c-top"], "wall")
    og.make_periodic_pair(["c-left"], ["c-right"])
    assert og.boundary_report().ok
    assert og.clear_boundary_type(["c-bottom"]) == ["wall"]
    assert og.unassigned_boundary_connectors() == ["c-bottom"]
    og.undo()
    assert og.boundary_report().ok
    og.redo()
    assert og.entity("wall").params["members"] == ({"connector": og.entity("c-top").id},)
    from opengrid.app.model_tree import patch_summary

    assert patch_summary(og.entity("wall"), og.store) == "wall, 1 connector"
    assert (
        patch_summary(og.entity("periodic-1a"), og.store) == "periodic ↔ periodic-1b, 1 connector"
    )
    assert (
        patch_summary(og.entity("periodic-1b"), og.store) == "periodic ↔ periodic-1a, 1 connector"
    )
    og.create_patch("lonely", "periodic", [])
    assert patch_summary(og.entity("lonely"), og.store) == "periodic, NOT PAIRED, 0 connectors"
    from opengrid.project.ogp import to_json

    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


@pytest.fixture
def window(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    win = MainWindow(off_screen=True)
    qtbot.addWidget(win)
    quad(win.session)
    yield win
    win.close()


def _ids(window, *names):
    return [window.session.entity(n).id for n in names]


def test_boundary_panel_follows_the_model(window, qtbot):
    dock = window.boundary_dock
    og = window.session
    qtbot.waitUntil(lambda: "4 of 4 boundaries have no type" in dock.status.text(), timeout=3000)
    rows = dock.rows("Boundaries")
    assert [row[0] for row in rows] == ["c-bottom", "c-right", "c-top", "c-left"]
    assert all(row[1] == "none" and row[3].startswith("no type yet") for row in rows)
    # Set type from the panel: a named patch; without a name, the existing patch of the type.
    window.set_selection(_ids(window, "c-bottom"))
    dock.patch_type.setCurrentText("wall")
    dock.name.setText("floor")
    dock.set_type.click()
    qtbot.waitUntil(lambda: dock.rows("Boundaries")[0][1] == "floor", timeout=3000)
    assert dock.rows("Boundaries")[0][2:] == ["wall", "ok"]
    assert window.patch_color_action.isChecked()
    dock.name.setText("")
    window.set_selection(_ids(window, "c-top"))
    dock.set_type.click()
    qtbot.waitUntil(lambda: dock.rows("Boundaries")[2][1] == "floor", timeout=3000)
    assert og.names("patch") == ["floor"]
    # Pair the sides: the pair row says the transform and that the nodes match.
    window.set_selection(_ids(window, "c-left", "c-right"))
    dock.pair.click()
    qtbot.waitUntil(lambda: "ready to export" in dock.status.text(), timeout=3000)
    (pair_row,) = dock.rows("Periodic pairs")
    assert pair_row[0] == "periodic-1a ↔ periodic-1b"
    assert pair_row[1] == "translation (2, 0, 0)"
    assert pair_row[3] == "match: 1 connector node by node"
    assert "All 4 boundaries have a type, 1 periodic pair(s) match" in dock.status.text()
    # Clicking a row selects the boundary; the pair row selects the source patch.
    dock._on_item_clicked(dock.tree.topLevelItem(0).child(1), 0)
    assert window.selected_names() == ["c-right"]
    dock._on_item_clicked(dock.tree.topLevelItem(1).child(0), 0)
    assert window.selected_names() == ["periodic-1a"]
    # Clear type: the boundary is unassigned again; an empty patch is listed as such.
    window.set_selection(_ids(window, "c-bottom"))
    dock.clear.click()
    qtbot.waitUntil(lambda: "1 of 4 boundaries have no type" in dock.status.text(), timeout=3000)
    assert dock.rows("Boundaries")[0][3].startswith("no type yet")
    og.create_patch("spare", "inlet", [])
    qtbot.waitUntil(lambda: dock.rows("Empty patches") != [], timeout=3000)
    assert dock.rows("Empty patches")[0][:3] == ["spare", "", "inlet"]
    og.undo()
    og.undo()
    qtbot.waitUntil(lambda: "ready to export" in dock.status.text(), timeout=3000)
    # A mismatched pair is said so, in the status and in the row.
    og.undo()  # the pairing
    og.modify("c-right", dimension=13)
    og.make_periodic_pair(["c-left"], ["c-right"])
    qtbot.waitUntil(lambda: "1 periodic pair(s) do not match" in dock.status.text(), timeout=3000)
    assert "same number of points" in dock.rows("Periodic pairs")[0][3]
    assert dock.rows("Boundaries")[1][3] == "ok"  # the connector itself is assigned
    # A 3D model: the panel waits for Refresh.
    og.undo()
    og.undo()
    og.make_periodic_pair(["c-left"], ["c-right"])
    assert window.extrude_model(thickness=0.1, ends="periodic")
    qtbot.waitUntil(lambda: "press Refresh" in dock.status.text(), timeout=3000)
    assert dock.stale
    dock.refresh.click()
    assert not dock.stale and "ready to export" in dock.status.text()
    assert [row[0] for row in dock.rows("Periodic pairs")] == [
        "periodic-1a ↔ periodic-1b",
        "front ↔ back",
    ]


def test_clicks_snap_to_the_reference_grid_in_2d(window, monkeypatch):
    # Offscreen the renderer is one pixel high, so every entity would be within snapping
    # reach: give the view a sane scale instead.
    monkeypatch.setattr(window.scene, "world_per_pixel", lambda: 0.001)
    point, target = window.snap_point((0.123, 0.456, 0.0))
    assert target is None and point == (0.123, 0.456, 0.0)
    window.set_grid(True)
    grid = window.scene.grid
    step = grid.step / grid.divisions
    point, target = window.snap_point((0.123, 0.456, 0.0))
    assert target is not None and target.label.startswith("grid")
    assert math.isclose(point[0], round(0.123 / step) * step, abs_tol=1e-12)
    assert math.isclose(point[1], round(0.456 / step) * step, abs_tol=1e-12)
    # An entity's point wins over the grid.
    point, target = window.snap_point((2.004, 0.0, 0.0))  # within 10 px * 0.001 of the end
    assert target.label == "bottom end" and point == (2.0, 0.0, 0.0)
    window.session.extrude_model(0.1)
    window._auto_navigation()
    assert window.navigation_mode == "3d"
    point, target = window.snap_point((0.123, 0.456, 0.0))
    assert target is None
