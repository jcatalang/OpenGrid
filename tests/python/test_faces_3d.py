"""Boundary types in 3D from the window: a click on a block or volume selects the face under
the pointer (a patch member), which then gets its type from its property panel, the Boundary
panel or the right-click menu; the Boundary panel's rows select their faces; colouring by
patch fills the faces."""

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")


@pytest.fixture
def live(qtbot):
    from test_gui_events import Live

    live = Live(qtbot)
    yield live
    live.close()


def _face_at(og, block, axis: int, value: float) -> str:
    """The name of the block face lying on the plane x[axis] = value."""
    import numpy as np

    from opengrid.model.blocks import face_array
    from opengrid.model.patches import BLOCK_FACES

    grid = np.asarray(og.block_grid(block))
    for face in BLOCK_FACES:
        if np.allclose(face_array(grid, face)[..., axis], value):
            return face
    raise AssertionError("no such face")


# Clicks: on faces the isometric view sees (x = 2, y = 1, z = 1), at points whose view ray
# then meets no curve (a curve seen through a face is picked first, by design).
def _box(og):
    lines = og.create_rectangle(center=(1.0, 0.5), width=2.0, height=1.0, name="r")
    domain = og.structured_domain_from_curves(lines, 9)
    return og.create_extruded_block(domain, "translate", vector=(0, 0, 1), dimension=5)


def test_a_block_face_becomes_an_inlet(live):
    window, og = live.window, live.og
    block = _box(og)
    block_id = og.entity(block).id
    live.widget.render()
    inlet, outlet = _face_at(og, block, 0, 2.0), _face_at(og, block, 1, 1.0)
    live.click((2.0, 0.7, 0.3))  # the x = 2 face
    face = {"block": block_id, "face": inlet}
    assert window.selected_faces == [face]
    assert window.scene._overlays.get("selected-faces") is not None  # filled, outlined
    panel = window.property_panel
    assert "face" in panel._title.text() and f"{block}.{inlet}" in panel._title.text()
    _kind, combo = panel._editors["boundary"]
    combo.setCurrentText("inlet")
    assert panel.apply()
    assert og.entity("inlet").params["members"] == (face,)
    # The right-click menu offers boundary types for a selected face.
    labels = [a.text().replace("&", "") for a in window.context_menu().actions()]
    assert "Boundary type" in labels
    # Two faces at once from the Boundary panel's Set type: the outlet and a wall.
    live.click((1.3, 1.0, 0.3))  # the y = 1 face
    assert window.selected_faces == [{"block": block_id, "face": outlet}]
    window._set_type_from_dock("outlet", "")
    assert og.entity("outlet").params["members"] == ({"block": block_id, "face": outlet},)
    # The Boundary panel lists every face; clicking a row selects that face.
    window.validate_boundary()
    rows = window.boundary_dock.rows("Boundaries")
    assigned = sorted(r[0] for r in rows if r[1] != "none")
    assert assigned == sorted([f"{block}.{inlet}", f"{block}.{outlet}"])
    top_face = _face_at(og, block, 2, 1.0)
    group = window.boundary_dock.tree.topLevelItem(0)
    rows = [group.child(i) for i in range(group.childCount())]
    row = next(r for r in rows if r.text(0) == f"{block}.{top_face}")
    window.boundary_dock._on_item_clicked(row, 0)
    assert window.selected_faces == [{"block": block_id, "face": top_face}]
    # Colour by patch: the faces filled with their patches' colours (unassigned: red).
    window.patch_color_action.setChecked(True)
    coloured = {m["face"]: c for m, c in window.scene.patch_faces if "face" in m}
    assert coloured[inlet] != coloured[outlet]
    assert coloured[top_face] == "#ff3b30"  # no type yet
    assert window.scene._overlays.get("patch-faces") is not None
    window.patch_color_action.setChecked(False)
    assert window.scene._overlays.get("patch-faces") is None


def test_extruded_volume_faces(live):
    window, og = live.window, live.og
    lines = og.create_rectangle(center=(1.0, 0.5), width=2.0, height=1.0, name="r")
    domain = og.unstructured_domain_from_curves(lines, 9)
    volume = og.create_extruded_volume(domain, "translate", vector=(0, 0, 1), dimension=3)
    live.widget.render()
    live.click((1.0, 0.5, 1.0))  # the top
    assert window.selected_faces == [{"volume": og.entity(volume).id, "face": "top"}]
    live.click((2.0, 0.7, 0.3))  # a side wall: the curve it was extruded from
    side = window.selected_faces[0]
    assert set(side) == {"connector"}
    assert "side wall" in window.property_panel._title.text()
    window._set_type_from_dock("outlet", "")
    assert og.entity("outlet").params["members"] == (side,)


def test_faces_pair_as_periodic(live, tmp_path):
    """The user's case (2026-10-21): two side faces of an extruded box paired, then the
    bottom (its base domain) with the top; the Boundary panel follows and the pairs export."""
    window, og = live.window, live.og
    block = _box(og)
    block_id = og.entity(block).id
    live.widget.render()
    sides = [_face_at(og, block, 1, 0.0), _face_at(og, block, 1, 1.0)]
    window.set_selection([], faces=[{"block": block_id, "face": f} for f in sides])
    assert window.pair_selection_as_periodic()
    base = og.entity(block).params["extrusion"]["domain"]
    window.set_selection([base], faces=[{"block": block_id, "face": _face_at(og, block, 2, 1.0)}])
    assert window.pair_selection_as_periodic()
    report = og.boundary_report()
    assert len(report.periodic) == 2 and all(
        type(match).__name__ == "FaceMatch" for _s, _t, match in report.periodic
    )
    dock = window.boundary_dock
    live.qtbot.waitUntil(lambda: len(dock.rows("Periodic pairs")) == 2, timeout=3000)
    assert not dock.stale
    others = [_face_at(og, block, 0, 0.0), _face_at(og, block, 0, 2.0)]
    og.set_boundary_type([{"block": block, "face": f} for f in others], "wall")
    summary = og.export_cgns(tmp_path / "box.cgns")
    assert summary["connections"] == 4  # two pairs, both ways


def test_face_identifiers_in_the_api():
    """Every face of a block or volume has an identifier the API takes wherever it takes
    boundary members; side walls are named as faces, not by the curve they rise from."""
    import opengrid
    from opengrid.model.faces import all_faces, face_identifier, face_member

    og = opengrid.Session()
    base = og.unstructured_domain_from_curves(og.create_rectangle(width=2, height=1), 9, name="b")
    og.create_extruded_volume(base, "translate", vector=(0, 0, 0.5), dimension=3, name="vol")
    square = og.structured_domain_from_curves(og.create_rectangle(center=(5, 0)), 9, name="sq")
    og.create_extruded_block(square, "translate", vector=(0, 0, 1), dimension=5, name="blk")
    identifiers = [face.identifier for face in all_faces(og.store)]
    assert identifiers == [
        "vol.bottom", "vol.top", "vol.side-1", "vol.side-2", "vol.side-3", "vol.side-4",
        "blk.imin", "blk.imax", "blk.jmin", "blk.jmax", "blk.kmin", "blk.kmax",
    ]  # fmt: skip
    for identifier in identifiers:
        assert face_identifier(og.store, face_member(og.store, identifier)) == identifier
    og.set_boundary_type(["vol.side-1", "blk.imin"], "inlet")
    og.make_periodic_pair(["blk.jmin"], ["blk.jmax"])
    assert [face_identifier(og.store, m) for m in og.entity("inlet").params["members"]] == [
        "vol.side-1",
        "blk.imin",
    ]
    assert "og.set_boundary_type(['vol.side-1', 'blk.imin'], 'inlet'" in og.journal_text()
    replayed = opengrid.replay(og.journal_text())
    assert replayed.entity("inlet").params == og.entity("inlet").params
    report = og.boundary_report()
    assert "vol.side-1" in report.boundary and "blk.jmax" in report.boundary


def test_the_tree_lists_faces_with_their_types(live):
    window, og = live.window, live.og
    block = _box(og)
    block_id = og.entity(block).id
    og.set_boundary_type([f"{block}.{_face_at(og, block, 0, 0.0)}"], "inlet")
    tree = window.model_tree
    rows = dict(tree.face_rows(block_id))
    assert set(rows) == {"imin", "imax", "jmin", "jmax", "kmin", "kmax"}
    inlet = _face_at(og, block, 0, 0.0)
    assert rows[inlet] == "face · inlet" and rows[_face_at(og, block, 0, 2.0)] == "face · no type"
    assert "Patches" not in tree.category_names()
    # A face picked in the tree: a hidden one too (here the bottom, under the block).
    item = tree._face_items[(("block", block_id), ("face", "kmin"))]
    item.setSelected(True)
    assert window.selected_faces == [{"block": block_id, "face": "kmin"}]
    assert "kmin" in window.property_panel._title.text()


def test_a_curve_type_is_removed_after_extruding(live):
    """A curve typed in 2D and then extruded is the volume's side wall: its type is removed by
    the curve's name (as the panel's "(none)" does) or by the face identifier, never added
    twice."""
    window, og = live.window, live.og
    lines = og.create_rectangle(width=2, height=1)
    base = og.unstructured_domain_from_curves(lines, 9, name="b")
    curve = og.entity(base).params["loops"][0][0]["connector"]
    name = og.store.get(curve).name
    og.set_boundary_type([name], "wall")
    og.create_extruded_volume(base, "translate", vector=(0, 0, 0.5), dimension=3, name="vol")
    og.add_to_patch("wall", [name])  # already there (as vol.side-1)
    assert len(og.entity("wall").params["members"]) == 1
    window.set_selection([curve])
    panel = window.property_panel
    panel._editors["boundary"][1].setCurrentText("(none)")
    panel.apply()
    assert panel._error.text() == ""
    assert og.entity("wall").params["members"] == ()
    og.undo()
    og.remove_from_patch("wall", ["vol.side-1"])
    assert og.entity("wall").params["members"] == ()
