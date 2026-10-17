"""Handles come with what moves the selection: a domain offers its curves' handles, anything
of a sketch the sketch's vertices, so the geometry can be dragged with the mesh in place and
remeshed (the user could not move the corners of a meshed sketch block)."""

import pytest

import opengrid
from opengrid.app.handles import handles_for


def test_handles_follow_the_selection_to_what_moves_it():
    og = opengrid.Session()
    og.create_sketch(
        [(0, 0), (2, 0), (1.8, 1), (0.2, 1)],
        [(0, 1, 2, 3)],
        spacing=0.1,
        unstructured=[0],
        walls={(0, 1): 0.01},
        name="s",
    )
    sketch_handles = [(h.entity, h.field, h.index) for h in handles_for(og, ["s"])]
    assert sketch_handles == [("s", "vertices", k) for k in range(4)]
    # The sketch's domain, connector and line all offer the sketch's vertices (nothing else).
    for name in ("s-d0", "s-c0-1", "s-e0-1"):
        assert [(h.entity, h.field, h.index) for h in handles_for(og, [name])] == sketch_handles
    assert len(handles_for(og, ["s-d0", "s-c0-1", "s"])) == 4  # once
    # A plain domain offers the handles of the curves under its connectors.
    for name, (a, b) in (
        ("b", ((0, 0), (1, 0))),
        ("r", ((1, 0), (1, 1))),
        ("t", ((1, 1), (0, 1))),
        ("l", ((0, 1), (0, 0))),
    ):
        og.create_line(a, b, name=name)
        og.create_connector(name, 7, name=f"c-{name}")
    og.create_domain([["c-b"], ["c-r"], ["c-t"], ["c-l"]], name="d")
    og.create_unstructured_domain(connectors=["c-b", "c-r", "c-t", "c-l"], name="u")
    assert [h.entity for h in handles_for(og, ["d"])] == ["b", "b", "r", "r", "t", "t", "l", "l"]
    assert [h.entity for h in handles_for(og, ["u"])] == ["b", "b", "r", "r", "t", "t", "l", "l"]
    assert [h.entity for h in handles_for(og, ["c-b"])] == ["b", "b"]


def test_dragging_a_corner_of_a_meshed_sketch_block(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    og.create_sketch(
        [(0, 0), (2, 0), (1.8, 1), (0.2, 1)],
        [(0, 1, 2, 3)],
        spacing=0.05,
        unstructured=[0],
        walls={(0, 1): 0.005},
        name="s",
    )
    before = og.unstructured_mesh("s-d0")
    assert before["layers"] > 0
    window.scene.world_per_pixel = lambda: 0.002
    window.set_selection([og.entity("s-d0").id])  # the domain, as a click in the mesh gives
    assert [h.entity for h in window.visible_handles()] == ["s"] * 4
    assert window.handle_press((1.8, 1.0, 0.0))
    for k in range(1, 6):
        window.handle_drag((1.8 + 0.04 * k, 1.0 + 0.02 * k, 0.0))
    assert window.handle_release((2.0, 1.1, 0.0)) == ["s"]
    assert og.entity("s").params["vertices"][2] == (2.0, 1.1, 0.0)
    after = og.unstructured_mesh("s-d0")
    assert after["layers"] > 0 and after is not before  # remeshed with the wall layers
    assert og.journal_text().splitlines()[-1].startswith("og.move_handle('s', 'vertices', 2")
    og.undo()
    assert og.entity("s").params["vertices"][2] == (1.8, 1.0, 0.0)
    window.close()
