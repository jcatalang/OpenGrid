"""Selecting what sits inside a box: a click on a curve seen through a surface picks the
curve, and View > Hide selected takes the box out of the way (the user could not select a
wing inside its far-field box)."""

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from PySide6.QtWidgets import QApplication
from test_gui_events import Live


def wing_in_a_box(og):
    og.create_box(size=(20, 20, 20), name="box")
    og.create_wing("0012", root_chord=2.0, tip_chord=1.0, span=5.0, name="wing")


def test_curves_are_picked_through_surfaces(qtbot):
    live = Live(qtbot)
    try:
        og = live.og
        wing_in_a_box(og)
        window = live.window
        assert window.navigation_mode == "3d"
        window.scene.view_preset("iso", negative=False)
        live.widget.render()
        QApplication.processEvents()
        # A point of the root's upper spline, seen through the box's front faces.
        curve = og.geometry("wing-root-upper")
        xyz = tuple(curve.eval(0.4))
        point = live.display(xyz)
        x = point.x()
        y = live.widget.height() - 1 - point.y()
        entity_id, _picked = window.scene.pick_display(x, y)
        assert entity_id is not None and og.store.get(entity_id).name == "wing-root-upper"
        # A click on the same spot selects it.
        live.click(xyz)
        assert window.selected_names() == ["wing-root-upper"]
        # Off the wing, the box face is what is there.
        live.click((-9.0, -9.0, 9.0))
        assert window.selected_names() and window.selected_names()[0].startswith("box-")
    finally:
        live.close()


def test_hide_selected_and_show_all(qtbot):
    live = Live(qtbot)
    try:
        og = live.og
        wing_in_a_box(og)
        window = live.window
        faces = [og.entity(f"box-{s}").id for s in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")]
        window.hide_selected()  # nothing selected: a hint, nothing hidden
        assert window.scene.hidden == set()
        window.set_selection(faces)
        window.hide_action.trigger()
        assert window.scene.hidden == set(faces) and window.selected_names() == []
        assert all(f not in window.scene._entity_actors for f in faces)
        assert og.entity("wing-upper").id in window.scene._entity_actors
        # Hidden entities are not picked: the click lands on nothing.
        live.click((-9.0, -9.0, 9.0))
        assert window.selected_names() == []
        # The model may change meanwhile: hidden stays hidden.
        og.create_sphere(name="ball")
        assert all(f not in window.scene._entity_actors for f in faces)
        window.show_all_action.trigger()
        assert window.scene.hidden == set()
        assert all(f in window.scene._entity_actors for f in faces)
        assert "H" in [a.shortcut().toString() for a in window.actions_with_single_keys()]
    finally:
        live.close()
