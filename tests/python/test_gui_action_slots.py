"""Menu actions reach their slots without Qt's "checked" flag: Mesh surfaces, Extrude 2D
model and File > New ran with their first optional parameter set to False (the user: Mesh
surfaces "does nothing")."""

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")

from PySide6.QtWidgets import QInputDialog, QMessageBox

from opengrid.app.main_window import MainWindow


def test_actions_triggered_from_the_menus(qtbot, monkeypatch):
    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    og.create_sphere(name="sp")
    window.set_selection([og.entity("sp-a").id, og.entity("sp-b").id])
    monkeypatch.setattr(QInputDialog, "getDouble", staticmethod(lambda *a, **k: (0.25, True)))
    window.mesh_surfaces_action.trigger()
    assert og.names("unstructured_domain") == ["sp-a-mesh", "sp-b-mesh"]
    # Extrude 2D model from its menu entry (both dialogs answered).
    og.new()
    window._refresh_all()
    og.create_rectangle((0, 0), 2, 1, name="r")
    for side in ("bottom", "right", "top", "left"):
        og.create_connector(f"r-{side}", 6, name=f"c-{side}")
    og.create_unstructured_domain(connectors=["c-bottom", "c-right", "c-top", "c-left"], name="d")
    monkeypatch.setattr(QInputDialog, "getItem", staticmethod(lambda *a, **k: ("periodic", True)))
    window.extrude_model_action.trigger()
    assert og.is_3d() and "front" in og.store
    # File > New asks before discarding changes.
    asked = []

    def question(*_args, **_kwargs):
        asked.append(True)
        return QMessageBox.StandardButton.Cancel

    monkeypatch.setattr(QMessageBox, "question", staticmethod(question))
    window.new_action.trigger()
    assert asked == [True] and "d" in og.store  # cancelled: nothing discarded
    window.close()
