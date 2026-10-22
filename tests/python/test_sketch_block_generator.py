"""Session.set_unstructured_generator on a sketch's block: the triangulation settings go to
the sketch's generator, the inflation to walls on the block's edges; the property panel of
such a domain works (the user could not change the smoothing, the algorithm or the element
type of a sketch block)."""

import pytest

import opengrid
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def sketch_with_unstructured_block(og):
    og.create_sketch(
        [(0, 0), (1, 0), (1, 1), (0, 1), (2, 0), (2, 1)],
        [(0, 1, 2, 3), (1, 4, 5, 2)],
        spacing=0.1,
        unstructured=[1],
        walls={(0, 1): 0.01},  # the structured block's wall spacing
        name="s",
    )
    return og


def test_settings_go_to_the_sketch():
    og = sketch_with_unstructured_block(opengrid.Session())
    assert og.entity("s-d1").kind == "unstructured_domain"
    og.set_unstructured_generator(
        "s-d1",
        min_angle=30.0,
        smoothing="laplacian",
        inflation={
            "connectors": ["s-c1-4"],
            "first_height": 0.005,
            "layers": 8,
            "growth_rate": 1.3,
        },
    )
    sketch = og.entity("s").params
    assert (
        sketch["generator"]["min_angle"] == 30.0 and sketch["generator"]["smoothing"] == "laplacian"
    )
    assert sketch["generator"]["inflation"] is None
    walls = {tuple(w["edge"]): w for w in sketch["walls"]}
    assert set(walls) == {(0, 1), (1, 4)}  # the structured block's wall is untouched
    assert walls[(1, 4)]["spacing"] == 0.005 and walls[(1, 4)]["layers"] == 8
    assert walls[(1, 4)]["growth"] == 1.3
    generator = og.entity("s-d1").params["generator"]
    assert generator["min_angle"] == 30.0
    assert generator["inflation"]["connectors"] == (og.entity("s-c1-4").id,)
    assert generator["inflation"]["first_height"] == 0.005
    assert og.unstructured_mesh("s-d1")["layers"] > 0
    assert "inflation" in og.journal_text().splitlines()[-1]
    assert "'s-c1-4'" in og.journal_text().splitlines()[-1]  # names, not ids
    # Dropping the inflation drops the wall; the shared edge's wall would be kept.
    og.set_unstructured_generator("s-d1", inflation=None)
    assert set(tuple(w["edge"]) for w in og.entity("s").params["walls"]) == {(0, 1)}
    assert og.entity("s-d1").params["generator"]["inflation"] is None
    og.undo()
    assert og.entity("s-d1").params["generator"]["inflation"] is not None
    # Refusals: a connector of another block, a total thickness.
    with pytest.raises(ValidationError, match=r"not an edge of block 1 of sketch 's'"):
        og.set_unstructured_generator(
            "s-d1", inflation={"connectors": ["s-c0-3"], "first_height": 0.005}
        )
    with pytest.raises(ValidationError, match="total thickness is not available"):
        og.set_unstructured_generator(
            "s-d1",
            inflation={"connectors": ["s-c1-4"], "first_height": 0.005, "total_thickness": 0.1},
        )
    # Layers past the core size and sliding layers go to the sketch's wall too.
    og.set_unstructured_generator(
        "s-d1",
        inflation={
            "connectors": ["s-c1-4"],
            "first_height": 0.005,
            "transition": False,
            "slide": True,
        },
    )
    wall = next(w for w in og.entity("s").params["walls"] if w["edge"] == (1, 4))
    assert not wall["transition"] and wall["slide"]
    inflation = og.entity("s-d1").params["generator"]["inflation"]
    assert not inflation["transition"] and inflation["slide"]
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_plain_domain_and_the_panel(qtbot):
    og = opengrid.Session()
    for name, (a, b) in (
        ("b", ((0, 0), (1, 0))),
        ("r", ((1, 0), (1, 1))),
        ("t", ((1, 1), (0, 1))),
        ("l", ((0, 1), (0, 0))),
    ):
        og.create_line(a, b, name=name)
        og.create_connector(name, 11, name=f"c-{name}")
    og.create_unstructured_domain(connectors=["c-b", "c-r", "c-t", "c-l"], name="d")
    og.set_unstructured_generator("d", {"min_angle": 28.0}, smoothing="laplacian")
    assert og.entity("d").params["generator"]["min_angle"] == 28.0
    assert og.entity("d").params["generator"]["smoothing"] == "laplacian"
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = sketch_with_unstructured_block(window.session)
    window.set_selection([og.entity("s-d1").id])
    panel = window.property_panel
    _kind, editor = panel._editors["generator"]
    editor.smoothing.setCurrentText("laplacian")
    editor.element_type.setCurrentIndex(1)
    chosen = editor.element_type.currentText()
    suggested = editor.suggested_first_height()
    editor.walls.setText("s-c1-4")  # first height left empty: the suggested one
    assert panel.apply(), panel._error.text()  # the panel is rebuilt afterwards
    sketch = og.entity("s").params
    assert sketch["generator"]["smoothing"] == "laplacian"
    assert sketch["generator"]["element_type"] == chosen
    walls = {tuple(w["edge"]): w for w in sketch["walls"]}
    assert walls[(1, 4)]["spacing"] == pytest.approx(suggested)
    assert og.entity("s-d1").params["generator"]["inflation"]["connectors"] == (
        og.entity("s-c1-4").id,
    )
    assert panel._error.text() == ""
    window.close()
