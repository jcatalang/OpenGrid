"""Messages the user sees: warnings and errors in the bar above the view (app/message_bar.py),
a plain Log tab (app/log_panel.py), and messages that say what happened (layers that did not
all grow and why, a structured domain's curve given layers)."""

import logging

import numpy as np
import pytest

import opengrid

LAYERS = {"first_height": 0.01, "growth_rate": 1.15, "layers": 20, "transition": True}


def sketch_box(og):
    """The user's case: a sketch block 0.4 x 0.2 at spacing 0.015."""
    og.create_sketch(
        vertices=((-0.2, 0.2, 0.0), (-0.2, 0.0, 0.0), (0.2, 0.0, 0.0), (0.2, 0.2, 0.0)),
        blocks=((0, 1, 2, 3),),
        spacing=0.015,
        name="sketch-1",
    )


@pytest.fixture
def window(qtbot):
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    yield window
    window.close()


def test_a_refusal_is_shown_above_the_view(window):
    og = window.session
    sketch_box(og)
    window.set_selection([og.entity("sketch-1-c1-2").id])
    assert window.inflation_on_selection(LAYERS) is None
    bar = window.message_bar
    assert bar.isVisibleTo(window) and bar._title.text() == "Not done"
    assert "a structured domain" in bar.lines()[0] and "Convert domain" in bar.lines()[0]
    text = window.log_panel.toPlainText()
    assert "not done: sketch-1-c1-2: bounds sketch-1-d0, a structured domain" in text
    assert "WARNING" not in text and "opengrid:" not in text  # no levels, no logger names
    bar.log_requested.emit()
    assert window.console_tabs.currentWidget() is window.log_panel


def test_layers_that_do_not_all_grow_say_why(window):
    og = window.session
    sketch_box(og)
    og.edit_sketch("sketch-1", unstructured=(0,))
    window.set_selection([og.entity("sketch-1-c1-2").id])
    window.message_bar.hide()
    assert window.inflation_on_selection(LAYERS) == ["sketch-1-d0"]
    mesh = og.unstructured_mesh("sketch-1-d0")
    assert mesh["layers"] == 4 and mesh["layer_stops"]["transition"] > 0
    (line,) = window.message_bar.lines()
    assert line.startswith("sketch-1-d0: 4 of the 20 inflation layers grow")
    assert "as thick as the cells beside it" in line and "about 0.015" in line
    assert window.message_bar._title.text() == "Warning"
    # A small first height: every layer grows, nothing to say about them (the sketch's edge
    # leaving the wall, clustered at 1e-4, is now stretched: that is shown).
    window.message_bar.hide()
    window.inflation_on_selection({**LAYERS, "first_height": 1e-4, "layers": 8})
    assert og.unstructured_mesh("sketch-1-d0")["layers"] == 8
    assert not any("inflation layers" in line for line in window.message_bar.lines())
    assert any("max stretching ratio" in line for line in window.message_bar.lines())


def test_the_log_counts_repeats_and_keeps_details_out(window):
    log = logging.getLogger("opengrid")
    log.warning("w: the same thing")
    log.warning("w: the same thing")
    log.warning("w: the same thing")
    log.debug("a detail for the log file only")
    text = window.log_panel.toPlainText()
    assert "warning: w: the same thing (3 times)" in text
    assert text.count("the same thing") == 1
    assert "a detail" not in text
    assert "Compute:" not in text  # (the compute summary is a detail)
    # The panel's "No changes" stays in the status bar.
    window.property_panel.message.emit("No changes", "status")
    assert "No changes" not in window.log_panel.toPlainText()
    assert window.statusBar().currentMessage() == "No changes"


def test_messages_shown_together_and_by_a_dialog():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from opengrid.app.message_bar import MAX_LINES, MessageBar

    QApplication.instance() or QApplication([])
    bar = MessageBar()
    for k in range(MAX_LINES + 3):
        bar.post(logging.WARNING, f"line {k}")
    bar.post(logging.ERROR, "line 0")  # (a repeat is not listed twice)
    assert len(bar.lines()) == MAX_LINES + 3
    assert "and 3 more" in bar._text.text() and bar._title.text() == "Not done"
    record = logging.LogRecord("opengrid", logging.ERROR, "", 0, "shown", (), None)
    record.on_screen = False
    bar.hide()
    bar.handler.emit(record)
    assert not bar.isVisible()


def test_sliding_layers_in_a_sketch_and_their_corner():
    """The user's case: walls on two sides of a sketch's unstructured block meeting in a
    concave corner, layers sliding into the other two sides."""
    og = opengrid.Session()
    sketch_box(og)
    og.edit_sketch("sketch-1", unstructured=(0,))
    for curve, growth, layers in (
        ("sketch-1-c1-2", 1.05, 20),
        ("sketch-1-c0-1", 1.15, 20),
        ("sketch-1-c1-2", 1.15, 10),
    ):
        og.set_inflation([curve], first_height=1e-3, growth_rate=growth, layers=layers, slide=True)
    inflation = og.entity("sketch-1-d0").params["generator"]["inflation"]
    assert inflation["slide"] and inflation["transition"]
    # The far sides got points at the layer heights (an explicit stack, then tanh).
    for side in ("sketch-1-c2-3", "sketch-1-c0-3"):
        assert og.entity(side).params["distribution"]["type"] == "explicit"
    mesh = og.unstructured_mesh("sketch-1-d0")
    assert mesh["layers"] == 10
    assert mesh["layer_stops"]["proximity"] == 0 and mesh["layer_stops"]["slope"] == 0
    quality = og.domain_quality("sketch-1-d0")
    assert np.min(quality["min_angle"]) > 20.0  # (3.4 degrees before: steps at the corner)
    replayed = opengrid.replay(og.journal_text())
    assert len(replayed.unstructured_mesh("sketch-1-d0")["quads"]) == len(mesh["quads"])


def test_wrapped_panel_text_is_not_cut(window, qtbot):
    """Every form of the panels lets its fields take the width available (macOS's default kept
    them at their preferred width, and word-wrapped labels lost their second line)."""
    from PySide6.QtWidgets import QFormLayout, QLabel

    from opengrid.app.inflation_dialog import InflationDialog

    og = window.session
    sketch_box(og)
    og.edit_sketch("sketch-1", unstructured=(0,))
    og.set_inflation(["sketch-1-c1-2"], first_height=0.01, growth_rate=1.15, layers=20)
    window.show()

    def cut():
        return [
            label.text()
            for label in window.property_panel.findChildren(QLabel)
            if label.isVisible()
            and label.wordWrap()
            and label.text()
            and label.height() + 1 < label.heightForWidth(label.width())
        ]

    for name in ("sketch-1-c1-2", "sketch-1-d0", "sketch-1"):
        window.set_selection([og.entity(name).id])
        forms = window.property_panel.findChildren(QFormLayout)
        assert forms
        for form in forms:
            assert form.fieldGrowthPolicy() == QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        qtbot.waitUntil(lambda: not cut(), timeout=3000)  # (once the layout has settled)
    dialog = InflationDialog(["c"])
    (form,) = dialog.findChildren(QFormLayout)
    assert form.fieldGrowthPolicy() == QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow


def test_sliding_layers_into_an_edge_between_two_walls():
    """The user's trapezoid: walls on three sides of a sketch block, the fourth (an outlet)
    between two of them: the layers of both walls run into it, at its two ends."""
    og = opengrid.Session()
    og.create_sketch(
        vertices=((-0.2, 0.15, 0.0), (-0.2, 0.0, 0.0), (0.2, 0.0, 0.0), (0.2, 0.2, 0.0)),
        blocks=((0, 1, 2, 3),),
        spacing=0.0144,
        name="sketch-1",
    )
    og.edit_sketch("sketch-1", unstructured=(0,))
    for curve in ("sketch-1-c1-2", "sketch-1-c0-1", "sketch-1-c0-3"):
        og.set_inflation([curve], first_height=1e-3, growth_rate=1.15, layers=10, slide=True)
    outlet = og.entity("sketch-1-c2-3").params
    assert outlet["distribution"]["type"] == "explicit" and len(outlet["breakpoints"]) == 2
    mesh = og.unstructured_mesh("sketch-1-d0")
    stops = mesh["layer_stops"]
    assert stops["slope"] == stops["proximity"] == stops["transition"] == 0
    assert np.min(og.domain_quality("sketch-1-d0")["min_angle"]) > 15.0  # (4 degrees before)
