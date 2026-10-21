"""Session.convert_domain: a domain changes kind on the same connectors, keeping its name."""

import pytest

import opengrid
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError


def rectangle(og, n=11, m=6):
    for a, b in (((0, 0), (2, 0)), ((2, 0), (2, 1)), ((2, 1), (0, 1)), ((0, 1), (0, 0))):
        og.create_line(a, b)
    for name, count in zip(og.names("line"), (n, m, n, m), strict=True):
        og.create_connector(name, count)
    return og.names("connector")


def test_convert_structured_to_unstructured_and_back(tmp_path):
    og = opengrid.Session()
    c = rectangle(og)
    og.create_domain([[c[0]], [c[1]], [c[2]], [c[3]]], name="d")
    og.create_zone("fluid", ["d"])
    og.create_patch("all", "wall", c)
    assert og.convert_domain("d") == "d"
    assert og.entity("d").kind == "unstructured_domain"
    mesh = og.unstructured_mesh("d")
    assert len(mesh["triangles"]) > 0 and len(mesh["points"]) >= 2 * 11 + 2 * 4
    assert og.store.get(og.entity("fluid").params["members"][0]["domain"]).name == "d"
    assert og.boundary_report().ok
    assert og.convert_domain("d") == "d"
    assert og.entity("d").kind == "domain"
    assert og.domain_dimensions("d") == (11, 6)
    assert og.store.get(og.entity("fluid").params["members"][0]["domain"]).name == "d"
    og.undo()
    assert og.entity("d").kind == "unstructured_domain"
    og.redo()
    assert og.entity("d").kind == "domain"
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
    og.save(tmp_path / "convert.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "convert.ogp")
    assert to_json(other.store) == to_json(og.store)


def test_convert_domain_refusals():
    og = opengrid.Session()
    c = rectangle(og)
    og.create_domain([[c[0]], [c[1]], [c[2]], [c[3]]], name="d")
    og.create_smoothing_group(["d"], name="g")
    with pytest.raises(ValidationError, match="smoothing group g"):
        og.convert_domain("d")
    og.delete("g")
    og.create_extruded_block("d", "translate", vector=(0, 0, 1), dimension=3, name="b")
    with pytest.raises(ValidationError, match="block b"):
        og.convert_domain("d")
    og.delete("b")
    og.create_circle((1, 0.5), 0.2, name="hole")
    og.create_connector("hole", 21, name="ch")
    og.convert_domain("d")
    og.modify(
        "d", loops=[[*og.entity("d").params["loops"][0]], [{"connector": "ch", "reversed": True}]]
    )
    with pytest.raises(ValidationError, match="has holes"):
        og.convert_domain("d")
    with pytest.raises(ValidationError, match="not a domain"):
        og.convert_domain(c[0])


def test_convert_a_sketch_domain_through_its_sketch():
    """A sketch's domain converts by its block joining (or leaving) the sketch's unstructured
    blocks (it was refused, 2026-10-21: the user's conversion "did not work")."""
    og = opengrid.Session()
    og.create_sketch([(3, 0), (4, 0), (4, 1), (3, 1)], [(0, 1, 2, 3)], spacing=0.25, name="sk")
    assert og.names("smoothing_group") == ["sk-group"]
    assert og.convert_domain("sk-d0") == "sk-d0"
    assert og.entity("sk-d0").kind == "unstructured_domain"
    assert og.entity("sk").params["unstructured"] == (0,)
    assert og.names("smoothing_group") == []  # unstructured blocks are not smoothed
    assert og.convert_domain("sk-d0") == "sk-d0"
    assert og.entity("sk-d0").kind == "domain" and og.entity("sk").params["unstructured"] == ()
    og.undo()
    assert og.entity("sk-d0").kind == "unstructured_domain"
    assert opengrid.replay(og.journal_text()).entity("sk-d0").kind == "unstructured_domain"


def test_convert_a_smoothed_domain_leaving_its_group():
    og = opengrid.Session()
    c = rectangle(og)
    og.create_domain([[c[0]], [c[1]], [c[2]], [c[3]]], name="d")
    og.create_smoothing_group(["d"], name="g")
    assert og.convert_domain("d", leave_groups=True) == "d"
    assert og.entity("d").kind == "unstructured_domain" and "g" not in og.store  # left empty
    og.undo()
    assert og.entity("d").kind == "domain" and "g" in og.store


def test_convert_in_the_window(qtbot, monkeypatch):
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    og.create_sketch([(0, 0), (1, 0), (1, 1), (0, 1)], [(0, 1, 2, 3)], spacing=0.25, name="sk")
    window.set_selection([og.entity("sk-d0").id])
    assert window.convert_selected_domain() == "sk-d0"
    assert og.entity("sk-d0").kind == "unstructured_domain"
    # A domain smoothed in a group: asked first; Cancel keeps it.
    lines = og.create_rectangle(center=(3, 0.5), width=2, height=1)
    og.structured_domain_from_curves(lines, 9, name="d")
    og.create_smoothing_group(["d"], name="g")
    window.set_selection([og.entity("d").id])
    asked = []
    monkeypatch.setattr(window, "_confirm_leave_groups", lambda d, g: asked.append((d, g)) or False)
    assert window.convert_selected_domain() is None and og.entity("d").kind == "domain"
    assert asked == [("d", ["g"])]
    monkeypatch.setattr(window, "_confirm_leave_groups", lambda d, g: True)
    assert window.convert_selected_domain() == "d"
    assert og.entity("d").kind == "unstructured_domain"
    window.close()
