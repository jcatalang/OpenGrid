import json

import pytest

import opengrid
from opengrid.project import from_json, load, save, to_json


def _rich_session() -> opengrid.Session:
    og = opengrid.Session()
    og.create_point((0.1, 0.2, 0.3))
    og.create_line((0, 0), (1, 0), name="chord")
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_polyline([(0, 0), (1e-7, 3.3), (2.0 / 3.0, 1.0)])
    og.create_spline([(0, 0), (1, 1), (2, 0)], method="clamped", end_tangent=(1, -1))
    og.create_nurbs(
        2, [0, 0, 0, 0.5, 1, 1, 1], [(0, 0), (1, 1), (2, 0), (3, 1)], weights=[1, 2, 1, 0.5]
    )
    og.delete("point-1")
    return og


def test_round_trip_is_byte_identical(tmp_path):
    og = _rich_session()
    first = tmp_path / "a.ogp"
    second = tmp_path / "b.ogp"
    og.save(first)
    other = opengrid.Session()
    other.load(first)
    other.save(second)
    assert first.read_bytes() == second.read_bytes()
    assert other.store.snapshot() == og.store.snapshot()
    assert other.store.next_id == og.store.next_id


def test_loaded_project_continues_ids_and_names(tmp_path):
    og = _rich_session()
    og.save(tmp_path / "p.ogp")
    other = opengrid.Session()
    other.load(tmp_path / "p.ogp")
    new = other.create_point((5, 5))
    assert new == "point-2"
    assert other.entity(new).id == og.store.next_id


def test_file_is_valid_sorted_json(tmp_path):
    og = _rich_session()
    text = to_json(og.store)
    data = json.loads(text)
    assert data["format"] == "opengrid-project"
    assert data["version"] == 1
    assert list(data) == sorted(data)
    assert [r["name"] for r in data["database"]] == og.names()


def test_load_after_session_clears_history(tmp_path):
    og = _rich_session()
    og.save(tmp_path / "p.ogp")
    og.create_point((1, 1))
    og.load(tmp_path / "p.ogp")
    assert not og.commands.can_undo
    assert "og.load(" in og.journal_text()
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.update(format="other"), "not an OpenGrid project"),
        (lambda d: d.update(version=99), "unsupported project version"),
        (lambda d: d["database"][0].update(type="blob"), "unknown type"),
        (lambda d: d["database"][1].update(id=d["database"][0]["id"]), "duplicate id"),
        (lambda d: d["database"][1].update(name=d["database"][0]["name"]), "duplicate name"),
        (lambda d: d["database"][0].update(p1=[0.0, 0.0, 0.0]), "coincide"),
        (lambda d: d.update(next_id=1), "next_id"),
        (lambda d: d.update(domains=[{"x": 1}]), "not supported"),
    ],
)
def test_invalid_projects_are_reported(mutate, message):
    og = _rich_session()
    data = json.loads(to_json(og.store))
    mutate(data)
    with pytest.raises(opengrid.ValidationError, match=message):
        from_json(json.dumps(data), source="bad.ogp")


def test_invalid_json_reports_line():
    with pytest.raises(opengrid.ValidationError, match=r"bad.ogp:2"):
        from_json('{\n "format": ,\n}', source="bad.ogp")


def test_module_level_save_load(tmp_path):
    og = _rich_session()
    save(og.store, tmp_path / "x.ogp")
    store = load(tmp_path / "x.ogp")
    assert store.snapshot() == og.store.snapshot()
