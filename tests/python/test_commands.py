"""Undo/redo of every command, rollback on failure, and journal replay."""

import math

import pytest
from helpers import naca0012_selig, write_selig

import opengrid
from opengrid.project import to_json


def _base_session() -> opengrid.Session:
    og = opengrid.Session()
    og.create_point((0.0, 0.0), name="origin")
    og.create_line((0.0, 0.0), (1.0, 0.0), name="chord")
    og.create_spline([(0, 0), (1, 1), (2, 0)], name="bump")
    return og


def _airfoil_file(tmp_path):
    path = tmp_path / "naca0012.dat"
    write_selig(path, naca0012_selig())
    return path


def _csv_file(tmp_path):
    path = tmp_path / "points.csv"
    path.write_text("x,y,z\n0,0,0\n1,0.5,0\n2,0,0.1\n")
    return path


COMMANDS = {
    "create_point": lambda og, tmp: og.create_point((3, 4, 5)),
    "create_line": lambda og, tmp: og.create_line((0, 1), (2, 3)),
    "create_arc": lambda og, tmp: og.create_arc((0, 0), 2.0, 0.0, 1.5),
    "create_arc_3pt": lambda og, tmp: og.create_arc_3pt((1, 0), (0, 1), (-1, 0)),
    "create_circle": lambda og, tmp: og.create_circle((0.5, 0), 20.0),
    "create_polyline": lambda og, tmp: og.create_polyline([(0, 0), (1, 0), (1, 1)]),
    "create_spline": lambda og, tmp: og.create_spline(
        [(0, 0), (1, 2), (3, 1)], method="catmull_rom"
    ),
    "create_spline_clamped": lambda og, tmp: og.create_spline(
        [(0, 0), (1, 2), (3, 1)], method="clamped", start_tangent=(0, 1)
    ),
    "create_nurbs": lambda og, tmp: og.create_nurbs(
        2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)]
    ),
    "modify": lambda og, tmp: og.modify("chord", p1=(2.0, 0.5)),
    "rename": lambda og, tmp: og.rename("bump", "hill"),
    "delete": lambda og, tmp: og.delete(["origin", "bump"]),
    "import_airfoil_spline": lambda og, tmp: og.import_airfoil_dat(_airfoil_file(tmp)),
    "import_airfoil_split": lambda og, tmp: og.import_airfoil_dat(_airfoil_file(tmp), mode="split"),
    "import_csv": lambda og, tmp: og.import_csv(_csv_file(tmp), mode="polyline"),
}


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_undo_redo_every_command(command, tmp_path):
    og = _base_session()
    before = og.store.snapshot()
    COMMANDS[command](og, tmp_path)
    after = og.store.snapshot()
    assert after != before
    og.undo()
    assert og.store.snapshot() == before
    og.redo()
    assert og.store.snapshot() == after
    # A second full cycle must behave identically.
    og.undo()
    assert og.store.snapshot() == before
    og.redo()
    assert og.store.snapshot() == after
    # Geometry is rebuilt correctly after redo.
    for entity in og.store:
        if entity.kind != "point":
            assert og.length(entity.name) > 0.0


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_journal_replay_reproduces_project(command, tmp_path):
    og = _base_session()
    COMMANDS[command](og, tmp_path)
    og.undo()
    og.redo()
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)


def test_failed_command_rolls_back():
    og = _base_session()
    before = og.store.snapshot()
    journal = og.journal_text()
    with pytest.raises(opengrid.ValidationError):
        og.modify("bump", points=[(0, 0)])
    with pytest.raises(opengrid.ValidationError, match="already in use"):
        og.rename("bump", "chord")
    with pytest.raises(opengrid.ValidationError, match="no such entity"):
        og.delete(["origin", "missing"])
    with pytest.raises(KeyError):
        og.modify("missing", p0=(0, 0))
    assert og.store.snapshot() == before
    assert og.journal_text() == journal
    og.undo()  # undoes the last *successful* command (create bump)
    assert "bump" not in og.store


def test_new_command_clears_redo():
    og = _base_session()
    og.undo()
    assert og.commands.can_redo
    og.create_point((1, 1))
    assert not og.commands.can_redo
    with pytest.raises(IndexError):
        og.redo()


def test_undo_restores_ids_and_names():
    og = _base_session()
    chord_id = og.entity("chord").id
    og.delete("chord")
    og.undo()
    assert og.entity("chord").id == chord_id
    # ids are never reused
    new = og.create_point((9, 9))
    assert og.entity(new).id > max(e.id for e in og.store if e.name != new)


def test_auto_names_are_unique_and_deterministic():
    og = opengrid.Session()
    names = [og.create_point((i, 0)) for i in range(3)]
    assert names == ["point-1", "point-2", "point-3"]
    og.delete("point-2")
    assert og.create_point((5, 5)) == "point-4"
    og.create_point((6, 6), name="point-9")
    assert og.create_point((7, 7)) == "point-10"


def test_invalid_names_rejected():
    og = opengrid.Session()
    for bad in ["", " a", "has space", "-lead", "x/y"]:
        with pytest.raises(opengrid.ValidationError):
            og.create_point((0, 0), name=bad)


def test_store_requires_transaction():
    og = opengrid.Session()
    with pytest.raises(RuntimeError, match="transaction"):
        og.store.add("point", {"xyz": (0.0, 0.0, 0.0)})


def test_listeners_notified():
    og = opengrid.Session()
    changes = []
    og.commands.add_listener(changes.append)
    og.create_point((0, 0))
    og.undo()
    og.redo()
    assert len(changes) == 3
    for change in changes:
        assert change.ids() == [1]
        assert change.before[1] is None
        assert change.after[1].name == "point-1"


def test_circle_canonical_params():
    og = opengrid.Session()
    name = og.create_circle((0, 0), 1.0)
    assert og.entity(name).params["end_angle"] == math.tau
