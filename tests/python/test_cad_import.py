"""STEP / IGES import (Phase 9a): opengrid._occt and Session.import_cad.

Fixtures in tests/data: block_hole.step / .igs (box [0, 2] x [0, 2] x [0, 1] minus a cylinder of
radius 0.4 about (1, 1, z), written by OCCT; the IGES file holds unsewn trimmed faces) and
sphere.step (a sphere: one face with two poles).
"""

from pathlib import Path

import numpy as np
import pytest

import opengrid
from opengrid.project import to_json
from opengrid.validation import ValidationError

pytest.importorskip("opengrid._occt")

DATA = Path(__file__).resolve().parents[1] / "data"
STEP = DATA / "block_hole.step"
IGES = DATA / "block_hole.igs"
SPHERE = DATA / "sphere.step"


def on_part(points) -> np.ndarray:
    """Distance of points to the part's surfaces (box planes and the hole's cylinder)."""
    p = np.atleast_2d(points)
    radial = np.abs(np.hypot(p[:, 0] - 1.0, p[:, 1] - 1.0) - 0.4)
    planes = np.abs(np.stack([p[:, 0], 2 - p[:, 0], p[:, 1], 2 - p[:, 1], p[:, 2], 1 - p[:, 2]]))
    return np.minimum(radial, planes.min(axis=0))


def test_read_cad_model():
    from opengrid import _occt

    model = _occt.read_cad(str(STEP))
    assert model["solids"] == 1
    assert len(model["faces"]) == 8  # the hole's periodic face is split in two
    assert all(edge["faces"] == 2 for edge in model["edges"])
    for edge in model["edges"]:
        poles = np.asarray(edge["curve"]["poles"])
        assert np.array_equal(poles[0], model["vertices"][edge["v0"]])
        assert np.array_equal(poles[-1], model["vertices"][edge["v1"]])
    with pytest.raises(ValueError):
        _occt.read_cad(str(STEP), tolerance=-1.0)


def test_import_creates_entities_on_the_cad_geometry():
    og = opengrid.Session()
    result = og.import_cad(STEP, spacing=0.1)
    assert len(result.curves) == len(result.connectors) == 18
    assert len(result.domains) == 8
    assert len(result.surfaces) == 2  # the two hole halves; planar faces need none
    assert not result.skipped_faces and not result.free_edges and not result.non_manifold_edges
    assert result.solids == 1
    for name in result.connectors:
        nodes = og.connector_nodes(name)
        assert on_part(nodes.points).max() < 1e-9
        spacing = np.linalg.norm(np.diff(nodes.points, axis=0), axis=1)
        assert spacing.max() <= 0.1 + 1e-9
    for name in result.domains:
        mesh = og.unstructured_mesh(name)
        assert on_part(mesh["points"]).max() < 1e-9
        assert len(mesh["triangles"]) > 0
    # One undo step removes the whole import; the journal replays it without the file.
    replayed = opengrid.replay(og.journal_text())
    assert to_json(replayed.store) == to_json(og.store)
    assert "imported from block_hole.step" in og.journal_text()
    og.undo()
    assert og.names() == []


def test_turning_angle_refines_curved_edges():
    og = opengrid.Session()
    result = og.import_cad(STEP, spacing=10.0, max_turn_angle=10.0)
    for name in result.connectors:
        nodes = og.connector_nodes(name).points
        tangent = np.diff(nodes, axis=0)
        tangent /= np.linalg.norm(tangent, axis=1)[:, None]
        if len(tangent) > 1:
            turn = np.degrees(np.arccos(np.clip(np.sum(tangent[1:] * tangent[:-1], 1), -1, 1)))
            assert turn.max() <= 10.0 + 1e-6
    half_circles = [n for n in result.connectors if len(og.connector_nodes(n).points) > 2]
    assert len(half_circles) == 4  # the hole's rims, split in two halves each
    assert all(len(og.connector_nodes(n).points) == 19 for n in half_circles)


def test_auto_domains_are_structured_where_four_sided():
    og = opengrid.Session()
    result = og.import_cad(STEP, domains="auto", spacing=0.15)
    # 4 sides of the box + 2 hole halves are four-sided; top and bottom have a hole.
    assert len(result.structured) == 6
    assert len(result.domains) == 8
    for name in result.structured:
        ni, nj = og.domain_dimensions(name)
        grid = og.domain_grid(name)
        assert grid.shape == (nj, ni, 3)
        assert on_part(grid.reshape(-1, 3)).max() < 1e-9
        assert "error" not in og.domain_status(name).lower()


def test_none_mode_and_option_checks():
    og = opengrid.Session()
    result = og.import_cad(STEP, domains="none", name="p")
    assert result.domains == [] and result.surfaces == []
    assert og.names("connector")[0] == "p-con-1"
    with pytest.raises(ValidationError):
        og.import_cad(STEP, domains="hex")
    with pytest.raises(ValidationError):
        og.import_cad(STEP, spacing=0.0)
    with pytest.raises(ValidationError):
        og.import_cad(DATA / "missing.step")
    with pytest.raises(ValidationError):
        og.import_cad(Path(__file__))  # unknown extension


def test_iges_faces_are_sewn_into_a_closed_shell():
    og = opengrid.Session()
    result = og.import_cad(IGES, spacing=0.15)
    assert any("sewing" in line for line in result.log)
    assert not result.free_edges
    assert len(result.domains) == 8
    loose = opengrid.Session()
    unsewn = loose.import_cad(IGES, spacing=0.15, sew=False)
    assert unsewn.free_edges
    assert "free edge" in unsewn.report()


def test_faces_with_poles_are_reported_not_meshed():
    og = opengrid.Session()
    result = og.import_cad(SPHERE)
    assert result.domains == []
    assert result.skipped_faces
    assert all("degenerate edge" in reason for _face, reason in result.skipped_faces)
    assert "no domain" in result.report()


def test_imported_solid_meshes_as_a_volume():
    pytest.importorskip("netgen.meshing")
    og = opengrid.Session()
    result = og.import_cad(STEP, spacing=0.15)
    og.create_volume(result.domains, name="solid")
    mesh = og.volume_mesh("solid")
    from opengrid.model.volumes import total_volume

    exact = 2.0 * 2.0 * 1.0 - np.pi * 0.4**2
    assert abs(total_volume(mesh) - exact) < 0.01 * exact
    assert mesh["quality"]["below_floor"] == 0


def test_gui_import_cad(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    result = window.import_cad(str(STEP), domains="auto")
    assert result is not None and len(result.structured) == 6
    assert window.import_cad(str(DATA / "missing.step"), domains="none") is None
    window.close()
