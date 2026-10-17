"""Session.mesh_surfaces: domains on the faces of a box, a cylinder (two halves and the
discs), a wing (ruled surfaces and flat caps) and an extruded surface, with shared edge
curves and connectors, so the shell is watertight for a volume (the user wanted to mesh the
cylinder inside its box)."""

import math

import pytest

import opengrid
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError

BOX_FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")


def test_box_faces_share_their_edges():
    og = opengrid.Session()
    og.create_box(size=(2, 2, 2), name="box")
    domains = og.mesh_surfaces([f"box-{s}" for s in BOX_FACES], spacing=0.5)
    assert domains == [f"box-{s}-mesh" for s in BOX_FACES]
    assert len(og.names("line")) == 12 and len(og.names("connector")) == 12  # one each, shared
    for name in domains:
        mesh = og.unstructured_mesh(name)
        assert len(mesh["triangles"]) >= 20
    assert og.entity("box-xmin-mesh").params["surface"] is None  # flat: planar domain
    og.undo()
    assert og.names("connector") == [] and og.names("line") == []
    og.redo()
    with pytest.raises(ValidationError, match="give a spacing or a number of points"):
        og.mesh_surfaces(["box-xmin"])
    with pytest.raises(ValidationError, match="positive length"):
        og.mesh_surfaces(["box-xmin"], spacing=-1)
    assert og.mesh_surfaces(["box-xmin"], points=4) == ["box-xmin-mesh-2"]  # connectors reused
    assert len(og.names("connector")) == 12
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_cylinder_in_a_box_becomes_a_volume():
    og = opengrid.Session()
    og.create_box(size=(4, 4, 4), name="box")
    parts = og.create_cylinder(name="cy")
    outer = og.mesh_surfaces([f"box-{s}" for s in BOX_FACES], spacing=0.6)
    curved = [n for n in parts if og.entity(n).kind == "revolved_surface"]
    inner = og.mesh_surfaces(curved, spacing=0.12)
    assert inner == ["cy-side-a-mesh", "cy-side-b-mesh", "cy-bottom-mesh", "cy-top-mesh"]
    assert og.entity("cy-side-a-mesh").params["surface"] == og.entity("cy-side-a").id
    assert og.entity("cy-top-mesh").params["surface"] is None  # a flat disc
    assert len(og.names("arc")) == 4  # the rims made by the shape, found, none added
    side = og.unstructured_mesh("cy-side-a-mesh")
    assert max(abs(math.hypot(p[0], p[1]) - 0.5) for p in side["points"]) < 1e-9
    top = og.unstructured_mesh("cy-top-mesh")
    assert all(abs(p[2] - 1.0) < 1e-12 for p in top["points"])
    # The disc's loop is the two rim halves the side meshes made: shared connectors.
    loop = og.entity("cy-top-mesh").params["loops"][0]
    assert {og.store.get(i["connector"]).name for i in loop} == {"c-cy-top-rim-a", "c-cy-top-rim-b"}
    volume = og.create_volume(outer + inner, name="flow")
    status = og.volume_status(volume)
    assert "2 shell(s)" in status and og.volume_mesh(volume)["tets"].shape[0] > 1000
    with pytest.raises(ValidationError, match="full turn"):
        og.create_revolved_surface("cy-seam-a", math.tau, (0, 0, 0), (0, 0, 1), name="whole")
        og.mesh_surfaces(["whole"], spacing=0.1)


def test_wing_extruded_surface_and_sphere():
    og = opengrid.Session()
    parts = og.create_wing("0012", 1.0, 0.5, 2.0, name="wg")
    surfaces = [n for n in parts if og.entity(n).kind == "ruled_surface"]
    domains = og.mesh_surfaces(surfaces, spacing=0.1)
    assert domains == ["wg-upper-mesh", "wg-lower-mesh", "wg-root-mesh", "wg-tip-mesh"]
    lines = og.names("line")
    assert len(lines) == 2  # the leading and trailing edges, shared by upper and lower
    for cap in ("wg-root-mesh", "wg-tip-mesh"):
        assert og.entity(cap).params["surface"] is None  # flat caps: planar domains
        loop = og.entity(cap).params["loops"][0]
        assert len(loop) == 2  # upper and lower section curves, the zero-length sides dropped
    upper = og.entity("wg-upper-mesh").params["loops"][0]
    assert len(upper) == 4
    volume = og.create_volume(domains, name="solid")  # the wing alone is a closed shell
    assert og.volume_mesh(volume)["tets"].shape[0] > 100
    # An extruded surface: its translated copy is made as a curve.
    og.create_arc((5, 0, 0), 1.0, 0.0, math.pi / 2, name="quarter")
    og.create_extruded_surface("quarter", (0, 0, 2), name="sheet")
    (sheet,) = og.mesh_surfaces(["sheet"], points=6)
    assert len(og.names("arc")) == 2 and og.entity(sheet).params["surface"] == og.entity("sheet").id
    og.create_arc((0, 0, 0), 1.0, 0.0, math.pi, normal=(0, 1, 0), x_axis=(0, 0, 1), name="m")
    og.create_revolved_surface("m", math.tau, (0, 0, 0), (0, 0, 1), name="whole-sphere")
    with pytest.raises(ValidationError, match="two halves"):
        og.mesh_surfaces(["whole-sphere"], spacing=0.1)
    with pytest.raises(ValidationError, match=r"not a surface|mesh_surfaces takes"):
        og.mesh_surfaces(["quarter"], spacing=0.1)
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def test_mesh_surfaces_from_the_gui(qtbot):
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    og = window.session
    assert window.mesh_selected_surfaces(spacing=0.5) is None  # nothing selected: a hint
    og.create_box(size=(2, 2, 2), name="box")
    window.set_selection([og.entity(f"box-{s}").id for s in BOX_FACES])
    names = window.mesh_selected_surfaces(spacing=0.5)
    assert names == [f"box-{s}-mesh" for s in BOX_FACES]
    assert window.selected_names() == names
    assert window.create_volume_on_selection() == "volume-1"
    window.close()


def test_sphere_meshes_with_fans_at_the_poles():
    """A sphere's halves are meshed on the surface: the pole nodes are opened into their
    degenerate side so the triangles fan around them, and the mesh is cleaned up where the
    parametrisation crowds (collapses, smoothing on the surface, flips)."""
    import numpy as np

    from opengrid.model.unstructured import _min_angles

    og = opengrid.Session()
    assert og.create_sphere(radius=1.0, name="sp") == [
        "sp-meridian-a",
        "sp-meridian-b",
        "sp-a",
        "sp-b",
    ]
    domains = og.mesh_surfaces(["sp-a", "sp-b"], spacing=0.15)
    assert domains == ["sp-a-mesh", "sp-b-mesh"]
    assert og.names("connector") == ["c-sp-meridian-b", "c-sp-meridian-a"]  # shared by both
    for name in domains:
        mesh = og.unstructured_mesh(name)
        points = np.asarray(mesh["points"])
        triangles = np.asarray(mesh["triangles"])
        assert np.max(np.abs(np.linalg.norm(points, axis=1) - 1.0)) < 1e-9  # on the sphere
        assert 350 < len(triangles) < 700  # about the spacing's worth, not crowded at the poles
        assert _min_angles(points, triangles).min() > 15.0
        edges = np.sort(
            np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]),
            axis=1,
        )
        _unique, counts = np.unique(edges, axis=0, return_counts=True)
        assert counts.max() <= 2  # manifold
        loop = og.entity(name).params["loops"][0]
        boundary = sum(
            len(og.connector_nodes(og.store.get(i["connector"]).name).points) - 1 for i in loop
        )
        assert np.allclose(points[:boundary, 1], 0.0, atol=1e-9)  # the meridian nodes come first
        north = int(np.argmin(np.linalg.norm(points - np.asarray([0.0, 0.0, 1.0]), axis=1)))
        assert np.sum(np.any(triangles == north, axis=1)) >= 3  # a fan at the pole
    og.create_box(size=(6, 6, 6), name="box")
    outer = og.mesh_surfaces([f"box-{s}" for s in BOX_FACES], spacing=0.8)
    volume = og.create_volume(outer + domains, name="flow")
    status = og.volume_status(volume)
    assert "2 shell(s)" in status and "miss the quality floor" not in status
    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
