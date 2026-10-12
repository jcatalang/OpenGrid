"""Phase 8c: quality optimisation of volume meshes (node moves and tetrahedron swaps)."""

import numpy as np
import pytest
from test_volumes import box, box_with_body

import opengrid
from opengrid import _core
from opengrid.model.volumes import cell_volumes, total_volume, volume_cells

pytest.importorskip("netgen.meshing")


def qualities(mesh):
    points = np.ascontiguousarray(mesh["points"])
    return {
        key: np.array([_core.cell_quality_3d(points, c.tolist()) for c in cells])
        for key, cells, _vtk in volume_cells(mesh)
    }


def test_optimisation_improves_poor_cells():
    og = box_with_body()
    plain = og.volume_mesh("fluid")
    assert plain["optimization"] is not None
    og.set_volume_generator("fluid", optimize=False)
    raw = og.volume_mesh("fluid")
    assert raw["optimization"] is None
    stats = plain["optimization"]
    # Same boundary nodes, same enclosed volume, every cell valid.
    n = raw["n_surface"]
    assert plain["n_surface"] == n
    np.testing.assert_array_equal(plain["points"][:n], raw["points"][:n])
    for c in raw["faces"]:
        for a, b in zip(raw["faces"][c], plain["faces"][c], strict=True):
            np.testing.assert_array_equal(a, b)
    assert total_volume(plain) == pytest.approx(raw["enclosed_volume"], rel=1e-12)
    for key, cells, _vtk in volume_cells(plain):
        assert np.all(cell_volumes(plain["points"], key, cells) > 0)
    # The worst cell is not worse, and fewer cells are below the threshold.
    q_raw = np.concatenate(list(qualities(raw).values()))
    q_opt = np.concatenate(list(qualities(plain).values()))
    assert stats["worst_before"] == pytest.approx(q_raw.min())
    assert stats["worst_after"] == pytest.approx(q_opt.min())
    assert stats["worst_after"] >= stats["worst_before"]
    assert stats["below_before"] == int(np.sum(q_raw < 0.3))
    assert stats["below_after"] == int(np.sum(q_opt < 0.3))
    if stats["below_before"]:
        assert stats["moved"] > 0
        assert stats["below_after"] < stats["below_before"]
    assert "optimisation" not in " ".join(plain["warnings"])


def test_pyramid_apexes_move_layers_stay():
    og = opengrid.Session()
    outer = box(og, (0, 0, 0), (3, 3, 3), 7, "o", "unstructured")
    inner = box(og, (1, 1, 1), (2, 2, 2), 6, "i", "structured")
    og.create_volume(
        outer + inner,
        generator={
            "max_size": 0.5,
            "inflation": {"faces": inner, "first_height": 0.02, "layers": 4},
        },
        name="fluid",
    )
    mesh = og.volume_mesh("fluid")
    og.set_volume_generator("fluid", optimize=False)
    raw = og.volume_mesh("fluid")
    # Layer cells and their nodes are untouched.
    for key in ("hexes", "prisms"):
        np.testing.assert_array_equal(mesh[key], raw[key])
    layer_nodes = np.unique(np.concatenate([raw["hexes"].reshape(-1), raw["prisms"].reshape(-1)]))
    np.testing.assert_array_equal(mesh["points"][layer_nodes], raw["points"][layer_nodes])
    n_layer_tets = raw["inflation"]["cells"]["tets"]
    np.testing.assert_array_equal(mesh["tets"][:n_layer_tets], raw["tets"][:n_layer_tets])
    assert len(mesh["pyramids"]) == len(raw["pyramids"])
    assert total_volume(mesh) == pytest.approx(raw["enclosed_volume"], rel=1e-12)


def test_optimize_flag_is_validated_and_journaled():
    og = box_with_body()
    with pytest.raises(opengrid.ValidationError, match="optimize"):
        og.set_volume_generator("fluid", optimize="yes")
    og.set_volume_generator("fluid", optimize=False)
    assert "'optimize': False" in og.journal_text()
    replayed = opengrid.replay(og.journal_text())
    assert replayed.entity("fluid").params["generator"]["optimize"] is False
