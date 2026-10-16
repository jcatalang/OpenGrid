"""Phase 11a: wall layers (the first layers off a wall marched normal to it and kept)."""

import numpy as np
import pytest
from test_templates import airfoil, channel, first_layer, worst

import opengrid
from opengrid.validation import ValidationError


def deviation(grid, k):
    """Deviation (degrees) from the normal of the grid line from row k to k + 1 off side
    i = 0, against the row's own tangent."""
    row, nxt = grid[:, k], grid[:, k + 1]
    t = row[2:] - row[:-2]
    v = nxt[1:-1] - row[1:-1]
    cos = np.abs(np.sum(t * v, axis=1)) / (np.linalg.norm(t, axis=1) * np.linalg.norm(v, axis=1))
    return np.degrees(np.arcsin(cos))


def test_cgrid_layers_keep_lines_normal_and_heights_geometric():
    og = opengrid.Session()
    og.create_template(
        "cgrid",
        body=airfoil(og),
        wall=1e-5,
        layers=30,
        growth=1.2,
        radial_points=(61, 31),
        name="c",
    )
    result = og.smooth_group("c-group")
    assert result.converged
    assert result.node_counts["layer"] > 0
    assert "in wall layers" in og.group_status("c-group")
    grid = og.domain_grid("c-d1")  # upper surface; i = 0 is the wall
    height, _ = first_layer(grid)
    assert np.allclose(height, 1e-5, rtol=1e-9)
    for k in (0, 8, 16, 28):
        dev = deviation(grid, k)
        assert dev.mean() < 0.1 and dev.max() < 0.5, (k, dev.mean(), dev.max())
        # Layer heights: spacing * growth^k along the line.
        step = np.linalg.norm(grid[2:-2, k + 1] - grid[2:-2, k], axis=1)
        assert np.allclose(step, 1e-5 * 1.2**k, rtol=1e-6)
    # Without layers the lines bend away from the normal within the same distance.
    og2 = opengrid.Session()
    og2.create_template("cgrid", body=airfoil(og2), wall=1e-5, radial_points=(61, 31), name="c")
    og2.smooth_group("c-group")
    assert deviation(og2.domain_grid("c-d1"), 16).mean() > 5.0
    assert worst(og) > 0.4


def test_3d_layers_on_an_extruded_bump_channel():
    og = opengrid.Session()
    og.create_template(
        "bump_channel",
        **channel(og, bump=True),
        blocks=(3, 1),
        wall=5e-3,
        length=1.0,
        three_d_blocks=True,
        layers=6,
    )
    result = og.smooth_group(og.names("smoothing_group")[0])
    assert result.converged
    assert result.node_counts["layer"] > 0
    assert worst(og) > 0.9


def _walls(og, sketch):
    return {w["edge"]: {"spacing": w["spacing"]} for w in og.entity(sketch).params["walls"]}


def test_layer_settings_are_validated_and_too_many_layers_reported():
    og = opengrid.Session()
    body = airfoil(og)
    with pytest.raises(ValidationError, match="whole number"):
        og.create_template("cgrid", body=body, wall=1e-3, layers=-1)
    with pytest.raises(ValidationError, match="at least 1"):
        og.create_template("cgrid", body=body, wall=1e-3, layers=5, growth=0.5)
    # Lines through the inner and outer rings have 61 nodes: 70 layers do not fit.
    og.create_template("cgrid", body=body, wall=1e-3, layers=70, name="c")
    with pytest.raises(ValueError, match="ends after"):
        og.smooth_group("c-group")
    # 40 layers of 1e-3 growing by 1.2 are thicker than the chord: they fold at the trailing
    # edge (reported, not clipped).
    og.edit_sketch("c", walls={k: {**w, "layers": 40} for k, w in _walls(og, "c").items()})
    with pytest.raises(ValueError, match="folds"):
        og.smooth_group("c-group")


def test_gui_wall_text_round_trips_layers():
    pytest.importorskip("PySide6")
    from opengrid.app.property_panel import format_group_walls, parse_group_walls

    og = opengrid.Session()
    og.create_line((0, 0), (1, 0), name="w")
    og.create_connector("w", dimension=11, name="cw")
    walls = parse_group_walls("cw 1e-5 30 1.15\ncw - 4")
    assert walls[0] == {"connector": "cw", "spacing": 1e-5, "layers": 30, "growth": 1.15}
    assert walls[1] == {"connector": "cw", "spacing": None, "layers": 4}
    named = [{**walls[0], "connector": og.entity("cw").id}]
    assert format_group_walls(named, og.store) == "cw 1e-05 30 1.15"
