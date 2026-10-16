"""Phase 10b: smoothing groups of blocks (multi-block hexahedral smoothing as one grid)."""

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pytest

import opengrid
from opengrid.model.blocks import FACE_AXIS, FACE_NAMES, face_array
from opengrid.project.ogp import to_json
from opengrid.validation import ValidationError

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _demo_09b():
    spec = importlib.util.spec_from_file_location("demo09b", EXAMPLES / "09b_topology_golden.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DEMO = _demo_09b()
R = DEMO.R


def ball():
    """Solid sphere of radius 0.5: core block + six shells (09b demo)."""
    og, _single = DEMO.ball_case()
    return og


def sphere_block(n=13):
    """The solid sphere as one six-face block (faces on the sphere; the 09b demo before its
    O-grid)."""
    og = opengrid.Session()
    og.create_arc(
        (0, 0, 0), R, -np.pi / 2, np.pi / 2, normal=(0, -1, 0), x_axis=(1, 0, 0), name="mer_z"
    )
    og.create_arc(
        (0, 0, 0), R, -np.pi / 2, np.pi / 2, normal=(0, 0, -1), x_axis=(0, 1, 0), name="mer_x"
    )
    og.create_revolved_surface(
        "mer_z", 2.5 * np.pi, axis=(0, 0, 1), start_angle=0.5 * np.pi, name="sphere_x"
    )
    og.create_revolved_surface("mer_z", 2 * np.pi, axis=(0, 0, 1), name="sphere_y")
    og.create_revolved_surface("mer_x", 2 * np.pi, axis=(1, 0, 0), name="sphere_z")
    signs = (-1, 1)
    corner = {
        (a, b, c): tuple(R * v / np.sqrt(3) for v in (a, b, c))
        for a in signs
        for b in signs
        for c in signs
    }
    arcs = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if key not in arcs:
            pa = np.array(corner[key[0]])
            pb = np.array(corner[key[1]])
            mid = R * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, n, name=f"c-{name}")
        return arcs[key]

    faces = []
    for axis, surface in ((0, "sphere_x"), (1, "sphere_y"), (2, "sphere_z")):
        for side in signs:
            a, b = [k for k in range(3) if k != axis]
            cs = []
            for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                c = [0, 0, 0]
                c[axis], c[a], c[b] = side, i, j
                cs.append(tuple(c))
            name = f"face-{'xyz'[axis]}{side}"
            og.create_structured_domain(
                connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                surface=surface,
                name=name,
            )
            faces.append(name)
    og.create_block(faces=faces, name="sphere")
    return og


def replays(og) -> bool:
    return to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)


def worst(og, names) -> float:
    return min(float(og.block_quality(n)["scaled_jacobian"].min()) for n in names)


def test_split_block_smoothed_as_a_group_matches_the_single_block():
    tolerance = 1e-12
    single = sphere_block()
    single.set_block_smoothing(
        "sphere", {"control": "laplace", "tolerance": tolerance, "max_iterations": 20000}
    )
    reference = single.block_grid("sphere")
    og = sphere_block()
    low, high, _face = og.split_block("sphere", "j", 5)
    group = og.create_smoothing_group(
        blocks=[low, high], slide=(), tolerance=tolerance, max_iterations=20000
    )
    result = og.smooth_group(group)
    assert result.converged
    assert result.node_counts["irregular"] == 0
    # Each part has its own (i, j, k) frame: compare the point sets.
    points = np.concatenate([og.block_grid(n).reshape(-1, 3) for n in (low, high)])
    points = np.unique(points + 0.0, axis=0)  # interface nodes are bitwise shared
    ref = reference.reshape(-1, 3)
    assert len(points) == len(ref)
    distance = np.sqrt(((ref[:, None, :] - points[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
    assert distance.max() < 1e-9


def test_ball_group_smooths_across_interfaces_and_slides_on_the_sphere():
    og = ball()
    names = og.names("block")
    before = worst(og, names)
    group = og.create_smoothing_group(blocks=names, name="g")
    result = og.smooth_group(group)
    assert result.converged
    assert result.node_counts["irregular"] > 0
    assert result.node_counts["sliding"] > 0
    assert set(result.roles.values()) == {"interface", "sliding"}
    after = worst(og, names)
    assert before < 0.6 and after >= 0.58
    assert result.max_kink < 3.0
    points = np.concatenate([og.block_grid(n).reshape(-1, 3) for n in names])
    radius = np.linalg.norm(points, axis=1)
    assert radius.max() <= R + 1e-12
    surface = radius > R * (1 - 1e-6)
    assert np.abs(radius[surface] - R).max() < 1e-9
    # Faces are shown where the group moved them.
    face = og.domain_grid("sphere-xmin")
    assert np.abs(np.linalg.norm(face, axis=-1) - R).max() < 1e-9
    assert "7 blocks" in og.group_status("g")
    assert og.blocks.evaluate(og.entity(names[0]).id)["group"] == og.entity("g").id


def test_block_group_is_thread_independent():
    before = opengrid.threads()
    digests = []
    try:
        for threads in (1, 4):
            opengrid.set_threads(threads)
            og = ball()
            names = og.names("block")
            og.smooth_group(og.create_smoothing_group(blocks=names))
            points = np.concatenate([og.block_grid(n) for n in names])
            digests.append(hashlib.sha256(points.tobytes()).hexdigest())
    finally:
        opengrid.set_threads(before)
    assert digests[0] == digests[1]


def test_block_group_walls_journal_project_and_edits(tmp_path):
    og = ball()
    names = og.names("block")
    group = og.create_smoothing_group(
        blocks=names, walls=[{"domain": "sphere-zmax", "spacing": 0.01}], name="g"
    )
    result = og.smooth_group(group)
    assert result.converged
    assert result.roles[og.entity("sphere-zmax").id] == "wall"
    # The first cell above the wall face has the requested height (away from its edges).
    shell = next(n for n in names if og.entity("sphere-zmax").id in og.entity(n).params["faces"])
    grid = og.block_grid(shell)
    for face in FACE_NAMES:
        a = face_array(grid, face)
        if np.allclose(np.linalg.norm(a, axis=-1), R) and a[..., 2].mean() > 0.2:
            break
    axis, side = FACE_AXIS[face]
    first = 0
    second = 1
    if side > 0:
        first = -1
        second = -2
    height = np.linalg.norm(
        np.take(grid, second, axis=2 - axis) - np.take(grid, first, axis=2 - axis), axis=-1
    )
    assert np.abs(height[3:-3, 3:-3] - 0.01).max() < 1e-4
    assert replays(og)
    path = tmp_path / "ball.ogp"
    og.save(path)
    loaded = opengrid.Session()
    loaded.load(path)
    assert to_json(loaded.store) == to_json(og.store)
    # Block edits keep the group up to date.
    single = sphere_block()
    single.create_smoothing_group(blocks=["sphere"], name="g")
    created = single.split_block("sphere", "k", 5)
    members = [single.store.get(b).name for b in single.entity("g").params["blocks"]]
    assert members == ["sphere", created[1]]
    single.undo()
    parts = single.ogrid_block("sphere", fraction=0.5, radial_points=5)
    members = [single.store.get(b).name for b in single.entity("g").params["blocks"]]
    assert sorted(members) == sorted(parts)
    assert single.smooth_group("g").converged
    assert replays(single)


def test_block_group_validation():
    og = ball()
    names = og.names("block")
    with pytest.raises(ValidationError, match="not both"):
        og.create_smoothing_group(domains=["sphere-xmin"], blocks=names)
    with pytest.raises(ValidationError, match="not a face of a member block"):
        og.create_smoothing_group(blocks=[names[0]], walls=["sphere-zmax"])
    og.create_extruded_block("sphere-xmin", "normal", distance=0.1, dimension=3, name="ext")
    with pytest.raises(ValidationError, match="extruded"):
        og.create_smoothing_group(blocks=["ext"])
