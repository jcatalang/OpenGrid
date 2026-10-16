"""Phase 11g: curved boundary elements - shape functions, the displacement extension, boundary
nodes on the geometry, convergence, thin layers, validity reports, exports."""

import inspect
import math

import numpy as np
import pytest
from test_high_order_elements import _read_msh
from test_periodic import channel

import opengrid
from opengrid import _core
from opengrid.io import lagrange
from opengrid.model import curving, shape


# ------------------------------------------------------------------ shape functions
@pytest.mark.parametrize("kind", lagrange.KINDS)
def test_shape_functions(kind):
    d = lagrange.DIM[kind]
    for p in range(1, 5):
        nodes = lagrange.as_array(lagrange.nodes("gmsh", kind, p))
        values, _ = shape.basis(kind, p, nodes)
        assert np.abs(values - np.eye(len(nodes))).max() < 1e-11  # Kronecker at the nodes
        points = shape.sample_points(kind, p)
        values, grads = shape.basis(kind, p, points)
        assert np.abs(values.sum(axis=1) - 1).max() < 1e-11
        assert np.abs(values @ nodes - points).max() < 1e-11  # linear maps reproduced
        assert np.abs(np.einsum("qnd,nc->qcd", grads, nodes) - np.eye(d)).max() < 1e-10
        x0 = points[len(points) // 2][None, :] * 0.9 + 0.1 * nodes.mean(axis=0)
        h = 1e-6
        for a in range(d):
            e = np.zeros(d)
            e[a] = h
            fd = (shape.basis(kind, p, x0 + e)[0] - shape.basis(kind, p, x0 - e)[0]) / (2 * h)
            assert np.abs(fd[0] - shape.basis(kind, p, x0)[1][0, :, a]).max() < 1e-6


def test_extend_displacement_binding():
    x = np.column_stack([np.linspace(0, 1, 11), np.zeros(11), np.zeros(11)])
    edges = np.column_stack([np.arange(10), np.arange(1, 11)]).astype(np.int64)
    fixed = np.zeros(11, dtype=np.uint8)
    fixed[[0, 10]] = 1
    shift = np.zeros((11, 3))
    shift[0, 1] = 1.0
    radius = np.zeros(11)
    radius[0] = 10.0
    out, solved, seeds, iterations = _core.extend_displacement(x, edges, fixed, shift, radius, 2.0)
    assert (solved, seeds, iterations) == (9, 1, 0)
    assert np.allclose(out[:, 1], 1 - x[:, 0])
    with pytest.raises(ValueError, match="out of range or degenerate"):
        _core.extend_displacement(x, np.array([[0, 0]], dtype=np.int64), fixed, shift, radius)


# ------------------------------------------------------------------ 2D
def ring(n: int, first: float = 0.002, radial: int = 0) -> opengrid.Session:
    og = opengrid.Session()
    og.create_circle((0, 0), 1.0, name="wall")
    og.create_circle((0, 0), 4.0, name="far")
    og.create_line((1, 0), (4, 0), name="cut")
    og.create_connector("wall", n + 1, name="c-wall")
    og.create_connector("far", n + 1, name="c-far")
    if radial == 0:
        radial = n // 4
    og.create_connector("cut", radial + 1, {"type": "tanh", "begin_spacing": first}, name="c-cut")
    og.create_structured_domain(connectors=["c-wall", "c-far", "c-cut"], name="ring")
    og.create_patch("p-wall", "wall", ["c-wall"])
    og.create_patch("p-far", "farfield", ["c-far"])
    return og


def edge_error(mesh, patch: str, radius: float) -> float:
    """Largest distance from the circle of the patch's element edges (sampled)."""
    lines = mesh.faces["line"][mesh.face_patch["line"] == mesh.patch_names.index(patch)]
    values, _ = shape.basis("line", mesh.order, np.linspace(0, 1, 41)[:, None])
    points = np.einsum("qn,enc->eqc", values, mesh.points[lines])
    return float(np.abs(np.linalg.norm(points[..., :2], axis=-1) - radius).max())


def test_curved_walls_converge_at_high_order():
    rates = {}
    for p in (1, 2, 3, 4):
        straight, curved = [], []
        for n in (16, 32):
            og = ring(n)
            flat = og.high_order(p, curved=False)
            mesh = og.high_order(p)
            straight.append(edge_error(flat, "p-wall", 1.0))
            curved.append(edge_error(mesh, "p-wall", 1.0))
            if p > 1:  # the far field is curved too
                assert edge_error(mesh, "p-far", 4.0) < 1e-2 * edge_error(flat, "p-far", 4.0)
            report = mesh.curving
            if p > 1:
                assert report["min_scaled_jacobian"] > report["min_scaled_jacobian_straight"] - 0.1
        rates[p] = math.log2(curved[0] / curved[1])
        assert math.log2(straight[0] / straight[1]) == pytest.approx(2.0, abs=0.1)
    # h^(p+1); even orders gain one more on a circle.
    assert rates[1] == pytest.approx(2.0, abs=0.1)
    assert rates[2] > 3.8 and rates[3] > 3.8 and rates[4] > 5.8


def test_thin_layers_follow_the_wall():
    """A 1e-5 first cell under a coarse wall: the chord sagitta is 2000 cells thick."""
    og = ring(16, first=1e-5, radial=60)
    mesh = og.high_order(4)
    report = mesh.curving
    assert report["max_displacement"] > 1000 * 1e-5
    assert report["min_scaled_jacobian"] > 0.7
    wall = curving.boundary_nodes(mesh)
    r = np.linalg.norm(mesh.points[wall, :2], axis=1)
    on_wall = r < 1.5
    assert np.abs(r[on_wall] - 1.0).max() < 1e-14


def test_airfoil_c_grid_without_patches():
    import test_templates as cases

    og = opengrid.Session()
    body = cases.airfoil(og)
    og.create_template("cgrid", body=body, wall=1e-5, layers=24, radial_points=(41, 25), name="c")
    og.smooth_group("c-group")
    assert not og.names("patch")
    straight = og.high_order(3, curved=False)
    mesh = og.high_order(3)
    assert mesh.faces == {}  # no patches: no boundary elements, but curved walls
    curves = [og.store.geometry(og.store.get(b).id) for b in body]
    moved = np.flatnonzero(np.linalg.norm(mesh.points - straight.points, axis=1) > 1e-12)
    moved = np.intersect1d(moved, curving.boundary_nodes(mesh))  # the extension moves others
    near = [
        i for i in moved if min(c.project(tuple(straight.points[i]))[-1] for c in curves) < 1e-3
    ]
    assert len(near) > 50
    assert max(min(c.project(tuple(mesh.points[i]))[-1] for c in curves) for i in near) < 1e-12
    # The curved elements include the far field's (scaled Jacobian 0.48 straight-sided).
    report = mesh.curving
    assert report["min_scaled_jacobian"] > report["min_scaled_jacobian_straight"] - 0.05


def test_straight_option_and_periodic_links():
    og = channel()
    og.pair_periodic("left", "right")
    mesh = og.high_order(3)
    assert mesh.curving["moved"] == 0  # straight walls: nothing to move
    assert np.array_equal(mesh.points, og.high_order(3, curved=False).points)
    ((_s, _t, _tr, targets, sources),) = mesh.periodic
    assert np.allclose(mesh.points[targets], mesh.points[sources] + [2.0, 0.0, 0.0])


def test_curved_gmsh_export(tmp_path):
    og = ring(16)
    path = tmp_path / "ring.msh"
    og.export_high_order(path, order=4)
    points, blocks, _names, _ = _read_msh(path)
    lines = next(
        rows for dim, _e, t, rows in blocks if dim == 1 and t == lagrange.gmsh_type("line", 4)
    )
    r = np.linalg.norm(points[lines, :2], axis=-1)
    on_wall = r.mean(axis=1) < 2.0
    assert np.abs(r[on_wall] - 1.0).max() < 1e-14
    straight = tmp_path / "straight.msh"
    og.export_high_order(straight, order=4, curved=False)
    points, blocks, _names, _ = _read_msh(straight)
    lines = next(
        rows for dim, _e, t, rows in blocks if dim == 1 and t == lagrange.gmsh_type("line", 4)
    )
    r = np.linalg.norm(points[lines, :2], axis=-1)
    assert np.abs(r[r.mean(axis=1) < 2.0] - 1.0).max() > 1e-3  # the chords


# ------------------------------------------------------------------ 3D
def _sphere_faces(kind: str, n: int):
    """The six sphere faces of the 09b ball (structured or unstructured domains on the
    sphere's revolved surfaces), without their block."""
    import test_groups3d as g3

    source = inspect.getsource(g3.sphere_block).replace(
        '    og.create_block(faces=faces, name="sphere")\n', ""
    )
    if kind == "unstructured":
        source = source.replace("og.create_structured_domain(", "og.create_unstructured_domain(")
    namespace = dict(g3.__dict__)
    exec(source, namespace)  # the test module's own helper, with one change
    og = namespace["sphere_block"](n)
    return og, g3.R


def sphere_in_box(kind: str):
    pytest.importorskip("netgen.meshing")
    from test_volumes import box

    og, radius = _sphere_faces(kind, 7)
    sphere = list(og.names("domain")) + list(og.names("unstructured_domain"))
    outer = box(og, (-2, -2, -2), (2, 2, 2), 5, "o", "structured")
    inflation = {"faces": sphere, "first_height": 1e-3, "growth_rate": 1.3, "layers": 4}
    og.create_volume(
        outer + sphere, generator={"max_size": 0.6, "inflation": inflation}, name="fluid"
    )
    og.create_patch("far", "farfield", outer)
    og.create_patch("sphere-wall", "wall", sphere)
    return og, radius


def sphere_error(mesh, radius: float) -> float:
    worst = 0.0
    for kind, conn in mesh.faces.items():
        values, _ = shape.basis(kind, mesh.order, shape.sample_points(kind, mesh.order))
        x = np.einsum("qn,enc->eqc", values, mesh.points[conn])
        on = np.linalg.norm(mesh.points[conn], axis=-1).max(axis=1) < 1.5 * radius
        if on.any():
            worst = max(worst, float(np.abs(np.linalg.norm(x[on], axis=-1) - radius).max()))
    return worst


@pytest.mark.parametrize("kind", ["structured", "unstructured"])
def test_sphere_in_a_box_with_wall_layers(kind):
    og, radius = sphere_in_box(kind)
    straight = og.high_order(3, curved=False)
    mesh = og.high_order(3)
    expected = {"tetra"}
    if kind == "structured":
        expected |= {"hexa", "pyramid"}
    else:
        expected |= {"prism"}
    assert expected <= set(mesh.cells)
    assert sphere_error(mesh, radius) < 1e-2 * sphere_error(straight, radius)  # 1/150 .. 1/400
    report = mesh.curving
    # Valid; the largest loss is a tetrahedron on top of the layer stack, whose bottom face
    # takes the wall's curvature (0.97 -> 0.64 on the triangulated sphere).
    assert report["min_scaled_jacobian"] > 0.5
    assert report["min_scaled_jacobian"] > report["min_scaled_jacobian_straight"] - 0.4


def test_pipe_walls_from_the_grid_and_planar_caps():
    """Pipe template: block faces without surfaces (grid interpolation); the end caps are
    planar and stay straight."""
    og = opengrid.Session()
    og.create_template(
        "pipe", boundary=[og.create_circle((0, 0), 0.5, name="c")], length=2.0, span_points=5
    )
    og.smooth_group(og.names("smoothing_group")[0])
    errors = []
    for p in (2, 3):
        straight = og.high_order(p, curved=False)
        mesh = og.high_order(p)
        b = curving.boundary_nodes(straight)
        caps = b[(straight.points[b, 2] < 1e-12) | (straight.points[b, 2] > 2 - 1e-12)]
        rim = caps[np.abs(np.linalg.norm(straight.points[caps, :2], axis=1) - 0.5) < 0.05]
        inside = np.setdiff1d(caps, rim)
        assert np.array_equal(mesh.points[inside], straight.points[inside])  # planar: kept
        assert np.abs(mesh.points[caps, 2] - straight.points[caps, 2]).max() < 1e-14  # in plane
        wall = b[np.abs(np.linalg.norm(straight.points[b, :2], axis=1) - 0.5) < 0.05]
        wall = np.setdiff1d(wall, caps)
        errors.append(
            (
                np.abs(np.linalg.norm(straight.points[wall, :2], axis=1) - 0.5).max(),
                np.abs(np.linalg.norm(mesh.points[wall, :2], axis=1) - 0.5).max(),
            )
        )
        assert mesh.curving["min_scaled_jacobian"] == pytest.approx(
            mesh.curving["min_scaled_jacobian_straight"], abs=0.01
        )
    assert errors[0][1] < 0.05 * errors[0][0] and errors[1][1] < 0.01 * errors[1][0]


def test_folds_are_reported():
    """One block wrapped onto a sphere has nearly flat corner cells; curving folds them, and
    the validator says which and where."""
    import test_groups3d as g3

    og = g3.sphere_block(7)
    og.create_patch("wall", "wall", list(og.names("domain")))
    with pytest.raises(opengrid.ValidationError, match="folds after curving"):
        og.high_order(4)
    assert og.high_order(4, curved=False).elements == 216
