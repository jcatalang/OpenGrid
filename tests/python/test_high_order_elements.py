"""Phase 11f: straight-sided high-order elements of every cell type, Gmsh / CGNS / VTKHDF / HDF5.

Node orders are checked against each format's own reference: Gmsh's element properties
(tests/python/data/gmsh_reference_nodes.json, generated with Gmsh outside the repository), VTK's
Lagrange cells (live), and CGNS through Gmsh's CGNS writer (cgns_reference_nodes.json) and the
SIDS node descriptions (cgns_sids_statements.json, for PENTA_40 / PENTA_75 that Gmsh does not
write). The exports are read back and every element's nodes are compared with the linear map of
its vertices at the format's reference coordinates.
"""

import json
import re
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

import h5py
import numpy as np
import pytest
from gmsh_tool import find_gmsh, gmsh_read
from test_periodic import channel
from test_unstructured import hybrid

import opengrid
from opengrid import _core
from opengrid.io import lagrange

DATA = Path(__file__).resolve().parent / "data"
GMSH_FAMILY = {
    "line": "Line",
    "triangle": "Triangle",
    "quad": "Quadrangle",
    "tetra": "Tetrahedron",
    "hexa": "Hexahedron",
    "prism": "Prism",
    "pyramid": "Pyramid",
}
CGNS_TYPE_IDS = {  # ElementType_t values (cgnslib.h)
    "BAR_2": 3, "BAR_3": 4, "BAR_4": 24, "BAR_5": 40,
    "TRI_3": 5, "TRI_6": 6, "TRI_10": 26, "TRI_15": 42,
    "QUAD_4": 7, "QUAD_9": 9, "QUAD_16": 28, "QUAD_25": 44,
    "TETRA_4": 10, "TETRA_10": 11, "TETRA_20": 30, "TETRA_35": 47,
    "PYRA_5": 12, "PYRA_14": 13, "PYRA_30": 33, "PYRA_55": 50,
    "PENTA_6": 14, "PENTA_18": 16, "PENTA_40": 36, "PENTA_75": 53,
    "HEXA_8": 17, "HEXA_27": 19, "HEXA_64": 39, "HEXA_125": 56,
}  # fmt: skip


def _exact(kind, coords, p) -> list:
    """Canonical reference coordinates as exact fractions (multiples of 1 / (2 p))."""
    d = lagrange.DIM[kind]
    den = 2 * p
    return [tuple(Fraction(round(float(c) * den), den) for c in x[:d]) for x in coords]


def _from_gmsh(kind, c) -> np.ndarray:
    c = np.array(c, dtype=float)
    if kind in ("line", "quad", "hexa"):
        return (c + 1) / 2
    if kind == "prism":
        c[:, 2] = (c[:, 2] + 1) / 2
    if kind == "pyramid":
        c[:, :2] = (c[:, :2] + 1) / 2
    return c


def _linear_map(kind, reference, vertices) -> np.ndarray:
    """Positions of reference points (canonical) in elements with the given vertices
    (elements, vertices, 3)."""
    w = np.array([[float(v) for v in lagrange.weights(kind, x)] for x in reference])
    return np.einsum("qk,ekd->eqd", w, vertices)


# ------------------------------------------------------------------ node orders
@pytest.mark.parametrize("kind", lagrange.KINDS)
def test_gmsh_node_orders(kind):
    ref = json.loads((DATA / "gmsh_reference_nodes.json").read_text())["nodes"]
    d = lagrange.DIM[kind]
    for p in range(1, 7):
        entry = ref[f"{GMSH_FAMILY[kind]}:{p}"]
        assert entry["type"] == lagrange.gmsh_type(kind, p)
        assert entry["nodes"] == lagrange.count(kind, p)
        coords = np.array(entry["coords"]).reshape(entry["nodes"], -1)[:, :d]
        assert list(lagrange.nodes("gmsh", kind, p)) == _exact(kind, _from_gmsh(kind, coords), p)


@pytest.mark.parametrize("kind", ["line", "triangle", "quad", "tetra", "hexa", "prism"])
def test_vtk_node_orders(kind):
    vtk = pytest.importorskip("vtk")
    classes = {
        "line": vtk.vtkLagrangeCurve,
        "triangle": vtk.vtkLagrangeTriangle,
        "quad": vtk.vtkLagrangeQuadrilateral,
        "tetra": vtk.vtkLagrangeTetra,
        "hexa": vtk.vtkLagrangeHexahedron,
        "prism": vtk.vtkLagrangeWedge,
    }
    for p in range(1, 7):
        n = lagrange.count(kind, p)
        cell = classes[kind]()
        cell.GetPointIds().SetNumberOfIds(n)
        cell.GetPoints().SetNumberOfPoints(n)
        for i in range(n):
            cell.GetPointIds().SetId(i, i)
        if kind == "quad":
            cell.SetOrder(p, p)
        elif kind == "hexa":
            cell.SetOrder(p, p, p)
        elif kind == "prism":
            cell.SetOrder(p, p, p, n)
        cell.Initialize()
        pc = cell.GetParametricCoords()
        coords = np.array([pc[k] for k in range(3 * n)]).reshape(n, 3)
        assert list(lagrange.nodes("vtk", kind, p)) == _exact(kind, coords, p), (kind, p)


def test_cgns_node_orders():
    ref = json.loads((DATA / "cgns_reference_nodes.json").read_text())["nodes"]
    checked = 0
    for kind in lagrange.KINDS:
        for p in range(1, lagrange.CGNS_MAX_ORDER + 1):
            entry = ref.get(f"{kind}:{p}")
            if entry is None:
                assert kind == "prism" and p > 2  # Gmsh writes no PENTA_40 / PENTA_75
                continue
            assert entry["cgns_type"] == CGNS_TYPE_IDS[lagrange.cgns_type(kind, p)]
            assert list(lagrange.nodes("cgns", kind, p)) == _exact(kind, entry["coords"], p)
            checked += 1
    assert checked == 26


def test_cgns_orders_match_the_sids_descriptions():
    """Each 'N is located at the centroid of sub-triangle A-B-C' names the node of the
    sub-entity of A, B, C nearest that centroid; 'at the intersection of the line A-B and the
    face C-D-E-F' the interior node nearest that point."""
    statements = json.loads((DATA / "cgns_sids_statements.json").read_text())
    count = 0
    for key, rows in [*statements["centroid"].items(), *statements["line_face"].items()]:
        kind, p = key.split(":")
        p = int(p)
        nodes = lagrange.nodes("cgns", kind, p)
        x = lagrange.as_array(nodes)
        support = [
            frozenset(k for k, w in enumerate(lagrange.weights(kind, n)) if w) for n in nodes
        ]
        line_face = key in statements["line_face"] and len(rows[0]) == 7
        for row in rows:
            if line_face:
                n, a, b, *face = row
                pa, pb = x[a - 1], x[b - 1]
                z = x[face[0] - 1][2]
                target = pa + (z - pa[2]) / (pb[2] - pa[2]) * (pb - pa)
                want = frozenset(range(len(lagrange.VERTICES[kind])))
            else:
                n, *tri = row
                target = x[np.array(tri) - 1].mean(axis=0)
                want = frozenset().union(*(support[k - 1] for k in tri))
            distance = np.linalg.norm(x - target, axis=1)
            distance[[s != want for s in support]] = np.inf
            order = np.argsort(distance, kind="stable")
            assert order[0] + 1 == n, (key, row)
            assert distance[order[1]] - distance[order[0]] > 1e-9, (key, row)
            count += 1
    assert count == 46


@pytest.mark.parametrize("fmt", ["gmsh", "vtk", "cgns"])
def test_permutations_are_bijections_keeping_vertices_first(fmt):
    for kind in lagrange.KINDS:
        if fmt == "vtk" and kind == "pyramid":
            with pytest.raises(ValueError, match="no order"):
                lagrange.permutation(fmt, kind, 2)
            continue
        top = 6
        if fmt == "cgns":
            top = 4
        for p in range(1, top + 1):
            perm = lagrange.permutation(fmt, kind, p)
            nv = len(lagrange.VERTICES[kind])
            assert sorted(perm.tolist()) == list(range(lagrange.count(kind, p)))
            assert perm[:nv].tolist() == list(range(nv))


# ------------------------------------------------------------------ meshes
def _check_mesh(mesh):
    """Every element's nodes are its linear map at the reference nodes; nodes are distinct
    and all used; boundary elements share the cells' nodes."""
    p = mesh.order
    used = np.zeros(len(mesh.points), dtype=bool)
    for table in (mesh.cells, mesh.faces):
        for kind, conn in table.items():
            nv = len(lagrange.VERTICES[kind])
            want = _linear_map(kind, lagrange.nodes("gmsh", kind, p), mesh.points[conn[:, :nv]])
            assert np.abs(want - mesh.points[conn]).max() < 1e-13
            if table is mesh.cells:
                used[conn.reshape(-1)] = True
    assert used.all()
    for conn in mesh.faces.values():
        assert used[conn.reshape(-1)].all()
    scale = np.abs(mesh.points).max()
    keys = np.round(mesh.points / (1e-9 * scale)).astype(np.int64)
    assert len(np.unique(keys, axis=0)) == len(mesh.points)


def _expected_nodes(mesh) -> int:
    """V + E (p - 1) + triangle faces (p - 1)(p - 2)/2 + quad faces (p - 1)^2 + interiors."""
    p = mesh.order
    edges, tri, quad = set(), set(), set()
    interior = 0
    for kind, conn in mesh.cells.items():
        nv = len(lagrange.VERTICES[kind])
        spec = lagrange.GMSH[kind]
        corners = conn[:, :nv].tolist()
        for row in corners:
            for a, b in spec["edges"]:
                edges.add(frozenset((row[a], row[b])))
            for face in spec.get("faces", ()):
                key = frozenset(row[v] for v in face)
                if len(face) == 3:
                    tri.add(key)
                else:
                    quad.add(key)
        total = lagrange.count(kind, p)
        boundary = nv + len(spec["edges"]) * (p - 1)
        for face in spec.get("faces", ()):
            if len(face) == 3:
                boundary += max(0, (p - 1) * (p - 2) // 2)
            else:
                boundary += (p - 1) ** 2
        interior += len(conn) * (total - boundary)
    return (
        mesh.vertices
        + len(edges) * (p - 1)
        + len(tri) * max(0, (p - 1) * (p - 2) // 2)
        + len(quad) * (p - 1) ** 2
        + interior
    )


def test_2d_hybrid_elements():
    og = hybrid()
    for p in (1, 2, 3, 5):
        mesh = og.high_order(p)
        assert set(mesh.cells) == {"triangle", "quad"}
        assert set(mesh.faces) == {"line"}
        assert mesh.vertices == len(og.high_order(1).points)
        assert len(mesh.points) == _expected_nodes(mesh)
        _check_mesh(mesh)
    with pytest.raises(ValueError, match="order must be"):
        og.high_order(7)


def mixed_3d():
    """Hexahedra (a block), tetrahedra and pyramids (a volume on its structured face)."""
    pytest.importorskip("netgen.meshing")
    from test_mixed3d import block_and_volume

    return block_and_volume()[0]


def layered_3d():
    """Tetrahedra, pyramids, prisms and hexahedra: a body in a box with hexahedral layers on
    its structured walls (prisms at the corners)."""
    pytest.importorskip("netgen.meshing")
    from test_inflation3d import body_in_box

    return body_in_box("structured", n=5, layers=3)[0]


def test_3d_elements_of_every_kind():
    og = layered_3d()
    for p in (2, 4):
        mesh = og.high_order(p)
        assert set(mesh.cells) == {"tetra", "pyramid", "prism", "hexa"}
        assert set(mesh.faces) == {"triangle", "quad"}
        assert len(mesh.points) == _expected_nodes(mesh)
        _check_mesh(mesh)


def _read_msh(path):
    """(points, [(dim, entity, type, (count, nodes) 0-based)], physical names, periodic)."""
    text = Path(path).read_text().split("\n")
    i = text.index("$Nodes")
    blocks, n, *_ = map(int, text[i + 1].split())
    assert blocks == 1
    tags = np.array(text[i + 3 : i + 3 + n], dtype=np.int64)
    assert (tags == np.arange(1, n + 1)).all()
    points = np.array([row.split() for row in text[i + 3 + n : i + 3 + 2 * n]], dtype=float)
    i = text.index("$Elements")
    count = int(text[i + 1].split()[0])
    k = i + 2
    out = []
    for _ in range(count):
        dim, entity, kind, m = map(int, text[k].split())
        rows = np.array([r.split()[1:] for r in text[k + 1 : k + 1 + m]], dtype=np.int64) - 1
        out.append((dim, entity, kind, rows))
        k += 1 + m
    i = text.index("$PhysicalNames")
    names = {}
    for row in text[i + 2 : i + 2 + int(text[i + 1])]:
        dim, tag, name = row.split(maxsplit=2)
        names[(int(dim), int(tag))] = name.strip('"')
    periodic = []
    if "$Periodic" in text:
        i = text.index("$Periodic")
        k = i + 2
        for _ in range(int(text[i + 1])):
            matrix = np.array(text[k + 1].split()[1:], dtype=float).reshape(4, 4)
            m = int(text[k + 2])
            pairs = np.array([r.split() for r in text[k + 3 : k + 3 + m]], dtype=np.int64) - 1
            periodic.append((text[k], matrix, pairs))
            k += 3 + m
    return points, out, names, periodic


def _gmsh_reference(kind, p) -> list:
    ref = json.loads((DATA / "gmsh_reference_nodes.json").read_text())["nodes"]
    entry = ref[f"{GMSH_FAMILY[kind]}:{p}"]
    coords = np.array(entry["coords"]).reshape(entry["nodes"], -1)[:, : lagrange.DIM[kind]]
    return _exact(kind, _from_gmsh(kind, coords), p)


def test_gmsh_export(tmp_path):
    og = mixed_3d()
    p = 3
    path = tmp_path / "mixed.msh"
    info = og.export_high_order(path, order=p)
    points, blocks, names, _ = _read_msh(path)
    kinds = {}
    for dim, entity, gmsh_type, rows in blocks:
        kind = next(k for k, types in lagrange.GMSH_TYPES.items() if gmsh_type in types)
        assert lagrange.GMSH_TYPES[kind].index(gmsh_type) + 1 == p
        assert (dim, entity) in names
        nv = len(lagrange.VERTICES[kind])
        want = _linear_map(kind, _gmsh_reference(kind, p), points[rows[:, :nv]])
        assert np.abs(want - points[rows]).max() < 1e-13
        label = kind
        if dim == 2:
            label = f"boundary {kind}"
        kinds[label] = kinds.get(label, 0) + len(rows)
    assert kinds == og.high_order(p).counts()
    assert set(kinds) == {"tetra", "pyramid", "hexa", "boundary triangle", "boundary quad"}
    assert info["elements"] == sum(kinds.values())
    if find_gmsh():
        log = gmsh_read(path)
        assert "Error" not in log and "Warning" not in log, log
        assert re.search(rf"\b{info['nodes']} nodes", log)


def _cgns_read(path):
    """{zone: (points (n, 3), {section: (type id, (count, nodes) 0-based, element range)},
    {bc name: (family, location, range)}, {connection: (donor, points, donor points,
    translation)})}."""
    out = {}
    with h5py.File(path, "r") as f:
        base = next(k for k in f if f[k].attrs.get("label") == b"CGNSBase_t")
        for zname, zone in f[base].items():
            if zone.attrs.get("label") != b"Zone_t":
                continue
            grid = zone["GridCoordinates"]
            axes = [a for a in ("CoordinateX", "CoordinateY", "CoordinateZ") if a in grid]
            xyz = np.zeros((len(grid[axes[0]][" data"]), 3))
            for k, a in enumerate(axes):
                xyz[:, k] = grid[a][" data"][()]
            sections = {}
            for sname, s in zone.items():
                if s.attrs.get("label") != b"Elements_t":
                    continue
                type_id = int(s[" data"][0])
                conn = s["ElementConnectivity"][" data"][()] - 1
                rng = s["ElementRange"][" data"][()]
                sections[sname] = (type_id, conn, rng)
            bcs = {}
            if "ZoneBC" in zone:
                for bname, bc in zone["ZoneBC"].items():
                    family = bc["FamilyName"][" data"][()].tobytes().decode()
                    location = bc["GridLocation"][" data"][()].tobytes().decode()
                    bcs[bname] = (family, location, bc["PointRange"][" data"][()])
            conns = {}
            if "ZoneGridConnectivity" in zone:
                for cname, c in zone["ZoneGridConnectivity"].items():
                    donor = c[" data"][()].tobytes().decode()
                    shift = None
                    if "GridConnectivityProperty" in c:
                        shift = c["GridConnectivityProperty/Periodic/Translation"][" data"][()]
                    conns[cname] = (
                        donor,
                        c["PointList"][" data"][()].reshape(-1) - 1,
                        c["PointListDonor"][" data"][()].reshape(-1) - 1,
                        shift,
                    )
            out[zname] = (xyz, sections, bcs, conns)
    return out


def _cgns_reference(kind, p) -> list:
    """The CGNS node order: Gmsh's CGNS writer, or the SIDS-checked order of PENTA_40/75."""
    ref = json.loads((DATA / "cgns_reference_nodes.json").read_text())["nodes"]
    entry = ref.get(f"{kind}:{p}")
    if entry is None:
        return list(lagrange.nodes("cgns", kind, p))
    return _exact(kind, entry["coords"], p)


def _check_cgns_elements(zones, p):
    names = {v: k for k, v in CGNS_TYPE_IDS.items()}
    kinds = {v: k for k, v in lagrange.CGNS_PREFIX.items()}
    for xyz, sections, _bcs, _conns in zones.values():
        for type_id, conn, rng in sections.values():
            prefix, n = names[type_id].split("_")
            kind = kinds[prefix]
            assert int(n) == lagrange.count(kind, p)
            rows = conn.reshape(-1, int(n))
            assert len(rows) == rng[1] - rng[0] + 1
            nv = len(lagrange.VERTICES[kind])
            want = _linear_map(kind, _cgns_reference(kind, p), xyz[rows[:, :nv]])
            assert np.abs(want - xyz[rows]).max() < 1e-13


def _cgnscheck(path):
    checker = shutil.which("cgnscheck")
    if checker is None:
        return
    result = subprocess.run(
        [checker, str(path)], capture_output=True, text=True, check=False, timeout=300
    )
    log = result.stdout + result.stderr
    assert result.returncode == 0, log
    assert "ERROR" not in log.upper() and "WARNING" not in log.upper(), log


@pytest.mark.skipif(not _core.cgns_available(), reason="built without CGNS")
def test_cgns_export_3d(tmp_path):
    og = layered_3d()
    for p in (2, 4):
        path = tmp_path / f"layers_p{p}.cgns"
        info = og.export_high_order(path, order=p)
        zones = _cgns_read(path)
        assert list(zones) == ["Elements"]
        xyz, sections, bcs, _ = zones["Elements"]
        types = {lagrange.CGNS_PREFIX[k] for k in ("tetra", "pyramid", "prism", "hexa")}
        type_names = {v: k for k, v in CGNS_TYPE_IDS.items()}
        present = {type_names[t].split("_")[0] for t, _c, _r in sections.values()}
        assert types | {"TRI", "QUAD"} == present
        assert {family for family, _loc, _r in bcs.values()} == {"far", "body"}
        assert {loc for _f, loc, _r in bcs.values()} == {"FaceCenter"}
        assert len(xyz) == info["nodes"]
        _check_cgns_elements(zones, p)
        _cgnscheck(path)


@pytest.mark.skipif(not _core.cgns_available(), reason="built without CGNS")
def test_cgns_export_2d_zones_and_periodic(tmp_path):
    og = channel()
    og.pair_periodic("left", "right")
    og.create_zone("bottom", ["lower"])
    og.create_zone("top", ["upper"])
    path = tmp_path / "channel.cgns"
    og.export_high_order(path, order=3)
    zones = _cgns_read(path)
    assert sorted(zones) == ["Elements_bottom", "Elements_top"]
    with h5py.File(path, "r") as f:
        assert f["OpenGrid"][" data"][()].tolist() == [2, 2]  # cell, physical dimension
        assert "CoordinateZ" not in f["OpenGrid/Elements_top/GridCoordinates"]
    _check_cgns_elements(zones, 3)
    for xyz, _sections, bcs, conns in zones.values():
        assert {loc for _f, loc, _r in bcs.values()} == {"EdgeCenter"}
        assert {family for family, _loc, _r in bcs.values()} == {"walls"}  # periodic: no BC
        for donor, here, there, shift in conns.values():
            donor_xyz = zones[donor][0]
            moved = donor_xyz[there]
            if shift is not None:
                moved = moved - np.array([*shift, 0.0])[: moved.shape[1]]
            assert np.abs(xyz[here] - moved).max() < 1e-12
    # Interface (shared nodes of the middle line at order 3: 21 vertices + 20 x 2) and
    # periodic links in both zones.
    shared = [c for z in zones.values() for name, c in z[3].items() if c[3] is None]
    assert {len(c[1]) for c in shared} == {21 + 20 * 2}
    periodic = [c for z in zones.values() for c in z[3].values() if c[3] is not None]
    assert len(periodic) == 4  # (left -> right and back) x (two zones)
    _cgnscheck(path)


@pytest.mark.skipif(not _core.cgns_available(), reason="built without CGNS")
def test_cgns_limits(tmp_path):
    og = hybrid()
    with pytest.raises(opengrid.ValidationError, match="up to order 4"):
        og.export_high_order(tmp_path / "x.cgns", order=5)
    with pytest.raises(ValueError, match="equispaced"):
        og.export_high_order(tmp_path / "x.msh", order=2, nodes="gll")
    with pytest.raises(ValueError, match="unknown high-order format"):
        og.export_high_order(tmp_path / "x.vtu", order=2)
    zone = ("z", [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0] * 3, [("s", "TRI_7", [0, 1, 2], "")], "")
    with pytest.raises(ValueError, match="unknown element type 'TRI_7'"):
        _core.write_cgns_elements(str(tmp_path / "y.cgns"), "B", 2, [zone], [], [])
    zone = ("z", [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0] * 3, [("s", "TRI_3", [0, 1, 2], "")], "")
    with pytest.raises(ValueError, match="cell dimension"):
        _core.write_cgns_elements(str(tmp_path / "y.cgns"), "B", 4, [zone], [], [])
    _core.write_cgns_elements(str(tmp_path / "y.cgns"), "B", 2, [zone], [], [])


def test_gmsh_periodic_links(tmp_path):
    og = channel()
    og.pair_periodic("left", "right")
    p = 4
    path = tmp_path / "channel.msh"
    info = og.export_high_order(path, order=p)
    assert info["periodic_links"] == 1
    points, _blocks, _names, periodic = _read_msh(path)
    ((header, matrix, pairs),) = periodic
    assert header.split()[0] == "1"
    # The right side: 21 vertices, 20 edges with p - 1 nodes each.
    assert len(pairs) == 21 + 20 * (p - 1)
    master = np.c_[points[pairs[:, 1]], np.ones(len(pairs))]
    assert np.abs((matrix @ master.T).T[:, :3] - points[pairs[:, 0]]).max() < 1e-12


def test_vtkhdf_export(tmp_path):
    vtk = pytest.importorskip("vtk")
    og = layered_3d()
    p = 3
    path = tmp_path / "layers.vtkhdf"
    info = og.export_high_order(path, order=p)
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfCells() == info["cells"]
    assert grid.GetNumberOfPoints() == info["points"]
    kinds = {68: "line", 69: "triangle", 70: "quad", 71: "tetra", 72: "hexa", 73: "prism"}
    rng = np.random.default_rng(0)
    seen = set()
    sub = vtk.reference(0)
    for c in range(0, grid.GetNumberOfCells(), 7):
        cell = grid.GetCell(c)
        t = cell.GetCellType()
        seen.add(t)
        n = cell.GetNumberOfPoints()
        vertices = np.array([grid.GetPoint(cell.GetPointId(k)) for k in range(n)])
        if t == 14:  # pyramid: linear (exact for a straight-sided one)
            assert n == 5
            continue
        kind = kinds[t]
        pc = rng.random(3) * 0.3
        weights = [0.0] * n
        x = [0.0, 0.0, 0.0]
        cell.EvaluateLocation(sub, list(pc), x, weights)
        nv = len(lagrange.VERTICES[kind])
        want = _linear_map(kind, [tuple(pc[: lagrange.DIM[kind]])], vertices[None, :nv])
        assert np.abs(np.array(x) - want[0, 0]).max() < 1e-12 * max(1.0, np.abs(vertices).max())
    assert seen == {71, 72, 73, 14}


def test_h5_node_sets(tmp_path):
    og = channel()  # structured quads only
    og.pair_periodic("left", "right")
    p = 3
    mesh = og.high_order(p)
    gll = mesh.element_nodes("quad", "gll")
    t = (np.asarray(_core.nodes_1d("gll", p + 1)) + 1) / 2
    reference = [(t[i], t[j]) for j in range(p + 1) for i in range(p + 1)]
    corners = mesh.points[mesh.cells["quad"][:, :4]]
    assert np.abs(gll - _linear_map("quad", reference, corners)).max() < 1e-14
    path = tmp_path / "channel_gll.h5"
    og.export_high_order(path, order=p, nodes="gll")
    with h5py.File(path) as f:
        assert f.attrs["node_set"] == "gll"
        assert np.allclose(f["quad/nodes"][()], gll)
        assert "points" not in f
    og.export_high_order(tmp_path / "e.h5", order=p)
    with h5py.File(tmp_path / "e.h5") as f:
        assert np.array_equal(f["quad/connectivity"][()], mesh.cells["quad"])
        assert "boundary/line/patch_id" in f
    with pytest.raises(ValueError, match="quadrilaterals and hexahedra"):
        hybrid().high_order(2).element_nodes("triangle", "gauss")


def test_gui_export(qtbot, tmp_path):
    from opengrid.app.main_window import MainWindow

    window = MainWindow(off_screen=True)
    qtbot.addWidget(window)
    window.session = channel()
    window.session.pair_periodic("left", "right")
    summary = window.export_high_order(str(tmp_path / "c.msh"), order=2)
    assert summary["order"] == 2 and summary["elements"] > 0
