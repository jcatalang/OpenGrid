"""Smoke test of an installed OpenGrid package: a structured grid, an unstructured domain, a
tetrahedral volume (Netgen) and the exports, without the GUI."""

import math
import tempfile
from pathlib import Path

import opengrid

og = opengrid.Session()
og.create_arc((0, 0, 0), 1.0, 0.0, math.pi / 2, name="inner")
og.create_arc((0, 0, 0), 2.0, 0.0, math.pi / 2, name="outer")
og.create_line((1, 0), (2, 0), name="bottom")
og.create_line((0, 2), (0, 1), name="left")
for curve in ("inner", "outer", "bottom", "left"):
    og.create_connector(curve, 21, name=f"c-{curve}")
og.create_domain([["c-bottom"], ["c-outer"], ["c-left"], [("c-inner", True)]], name="d")
og.set_domain_generator("d", elliptic={"walls": [{"edge": 0, "spacing": 0.01}]})
assert og.smooth_domain("d").converged
og.create_patch("p-all", "wall", ["c-inner", "c-outer", "c-bottom", "c-left"])

# Unit cube: twelve connectors, six unstructured faces, tetrahedra from Netgen.
cube = opengrid.Session()
vertices = {}
for c in ((i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)):
    vertices[c] = tuple(float(v) for v in c)
edges = {}


def connector(a, b):
    key = tuple(sorted((a, b)))
    if key not in edges:
        name = f"e{len(edges)}"
        cube.create_line(vertices[key[0]], vertices[key[1]], name=name)
        edges[key] = cube.create_connector(name, 7, name=f"c-{name}")
    return edges[key]


faces = []
for axis in range(3):
    for side in (0, 1):
        a, b = [x for x in range(3) if x != axis]
        loop = []
        for i, j in ((0, 0), (1, 0), (1, 1), (0, 1)):
            c = [0, 0, 0]
            c[axis] = side
            c[a] = i
            c[b] = j
            loop.append(tuple(c))
        names = [connector(loop[k], loop[(k + 1) % 4]) for k in range(4)]
        faces.append(cube.create_unstructured_domain(connectors=names, name=f"f{axis}{side}"))
cube.create_volume(faces, generator={"max_size": 0.2}, name="v")
cube.create_patch("p-walls", "wall", faces)
volume = cube.volume_quality("v")["volume"].sum()
assert abs(volume - 1.0) < 1e-12, volume

with tempfile.TemporaryDirectory() as tmp:
    og.export_vtkhdf(str(Path(tmp) / "d.vtkhdf"))
    og.export_cgns(str(Path(tmp) / "d.cgns"))
    cube.export_vtkhdf(str(Path(tmp) / "v.vtkhdf"))
print("OpenGrid smoke test passed:", og.domain_status("d"), "|", cube.volume_status("v"))
