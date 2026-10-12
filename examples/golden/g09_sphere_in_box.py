"""Golden case 9: sphere (radius 0.5) in a box: prism layers (first height 2e-3, growth 1.2, up
to 12 layers) on a cubed-sphere surface (six unstructured patches on revolved surfaces) and
Netgen tetrahedra in the core (needs the netgen package).

Checks: prisms and tetrahedra present, no negative volumes, filled volume equal to the volume
enclosed by the boundary to 1e-9, no cell below the quality floor, patches complete.
"""

import itertools
import math

import numpy as np

import opengrid
from opengrid.model.volumes import total_volume

R = 0.5


def box_faces(og, lo, hi, n, tag):
    corners = {
        c: tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3))
        for c in itertools.product((0, 1), repeat=3)
    }
    edges = {}

    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"{tag}{len(edges)}"
            og.create_line(corners[key[0]], corners[key[1]], name=name)
            edges[key] = og.create_connector(name, n, name=f"c{name}")
        return edges[key], key != (a, b)

    faces = {}
    for axis in range(3):
        for side in (0, 1):
            u, v = [x for x in range(3) if x != axis]
            loop = []
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                c = [0, 0, 0]
                c[axis], c[u], c[v] = side, a, b
                loop.append(tuple(c))
            name = f"{tag}_{'xyz'[axis]}{('min', 'max')[side]}"
            og.create_domain([[edge(loop[k], loop[(k + 1) % 4])] for k in range(4)], name=name)
            faces[name] = name
    return list(faces)


def sphere_patches(og, n):
    og.create_arc(
        (0, 0, 0), R, -math.pi / 2, math.pi / 2, normal=(0, -1, 0), x_axis=(1, 0, 0), name="mer_z"
    )
    og.create_arc(
        (0, 0, 0), R, -math.pi / 2, math.pi / 2, normal=(0, 0, -1), x_axis=(0, 1, 0), name="mer_x"
    )
    og.create_revolved_surface(
        "mer_z", 2.5 * math.pi, axis=(0, 0, 1), start_angle=0.5 * math.pi, name="sphere_x"
    )
    og.create_revolved_surface("mer_z", 2 * math.pi, axis=(0, 0, 1), name="sphere_y")
    og.create_revolved_surface("mer_x", 2 * math.pi, axis=(1, 0, 0), name="sphere_z")
    corner = {
        c: tuple(R * v / math.sqrt(3) for v in c) for c in itertools.product((-1, 1), repeat=3)
    }
    arcs = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if key not in arcs:
            pa, pb = np.array(corner[key[0]]), np.array(corner[key[1]])
            mid = R * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, n, name=f"c{name}")
        return arcs[key]

    patches = []
    for axis, surface in ((0, "sphere_x"), (1, "sphere_y"), (2, "sphere_z")):
        for side in (-1, 1):
            a, b = [k for k in range(3) if k != axis]
            cs = []
            for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                c = [0, 0, 0]
                c[axis], c[a], c[b] = side, i, j
                cs.append(tuple(c))
            name = f"sphere_{'xyz'[axis]}{side}"
            og.create_unstructured_domain(
                connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                surface=surface,
                name=name,
            )
            patches.append(name)
    return patches


def build() -> opengrid.Session:
    og = opengrid.Session()
    patches = sphere_patches(og, 13)
    outer = box_faces(og, (-2.5, -2.5, -2.5), (5.0, 2.5, 2.5), 7, "far")
    og.create_volume(
        outer + patches,
        generator={
            "max_size": 0.7,
            "inflation": {"faces": patches, "first_height": 2e-3, "growth_rate": 1.2, "layers": 12},
        },
        name="fluid",
    )
    og.create_patch("sphere", "wall", patches)
    og.create_patch("far", "farfield", outer)
    return og


def check(og) -> dict:
    mesh = og.volume_mesh("fluid")
    assert len(mesh["prisms"]) > 0 and len(mesh["tets"]) > 0
    filled = total_volume(mesh)
    assert abs(filled - mesh["enclosed_volume"]) < 1e-9 * filled
    assert mesh["quality"]["below_floor"] == 0
    assert og.volume_quality("fluid")["scaled_jacobian"].min() > 0
    assert og.boundary_report().ok
    return {
        "tets": len(mesh["tets"]),
        "prisms": len(mesh["prisms"]),
        "pyramids": len(mesh["pyramids"]),
        "volume": filled,
        "worst_jacobian": mesh["quality"]["worst_jacobian"],
    }


if __name__ == "__main__":
    print(check(build()))
