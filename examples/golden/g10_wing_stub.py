"""Golden case 10: swept NACA 0012 wing stub on a symmetry plane (root chord 1, tip 0.6, span
1.5): structured upper / lower surfaces on ruled surfaces (hexahedral layers), a planar tip
cap (prism layers), the symmetry plane re-meshed around the sliding layers, Netgen tetrahedra
in the core (needs the netgen package).

Checks: hexahedra, prisms and tetrahedra present, filled volume equal to the enclosed volume to
1e-9, no cell below the quality floor, patches complete.
"""

import math

import opengrid
from opengrid.model.volumes import total_volume

SPAN, TIP, SWEEP = 1.5, 0.6, 0.3


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def build() -> opengrid.Session:
    og = opengrid.Session()
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 80)) for i in range(81)]
    for z, chord, x0, tag in ((0.0, 1.0, 0.0, "root"), (SPAN, TIP, SWEEP, "tip")):
        for sign, side in ((1, "upper"), (-1, "lower")):
            og.create_spline(
                [(x0 + chord * x, sign * chord * naca0012(x), z) for x in xs],
                "clamped",
                start_tangent=(0, sign, 0),
                name=f"{tag}_{side}",
            )
            spacing = {"type": "tanh", "begin_spacing": 0.006 * chord, "end_spacing": 0.01 * chord}
            og.create_connector(f"{tag}_{side}", 31, spacing, name=f"c_{tag}_{side}")
    og.create_line((0, 0, 0), (SWEEP, 0, SPAN), name="le")
    og.create_line((1, 0, 0), (SWEEP + TIP, 0, SPAN), name="te")
    og.create_connector("le", 41, name="c_le")
    og.create_connector("te", 41, name="c_te")
    for side in ("upper", "lower"):
        og.create_ruled_surface(f"root_{side}", f"tip_{side}", name=f"surface_{side}")
        og.create_domain(
            [[f"c_root_{side}"], ["c_te"], [(f"c_tip_{side}", True)], [("c_le", True)]],
            surface=f"surface_{side}",
            name=f"wing_{side}",
        )
    og.create_unstructured_domain(connectors=["c_tip_upper", "c_tip_lower"], name="wing_tip")
    # Farfield box [-2, 4] x [-2.5, 2.5] x [0, 3]; the symmetry plane z = 0 holds the root.
    lo, hi = (-2.0, -2.5, 0.0), (4.0, 2.5, 3.0)
    corners = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                corners[(i, j, k)] = (
                    lo[0] + i * (hi[0] - lo[0]),
                    lo[1] + j * (hi[1] - lo[1]),
                    lo[2] + k * (hi[2] - lo[2]),
                )
    edges = {}

    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"far{len(edges)}"
            og.create_line(corners[key[0]], corners[key[1]], name=name)
            edges[key] = og.create_connector(name, 7, name=f"c{name}")
        return edges[key]

    faces = {}
    for axis in range(3):
        for side in (0, 1):
            u, v = [x for x in range(3) if x != axis]
            loop = []
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                c = [0, 0, 0]
                c[axis], c[u], c[v] = side, a, b
                loop.append(tuple(c))
            name = f"far_{'xyz'[axis]}{('min', 'max')[side]}"
            connectors = [edge(loop[m], loop[(m + 1) % 4]) for m in range(4)]
            if name == "far_zmin":
                og.create_unstructured_domain(
                    connectors=[*connectors, "c_root_upper", "c_root_lower"], name=name
                )
            else:
                og.create_unstructured_domain(connectors=connectors, name=name)
            faces[name] = name
    walls = ["wing_upper", "wing_lower", "wing_tip"]
    og.create_volume(
        list(faces) + walls,
        generator={
            "max_size": 0.6,
            "inflation": {"faces": walls, "first_height": 2e-3, "growth_rate": 1.2, "layers": 12},
        },
        name="fluid",
    )
    og.create_patch("wing", "wall", walls)
    og.create_patch("symmetry", "symmetry", ["far_zmin"])
    og.create_patch("far", "farfield", [f for f in faces if f != "far_zmin"])
    return og


def check(og) -> dict:
    mesh = og.volume_mesh("fluid")
    for key in ("tets", "prisms", "hexes"):
        assert len(mesh[key]) > 0, key
    filled = total_volume(mesh)
    assert abs(filled - mesh["enclosed_volume"]) < 1e-9 * filled
    assert mesh["quality"]["below_floor"] == 0
    assert og.boundary_report().ok
    return {
        "tets": len(mesh["tets"]),
        "prisms": len(mesh["prisms"]),
        "hexes": len(mesh["hexes"]),
        "max_layers": int(mesh["inflation"]["max_layers"]),
        "worst_jacobian": mesh["quality"]["worst_jacobian"],
    }


if __name__ == "__main__":
    print(check(build()))
