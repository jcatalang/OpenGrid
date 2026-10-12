"""Golden case 7: 3D rectangular duct 4 x 1 x 1, one six-face structured block with tanh wall
clustering on the cross-section, and a version extruded from its inlet face.

Checks: volume 4 to 1e-12, all hexahedra perfect (scaled Jacobian 1), both builds give the same
points, patches complete (inlet, outlet, walls).
"""

import itertools

import numpy as np

import opengrid

L, H = 4.0, 1.0


def build() -> opengrid.Session:
    og = opengrid.Session()
    size = (L, H, H)
    corner = {
        c: tuple(c[k] * size[k] for k in range(3)) for c in itertools.product((0, 1), repeat=3)
    }
    walls = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.01}
    connectors = {}
    for c in corner:
        for axis in range(3):
            if c[axis] == 0:
                d = list(c)
                d[axis] = 1
                d = tuple(d)
                name = f"e{''.join(map(str, c))}{'xyz'[axis]}"
                og.create_line(corner[c], corner[d], name=name)
                if axis == 0:
                    connectors[(c, d)] = og.create_connector(name, 41, name=f"c-{name}")
                else:
                    connectors[(c, d)] = og.create_connector(name, 21, walls, name=f"c-{name}")

    def edge(a, b):
        if (a, b) in connectors:
            return [connectors[(a, b)]]
        return [(connectors[(b, a)], True)]

    faces = {}
    for axis in range(3):
        for side in (0, 1):
            u, v = [x for x in range(3) if x != axis]
            loop = []
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                c = [0, 0, 0]
                c[axis], c[u], c[v] = side, a, b
                loop.append(tuple(c))
            name = f"{'xyz'[axis]}{('min', 'max')[side]}"
            og.create_domain([edge(loop[k], loop[(k + 1) % 4]) for k in range(4)], name=name)
            faces[name] = name
    og.create_block(faces=list(faces), name="duct")
    og.create_extruded_block("xmin", "translate", vector=(L, 0, 0), dimension=41, name="extruded")
    og.create_patch("inlet", "inlet", ["xmin"])
    og.create_patch("outlet", "outlet", ["xmax"])
    og.create_patch("walls", "wall", ["ymin", "ymax", "zmin", "zmax"])
    return og


def check(og) -> dict:
    q = og.block_quality("duct")
    volume = float(q["volume"].sum())
    assert abs(volume - L * H * H) < 1e-12
    assert q["scaled_jacobian"].min() > 1 - 1e-12
    a = np.unique(np.round(og.block_grid("duct").reshape(-1, 3), 12), axis=0)
    b = np.unique(np.round(og.block_grid("extruded").reshape(-1, 3), 12), axis=0)
    assert np.array_equal(a, b)
    og.delete(["extruded"])
    assert og.boundary_report().ok
    return {"volume": volume, "cells": int(q["volume"].size)}


if __name__ == "__main__":
    print(check(build()))
