"""Golden case 5: backward-facing step, three structured blocks (inlet channel, upper and lower
blocks behind the step) with 1-to-1 interfaces and wall clustering.

Checks: total area 45, the two block interfaces are detected, all boundaries are in patches,
all cells valid.
"""

import opengrid

WALL = 2e-3


def build() -> opengrid.Session:
    og = opengrid.Session()
    lines = {
        "in-bottom": ((-5, 1), (0, 1)),
        "in-top": ((-5, 2), (0, 2)),
        "inlet": ((-5, 1), (-5, 2)),
        "interface": ((0, 1), (0, 2)),
        "step": ((0, 0), (0, 1)),
        "lo-bottom": ((0, 0), (20, 0)),
        "mid": ((0, 1), (20, 1)),
        "up-top": ((0, 2), (20, 2)),
        "out-lo": ((20, 0), (20, 1)),
        "out-up": ((20, 1), (20, 2)),
    }
    for name, (a, b) in lines.items():
        og.create_line(a, b, name=name)
    walls = {"type": "tanh", "begin_spacing": WALL, "end_spacing": WALL}
    toward = {"type": "tanh", "begin_spacing": 0.2, "end_spacing": 0.01}
    behind = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.5}
    dims = {"in-bottom": (31, toward), "in-top": (31, toward), "lo-bottom": (101, behind)}
    dims.update({"mid": (101, behind), "up-top": (101, behind)})
    for name in lines:
        n, dist = dims.get(name, (41, walls))
        og.create_connector(name, n, dist, name=f"c-{name}")
    og.create_domain(
        [["c-in-bottom"], ["c-interface"], [("c-in-top", True)], [("c-inlet", True)]],
        name="inlet-block",
    )
    og.create_domain(
        [["c-lo-bottom"], ["c-out-lo"], [("c-mid", True)], [("c-step", True)]],
        name="lower-block",
    )
    og.create_domain(
        [["c-mid"], ["c-out-up"], [("c-up-top", True)], [("c-interface", True)]],
        name="upper-block",
    )
    og.create_patch("inflow", "inlet", ["c-inlet"])
    og.create_patch("outflow", "outlet", ["c-out-lo", "c-out-up"])
    og.create_patch(
        "walls", "wall", ["c-in-bottom", "c-in-top", "c-step", "c-lo-bottom", "c-up-top"]
    )
    return og


def check(og) -> dict:
    area = 0.0
    for name in ("inlet-block", "lower-block", "upper-block"):
        q = og.domain_quality(name)
        assert q["scaled_jacobian"].min() > 0
        area += float(q["area"].sum())
    assert abs(area - 45.0) < 1e-10
    report = og.boundary_report()
    assert report.ok, report.issues
    names = sorted(text.split(" ")[0] for text in report.interfaces)
    assert names == ["c-interface", "c-mid"]
    return {"area": area, "interfaces": sorted(report.interfaces)}


if __name__ == "__main__":
    print(check(build()))
