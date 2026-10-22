"""Inflation layers set from the curves, tapering or running into the neighbouring curves.

A channel 4 long and 1 high: a wall at the bottom, an outlet at x = 4 and an inlet leaning at
60 degrees to the wall. ``set_inflation`` gives the wall 10 layers (first height 1e-3, growth
1.25):
  * as before, the layers taper towards the inlet and the outlet in steps (the default);
  * with ``slide=True`` they run into the inlet and the outlet along their own points, which
    set_inflation places at the layer heights near the wall (on the leaning inlet 1 / sin 60
    farther apart along it, so the height across the layers is the layer height);
  * a sharper wedge (30 degrees) is left tapering: layers squeezed into it would make thin,
    distorted cells (layers slide where the wall and its neighbour meet between 45 and 135
    degrees).

Writes 12c_sliding_layers.png (the inlet corner of each case), the projects and journals.

    python examples/12c_sliding_layers.py
"""

import math
import sys
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
LAYERS = {"first_height": 1e-3, "growth_rate": 1.25, "layers": 10, "transition": False}


def channel(inlet_deg: float) -> opengrid.Session:
    """The channel with its curves, points and the unstructured domain "flow"."""
    og = opengrid.Session()
    lean = 1.0 / math.tan(math.radians(inlet_deg))
    og.create_line((0, 0), (4, 0), name="bottom")
    og.create_line((4, 0), (4, 1), name="right")
    og.create_line((4, 1), (lean, 1), name="ceiling")
    og.create_line((lean, 1), (0, 0), name="left")
    for curve, points, name in (
        ("bottom", 81, "wall"),
        ("right", 21, "outlet"),
        ("ceiling", 61, "top"),
        ("left", 21, "inlet"),
    ):
        og.create_connector(curve, points, name=name)
    og.create_unstructured_domain(connectors=["wall", "outlet", "top", "inlet"], name="flow")
    og.set_boundary_type(["wall"], "wall")
    og.set_boundary_type(["inlet"], "inlet")
    og.set_boundary_type(["outlet"], "outlet")
    og.set_boundary_type(["top"], "symmetry")
    return og


def summary(og: opengrid.Session, label: str) -> dict:
    mesh = og.unstructured_mesh("flow")
    quality = og.domain_quality("flow")
    worst = float(np.min(quality["scaled_jacobian"]))
    print(
        f"{label}: {len(mesh['quads'])} layer quads, {len(mesh['triangles'])} triangles, "
        f"worst scaled Jacobian {worst:.3f}, worst angle {np.min(quality['min_angle']):.1f} deg"
    )
    return mesh


def polydata(mesh):
    import pyvista as pv

    faces = []
    for cells in (mesh["triangles"], mesh["quads"]):
        cells = np.asarray(cells, dtype=np.int64)
        if len(cells):
            faces.append(np.hstack([np.full((len(cells), 1), cells.shape[1]), cells]).ravel())
    return pv.PolyData(np.asarray(mesh["points"], dtype=float), faces=np.concatenate(faces))


def render(cases, path: Path) -> None:
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2400, 900), border=False)
    for k, (mesh, title, centre, half) in enumerate(cases):
        plotter.subplot(0, k)
        plotter.set_background("#1f2329")
        plotter.add_mesh(polydata(mesh), color="#e8c547", show_edges=True, edge_color="#30363d")
        plotter.add_text(title, font_size=11, color="white")
        plotter.camera_position = [(*centre, 10.0), (*centre, 0.0), (0, 1, 0)]
        plotter.enable_parallel_projection()
        plotter.camera.parallel_scale = half
    plotter.screenshot(str(path))
    plotter.close()
    from PIL import Image

    image = Image.open(path).convert("RGB").quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    image.save(path, optimize=True)


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    tapered = channel(60.0)
    tapered.set_inflation(["wall"], **LAYERS)
    taper_mesh = summary(tapered, "60 degrees, tapering")
    slid = channel(60.0)
    slid.set_inflation(["wall"], slide=True, **LAYERS)
    slide_mesh = summary(slid, "60 degrees, sliding")
    assert slid.last_warnings == []
    assert len(slide_mesh["quads"]) == 80 * LAYERS["layers"]  # every layer, end to end
    wedge = channel(30.0)
    wedge.set_inflation(["wall"], slide=True, **LAYERS)
    wedge_mesh = summary(wedge, "30 degrees, slide asked: tapering")
    print("  ", wedge.last_warnings[0])
    assert slid.boundary_report().ok
    for og, name in ((slid, "12c_sliding_layers"), (wedge, "12c_wedge_30deg")):
        og.save(OUTPUT / f"{name}.ogp")
        og.save_journal(OUTPUT / f"{name}_journal.py")
    path = OUTPUT / "12c_sliding_layers.png"
    render(
        [
            (taper_mesh, "inlet at 60 deg: layers taper (the default)", (0.09, 0.04), 0.06),
            (slide_mesh, "slide=True: the layers run into the inlet", (0.045, 0.03), 0.045),
            (wedge_mesh, "inlet at 30 deg: too sharp a wedge, they taper", (0.2, 0.06), 0.08),
        ],
        path,
    )
    print(f"wrote {path.relative_to(HERE.parent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
