"""Phase 6 demo: boundary-layer inflation (hybrid quad layers + unstructured core).

Case 1 - "naca_hybrid": NACA 0012 (sharp trailing edge) in a circular farfield (R = 20 c).
         Layers from the airfoil with first height 1e-5 c and growth 1.2, stopping where the
         cells become isotropic; the core is triangles (advancing front) and, as a variant,
         all-quad.
Case 2 - "two_element": main element plus a 30 % chord flap deflected 20 degrees with a small
         gap; the layers of both elements stop short of each other in the gap (proximity).

Exports CGNS (cgnscheck when available), Gmsh 4.1 (read back with meshio and Gmsh) and VTKHDF,
writes 06_inflation.png and the projects/journals.

    python examples/06_inflation.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
CORE = {"algorithm": "advancing_front", "growth_rate": 1.15, "smoothing": "angle_based"}


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def element(og, prefix, n, scale=1.0, dx=0.0, dy=0.0, angle=0.0, spacing=(2e-3, 1e-3)):
    """NACA 0012 surfaces (rotated nose-down by ``angle`` degrees about the leading edge)."""
    c = math.cos(math.radians(angle))
    s = math.sin(math.radians(angle))

    def place(x, y):
        return (dx + scale * (c * x + s * y), dy + scale * (-s * x + c * y))

    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    tangent = place(0, 1)
    up = (tangent[0] - dx, tangent[1] - dy)
    og.create_spline(
        [place(x, naca0012(x)) for x in xs], "clamped", start_tangent=up, name=f"{prefix}upper"
    )
    og.create_spline(
        [place(x, -naca0012(x)) for x in xs],
        "clamped",
        start_tangent=(-up[0], -up[1]),
        name=f"{prefix}lower",
    )
    distribution = {
        "type": "tanh",
        "begin_spacing": spacing[0] * scale,
        "end_spacing": spacing[1] * scale,
    }
    og.create_connector(f"{prefix}upper", n, distribution, name=f"c-{prefix}upper")
    og.create_connector(f"{prefix}lower", n, distribution, name=f"c-{prefix}lower")
    return [f"c-{prefix}upper", f"c-{prefix}lower"]


def naca_case(element_type: str) -> opengrid.Session:
    og = opengrid.Session()
    walls = element(og, "", 101)
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_unstructured_domain(
        connectors=[*walls, "c-farfield"],
        generator={
            **CORE,
            "element_type": element_type,
            "inflation": {
                "connectors": walls,
                "first_height": 1e-5,
                "growth_rate": 1.2,
                "layers": 60,
            },
        },
        name="flow",
    )
    og.create_patch("airfoil", "wall", walls)
    og.create_patch("far", "farfield", ["c-farfield"])
    return og


def two_element_case() -> opengrid.Session:
    og = opengrid.Session()
    main = element(og, "main-", 101)
    flap = element(og, "flap-", 61, scale=0.3, dx=0.97, dy=-0.045, angle=20.0)
    og.create_circle((0.6, 0.0), 20.0, name="farfield")
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_unstructured_domain(
        connectors=[*main, *flap, "c-farfield"],
        generator={
            **CORE,
            "inflation": {
                "connectors": main + flap,
                "first_height": 1e-5,
                "growth_rate": 1.2,
                "layers": 60,
            },
        },
        name="flow",
    )
    og.create_patch("main", "wall", main)
    og.create_patch("flap", "wall", flap)
    og.create_patch("far", "farfield", ["c-farfield"])
    return og


def summary(og: opengrid.Session, label: str) -> None:
    start = time.perf_counter()
    mesh = og.unstructured_mesh("flow")
    seconds = time.perf_counter() - start
    q = og.domain_quality("flow")
    print(
        f"  {label}: {og.domain_status('flow')} ({seconds:.2f} s); scaled Jacobian min "
        f"{q['scaled_jacobian'].min():.3f} mean {q['scaled_jacobian'].mean():.3f}; "
        f"warnings {mesh['warnings']}"
    )


def export(og: opengrid.Session, label: str) -> None:
    result = og.boundary_report()
    assert result.ok, result.issues
    if _core.cgns_available():
        cgns = OUTPUT / f"06_{label}.cgns"
        print(f"    CGNS {og.export_cgns(cgns)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run(
                [checker, str(cgns)], capture_output=True, text=True, check=False
            ).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"    cgnscheck: {status}")
    msh = OUTPUT / f"06_{label}.msh"
    info = og.export_gmsh(msh)
    import meshio

    mesh = meshio.read(msh)
    print(f"    Gmsh {info}; meshio cells { {b.type: len(b.data) for b in mesh.cells} }")
    sys.path.insert(0, str(HERE.parent / "tests" / "python"))
    from gmsh_tool import find_gmsh, gmsh_read

    if find_gmsh():
        log = gmsh_read(msh)
        status = "clean"
        if "Error" in log or "Warning" in log:
            status = "PROBLEMS:\n" + log
        print(f"    Gmsh executable read: {status}")
    print(f"    VTKHDF {og.export_vtkhdf(OUTPUT / f'06_{label}.vtkhdf')}")
    og.save(OUTPUT / f"06_{label}.ogp")
    og.save_journal(OUTPUT / f"06_{label}_journal.py")
    assert to_json(
        opengrid.replay((OUTPUT / f"06_{label}_journal.py").read_text()).store
    ) == to_json(og.store)


def first_layer_check(og: opengrid.Session, walls) -> float:
    """Largest relative error of the first-layer height (distance to the wall polygon)."""
    mesh = og.unstructured_mesh("flow")
    segments = []
    for name in walls:
        nodes = og.connector_nodes(name).points[:, :2]
        segments.append((nodes[:-1], nodes[1:]))
    a = np.concatenate([s[0] for s in segments])
    b = np.concatenate([s[1] for s in segments])
    n_boundary = sum(len(og.connector_nodes(c).points) - 1 for c in [*walls, "c-farfield"])
    candidates = mesh["points"][n_boundary:, :2]
    ab = b - a
    distance = np.full(len(candidates), np.inf)
    for k in range(len(a)):
        t = np.clip(((candidates - a[k]) @ ab[k]) / (ab[k] @ ab[k]), 0.0, 1.0)
        d = np.linalg.norm(candidates - (a[k] + t[:, None] * ab[k]), axis=1)
        distance = np.minimum(distance, d)
    first = distance[distance < 1.5e-5]
    return float(np.max(np.abs(first / 1e-5 - 1.0))), len(first)


def render(cases, path):
    import pyvista as pv

    from opengrid.app.scene import triangle_mesh

    views = [
        (cases["naca"], (0.5, 0.0), 0.7, "NACA 0012: layers + advancing-front core"),
        (cases["naca"], (1.0, 0.0), 0.02, "trailing-edge fan"),
        (cases["naca"], (0.0, 0.0), 0.004, "leading edge (first cell 1e-5 c)"),
        (cases["naca_quads"], (0.5, 0.0), 0.7, "all-quad core"),
        (cases["two"], (0.95, -0.05), 0.3, "two-element airfoil"),
        (cases["two"], (0.985, -0.03), 0.03, "gap: proximity termination"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(2400, 1300), border=False)
    for k, (og, center, scale, title) in enumerate(views):
        plotter.subplot(k // 3, k % 3)
        plotter.set_background("#1f2329")
        m = og.unstructured_mesh("flow")
        plotter.add_mesh(
            triangle_mesh(m["points"], m["triangles"], m["quads"]),
            style="wireframe",
            color="#9fb3c8",
            line_width=1,
        )
        plotter.add_text(title, font_size=11, color="white")
        plotter.enable_parallel_projection()
        plotter.camera.focal_point = (center[0], center[1], 0.0)
        plotter.camera.position = (center[0], center[1], 10.0)
        plotter.camera.up = (0, 1, 0)
        plotter.camera.parallel_scale = scale
    plotter.screenshot(str(path))
    plotter.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    cases = {}
    print("naca_hybrid (first height 1e-5 c, growth 1.2):")
    cases["naca"] = naca_case("triangles")
    summary(cases["naca"], "triangle core")
    error, count = first_layer_check(cases["naca"], ["c-upper", "c-lower"])
    print(f"    first layer: {count} nodes, max relative height error {error:.2e}")
    export(cases["naca"], "naca_hybrid")
    cases["naca_quads"] = naca_case("quads")
    summary(cases["naca_quads"], "all-quad core")
    export(cases["naca_quads"], "naca_quads")
    print("two_element:")
    cases["two"] = two_element_case()
    summary(cases["two"], "main + flap")
    export(cases["two"], "two_element")
    render(cases, OUTPUT / "06_inflation.png")
    if args.gui:
        return run_gui(OUTPUT / "06_two_element.ogp", args.screenshot)
    return 0


def run_gui(project, screenshot):
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow()
    window.open_project(str(project))
    window.show()
    window.set_selection([window.session.entity("flow").id])
    window.quality_dock.metric.setCurrentIndex(0)  # scaled Jacobian
    window.quality_action.setChecked(True)
    camera = window.viewport.camera
    camera.focal_point = (0.97, -0.04, 0.0)
    camera.position = (0.97, -0.04, 10.0)
    camera.parallel_scale = 0.12
    if screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(2000, save_and_quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
