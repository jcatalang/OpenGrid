"""Phase 5a demo: unstructured triangle meshes around a NACA 0012.

Case 1 - "triangles": the airfoil (LE/TE-clustered connectors) inside a circular farfield
         (R = 20 c), triangulated with growth rate 1.15 and minimum angle 28 degrees.
Case 2 - "hybrid": a structured O-ring around the airfoil out to an ellipse (exact rational
         NURBS), elliptically smoothed with first-cell height 1e-4 c, and triangles from the
         ellipse to the farfield. The ellipse connector is shared, so the mesh is conforming.

Exports Gmsh 4.1 (read back with meshio and, if available, Gmsh) and VTKHDF, prints quality
statistics, writes 05a_unstructured.png and the projects/journals.

    python examples/05a_unstructured.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import sys
import time
from collections import Counter
from pathlib import Path

import opengrid
from opengrid.project import to_json

OUTPUT = Path(__file__).resolve().parent / "output"
N_SURFACE = 101


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def airfoil(og):
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    og.create_spline([(x, naca0012(x)) for x in xs], "clamped", start_tangent=(0, 1), name="upper")
    og.create_spline(
        [(x, -naca0012(x)) for x in xs], "clamped", start_tangent=(0, -1), name="lower"
    )
    surface = {"type": "tanh", "begin_spacing": 2e-3, "end_spacing": 1e-3}
    og.create_connector("upper", N_SURFACE, surface, name="c-upper")
    og.create_connector("lower", N_SURFACE, surface, name="c-lower")


def triangles_case() -> opengrid.Session:
    og = opengrid.Session()
    airfoil(og)
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_unstructured_domain(
        connectors=["c-upper", "c-lower", "c-farfield"],
        generator={"growth_rate": 1.15, "min_angle": 28.0},
        name="triangles",
    )
    og.create_patch("airfoil", "wall", ["c-upper", "c-lower"])
    og.create_patch("far", "farfield", ["c-farfield"])
    return og


def hybrid_case() -> opengrid.Session:
    og = opengrid.Session()
    airfoil(og)
    # Ellipse x = 0.5 + 0.75 cos t, y = 0.3 sin t as the affine image of the NURBS circle.
    w = math.sqrt(2.0) / 2.0
    unit = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)]
    control = [(0.5 + 0.75 * x, 0.3 * y) for x, y in unit]
    og.create_nurbs(
        2,
        [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1],
        control,
        [1, w, 1, w, 1, w, 1, w, 1],
        name="ellipse",
    )
    og.create_line((1.0, 0.0), (1.25, 0.0), name="cut")
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_connector("ellipse", 2 * N_SURFACE - 1, name="c-ellipse")
    og.create_connector(
        "cut", 41, {"type": "tanh", "begin_spacing": 1e-4, "end_spacing": 0.02}, name="c-cut"
    )
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_domain(
        [[("c-lower", True), "c-upper"], ["c-cut"], ["c-ellipse"], [("c-cut", True)]],
        generator={
            "method": "hermite_tfi",
            "elliptic": {
                "control": "thomas_middlecoff",
                "walls": [{"edge": 0, "spacing": 1e-4}],
                "max_iterations": 20000,
            },
        },
        name="ring",
    )
    og.create_unstructured_domain(
        connectors=["c-ellipse", "c-farfield"],
        generator={"growth_rate": 1.15, "min_angle": 28.0},
        name="outside",
    )
    og.create_patch("airfoil", "wall", ["c-upper", "c-lower"])
    og.create_patch("far", "farfield", ["c-farfield"])
    og.create_zone("fluid", ["ring", "outside"])
    return og


def report(og: opengrid.Session, label: str) -> None:
    for name in og.names("domain") + og.names("unstructured_domain"):
        start = time.perf_counter()
        if og.entity(name).kind == "domain":
            result = og.smooth_domain(name)
            print(f"  {name}: {og.domain_status(name)} ({time.perf_counter() - start:.1f} s)")
            if result.warnings:
                print(f"    warnings: {result.warnings}")
        else:
            og.unstructured_mesh(name)
            print(f"  {name}: {og.domain_status(name)} ({time.perf_counter() - start:.2f} s)")
        q = og.domain_quality(name)
        print(
            f"    min angle {q['min_angle'].min():.2f} deg, "
            f"max skewness {q['skewness'].max():.3f}, "
            f"min scaled Jacobian {q['scaled_jacobian'].min():.3f}, "
            f"max area ratio {q['area_ratio'].max():.2f}"
        )
    result = og.boundary_report()
    assert result.ok, result.issues
    print(f"  interfaces: {', '.join(result.interfaces) or 'none'}")
    msh = OUTPUT / f"05a_{label}.msh"
    summary = og.export_gmsh(msh)
    import meshio

    mesh = meshio.read(msh)
    types = Counter()
    for block in mesh.cells:
        types[block.type] += len(block.data)
    print(f"  Gmsh {summary}; meshio cells {dict(types)}")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests" / "python"))
    from gmsh_tool import find_gmsh, gmsh_read

    if find_gmsh():
        log = gmsh_read(msh)
        status = "clean"
        if "Error" in log or "Warning" in log:
            status = "PROBLEMS:\n" + log
        print(f"  Gmsh executable read: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'05a_{label}.vtkhdf')}")
    og.save(OUTPUT / f"05a_{label}.ogp")
    og.save_journal(OUTPUT / f"05a_{label}_journal.py")
    assert to_json(
        opengrid.replay((OUTPUT / f"05a_{label}_journal.py").read_text()).store
    ) == to_json(og.store)


def render(cases, path):
    import pyvista as pv

    from opengrid.app.scene import quad_mesh, triangle_mesh

    views = [
        (cases[0], (0.5, 0.0), 1.2, "triangles"),
        (cases[0], (0.0, 0.0), 0.06, "triangles: LE"),
        (cases[1], (0.5, 0.0), 1.2, "hybrid"),
        (cases[1], (1.1, 0.0), 0.25, "hybrid: TE / interface"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(1800, 1200), border=False)
    for k, (og, center, scale, title) in enumerate(views):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        for name in og.names("domain"):
            plotter.add_mesh(
                quad_mesh(og.domain_grid(name)), style="wireframe", color="#ffb86b", line_width=1
            )
        for name in og.names("unstructured_domain"):
            m = og.unstructured_mesh(name)
            plotter.add_mesh(
                triangle_mesh(m["points"], m["triangles"]),
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
    cases = []
    for label, build in (("triangles", triangles_case), ("hybrid", hybrid_case)):
        print(f"{label}:")
        og = build()
        report(og, label)
        cases.append(og)
    render(cases, OUTPUT / "05a_unstructured.png")
    if args.gui:
        return run_gui(OUTPUT / "05a_hybrid.ogp", args.screenshot)
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
    window.session.smooth_domain("ring")
    window.scene.rebuild()
    window.set_selection([window.session.entity("outside").id])
    window.quality_dock.metric.setCurrentIndex(3)  # min angle
    window.quality_action.setChecked(True)
    camera = window.viewport.camera
    camera.focal_point = (0.6, 0.0, 0.0)
    camera.position = (0.6, 0.0, 10.0)
    camera.parallel_scale = 0.7
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
