"""Phase 5c demo: quad-dominant and all-quad unstructured meshes, smoothing, CGNS export.

Case 1 - the NACA 0012 of the Phase 5a demo in a circular farfield (R = 20 c), growth 1.15,
         meshed as triangles, quad-dominant and all-quad (advancing-front triangulation,
         angle-based smoothing); quality statistics side by side.
Case 2 - "quad_hybrid": the structured O-ring of the Phase 5a demo with an all-quad region
         outside it, exported to CGNS (checked with cgnscheck when available), Gmsh 4.1 and
         VTKHDF.

Writes 05c_quads.png and the projects/journals.

    python examples/05c_quads.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib.util
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
GENERATOR = {
    "algorithm": "advancing_front",
    "growth_rate": 1.15,
    "smoothing": "angle_based",
    "smoothing_iterations": 5,
}


def phase5a():
    spec = importlib.util.spec_from_file_location("demo5a", HERE / "05a_unstructured.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def airfoil_case(demo5a, element_type: str) -> opengrid.Session:
    og = opengrid.Session()
    demo5a.airfoil(og)
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_unstructured_domain(
        connectors=["c-upper", "c-lower", "c-farfield"],
        generator={**GENERATOR, "element_type": element_type},
        name="mesh",
    )
    og.create_patch("airfoil", "wall", ["c-upper", "c-lower"])
    og.create_patch("far", "farfield", ["c-farfield"])
    return og


def hybrid_case(demo5a) -> opengrid.Session:
    og = demo5a.hybrid_case()
    og.set_unstructured_generator("outside", **GENERATOR, element_type="quads")
    return og


def statistics(og: opengrid.Session, name: str) -> dict:
    start = time.perf_counter()
    mesh = og.unstructured_mesh(name)
    seconds = time.perf_counter() - start
    q = og.domain_quality(name)
    return {
        "triangles": len(mesh["triangles"]),
        "quads": len(mesh["quads"]),
        "seconds": seconds,
        "jacobian": q["scaled_jacobian"].ravel(),
        "skewness": q["skewness"].ravel(),
        "warnings": list(mesh["warnings"]),
    }


def describe(label, s) -> str:
    return (
        f"  {label:14s} {s['triangles']:6d} tris {s['quads']:6d} quads  "
        f"scaled Jacobian min {s['jacobian'].min():.3f} mean {s['jacobian'].mean():.3f}  "
        f"max skewness {s['skewness'].max():.3f}  {s['seconds']:.2f} s  "
        f"warnings {s['warnings']}"
    )


def export(og: opengrid.Session, label: str) -> None:
    result = og.boundary_report()
    assert result.ok, result.issues
    print(f"  interfaces: {', '.join(result.interfaces) or 'none'}")
    if _core.cgns_available():
        cgns = OUTPUT / f"05c_{label}.cgns"
        print(f"  CGNS {og.export_cgns(cgns)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run(
                [checker, str(cgns)], capture_output=True, text=True, check=False
            ).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    msh = OUTPUT / f"05c_{label}.msh"
    summary = og.export_gmsh(msh)
    import meshio

    mesh = meshio.read(msh)
    print(f"  Gmsh {summary}; meshio cells { {b.type: len(b.data) for b in mesh.cells} }")
    sys.path.insert(0, str(HERE.parent / "tests" / "python"))
    from gmsh_tool import find_gmsh, gmsh_read

    if find_gmsh():
        log = gmsh_read(msh)
        status = "clean"
        if "Error" in log or "Warning" in log:
            status = "PROBLEMS:\n" + log
        print(f"  Gmsh executable read: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'05c_{label}.vtkhdf')}")
    og.save(OUTPUT / f"05c_{label}.ogp")
    og.save_journal(OUTPUT / f"05c_{label}_journal.py")
    assert to_json(
        opengrid.replay((OUTPUT / f"05c_{label}_journal.py").read_text()).store
    ) == to_json(og.store)


def render(cases, hybrid, stats, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyvista as pv

    from opengrid.app.scene import quad_mesh, triangle_mesh

    views = [
        (cases["triangles"], (0.5, 0.0), 1.0, "triangles"),
        (cases["quad_dominant"], (0.5, 0.0), 1.0, "quad-dominant"),
        (cases["quads"], (0.5, 0.0), 1.0, "all-quad"),
        (cases["quads"], (0.0, 0.0), 0.05, "all-quad: LE"),
        (hybrid, (0.6, 0.0), 0.9, "hybrid: O-ring + all-quad"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(2400, 1300), border=False)
    for k, (og, center, scale, title) in enumerate(views):
        plotter.subplot(k // 3, k % 3)
        plotter.set_background("#1f2329")
        for name in og.names("domain"):
            plotter.add_mesh(
                quad_mesh(og.domain_grid(name)), style="wireframe", color="#ffb86b", line_width=1
            )
        for name in og.names("unstructured_domain"):
            m = og.unstructured_mesh(name)
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
    image = plotter.screenshot(return_img=True)
    plotter.close()

    figure = plt.figure(figsize=(24, 13), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1))
    axes.imshow(image)
    axes.axis("off")
    histogram = figure.add_axes((0.70, 0.07, 0.27, 0.36))
    bins = np.linspace(0, 1, 41)
    for label, color in (
        ("triangles", "#4c8bf5"),
        ("quad_dominant", "#f5a623"),
        ("quads", "#7ed321"),
    ):
        s = stats[label]
        histogram.hist(
            s["jacobian"],
            bins=bins,
            alpha=0.55,
            color=color,
            label=f"{label}: {s['triangles']} tris + {s['quads']} quads",
        )
    histogram.set_xlabel("scaled Jacobian")
    histogram.set_ylabel("cells")
    histogram.legend(fontsize=9)
    figure.savefig(path)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    demo5a = phase5a()

    print("NACA 0012, R = 20 c, advancing front + angle-based smoothing:")
    cases = {}
    stats = {}
    for element_type in ("triangles", "quad_dominant", "quads"):
        og = airfoil_case(demo5a, element_type)
        cases[element_type] = og
        stats[element_type] = statistics(og, "mesh")
        print(describe(element_type, stats[element_type]))
    export(cases["quads"], "quads")

    print("quad_hybrid:")
    hybrid = hybrid_case(demo5a)
    start = time.perf_counter()
    hybrid.smooth_domain("ring")
    print(f"  ring: {hybrid.domain_status('ring')} ({time.perf_counter() - start:.1f} s)")
    print(describe("outside", statistics(hybrid, "outside")))
    export(hybrid, "quad_hybrid")

    render(cases, hybrid, stats, OUTPUT / "05c_quads.png")
    if args.gui:
        return run_gui(OUTPUT / "05c_quad_hybrid.ogp", args.screenshot)
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
    window.quality_dock.metric.setCurrentIndex(0)  # scaled Jacobian
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
