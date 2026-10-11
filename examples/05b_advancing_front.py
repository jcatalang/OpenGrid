"""Phase 5b demo: advancing-front triangulation, compared with the constrained Delaunay mesher.

Case 1 - "af_triangles": the NACA 0012 of the Phase 5a demo (LE/TE-clustered connectors) in a
         circular farfield (R = 20 c), growth rate 1.15, meshed with both generators
         (algorithm "delaunay" and "advancing_front"); quality statistics side by side.
Case 2 - "af_hybrid": the structured O-ring of the Phase 5a demo with the outer region meshed by
         the advancing front; the shared ellipse connector keeps the mesh conforming.

Exports Gmsh 4.1 (read back with meshio and, if available, Gmsh) and VTKHDF, writes
05b_advancing_front.png and the projects/journals.

    python examples/05b_advancing_front.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def phase5a():
    spec = importlib.util.spec_from_file_location("demo5a", HERE / "05a_unstructured.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def triangles_case(demo5a, algorithm: str) -> opengrid.Session:
    og = opengrid.Session()
    demo5a.airfoil(og)
    og.create_circle((0.5, 0.0), 20.0, name="farfield")
    og.create_connector("farfield", 121, name="c-farfield")
    og.create_unstructured_domain(
        connectors=["c-upper", "c-lower", "c-farfield"],
        generator={"algorithm": algorithm, "growth_rate": 1.15, "min_angle": 28.0},
        name="triangles",
    )
    og.create_patch("airfoil", "wall", ["c-upper", "c-lower"])
    og.create_patch("far", "farfield", ["c-farfield"])
    return og


def hybrid_case(demo5a) -> opengrid.Session:
    og = demo5a.hybrid_case()
    og.set_unstructured_generator("outside", algorithm="advancing_front")
    return og


def statistics(og: opengrid.Session, name: str) -> dict:
    start = time.perf_counter()
    mesh = og.unstructured_mesh(name)
    seconds = time.perf_counter() - start
    angles = og.domain_quality(name)["min_angle"].ravel()
    return {
        "triangles": len(mesh["triangles"]),
        "points": len(mesh["points"]),
        "seconds": seconds,
        "min": float(angles.min()),
        "mean": float(angles.mean()),
        "below30": float(100.0 * np.mean(angles < 30.0)),
        "angles": angles,
        "warnings": list(mesh["warnings"]),
    }


def export(og: opengrid.Session, label: str) -> None:
    result = og.boundary_report()
    assert result.ok, result.issues
    print(f"  interfaces: {', '.join(result.interfaces) or 'none'}")
    msh = OUTPUT / f"05b_{label}.msh"
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
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'05b_{label}.vtkhdf')}")
    og.save(OUTPUT / f"05b_{label}.ogp")
    og.save_journal(OUTPUT / f"05b_{label}_journal.py")
    assert to_json(
        opengrid.replay((OUTPUT / f"05b_{label}_journal.py").read_text()).store
    ) == to_json(og.store)


def render(delaunay, front, hybrid, stats, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyvista as pv

    from opengrid.app.scene import quad_mesh, triangle_mesh

    views = [
        (delaunay, (0.5, 0.0), 1.2, "Delaunay"),
        (front, (0.5, 0.0), 1.2, "advancing front"),
        (hybrid, (0.6, 0.0), 0.9, "hybrid: O-ring + advancing front"),
        (delaunay, (0.0, 0.0), 0.06, "Delaunay: LE"),
        (front, (0.0, 0.0), 0.06, "advancing front: LE"),
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
    image = plotter.screenshot(return_img=True)
    plotter.close()

    figure = plt.figure(figsize=(24, 13), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1))
    axes.imshow(image)
    axes.axis("off")
    histogram = figure.add_axes((0.70, 0.07, 0.27, 0.36))
    bins = np.linspace(20, 60, 41)
    for label, color in (("delaunay", "#4c8bf5"), ("advancing_front", "#f5a623")):
        s = stats[label]
        histogram.hist(
            s["angles"],
            bins=bins,
            alpha=0.6,
            color=color,
            label=f"{label}: {s['triangles']} tris, min {s['min']:.1f}, mean {s['mean']:.1f}",
        )
    histogram.set_xlabel("min included angle [deg]")
    histogram.set_ylabel("triangles")
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

    print("af_triangles (NACA 0012, R = 20 c, growth 1.15):")
    sessions = {}
    stats = {}
    for algorithm in ("delaunay", "advancing_front"):
        og = triangles_case(demo5a, algorithm)
        sessions[algorithm] = og
        s = statistics(og, "triangles")
        stats[algorithm] = s
        print(
            f"  {algorithm:16s} {s['triangles']:6d} tris {s['points']:6d} pts "
            f"min angle {s['min']:.2f}  mean {s['mean']:.2f}  <30 deg {s['below30']:.2f}%  "
            f"{s['seconds']:.2f} s  warnings {s['warnings']}"
        )
    export(sessions["advancing_front"], "af_triangles")

    print("af_hybrid:")
    hybrid = hybrid_case(demo5a)
    start = time.perf_counter()
    hybrid.smooth_domain("ring")
    print(f"  ring: {hybrid.domain_status('ring')} ({time.perf_counter() - start:.1f} s)")
    s = statistics(hybrid, "outside")
    print(
        f"  outside (advancing front): {s['triangles']} tris, min angle {s['min']:.2f}, "
        f"mean {s['mean']:.2f}, warnings {s['warnings']}"
    )
    export(hybrid, "af_hybrid")

    render(
        sessions["delaunay"],
        sessions["advancing_front"],
        hybrid,
        stats,
        OUTPUT / "05b_advancing_front.png",
    )
    if args.gui:
        return run_gui(OUTPUT / "05b_af_triangles.ogp", args.screenshot)
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
    window.set_selection([window.session.entity("triangles").id])
    window.quality_dock.metric.setCurrentIndex(3)  # min angle
    window.quality_action.setChecked(True)
    camera = window.viewport.camera
    camera.focal_point = (0.5, 0.0, 0.0)
    camera.position = (0.5, 0.0, 10.0)
    camera.parallel_scale = 0.8
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
