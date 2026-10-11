"""Phase 7a demo: structured 3D blocks.

Case 1 - "sector_duct": annular sector duct r in [1, 2], theta in [0, 60] degrees, z in [0, 1],
         made of two blocks (theta 0..30 and 30..60) built from six faces each. The cylinder
         faces lie on extruded surfaces (parametric TFI); the blocks share the face at 30 degrees
         (1-to-1 interface). Volume compared with the analytic value.
Case 2 - "wing_section": the elliptically smoothed O-ring around a NACA 0012 of the Phase 5a
         demo, extruded over a span of 2 c (21 stations clustered towards the root); the wake
         cut becomes a block-internal 1-to-1 interface.

Exports CGNS 3D (checked with cgnscheck when available), VTKHDF (hexahedra + boundary quads),
Gmsh 4.1 and Plot3D 3D; writes 07a_blocks.png and the projects/journals.

    python examples/07a_blocks.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib.util
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
T = math.pi / 3


def polar(r, t, z):
    return (r * math.cos(t), r * math.sin(t), z)


def sector_block(og, t0, t1, n_theta, tag, shared=None):
    """Six-face block of the annular sector theta in [t0, t1]; ``shared`` is the theta = t0
    side (connectors and domain) of the previous block."""
    for z, s in ((0.0, "0"), (1.0, "1")):
        og.create_arc((0, 0, z), 1.0, t0, t1, name=f"inner{s}-{tag}")
        og.create_arc((0, 0, z), 2.0, t0, t1, name=f"outer{s}-{tag}")
        og.create_connector(f"inner{s}-{tag}", n_theta, name=f"c-inner{s}-{tag}")
        og.create_connector(f"outer{s}-{tag}", n_theta, name=f"c-outer{s}-{tag}")
    sides = {}
    radial = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.01}
    axial = {"type": "tanh", "begin_spacing": 0.02, "end_spacing": 0.02}
    for t, which in ((t0, "a"), (t1, "b")):
        if which == "a" and shared is not None:
            sides[which] = shared
            continue
        names = {}
        label = f"{round(math.degrees(t))}"
        for z, s in ((0.0, "0"), (1.0, "1")):
            og.create_line(polar(1, t, z), polar(2, t, z), name=f"radial{s}-{label}")
            names[f"r{s}"] = og.create_connector(
                f"radial{s}-{label}", 21, radial, name=f"c-radial{s}-{label}"
            )
        for r, s in ((1, "i"), (2, "o")):
            og.create_line(polar(r, t, 0), polar(r, t, 1), name=f"axial{s}-{label}")
            names[f"v{s}"] = og.create_connector(
                f"axial{s}-{label}", 21, axial, name=f"c-axial{s}-{label}"
            )
        names["domain"] = og.create_domain(
            [[names["r0"]], [names["vo"]], [(names["r1"], True)], [(names["vi"], True)]],
            name=f"side-{label}",
        )
        sides[which] = names
    a, b = sides["a"], sides["b"]
    og.create_extruded_surface(f"inner0-{tag}", (0, 0, 1), name=f"cyl-inner-{tag}")
    og.create_extruded_surface(f"outer0-{tag}", (0, 0, 1), name=f"cyl-outer-{tag}")
    og.create_domain(
        [[a["r0"]], [f"c-outer0-{tag}"], [(b["r0"], True)], [(f"c-inner0-{tag}", True)]],
        name=f"bottom-{tag}",
    )
    og.create_domain(
        [[a["r1"]], [f"c-outer1-{tag}"], [(b["r1"], True)], [(f"c-inner1-{tag}", True)]],
        name=f"top-{tag}",
    )
    og.create_domain(
        [[f"c-inner0-{tag}"], [b["vi"]], [(f"c-inner1-{tag}", True)], [(a["vi"], True)]],
        name=f"wall-inner-{tag}",
        surface=f"cyl-inner-{tag}",
    )
    og.create_domain(
        [[f"c-outer0-{tag}"], [b["vo"]], [(f"c-outer1-{tag}", True)], [(a["vo"], True)]],
        name=f"wall-outer-{tag}",
        surface=f"cyl-outer-{tag}",
    )
    og.create_block(
        faces=[
            f"bottom-{tag}",
            f"top-{tag}",
            f"wall-inner-{tag}",
            f"wall-outer-{tag}",
            a["domain"],
            b["domain"],
        ],
        name=f"duct-{tag}",
    )
    return b


def sector_duct() -> opengrid.Session:
    og = opengrid.Session()
    side = sector_block(og, 0.0, T / 2, 41, "a")
    sector_block(og, T / 2, T, 41, "b", shared=side)
    og.create_patch(
        "walls", "wall", ["wall-inner-a", "wall-outer-a", "wall-inner-b", "wall-outer-b"]
    )
    og.create_patch("inlet", "inlet", ["bottom-a", "bottom-b"])
    og.create_patch("outlet", "outlet", ["top-a", "top-b"])
    og.create_patch("sides", "symmetry", ["side-0", "side-60"])
    og.create_zone("duct", ["duct-a", "duct-b"])
    return og


def wing_section() -> opengrid.Session:
    spec = importlib.util.spec_from_file_location("demo5a", HERE / "05a_unstructured.py")
    demo5a = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo5a)
    og = demo5a.hybrid_case()
    og.delete(["outside", "far", "fluid"])
    og.delete(["airfoil"])
    og.smooth_domain("ring")
    og.create_extruded_block(
        "ring",
        "translate",
        vector=(0.0, 0.0, 2.0),
        dimension=21,
        distribution={"type": "tanh", "begin_spacing": 0.02},
        name="wing",
    )
    og.create_patch("airfoil", "wall", [{"block": "wing", "face": "jmin"}])
    og.create_patch("outer", "farfield", [{"block": "wing", "face": "jmax"}])
    og.create_patch("root", "symmetry", ["ring"])
    og.create_patch("tip", "symmetry", [{"block": "wing", "face": "kmax"}])
    return og


def report(og: opengrid.Session, label: str) -> None:
    for name in og.names("block"):
        start = time.perf_counter()
        status = og.block_status(name)
        q = og.block_quality(name)
        print(
            f"  {name}: {status} ({time.perf_counter() - start:.2f} s); scaled Jacobian min "
            f"{q['scaled_jacobian'].min():.3f}, max skewness {q['skewness'].max():.3f}, "
            f"max aspect ratio {q['aspect_ratio'].max():.1f}, volume {q['volume'].sum():.6f}"
        )
    result = og.boundary_report()
    assert result.ok, result.issues
    print(f"  interfaces: {', '.join(result.interfaces)}")
    if _core.cgns_available():
        cgns = OUTPUT / f"07a_{label}.cgns"
        print(f"  CGNS {og.export_cgns(cgns)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(cgns)], capture_output=True, text=True).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'07a_{label}.vtkhdf')}")
    msh = OUTPUT / f"07a_{label}.msh"
    print(f"  Gmsh {og.export_gmsh(msh)}")
    sys.path.insert(0, str(HERE.parent / "tests" / "python"))
    from gmsh_tool import find_gmsh, gmsh_read

    if find_gmsh():
        log = gmsh_read(msh)
        status = "clean"
        if "Error" in log or "Warning" in log:
            status = "PROBLEMS:\n" + log
        print(f"  Gmsh executable read: {status}")
    og.export_plot3d(OUTPUT / f"07a_{label}.xyz")
    og.save(OUTPUT / f"07a_{label}.ogp")
    og.save_journal(OUTPUT / f"07a_{label}_journal.py")
    replayed = opengrid.replay((OUTPUT / f"07a_{label}_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)


def structured(grid):
    import pyvista as pv

    return pv.StructuredGrid(
        np.ascontiguousarray(grid[..., 0].transpose(2, 1, 0)),
        np.ascontiguousarray(grid[..., 1].transpose(2, 1, 0)),
        np.ascontiguousarray(grid[..., 2].transpose(2, 1, 0)),
    )


def render(duct, wing, path):
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2400, 900), border=False)
    colours = ("#7fd1b9", "#f5a623")
    # Duct: blocks coloured, shrunk cells cut at the mid plane z = 0.5.
    plotter.subplot(0, 0)
    for name, colour in zip(duct.names("block"), colours, strict=True):
        mesh = structured(duct.block_grid(name))
        plotter.add_mesh(
            mesh.extract_surface(algorithm="dataset_surface"),
            color=colour,
            show_edges=True,
            edge_color="#30363d",
        )
    plotter.add_text("annular sector duct (2 blocks)", font_size=11, color="white")
    plotter.subplot(0, 1)
    for name in duct.names("block"):
        mesh = structured(duct.block_grid(name))
        q = duct.block_quality(name)["scaled_jacobian"].reshape(-1)
        mesh.cell_data["scaled Jacobian"] = q
        cut = mesh.clip(normal=(0, 0, 1), origin=(0, 0, 0.5), crinkle=True).shrink(0.85)
        plotter.add_mesh(
            cut, scalars="scaled Jacobian", clim=(0.9, 1.0), cmap="viridis", show_edges=False
        )
    plotter.add_text("cut z = 0.5, shrunk cells", font_size=11, color="white")
    plotter.subplot(0, 2)
    mesh = structured(wing.block_grid("wing"))
    plotter.add_mesh(
        mesh.extract_surface(algorithm="dataset_surface"),
        color="#9fb3c8",
        show_edges=True,
        edge_color="#30363d",
        opacity=0.35,
    )
    sliced = mesh.clip(normal=(0, 0, 1), origin=(0, 0, 0.05), crinkle=True)
    plotter.add_mesh(sliced, color="#ffb86b", show_edges=True, edge_color="#30363d")
    plotter.add_text("extruded wing section", font_size=11, color="white")
    for k in range(3):
        plotter.subplot(0, k)
        plotter.set_background("#1f2329")
        plotter.camera_position = "iso"
        plotter.reset_camera()
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
    print("sector_duct:")
    duct = sector_duct()
    report(duct, "sector_duct")
    volume = sum(duct.block_quality(n)["volume"].sum() for n in duct.names("block"))
    analytic = 0.5 * T * (4.0 - 1.0) * 1.0
    print(f"  volume {volume:.8f}, analytic {analytic:.8f}, relative {volume / analytic - 1:.2e}")
    print("wing_section:")
    wing = wing_section()
    report(wing, "wing_section")
    render(duct, wing, OUTPUT / "07a_blocks.png")
    if args.gui:
        return run_gui(OUTPUT / "07a_sector_duct.ogp", args.screenshot)
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
    window.set_selection([window.session.entity("duct-a").id])
    window.quality_dock.metric.setCurrentIndex(0)  # scaled Jacobian
    window.quality_action.setChecked(True)
    window.shrink_action.setChecked(True)
    window.cut_action.setChecked(True)
    window.viewport.camera_position = "iso"
    window.viewport.reset_camera()
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
