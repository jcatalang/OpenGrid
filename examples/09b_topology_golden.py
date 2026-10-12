"""Phase 9b demo: topology helpers for multi-block hexahedral meshes, and the golden cases.

  * "pipe": a circular pipe (radius 0.5, length 3). The cross-section as one structured domain
    has corner cells with a scaled Jacobian of about 0.1; ogrid_domain turns it into an O-H
    (butterfly) layout of five domains with a 1e-3 wall spacing, and each is extruded into a
    block: an all-hexahedral pipe of five blocks with 1-to-1 interfaces.
  * "ball": a solid sphere (radius 0.5) from six structured faces on the sphere. One block has
    inverted cells at the corners; ogrid_block gives a core block and six shell blocks, all
    hexahedral and valid.
  * "duct": a block split at index 20 of its long axis into two blocks and joined back (same
    points).
  * The ten golden cases of SPEC.md (examples/golden) are built from their scripts and checked;
    the results go to output/09b_golden_report.txt.
Exports CGNS (structured multi-block with 1-to-1 connectivity; cgnscheck when available) and
VTKHDF for the pipe and the ball; writes 09b_topology.png, projects and journals.

    python examples/09b_topology_golden.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib
import importlib.util
import itertools
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
GOLDEN = HERE / "golden"
R = 0.5


def circle_domain(og, n=33, name="pipe"):
    for k in range(4):
        a0 = -3 * math.pi / 4 + k * math.pi / 2
        og.create_arc((0, 0, 0), R, a0, a0 + math.pi / 2, name=f"arc{k}")
        og.create_connector(f"arc{k}", n, name=f"c-arc{k}")
    return og.create_domain([[f"c-arc{k}"] for k in range(4)], name=name)


def pipe_case() -> tuple[opengrid.Session, float]:
    og = opengrid.Session()
    circle_domain(og)
    single = float(og.domain_quality("pipe")["scaled_jacobian"].min())
    parts = og.ogrid_domain("pipe", fraction=0.6, radial_points=13, radial_spacing=1e-3)
    axial = {"type": "tanh", "begin_spacing": 0.02, "end_spacing": 0.02}
    blocks = []
    for name in parts:
        blocks.append(
            og.create_extruded_block(
                name,
                "translate",
                vector=(0, 0, 3.0),
                dimension=61,
                distribution=axial,
                name=f"block-{name}",
            )
        )
    og.create_patch("inlet", "inlet", parts)
    og.create_patch("outlet", "outlet", [{"block": b, "face": "kmax"} for b in blocks])
    walls = [{"block": f"block-pipe-{k}", "face": "jmin"} for k in range(4)]
    og.create_patch("wall", "wall", walls)
    return og, single


def ball_case() -> tuple[opengrid.Session, float]:
    og = opengrid.Session()
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
            pa = np.array(corner[key[0]])
            pb = np.array(corner[key[1]])
            mid = R * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, 21, name=f"c-{name}")
        return arcs[key]

    faces = []
    for axis, surface in ((0, "sphere_x"), (1, "sphere_y"), (2, "sphere_z")):
        for side in (-1, 1):
            a, b = [k for k in range(3) if k != axis]
            cs = []
            for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                c = [0, 0, 0]
                c[axis], c[a], c[b] = side, i, j
                cs.append(tuple(c))
            name = f"sphere-{'xyz'[axis]}{('min', 'max')[side > 0]}"
            og.create_structured_domain(
                connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                surface=surface,
                name=name,
            )
            faces.append(name)
    og.create_block(faces=faces, name="ball")
    single = float(og.block_quality("ball")["scaled_jacobian"].min())
    og.ogrid_block("ball", fraction=0.55, radial_points=11)
    og.create_patch("surface", "wall", faces)
    return og, single


def duct_case() -> opengrid.Session:
    spec = importlib.util.spec_from_file_location("g07", GOLDEN / "g07_duct.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    og = module.build()
    og.delete(["extruded"])
    before = og.block_grid("duct").copy()
    sizes = before.shape[2::-1]  # (ni, nj, nk)
    axis = "ijk"[int(np.argmax(sizes))]
    low, high, middle = og.split_block("duct", axis, 20)
    print(
        f"duct: split at {axis} = 20 into {low} {og.block_grid(low).shape[2::-1]} and {high} "
        f"{og.block_grid(high).shape[2::-1]} (new face {middle}); boundary "
        f"{boundary_state(og)}"
    )
    og.join_blocks(low, high)
    same = np.array_equal(
        np.unique(np.round(before.reshape(-1, 3), 12), axis=0),
        np.unique(np.round(og.block_grid("duct").reshape(-1, 3), 12), axis=0),
    )
    print(f"duct: joined back, same points: {same}")
    return og


def boundary_state(og) -> str:
    if og.boundary_report().ok:
        return "complete"
    return "INCOMPLETE"


def block_summary(og, label, single) -> None:
    worst = 1.0
    cells = 0
    volume = 0.0
    for name in og.names("block"):
        q = og.block_quality(name)
        worst = min(worst, float(q["scaled_jacobian"].min()))
        cells += q["volume"].size
        volume += float(q["volume"].sum())
    report = og.boundary_report()
    print(
        f"{label}: {len(og.names('block'))} blocks, {cells} hexahedra, volume {volume:.6f}; "
        f"min scaled Jacobian {worst:.3f} (one block: {single:.3f}); "
        f"{len(report.interfaces)} interfaces, boundary {boundary_state(og)}"
    )


def export(og, stem) -> None:
    if _core.cgns_available():
        path = OUTPUT / f"{stem}.cgns"
        print(f"  CGNS {og.export_cgns(path)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'{stem}.vtkhdf')}")
    og.save(OUTPUT / f"{stem}.ogp")
    og.save_journal(OUTPUT / f"{stem}_journal.py")
    replayed = opengrid.replay((OUTPUT / f"{stem}_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)


def golden_report(path) -> None:
    lines = ["Golden validation cases (examples/golden), built from their scripts and checked.", ""]
    for script in sorted(GOLDEN.glob("g[0-9][0-9]_*.py")):
        spec = importlib.util.spec_from_file_location(script.stem, script)
        module = importlib.util.module_from_spec(spec)
        sys.modules[script.stem] = module
        spec.loader.exec_module(module)
        start = time.perf_counter()
        og = module.build()
        result = module.check(og)
        elapsed = time.perf_counter() - start
        replayed = to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
        summary = ", ".join(f"{k}={_short(v)}" for k, v in result.items())
        line = f"{script.stem}: ok ({elapsed:.1f} s, journal replays: {replayed}); {summary}"
        lines.append(line)
        print(f"  {line}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _short(value) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def render(pipe, ball, single_pipe, path) -> None:
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)

    def domain_surface(grid, scalars):
        surface = pv.StructuredGrid(grid[..., 0], grid[..., 1], grid[..., 2])
        surface.cell_data["sj"] = scalars.ravel()
        return surface

    def show(index, title):
        plotter.subplot(*index)
        plotter.set_background("#1f2329")
        plotter.add_text(title, font_size=11, color="white")

    # (a) one block: corner cells collapse.
    show((0, 0), f"pipe section, one block: min scaled Jacobian {single_pipe:.2f}")
    og = opengrid.Session()
    circle_domain(og)
    q = og.domain_quality("pipe")["scaled_jacobian"]
    plotter.add_mesh(
        domain_surface(og.domain_grid("pipe"), q),
        scalars="sj",
        clim=[0, 1],
        cmap="viridis",
        show_edges=True,
        edge_color="#20303f",
        show_scalar_bar=False,
    )
    plotter.view_xy()
    # (b) butterfly section (base faces of the pipe blocks).
    worst = min(
        float(pipe.domain_quality(n)["scaled_jacobian"].min()) for n in pipe.names("domain")
    )
    show((0, 1), f"butterfly (O-H): min scaled Jacobian {worst:.2f}")
    for name in pipe.names("domain"):
        q = pipe.domain_quality(name)["scaled_jacobian"]
        plotter.add_mesh(
            domain_surface(pipe.domain_grid(name), q),
            scalars="sj",
            clim=[0, 1],
            cmap="viridis",
            show_edges=True,
            edge_color="#20303f",
            show_scalar_bar=False,
        )
    plotter.view_xy()
    # (c) pipe blocks in 3D.
    show((1, 0), "all-hexahedral pipe: 5 blocks")
    colours = ["#6fa8dc", "#f5a623", "#7fd1b9", "#e8c547", "#d95f5f", "#b07fd1", "#9fc5e8"]
    for k, name in enumerate(pipe.names("block")):
        g = pipe.block_grid(name)
        grid = pv.StructuredGrid(g[..., 0], g[..., 1], g[..., 2])
        plotter.add_mesh(
            grid.extract_surface(), color=colours[k], show_edges=True, edge_color="#20303f"
        )
    plotter.camera_position = [(4.5, -4.5, 4.5), (0, 0, 1.5), (0, 0, 1)]
    # (d) ball: cut z = 0 through the seven blocks.
    show((1, 1), "solid sphere: core + 6 shell blocks, cut z = 0")
    for k, name in enumerate(ball.names("block")):
        g = ball.block_grid(name)
        grid = pv.StructuredGrid(g[..., 0], g[..., 1], g[..., 2])
        cut = grid.slice(normal=(0, 0, 1), origin=(0, 0, 1e-9))
        if cut.n_cells:
            plotter.add_mesh(cut, color=colours[k], show_edges=True, edge_color="#20303f")
    plotter.view_xy()
    plotter.screenshot(str(path))
    plotter.close()


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
    window.set_selection([window.session.entity("ball-core").id])
    window.quality_action.setChecked(True)
    window.viewport.camera_position = [(1.6, -1.4, 1.2), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    if screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(3000, save_and_quit)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    pipe, single_pipe = pipe_case()
    block_summary(pipe, "pipe", single_pipe)
    exact = math.pi * R**2 * 3.0
    print(f"  exact pipe volume {exact:.6f} (the faceted section is slightly smaller)")
    export(pipe, "09b_pipe")
    ball, single_ball = ball_case()
    block_summary(ball, "ball", single_ball)
    export(ball, "09b_ball")
    duct_case()
    print("golden cases:")
    golden_report(OUTPUT / "09b_golden_report.txt")
    render(pipe, ball, single_pipe, OUTPUT / "09b_topology.png")
    if args.gui:
        return run_gui(OUTPUT / "09b_ball.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
