"""Phase 10a demo: multi-block structured grids smoothed as one grid (smoothing groups).

  * "cylinder": a cylinder (radius 0.5) in a channel (6 x 3), four structured domains between
    the cylinder and the channel walls, joined along the diagonals. Smoothed one domain at a
    time, the grid lines turn sharply on the diagonals (the interfaces stay where the
    connectors put them). As a smoothing group, the diagonals move, the lines cross them
    smoothly, the cylinder keeps its first-cell height (wall control) and the channel walls
    slide (the grid meets them at right angles).
  * "pipe": a circle as an O-grid (butterfly): a core and four domains around it, with four
    irregular points (three domains around a node). As a group, the core corners and the
    nodes on the circle move; the result does not depend on the initial core size.
Writes 10a_multiblock.png (per-domain vs group smoothing, coloured by scaled Jacobian),
10a_multiblock_report.txt, projects and journals.

    python examples/10a_multiblock_smoothing.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.model.groups import interface_kinks

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
R = 0.5
WALL_SPACING = 0.01


def cylinder_case(n_arc=33, n_radial=41):
    """Four domains between the cylinder and the channel; returns (session, domains, walls)."""
    og = opengrid.Session()
    corners = [(-3.0, -1.5), (3.0, -1.5), (3.0, 1.5), (-3.0, 1.5)]
    # Circle points towards the box corners (45 degree directions keep the cells balanced).
    angles = [math.atan2(y, x) for x, y in corners]
    for k in range(4):
        a0 = angles[k]
        a1 = angles[(k + 1) % 4]
        if a1 < a0:
            a1 += 2 * math.pi
        og.create_arc((0, 0, 0), R, a0, a1, name=f"arc{k}")
        og.create_connector(f"arc{k}", n_arc, name=f"c-arc{k}")
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"side{k}")
        og.create_connector(f"side{k}", n_arc, name=f"c-side{k}")
        p = (R * math.cos(a0), R * math.sin(a0))
        og.create_line(p, corners[k], name=f"diag{k}")
        og.create_connector(
            f"diag{k}",
            n_radial,
            {"type": "tanh", "begin_spacing": WALL_SPACING, "end_spacing": 0.15},
            name=f"c-diag{k}",
        )
    domains = []
    for k in range(4):
        # j = 0 on the cylinder, counter-clockwise loop: the arc clockwise (reversed), this
        # diagonal outwards, the side, the next diagonal back to the cylinder.
        name = og.create_domain(
            [
                [{"connector": f"c-arc{k}", "reversed": True}],
                [f"c-diag{k}"],
                [f"c-side{k}"],
                [{"connector": f"c-diag{(k + 1) % 4}", "reversed": True}],
            ],
            name=f"cyl{k}",
        )
        domains.append(name)
    return og, domains, [f"c-arc{k}" for k in range(4)]


def pipe_case(fraction=0.55):
    og = opengrid.Session()
    for k in range(4):
        a0 = -3 * math.pi / 4 + k * math.pi / 2
        og.create_arc((0, 0, 0), R, a0, a0 + math.pi / 2, name=f"arc{k}")
        og.create_connector(f"arc{k}", 25, name=f"c-arc{k}")
    og.create_domain([[f"c-arc{k}"] for k in range(4)], name="pipe")
    names = og.ogrid_domain("pipe", fraction=fraction, radial_points=13, radial_spacing=0.01)
    return og, names


def worst_jacobian(og, names) -> float:
    return min(float(og.domain_quality(n)["scaled_jacobian"].min()) for n in names)


def wall_angle(og, names) -> float:
    """Largest deviation from 90 degrees of the first grid line at j = 0 (cylinder domains)."""
    worst = 0.0
    for name in names:
        g = og.domain_grid(name)[..., :2]
        tangent = g[0, 2:] - g[0, :-2]
        normal = g[1, 1:-1] - g[0, 1:-1]
        cosine = np.abs(np.sum(tangent * normal, axis=1)) / (
            np.linalg.norm(tangent, axis=1) * np.linalg.norm(normal, axis=1)
        )
        worst = max(worst, float(np.degrees(np.arcsin(np.clip(cosine, 0, 1))).max()))
    return worst


def first_height(og, names) -> tuple[float, float]:
    heights = np.concatenate(
        [
            np.linalg.norm(og.domain_grid(n)[1, 1:-1, :2] - og.domain_grid(n)[0, 1:-1, :2], axis=1)
            for n in names
        ]
    )
    return float(heights.min()), float(heights.max())


def run_cases(report: list[str]):
    states = {}
    # ---------------------------------------------------------------- cylinder
    og, domains, walls = cylinder_case()
    for name in domains:
        og.set_domain_generator(
            name, elliptic={"control": "laplace", "walls": [{"edge": 0, "spacing": WALL_SPACING}]}
        )
    t = time.perf_counter()
    for name in domains:
        og.smooth_domain(name)
    per_domain_time = time.perf_counter() - t
    per_domain = {og.entity(n).id: og.domain_grid(n).copy() for n in domains}
    per_domain_sj = worst_jacobian(og, domains)
    per_domain_angle = wall_angle(og, domains)
    group = og.create_smoothing_group(
        domains, walls=[{"connector": w, "spacing": WALL_SPACING} for w in walls], name="cylinder"
    )
    t = time.perf_counter()
    result = og.smooth_group(group)
    group_time = time.perf_counter() - t
    params = og.entity(group).params
    kink_before = max(interface_kinks(og.store, params, per_domain, result.singular).values())
    low, high = first_height(og, domains)
    report += [
        "cylinder in a channel: 4 domains, "
        f"{sum(np.prod(og.domain_dimensions(n)) for n in domains)} points",
        "  per-domain smoothing (Laplace + wall control, interfaces fixed): "
        f"{per_domain_time:.2f} s",
        f"    interface kink {kink_before:.1f} deg, min scaled Jacobian {per_domain_sj:.3f}, "
        f"wall angle error {per_domain_angle:.2f} deg",
        f"  smoothing group: {og.group_status(group)}, {group_time:.2f} s",
        f"    interface kink {result.max_kink:.2f} deg, min scaled Jacobian "
        f"{worst_jacobian(og, domains):.3f}, wall angle error {wall_angle(og, domains):.2f} deg",
        f"    first-cell height on the cylinder {low:.5f} .. {high:.5f} (requested {WALL_SPACING})",
        f"    roles: {', '.join(f'{og.store.get(c).name} {r}' for c, r in result.roles.items())}",
    ]
    states["cylinder"] = (og, domains, per_domain, per_domain_sj, kink_before, result)
    og.save(OUTPUT / "10a_cylinder.ogp")
    og.save_journal(OUTPUT / "10a_cylinder_journal.py")
    # -------------------------------------------------------------------- pipe
    og, names = pipe_case()
    for name in names:
        og.set_domain_generator(name, elliptic={"control": "laplace"})
        og.smooth_domain(name)
    per_domain = {og.entity(n).id: og.domain_grid(n).copy() for n in names}
    per_domain_sj = worst_jacobian(og, names)
    group = og.create_smoothing_group(names, name="pipe-group")
    result = og.smooth_group(group)
    params = og.entity(group).params
    kink_before = max(interface_kinks(og.store, params, per_domain, result.singular).values())
    other, other_names = pipe_case(fraction=0.35)
    other_group = other.create_smoothing_group(other_names)
    other.smooth_group(other_group)
    difference = max(
        float(np.abs(og.domain_grid(a) - other.domain_grid(b)).max())
        for a, b in zip(names, other_names, strict=True)
    )
    points = np.concatenate([og.domain_grid(n).reshape(-1, 3) for n in names])
    radius = np.linalg.norm(points[:, :2], axis=1)
    on_circle = np.abs(radius[radius > 0.999 * R] - R).max()
    report += [
        "",
        "pipe (O-grid, core + 4): ",
        f"  per-domain smoothing (Laplace, interfaces fixed): interface kink {kink_before:.1f} "
        f"deg, min scaled Jacobian {per_domain_sj:.3f}",
        f"  smoothing group: {og.group_status(group)}",
        f"    interface kink {result.max_kink:.2f} deg (away from the 4 irregular points), "
        f"min scaled Jacobian {worst_jacobian(og, names):.3f}",
        f"    sliding nodes on the circle within {on_circle:.1e}",
        f"    initial core size 0.55 vs 0.35: grids differ by {difference:.1e} (solver tolerance)",
    ]
    states["pipe"] = (og, names, per_domain, per_domain_sj, kink_before, result)
    og.save(OUTPUT / "10a_pipe.ogp")
    og.save_journal(OUTPUT / "10a_pipe_journal.py")
    return states


def render(states, path: Path) -> None:
    import pyvista as pv

    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1500), border=False)

    def add(og, names, grids, row, column, title):
        plotter.subplot(row, column)
        plotter.set_background("#1f2329")
        plotter.add_text(title, font_size=10, color="white")
        for name in names:
            grid = grids[og.entity(name).id]
            q = opengrid._core.quality_2d(np.ascontiguousarray(grid))["scaled_jacobian"]
            # Transposed: VTK orders points (and cells) with i fastest, as the quality arrays.
            surface = pv.StructuredGrid(grid[..., 0].T, grid[..., 1].T, grid[..., 2].T)
            surface.cell_data["sj"] = q.ravel()
            plotter.add_mesh(
                surface,
                scalars="sj",
                clim=[0, 1],
                cmap="viridis",
                show_edges=True,
                edge_color="#20303f",
                show_scalar_bar=False,
            )
        plotter.view_xy()

    for row, key in enumerate(("cylinder", "pipe")):
        og, names, per_domain, sj, kink, result = states[key]
        title = f"{key}: per domain - kink {kink:.0f} deg, min SJ {sj:.2f}"
        add(og, names, per_domain, row, 0, title)
        after = worst_jacobian(og, names)
        add(
            og,
            names,
            result.grids,
            row,
            1,
            f"{key}: smoothing group - kink {result.max_kink:.1f} deg, min SJ {after:.2f}",
        )
    plotter.subplot(0, 0)
    plotter.camera.zoom(1.6)
    plotter.subplot(0, 1)
    plotter.camera.zoom(1.6)
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
    group = window.session.entity("cylinder").id
    window.set_selection([group])
    window.quality_action.setChecked(True)
    window.run_smoother("cylinder")
    if screenshot is not None:

        def save_and_quit():
            window.smoother_dock.wait()
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
    report: list[str] = []
    states = run_cases(report)
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "10a_multiblock_report.txt").write_text(text)
    render(states, OUTPUT / "10a_multiblock.png")
    if args.gui:
        return run_gui(OUTPUT / "10a_cylinder.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
