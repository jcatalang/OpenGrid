"""Phase 10b demo: multi-block hexahedral grids smoothed as one grid (smoothing groups of blocks).

  * "box": a sphere (radius 0.5) in a box (side 3) - six blocks between the six patches of a
    cubed sphere and the six box faces, joined along planar side faces; the eight diagonal
    edges from the sphere's patch corners to the box corners are each shared by three blocks
    (irregular edges). As a group: the side faces and diagonals move, the sphere is a wall
    with a first-cell height of 0.01, the box faces slide on their planes and the box edges
    along their lines. Exported to CGNS (multi-block, 1-to-1; cgnscheck when available).
  * "ball": the solid sphere as a core block and six shells (09b demo). As a group the core
    corners and the edges between shells move and the sphere's nodes slide on it.
  * A 1M-cell ball (core + six shells) timed as a group.
Writes 10b_multiblock_3d.png (cut z = 0, per-block vs group, coloured by scaled Jacobian),
10b_multiblock_3d_report.txt, projects and journals.

    python examples/10b_multiblock_3d.py [--quick] [--gui [--screenshot out.png]]
"""

import argparse
import itertools
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
R = 0.5
BOX = 1.5
WALL_SPACING = 0.01


def sphere_surfaces(og):
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
    return {0: "sphere_x", 1: "sphere_y", 2: "sphere_z"}


def face_corners(axis, side):
    a, b = [k for k in range(3) if k != axis]
    out = []
    for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        c = [0, 0, 0]
        c[axis], c[a], c[b] = side, i, j
        out.append(tuple(c))
    return out


def box_case(n=21, n_radial=25):
    """Sphere in a box: six blocks; returns (session, blocks, sphere faces, box faces)."""
    og = opengrid.Session()
    surfaces = sphere_surfaces(og)
    signs = itertools.product((-1, 1), repeat=3)
    inner = {c: tuple(R * v / math.sqrt(3) for v in c) for c in signs}
    outer = {c: tuple(BOX * v for v in c) for c in inner}
    arcs = {}
    edges = {}
    diagonals = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if key not in arcs:
            pa = np.array(inner[key[0]])
            pb = np.array(inner[key[1]])
            mid = R * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, n, name=f"c-{name}")
        return arcs[key]

    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"edge{len(edges)}"
            og.create_line(outer[key[0]], outer[key[1]], name=name)
            edges[key] = og.create_connector(name, n, name=f"c-{name}")
        return edges[key]

    def diagonal(c):
        if c not in diagonals:
            name = f"diag{len(diagonals)}"
            og.create_line(inner[c], outer[c], name=name)
            diagonals[c] = og.create_connector(
                name,
                n_radial,
                {"type": "tanh", "begin_spacing": WALL_SPACING, "end_spacing": 0.15},
                name=f"c-{name}",
            )
        return diagonals[c]

    sides = {}

    def side_face(a, b):
        key = tuple(sorted((a, b)))
        if key not in sides:
            name = f"side{len(sides)}"
            og.create_structured_domain(
                connectors=[arc(a, b), diagonal(a), edge(a, b), diagonal(b)], name=name
            )
            sides[key] = name
        return sides[key]

    blocks = []
    sphere_faces = []
    box_faces = []
    for axis in range(3):
        for side in (-1, 1):
            cs = face_corners(axis, side)
            label = f"{'xyz'[axis]}{('min', 'max')[side > 0]}"
            sphere = og.create_structured_domain(
                connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                surface=surfaces[axis],
                name=f"sphere-{label}",
            )
            box = og.create_structured_domain(
                connectors=[edge(cs[k], cs[(k + 1) % 4]) for k in range(4)], name=f"box-{label}"
            )
            faces = [sphere, box] + [side_face(cs[k], cs[(k + 1) % 4]) for k in range(4)]
            blocks.append(og.create_block(faces=faces, name=f"block-{label}"))
            sphere_faces.append(sphere)
            box_faces.append(box)
    og.create_patch("sphere", "wall", sphere_faces)
    og.create_patch("farfield", "farfield", box_faces)
    return og, blocks, sphere_faces, box_faces


def ball_case(n=21, radial=11):
    og = opengrid.Session()
    surfaces = sphere_surfaces(og)
    signs = itertools.product((-1, 1), repeat=3)
    corner = {c: tuple(R * v / math.sqrt(3) for v in c) for c in signs}
    arcs = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if key not in arcs:
            pa = np.array(corner[key[0]])
            pb = np.array(corner[key[1]])
            mid = R * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, n, name=f"c-{name}")
        return arcs[key]

    faces = []
    for axis in range(3):
        for side in (-1, 1):
            cs = face_corners(axis, side)
            faces.append(
                og.create_structured_domain(
                    connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                    surface=surfaces[axis],
                    name=f"sphere-{'xyz'[axis]}{('min', 'max')[side > 0]}",
                )
            )
    og.create_block(faces=faces, name="ball")
    names = og.ogrid_block("ball", fraction=0.55, radial_points=radial)
    return og, names


def worst(og, names) -> float:
    return min(float(og.block_quality(n)["scaled_jacobian"].min()) for n in names)


def cells(og, names) -> int:
    return int(sum(np.prod(np.array(og.block_grid(n).shape[:3]) - 1) for n in names))


def first_heights(og, names, faces):
    """First-cell heights above the given faces (domains), away from their edges."""
    from opengrid.model.blocks import FACE_AXIS, FACE_NAMES

    heights = []
    face_ids = {og.entity(f).id for f in faces}
    for name in names:
        result = og.blocks.evaluate(og.entity(name).id)
        grid = result["grid"]
        for face in FACE_NAMES:
            if result["faces"][face] not in face_ids:
                continue
            axis, side = FACE_AXIS[face]
            first, second = 0, 1
            if side > 0:
                first, second = -1, -2
            h = np.linalg.norm(
                np.take(grid, second, axis=2 - axis) - np.take(grid, first, axis=2 - axis), axis=-1
            )
            heights.append(h[2:-2, 2:-2].ravel())
    values = np.concatenate(heights)
    return float(values.min()), float(values.max())


def export_box(og, report: list[str]) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "box.cgns"
        og.export_cgns(path)
        size = path.stat().st_size
        checker = shutil.which("cgnscheck")
        verdict = "cgnscheck not found"
        if checker is not None:
            run = subprocess.run([checker, str(path)], capture_output=True, text=True)
            errors = [line for line in run.stdout.splitlines() if "ERROR" in line.upper()]
            verdict = f"cgnscheck: {len(errors)} error(s)"
        report.append(f"  CGNS export: {size / 1e6:.1f} MB, {verdict}")


def run_cases(report: list[str], quick: bool):
    states = {}
    # ------------------------------------------------------------------- sphere in a box
    og, blocks, sphere_faces, box_faces = box_case()
    per_block = {og.entity(b).id: og.block_grid(b).copy() for b in blocks}
    before = worst(og, blocks)
    group = og.create_smoothing_group(
        blocks=blocks,
        walls=[{"domain": f, "spacing": WALL_SPACING} for f in sphere_faces],
        name="box-group",
    )
    t = time.perf_counter()
    result = og.smooth_group(group)
    elapsed = time.perf_counter() - t
    low, high = first_heights(og, blocks, sphere_faces)
    planes = max(
        float(np.abs(np.abs(og.domain_grid(f)).max(axis=(0, 1)) - BOX).min()) for f in box_faces
    )
    roles = {}
    for d, role in result.roles.items():
        roles.setdefault(role, []).append(og.store.get(d).name)
    report += [
        f"sphere in a box: 6 blocks, {cells(og, blocks)} cells",
        f"  per block (TFI, faces fixed): min scaled Jacobian {before:.3f}",
        f"  smoothing group: {og.group_status(group)}, {elapsed:.2f} s",
        f"    min scaled Jacobian {worst(og, blocks):.3f}, "
        f"interface kink {result.max_kink:.2f} deg",
        f"    first-cell height on the sphere {low:.5f} .. {high:.5f} (requested {WALL_SPACING})",
        f"    box faces stay on their planes within {planes:.1e}",
        "    roles: " + "; ".join(f"{role}: {len(names)}" for role, names in sorted(roles.items())),
    ]
    report.append(
        "    (the kink is inherited from the sphere's own surface grid, fixed as a wall: its grid"
        " lines turn by about 36 deg where the cubed-sphere patches meet, and wall control"
        " copies the wall grid into the first layers)"
    )
    free, free_blocks, _sf, _bf = box_case()
    free_group = free.create_smoothing_group(blocks=free_blocks)
    free_result = free.smooth_group(free_group)
    report.append(
        f"  without wall control (sphere sliding): interface kink {free_result.max_kink:.2f} deg,"
        f" min scaled Jacobian {worst(free, free_blocks):.3f}"
    )
    report_ok = og.boundary_report().ok
    report.append(f"  boundary report ok: {report_ok}")
    if report_ok:
        export_box(og, report)
    states["box"] = (og, blocks, per_block, before, result)
    og.save(OUTPUT / "10b_box.ogp")
    og.save_journal(OUTPUT / "10b_box_journal.py")
    # ------------------------------------------------------------------------- ball
    og, names = ball_case()
    per_block = {og.entity(b).id: og.block_grid(b).copy() for b in names}
    before = worst(og, names)
    group = og.create_smoothing_group(blocks=names, name="ball-group")
    t = time.perf_counter()
    result = og.smooth_group(group)
    elapsed = time.perf_counter() - t
    points = np.concatenate([og.block_grid(n).reshape(-1, 3) for n in names])
    radius = np.linalg.norm(points, axis=1)
    on_sphere = float(np.abs(radius[radius > R * (1 - 1e-6)] - R).max())
    report += [
        "",
        f"ball (core + 6 shells): {cells(og, names)} cells",
        f"  per block (TFI, faces fixed): min scaled Jacobian {before:.3f}",
        f"  smoothing group: {og.group_status(group)}, {elapsed:.2f} s",
        f"    min scaled Jacobian {worst(og, names):.3f}, interface kink {result.max_kink:.2f} deg"
        " (away from irregular edges)",
        f"    sliding nodes on the sphere within {on_sphere:.1e}",
    ]
    states["ball"] = (og, names, per_block, before, result)
    og.save(OUTPUT / "10b_ball.ogp")
    og.save_journal(OUTPUT / "10b_ball_journal.py")
    # ------------------------------------------------------------- 1M-cell timing
    if not quick:
        big, big_names = ball_case(n=65, radial=33)
        group = big.create_smoothing_group(blocks=big_names)
        t = time.perf_counter()
        big.smooth_group(group)
        elapsed = time.perf_counter() - t
        report += [
            "",
            f"1M-cell ball: {cells(big, big_names)} cells, {opengrid.threads()} threads: "
            f"{big.group_status(group)}, {elapsed:.1f} s",
        ]
    return states


def render(states, path: Path) -> None:
    import pyvista as pv

    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 2000), border=False)
    for row, key in enumerate(("box", "ball")):
        og, names, per_block, before, result = states[key]
        for column, (grids, title) in enumerate(
            (
                (per_block, f"{key}: per block - min SJ {before:.2f}"),
                (
                    result.grids,
                    f"{key}: smoothing group - min SJ {worst(og, names):.2f}, "
                    f"kink {result.max_kink:.1f} deg",
                ),
            )
        ):
            plotter.subplot(row, column)
            plotter.set_background("#1f2329")
            plotter.add_text(title, font_size=10, color="white")
            for name in names:
                g = grids[og.entity(name).id]
                q = opengrid._core.quality_3d(np.ascontiguousarray(g))["scaled_jacobian"]
                mesh = pv.StructuredGrid(
                    np.ascontiguousarray(g[..., 0].transpose(2, 1, 0)),
                    np.ascontiguousarray(g[..., 1].transpose(2, 1, 0)),
                    np.ascontiguousarray(g[..., 2].transpose(2, 1, 0)),
                )
                mesh.cell_data["sj"] = q.reshape(-1)
                cut = mesh.clip(normal=(0, 0, 1), origin=(0, 0, 1e-9), crinkle=True)
                if cut.n_cells == 0:
                    continue
                plotter.add_mesh(
                    cut,
                    scalars="sj",
                    clim=[0, 1],
                    cmap="viridis",
                    show_edges=True,
                    edge_color="#20303f",
                    show_scalar_bar=False,
                )
            plotter.camera_position = [(0.0, -0.4, 5.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
            if key == "ball":
                plotter.camera_position = [(0.0, -0.15, 1.6), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
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
    window.set_selection([window.session.entity("box-group").id])
    window.quality_action.setChecked(True)
    window.cut_action.setChecked(True)
    window.run_smoother("box-group")
    window.viewport.camera_position = [(2.4, -4.2, 3.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    if screenshot is not None:

        def save_and_quit():
            window.smoother_dock.wait()
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(5000, save_and_quit)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--quick", action="store_true", help="skip the 1M-cell timing")
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    report: list[str] = []
    states = run_cases(report, args.quick or args.gui)
    text = "\n".join(report) + "\n"
    print(text, end="")
    if not (args.quick or args.gui):
        (OUTPUT / "10b_multiblock_3d_report.txt").write_text(text)
        render(states, OUTPUT / "10b_multiblock_3d.png")
    if args.gui:
        return run_gui(OUTPUT / "10b_box.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
