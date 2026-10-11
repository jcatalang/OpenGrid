"""Phase 4 demo: patches, zones and solver-ready export.

Two cases (geometry from the Phase 3 demos):
  * backward-facing step (3 blocks, 2 interfaces): patches walls / inflow / outflow,
    zones "channel" and "expansion";
  * NACA 0012 C-grid (elliptic, wall spacing 1e-4 c; the wake is a cut within the block):
    patches airfoil (wall), farfield, outflow.
Each case is validated (every boundary connector in exactly one patch, interfaces detected
automatically) and exported to CGNS (checked with cgnscheck when available), Gmsh 4.1 (read
back with meshio) and VTKHDF (read back with VTK). Prints patch face counts per format.

    python examples/04_solver_ready.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

import opengrid
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def step_case() -> opengrid.Session:
    og = opengrid.Session()
    _load("03_backward_step").build(og)
    og.create_patch(
        "walls", "wall", ["c-in-bottom", "c-in-top", "c-step", "c-lo-bottom", "c-up-top"]
    )
    og.create_patch("inflow", "inlet", ["c-inlet"])
    og.create_patch("outflow", "outlet", ["c-out-lo", "c-out-up"])
    og.create_zone("channel", ["inlet-block"])
    og.create_zone("expansion", ["lower-block", "upper-block"])
    return og


def airfoil_case() -> opengrid.Session:
    og = opengrid.Session()
    _load("03_airfoil_grids").c_grid(og)
    og.create_patch("airfoil", "wall", ["c-c-upper", "c-c-lower"])
    og.create_patch("farfield", "farfield", ["c-c-top", "c-c-front", "c-c-bottom"])
    og.create_patch("outflow", "outlet", ["c-c-outflow-top", "c-c-outflow-bottom"])
    og.create_zone("fluid", ["c-grid"])
    return og


def export(og: opengrid.Session, name: str) -> None:
    report = og.boundary_report()
    print(
        f"{name}: {len(report.boundary)} boundary connectors, "
        f"interfaces: {', '.join(report.interfaces)}"
    )
    assert report.ok, report.issues
    expected = {}
    for patch in og.store.entities("patch"):
        expected[patch.name] = sum(
            og.store.get(m["connector"]).params["dimension"] - 1 for m in patch.params["members"]
        )

    cgns = OUTPUT / f"04_{name}.cgns"
    print(f"  CGNS   {og.export_cgns(cgns)}")
    checker = shutil.which("cgnscheck")
    if checker is not None:
        out = subprocess.run(
            [checker, str(cgns)], capture_output=True, text=True, check=False
        ).stdout
        status = "clean"
        if "warning" in out.lower() or "error" in out.lower():
            status = "PROBLEMS:\n" + out
        print(f"         cgnscheck: {status}")

    msh = OUTPUT / f"04_{name}.msh"
    print(f"  Gmsh   {og.export_gmsh(msh)}")
    import meshio

    mesh = meshio.read(msh)
    lines = {}
    for block, tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"], strict=True):
        if block.type == "line":
            for tag in np.unique(tags):
                group = next(n for n, (t, d) in mesh.field_data.items() if t == tag and d == 1)
                lines[group] = lines.get(group, 0) + int(np.count_nonzero(tags == tag))
    print(f"         meshio physical groups: {sorted(mesh.field_data)}; boundary faces {lines}")
    assert lines == expected

    vtkhdf = OUTPUT / f"04_{name}.vtkhdf"
    summary = og.export_vtkhdf(vtkhdf)
    import vtk

    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(vtkhdf))
    reader.Update()
    pdc = reader.GetOutputDataObject(0)
    counts = {}
    for i in range(pdc.GetNumberOfPartitionedDataSets()):
        block = pdc.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
        counts[block] = pdc.GetPartitionedDataSet(i).GetPartition(0).GetNumberOfCells()
    print(f"  VTKHDF {summary['blocks']} read back by VTK: {counts}")
    assert {k: v for k, v in counts.items() if k in expected} == expected

    og.save(OUTPUT / f"04_{name}.ogp")
    og.save_journal(OUTPUT / f"04_{name}_journal.py")
    assert to_json(
        opengrid.replay((OUTPUT / f"04_{name}_journal.py").read_text()).store
    ) == to_json(og.store)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    export(step_case(), "backward_step")
    export(airfoil_case(), "airfoil_cgrid")
    if args.gui:
        return run_gui(OUTPUT / "04_backward_step.ogp", args.screenshot)
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
    window.patch_color_action.setChecked(True)
    window.validate_boundary()
    window.boundary_dock_widget.raise_()
    window.set_selection([window.session.entity("walls").id])
    camera = window.viewport.camera  # close-up of the step corner and the interfaces
    camera.focal_point = (0.5, 1.0, 0.0)
    camera.position = (0.5, 1.0, 10.0)
    camera.parallel_scale = 1.3
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
