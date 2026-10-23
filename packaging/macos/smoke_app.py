"""Smoke test of the macOS app, run with its own Python and an empty environment (env -i,
python -I; packaging/macos/test_app.sh): everything it loads must come from the app. Imports
(the standard library's readline and curses included), 2D layers and exports, 3D Netgen volumes,
the window offscreen with Qt's plugins from the app. Prints PASS at the end."""

import os
import sys
import tempfile

out = tempfile.mkdtemp()
print("python", sys.executable, sys.prefix)
import curses  # noqa: E402, F401 (libncurses -> libtinfo: relative links)
import readline  # noqa: E402, F401

import numpy as np  # noqa: E402

import opengrid  # noqa: E402

print("opengrid", opengrid.__version__, "core", opengrid.core_version(), opengrid.__file__)
from importlib.metadata import metadata  # noqa: E402

# The OpenGrid just built, not a cached older package of the same version.
assert metadata("opengrid")["License-Expression"] == "MIT", "a stale OpenGrid package"
assert sys.prefix in opengrid.__file__

import netgen  # noqa: E402, F401

from opengrid import _occt  # noqa: E402, F401

# 2D with sliding layers, exports.
og = opengrid.Session()
og.create_line((0, 0), (4, 0), name="b")
og.create_line((4, 0), (4, 1), name="r")
og.create_line((4, 1), (0, 1), name="t")
og.create_line((0, 1), (0, 0), name="l")
for curve, n, name in (
    ("b", 81, "wall"),
    ("r", 21, "outlet"),
    ("t", 61, "top"),
    ("l", 21, "inlet"),
):
    og.create_connector(curve, n, name=name)
og.create_unstructured_domain(connectors=["wall", "outlet", "top", "inlet"], name="flow")
og.set_inflation(["wall"], first_height=1e-3, growth_rate=1.25, layers=8, slide=True)
assert len(og.unstructured_mesh("flow")["quads"]) == 640
og.set_boundary_type(["wall"], "wall")
og.set_boundary_type(["inlet", "outlet", "top"], "farfield")
og.export_cgns(os.path.join(out, "flow.cgns"))  # (HDF5 / CGNS from the app)
og.export_gmsh(os.path.join(out, "flow.msh"))
print("2D ok", sorted(os.listdir(out)))

# 3D: a sphere in a box, Netgen tetrahedra with prism layers.
og = opengrid.Session()
parts = og.create_sphere(radius=0.5, name="sphere")
og.create_box(center=(0.5, 0.0, 0.0), size=(4.0, 3.0, 3.0), name="box")
body = og.mesh_surfaces([n for n in parts if og.entity(n).kind == "revolved_surface"], spacing=0.12)
outer = og.mesh_surfaces(
    [f"box-{s}" for s in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")], spacing=0.6
)
volume = og.create_volume(outer + body, generator={"max_size": 0.6}, name="fluid")
og.set_inflation(body, first_height=5e-3, layers=3)
mesh = og.volume_mesh(volume)
print("3D ok", {k: len(v) for k, v in mesh.items() if k in ("tetrahedra", "prisms")})
assert len(mesh["prisms"]) > 0

# The window, offscreen; Qt's plugins come from the bundle.
os.environ["QT_QPA_PLATFORM"] = "offscreen"
from PySide6.QtCore import QLibraryInfo  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

plugins = QLibraryInfo.path(QLibraryInfo.LibraryPath.PluginsPath)
print("Qt plugins", plugins)
assert plugins.startswith(sys.prefix) or os.path.realpath(plugins).startswith(
    os.path.realpath(sys.prefix)
), plugins
cocoa = os.path.join(plugins, "platforms", "libqcocoa.dylib")
assert os.path.exists(cocoa), cocoa
app = QApplication.instance() or QApplication([])
from opengrid.app.main_window import MainWindow  # noqa: E402

window = MainWindow(off_screen=True)
window.session.create_line((0, 0), (1, 0))
window.close()
import matplotlib  # noqa: E402, F401
import qtconsole  # noqa: E402, F401

print("window ok")
print("np", np.__version__)
print("PASS")
