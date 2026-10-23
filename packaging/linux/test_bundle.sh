#!/usr/bin/env bash
# Check the self-contained bundle in a clean container: no conda, no compiler, only what a
# desktop install has (Mesa, fonts) plus a virtual display. Extracts the archive, then: the
# imports (OpenGrid, its C++ core, OpenCASCADE, Netgen), a 2D case, a 3D volume case, the CAD
# import, the application window (a screenshot), a moved folder and install.sh.
#
#   packaging/linux/test_bundle.sh [archive] [image]   # default: the newest archive in
#                                                      # packaging/dist/linux, openSUSE Leap 15.6
#   packaging/linux/test_bundle.sh "" rockylinux:8      # glibc 2.28, the bundle's minimum
#   packaging/linux/test_bundle.sh "" registry.suse.com/suse/sle15:15.3   # SLES 15 SP3
# Images whose repositories have no desktop packages (SLES base images) get the headless
# checks only: imports, meshing, exports.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/packaging/dist/linux"
ARCHIVE="${1:-}"
if [ -z "$ARCHIVE" ]; then
  ARCHIVE="$(ls -t "$OUT"/OpenGrid-*-linux-x86_64.tar.gz | head -1)"
fi
IMAGE="${2:-opensuse/leap:15.6}"
PLATFORM="${OG_PLATFORM:-linux/amd64}"
cat > "$OUT/test_inside.sh" <<'INSIDE'
set -euo pipefail
step() { echo "== $(date +%H:%M:%S) $*"; }
step "$(. /etc/os-release && echo "$PRETTY_NAME"), $(uname -m), glibc $(ldd --version | head -1 | awk '{print $NF}')"
# (Meshes depend slightly on the CPU: Netgen amplifies rounding differences of CPU-specific
# kernels. The machine is recorded with every result.)
echo "   CPU:$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2), $(nproc) cores"
step "desktop packages only (Mesa, fonts, X11 helpers) and a virtual display"
if command -v zypper >/dev/null; then
  zypper --non-interactive --quiet refresh || echo "   (repositories not reachable)"
  if ! command -v tar >/dev/null || ! command -v gzip >/dev/null; then
    zypper --non-interactive --quiet install --no-recommends tar gzip which
  fi
  zypper --non-interactive --quiet install --no-recommends xauth Mesa-libGL1 Mesa-libEGL1 \
    Mesa-dri libxkbcommon-x11-0 fontconfig dejavu-fonts xvfb-run \
    || echo "   no desktop packages in this image's repositories: headless checks only"
elif command -v dnf >/dev/null; then
  dnf install -y -q tar gzip which mesa-libGL mesa-libEGL mesa-dri-drivers libxkbcommon-x11 \
    fontconfig dejavu-sans-fonts xorg-x11-server-Xvfb xorg-x11-xauth >/dev/null
else
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq && apt-get install -y -qq libgl1 libegl1 libxkbcommon-x11-0 fontconfig \
    fonts-dejavu-core xvfb xauth >/dev/null
fi
for tool in conda gcc g++ c++ cmake; do
  if command -v "$tool" >/dev/null; then echo "   unexpected: $tool is installed"; exit 1; fi
done
echo "   no conda, no compiler, no cmake on this machine"
step "extract"
mkdir -p /home/user && cd /home/user
tar -xzf /bundle.tar.gz
B="$(ls -d /home/user/OpenGrid-*-linux-x86_64)"
cd "$B"
step "licences (OpenGrid's and every bundled package's)"
grep -q "MIT License" LICENSE
grep -q "^qt6-main " licenses/README.txt && test -d licenses/qt6-main && test -d licenses/netgen
echo "   $(find licenses -mindepth 1 -maxdepth 1 -type d | wc -l) packages' licence texts"
step "first start and the imports (the OpenGrid just built: MIT-licensed, not a cached older build)"
./opengrid-python -c "from importlib.metadata import metadata; assert metadata('opengrid')['License-Expression'] == 'MIT', 'a stale OpenGrid package'"
./opengrid-python -c "import opengrid, opengrid._core as c, opengrid._occt, netgen.meshing; print('   OpenGrid', opengrid.__version__, '- core', opengrid.core_version(), '- CGNS', c.cgns_available(), '- OpenMP', c.openmp_enabled(), '- threads', opengrid.threads())"
step "headless: a 2D case extruded, CGNS and Gmsh (no display, no OpenGL)"
./opengrid-python - <<'PY'
import opengrid
og = opengrid.Session()
box = og.create_rectangle(width=2.0, height=1.0, name="box")
body = og.create_ellipse(center=(0.0, 0.0), a=0.3, b=0.15, name="body")
og.unstructured_domain_from_curves([*box, body], 31, name="flow")
walls = og.names("connector")
og.set_boundary_type(walls[4:], "wall")
og.set_boundary_type(walls[:4], "farfield")
og.extrude_model(0.05)
print("   ", og.export_cgns("/tmp/headless.cgns"), og.export_gmsh("/tmp/headless.msh"))
PY
if ! command -v xvfb-run >/dev/null; then
  step "done (headless checks only)"
  exit 0
fi
./opengrid-python -c "import PySide6, vtk; print('   Qt and VTK load')"
step "2D: airfoil cases (shapes, inflation, periodic pair, extrusion, CGNS / Gmsh)"
./opengrid-python examples/12a_airfoil_cases.py | tail -3
step "3D: bodies in a box (surface meshes, Netgen volumes with prism layers)"
./opengrid-python examples/12b_bodies_in_a_box.py | tail -3
step "CAD import (OpenCASCADE: STEP, IGES)"
./opengrid-python examples/09a_cad_import.py | tail -3
step "the application window (virtual display, screenshot)"
xvfb-run -a -s "-screen 0 1600x1000x24" ./opengrid-python examples/00_launch.py --screenshot "/out/bundle_window_$OG_TAG.png"
step "a project opened in the application (drawn, coloured by quality, screenshot)"
cat > /tmp/open_project.py <<'PY'
import os
import sys
from PySide6.QtWidgets import QApplication
app = QApplication(sys.argv)
from opengrid.app.main_window import MainWindow
window = MainWindow()
window.resize(1600, 1000)
window.show()
assert window.open_project("examples/output/12a_far_field.ogp"), "the project did not open"
window.set_selection([window.session.entity("flow").id])  # the quality panel follows
for _ in range(10):  # each pass may redraw (software OpenGL is slow under emulation)
    app.processEvents()
window.grab().save(f"/out/bundle_project_{os.environ['OG_TAG']}.png")
print("   project opened:", len(window.session.store), "entities; screenshot saved")
PY
xvfb-run -a -s "-screen 0 1600x1000x24" ./opengrid-python /tmp/open_project.py
step "a moved folder"
mv "$B" /home/user/moved && cd /home/user/moved
./opengrid-python -c "import opengrid, opengrid._core; print('   imports fine after the move')"
mv /home/user/moved "$B" && cd "$B"
step "install.sh (menu entry, ~/.local/bin links)"
HOME=/home/user ./install.sh
HOME=/home/user /home/user/.local/bin/opengrid-python -c "import opengrid; print('   ~/.local/bin/opengrid-python works')"
ls /home/user/.local/share/applications/opengrid.desktop
step "done"
INSIDE
TAG="$(echo "$IMAGE" | tr '/:' '__')"
docker rm -f "opengrid-bundle-test-$TAG" >/dev/null 2>&1 || true
docker run --rm --name "opengrid-bundle-test-$TAG" --platform "$PLATFORM" -e OG_TAG="$TAG" \
  -v "$ARCHIVE:/bundle.tar.gz:ro" -v "$OUT/test_inside.sh:/inside.sh:ro" -v "$OUT:/out" \
  "$IMAGE" bash -c 'bash /inside.sh && echo "OPENGRID BUNDLE TEST: PASS" || { echo "OPENGRID BUNDLE TEST: FAIL"; exit 1; }' \
  2>&1 | tee "$OUT/test-$TAG.log"
