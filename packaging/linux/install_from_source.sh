#!/usr/bin/env bash
# Build and install OpenGrid from this checkout on Linux, with everything it needs, for this user
# (no root, nothing from the system but glibc and OpenGL):
#   1. Miniforge into a private folder ($OPENGRID_HOME, default ~/.local/share/opengrid/miniforge),
#      unless it is there already; no system or shared conda is used or changed;
#   2. the environment "opengrid" from environment.yml plus conda-forge's GCC (cxx-compiler);
#   3. OpenGrid built and installed in it, editable (after a git pull, run this script again);
#   4. launchers ~/.local/bin/opengrid and opengrid-python, and a menu entry.
#
#   packaging/linux/install_from_source.sh          # install or update
#   packaging/linux/install_from_source.sh --test   # ... and run the C++ and Python tests
#
# Needs internet for the downloads (about 1.5 GB; 3 GB on disk) and 10-20 minutes.
set -euo pipefail
ROOT="$(cd "$(dirname "$(readlink -f "$0")")/../.." && pwd)"
PREFIX="${OPENGRID_HOME:-$HOME/.local/share/opengrid/miniforge}"
ENV_NAME=opengrid
JOBS="${OG_JOBS:-$(nproc)}"
step() { echo "== $*"; }

if [ "$(uname -s)" != "Linux" ]; then
  echo "install_from_source.sh: Linux only (macOS: see README.md)" >&2
  exit 1
fi
if [ ! -x "$PREFIX/bin/conda" ]; then
  step "Miniforge into $PREFIX"
  mkdir -p "$(dirname "$PREFIX")"
  installer="$(mktemp -d)/miniforge.sh"
  curl -fsSL -o "$installer" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash "$installer" -b -p "$PREFIX" >/dev/null
  rm -f "$installer"
fi
# Only this prefix's conda, whatever the user's shell has.
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL PYTHONHOME PYTHONPATH
source "$PREFIX/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  step "updating the environment $ENV_NAME from environment.yml"
  conda env update -q -n "$ENV_NAME" -f "$ROOT/environment.yml"
else
  step "creating the environment $ENV_NAME from environment.yml"
  conda env create -q -n "$ENV_NAME" -f "$ROOT/environment.yml"
fi
conda install -y -q -n "$ENV_NAME" -c conda-forge --override-channels cxx-compiler
conda activate "$ENV_NAME"
# The prefixed compiler finds its own archiver and linker (no binutils needed on the system).
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++"
if [ "$(uname -m)" = "aarch64" ]; then
  export CXX="$CONDA_PREFIX/bin/aarch64-conda-linux-gnu-c++"
fi
step "building OpenGrid ($("$CXX" --version | head -1))"
cd "$ROOT"
CMAKE_BUILD_PARALLEL_LEVEL="$JOBS" pip install --no-build-isolation -q -e .
python -c "import opengrid; print('   OpenGrid', opengrid.__version__, '- core', opengrid.core_version())"

step "launchers"
BIN="$HOME/.local/bin"
mkdir -p "$BIN"
for name in opengrid opengrid-python; do
  run="python"
  if [ "$name" = "opengrid" ]; then
    run="python -m opengrid"
  fi
  cat > "$BIN/$name" <<LAUNCHER
#!/bin/bash
# OpenGrid from $ROOT, environment $ENV_NAME in $PREFIX (install_from_source.sh)
unset PYTHONHOME PYTHONPATH QT_PLUGIN_PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL
source "$PREFIX/etc/profile.d/conda.sh"
conda activate $ENV_NAME
exec $run "\$@"
LAUNCHER
  chmod +x "$BIN/$name"
done
APPS="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
ICONS="${XDG_DATA_HOME:-$HOME/.local/share}/icons/hicolor/256x256/apps"
mkdir -p "$APPS" "$ICONS"
cp "$ROOT/docs/images/opengrid_icon.png" "$ICONS/opengrid.png"
cat > "$APPS/opengrid.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=OpenGrid
Comment=Structured, unstructured and hybrid mesh generator
Exec=$BIN/opengrid %F
Icon=opengrid
Terminal=false
Categories=Science;Engineering;Education;
DESKTOP
echo "   $BIN/opengrid, $BIN/opengrid-python and the OpenGrid menu entry"
case ":$PATH:" in *":$BIN:"*) ;; *) echo "   note: add $BIN to your PATH to type 'opengrid'";; esac

if [ "${1:-}" = "--test" ]; then
  step "C++ tests"
  cmake --preset release >/dev/null
  cmake --build --preset release -j "$JOBS" >/dev/null
  ctest --preset release -j "$JOBS" | tail -3
  step "Python tests"
  runner=()
  if [ -z "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ] && command -v xvfb-run >/dev/null; then
    runner=(xvfb-run -a)
  fi
  "${runner[@]}" python -m pytest -q tests/python
fi
step "done: run 'opengrid' (or $BIN/opengrid)"
