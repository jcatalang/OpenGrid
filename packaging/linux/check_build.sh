#!/usr/bin/env bash
# Build OpenGrid and run its tests in a clean Linux container, the way a user of that system
# would (docs/INSTALL_LINUX.md, by hand): Miniforge, the conda-forge environment of
# environment.yml plus conda-forge's GCC (cxx-compiler: the distribution's own compiler does not
# matter), the C++ build with warnings as errors, ctest, the Python package and pytest under a
# virtual X display.
#
#   packaging/linux/check_build.sh                      # openSUSE Leap 15.6 on x86_64
#   packaging/linux/check_build.sh ubuntu:24.04         # another distribution
#   OG_JOBS=6 packaging/linux/check_build.sh            # parallel compile jobs (default 2)
#   OG_PYTEST="-k periodic" packaging/linux/check_build.sh   # a subset of the Python tests
#   OG_DETACH=1 packaging/linux/check_build.sh          # run in the background (docker logs -f)
#
# Needs Docker. On Apple silicon the x86_64 image runs emulated (correct, but several times
# slower). The source is the working tree (tracked and new files, not build/ or .git);
# Miniforge and its package cache persist in a Docker volume, so later runs skip the downloads.
# The log is build/linux-check/<image>.log; its last line says PASS or FAIL.
set -euo pipefail
IMAGE="${1:-opensuse/leap:15.6}"
PLATFORM="${OG_PLATFORM:-linux/amd64}"
JOBS="${OG_JOBS:-2}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/build/linux-check"
TAG="$(echo "$IMAGE" | tr '/:' '__')"
mkdir -p "$OUT"
cd "$ROOT"
# COPYFILE_DISABLE: macOS tar would add ._name files (extended attributes) to the archive.
export COPYFILE_DISABLE=1
git ls-files -z --cached --others --exclude-standard \
  | while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done \
  | tar --null -T - -czf "$OUT/source.tar.gz"
cat > "$OUT/inside.sh" <<'INSIDE'
set -euo pipefail
step() { echo "== $(date +%H:%M:%S) $*"; }
. /etc/os-release
step "$PRETTY_NAME, $(uname -m), $(nproc) CPUs"
step "system packages (OpenGL, virtual display, fonts)"
if command -v zypper >/dev/null; then
  zypper --non-interactive --quiet refresh
  zypper --non-interactive --quiet install --no-recommends curl tar gzip which fontconfig
  for p in xvfb-run Mesa-libGL1 Mesa-libEGL1 Mesa-dri libxkbcommon-x11-0 dejavu-fonts; do
    zypper --non-interactive --quiet install --no-recommends "$p" || echo "   (no package $p)"
  done
elif command -v apt-get >/dev/null; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq
  apt-get install -y -qq curl ca-certificates xvfb xauth libgl1 libegl1 libxkbcommon-x11-0 \
    fontconfig fonts-dejavu-core >/dev/null
else
  echo "no zypper or apt-get: install curl, Xvfb and Mesa yourself"; exit 1
fi
CONDA=/opt/miniforge
if [ ! -x "$CONDA/bin/conda" ]; then
  step "Miniforge"
  curl -fsSL -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash /tmp/miniforge.sh -b -u -p "$CONDA" >/dev/null
fi
source "$CONDA/etc/profile.d/conda.sh"
mkdir -p /work && tar -xzf /source.tar.gz -C /work && cd /work
step "conda environment from environment.yml, plus cxx-compiler"
conda env remove -y -q -n og >/dev/null 2>&1 || true
conda env create -q -n og -f environment.yml
conda install -y -q -n og -c conda-forge --override-channels cxx-compiler
conda activate og
# conda-forge's GCC under its prefixed name: it finds its own archiver and linker (docs/INSTALL_LINUX.md).
export CXX="$CONDA_PREFIX/bin/$(uname -m)-conda-linux-gnu-c++"
step "compiler: $("$CXX" --version | head -1); cmake $(cmake --version | head -1 | awk '{print $3}')"
step "C++ build (warnings as errors) and ctest"
cmake --preset ci > /tmp/configure.log || { cat /tmp/configure.log; exit 1; }
grep -E "OpenGrid:|OpenCASCADE|CGNS" /tmp/configure.log || true
cmake --build --preset ci -j "$OG_JOBS"
ctest --preset ci -j "$OG_JOBS"
step "Python package (reusing the C++ build)"
CMAKE_BUILD_PARALLEL_LEVEL="$OG_JOBS" pip install --no-build-isolation -q -C build-dir=build/ci .
python -c "import opengrid, opengrid._core; print('opengrid', opengrid.__version__, 'core', opengrid.core_version())"
python -c "import opengrid._occt" && echo "   STEP / IGES import: yes" || echo "   STEP / IGES import: no"
step "pytest (virtual display)"
xvfb-run -a -s "-screen 0 1600x1200x24" python -m pytest -q -p no:cacheprovider ${OG_PYTEST:-}
step "done"
INSIDE
NAME="opengrid-linux-check-$TAG"
docker rm -f "$NAME" >/dev/null 2>&1 || true
ARGS=(--name "$NAME" --platform "$PLATFORM"
  -e OG_JOBS="$JOBS" -e OG_PYTEST="${OG_PYTEST:-}"
  -v "opengrid-linux-conda-${PLATFORM//\//-}:/opt/miniforge"
  -v "$OUT/source.tar.gz:/source.tar.gz:ro" -v "$OUT/inside.sh:/inside.sh:ro"
  "$IMAGE" bash -c 'bash /inside.sh && echo "OPENGRID LINUX CHECK: PASS" || { echo "OPENGRID LINUX CHECK: FAIL"; exit 1; }')
if [ -n "${OG_DETACH:-}" ]; then
  docker run -d "${ARGS[@]}" >/dev/null
  echo "running in the background: docker logs -f $NAME"
else
  docker run --rm "${ARGS[@]}" 2>&1 | tee "$OUT/$TAG.log"
fi
