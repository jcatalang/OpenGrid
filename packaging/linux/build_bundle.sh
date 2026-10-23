#!/usr/bin/env bash
# Build the self-contained Linux bundle OpenGrid-<version>-linux-x86_64.tar.gz: Python, Qt, VTK,
# Netgen, OpenCASCADE, HDF5 / CGNS and the compiled OpenGrid in one relocatable environment
# (conda-pack), with launchers. The target machine needs nothing else: no conda, no compiler,
# no root (only glibc >= 2.28 and an OpenGL driver for the window, which every desktop has).
#
#   packaging/linux/build_bundle.sh            # -> packaging/dist/linux/ (in Docker)
#   OG_DETACH=1 packaging/linux/build_bundle.sh   # in the background (docker logs -f ...)
#   OG_NATIVE=1 packaging/linux/build_bundle.sh   # directly on this Linux x86_64 machine
#
# The conda package from packaging/conda/recipe.yaml (rattler-build, conda-forge's GCC), a
# fresh runtime environment with it, conda-pack, the launchers; in a Docker container
# (openSUSE Leap 15.6 x86_64; emulated, hence slow, on Apple silicon) or, with OG_NATIVE, on
# the machine itself (a Linux x86_64 with internet: the GitHub workflow
# .github/workflows/linux-bundle.yml does this). Then check it with
# packaging/linux/test_bundle.sh in clean containers.
set -euo pipefail
IMAGE="${OG_IMAGE:-opensuse/leap:15.6}"
PLATFORM="${OG_PLATFORM:-linux/amd64}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT/pyproject.toml")"
OUT="$ROOT/packaging/dist/linux"
mkdir -p "$OUT" "$ROOT/build/linux-check"
cd "$ROOT"
# COPYFILE_DISABLE: macOS tar would add ._name files (extended attributes) to the archive.
export COPYFILE_DISABLE=1
git ls-files -z --cached --others --exclude-standard \
  | while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done \
  | tar --null -T - -czf "$ROOT/build/linux-check/source.tar.gz"
cat > "$ROOT/build/linux-check/bundle_inside.sh" <<'INSIDE'
set -euo pipefail
step() { echo "== $(date +%H:%M:%S) $*"; }
step "$(. /etc/os-release && echo "$PRETTY_NAME"), $(uname -m)"
command -v tar >/dev/null && command -v gzip >/dev/null || {
  zypper --non-interactive --quiet refresh && zypper --non-interactive --quiet install tar gzip curl
}
# Paths: the container's, or (OG_NATIVE) folders of the checkout given by the caller.
CONDA="${OG_CONDA:-/opt/miniforge}"
WORK="${OG_WORK:-/work}"
BUILD="${OG_BUILD:-/build}"
OUTDIR="${OG_OUT:-/out}"
if [ ! -x "$CONDA/bin/conda" ]; then
  step "Miniforge"
  curl -fsSL -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash /tmp/miniforge.sh -b -u -p "$CONDA" >/dev/null
fi
source "$CONDA/etc/profile.d/conda.sh"
if [ ! -x "$CONDA/envs/pack/bin/rattler-build" ]; then
  step "tools: rattler-build, conda-pack"
  conda create -y -q -n pack -c conda-forge --override-channels rattler-build conda-pack
fi
mkdir -p "$WORK" "$BUILD"
if [ -z "${OG_NATIVE:-}" ]; then
  tar -xzf /source.tar.gz -C "$WORK"
fi
rm -rf "$BUILD/channel" "$BUILD/env" "$BUILD/env.tar.gz" "$BUILD"/OpenGrid-*-linux-x86_64
step "conda package (rattler-build; its own build and test environments)"
"$CONDA/envs/pack/bin/rattler-build" build --recipe "$WORK/packaging/conda/recipe.yaml" \
  -m "$WORK/packaging/conda/variants.yaml" -c conda-forge --output-dir "$BUILD/channel" \
  > "$BUILD/rattler.log" 2>&1 || { tail -60 "$BUILD/rattler.log"; exit 1; }
ls "$BUILD/channel/linux-64/"
# The intermediate environments go as soon as they are done (a CI runner has ~14 GB free).
rm -rf "$BUILD/channel/bld" "$BUILD/channel/src_cache"
step "runtime environment"
conda create -y -q -p "$BUILD/env" --override-channels -c "file://$BUILD/channel" -c conda-forge \
  "opengrid=$OG_VERSION" python=3.12
step "licences of the bundled packages (from the package cache: a packed environment has none)"
"$CONDA/envs/pack/bin/python" "$WORK/packaging/collect_licenses.py" "$BUILD/env" "$BUILD/licenses"
step "conda-pack"
"$CONDA/envs/pack/bin/conda-pack" -p "$BUILD/env" -o "$BUILD/env.tar.gz" --quiet
rm -rf "$BUILD/env"
NAME="OpenGrid-$OG_VERSION-linux-x86_64"
B="$BUILD/$NAME"
mkdir -p "$B/env"
tar -xzf "$BUILD/env.tar.gz" -C "$B/env"
rm -f "$BUILD/env.tar.gz"
step "launchers, examples, README"
cp -r "$WORK/packaging/linux/bundle/." "$B/"
cp "$WORK/LICENSE" "$B/LICENSE"
rm -rf "$B/licenses" && mv "$BUILD/licenses" "$B/licenses"
chmod +x "$B/opengrid" "$B/opengrid-python" "$B/install.sh"
cp "$WORK/docs/images/opengrid_icon.png" "$B/opengrid.png"
mkdir -p "$B/examples" "$B/tests"
(cd "$WORK/examples" && find . -name "*.py" -not -path "./output/*" -exec cp --parents {} "$B/examples/" \;)
cp "$WORK/examples/README.md" "$B/examples/"
mkdir -p "$B/examples/output"
cp -r "$WORK/tests/data" "$B/tests/"
mkdir -p "$B/tests/python" && cp "$WORK"/tests/python/*.py "$B/tests/python/"
echo "$OG_VERSION" > "$B/VERSION"
step "the oldest glibc the bundle runs on"
# The newest GLIBC_2.x symbol version any binary needs (conda's binutils, a Python dependency
# OpenGrid never runs, left out); the launcher compares it with the machine's.
"$B/env/bin/python" - "$B/env" > "$B/GLIBC_MIN" <<'SCAN'
import os, re, sys
pattern = re.compile(rb"GLIBC_2\.(\d+)(?:\.(\d+))?")
best, where = (2, 17, 0), "glibc"
for folder, _dirs, files in os.walk(sys.argv[1]):
    for name in files:
        path = os.path.join(folder, name)
        if "conda-linux-gnu-" in name or os.path.islink(path):
            continue
        try:
            with open(path, "rb") as handle:
                if handle.read(4) != b"\x7fELF":
                    continue
                data = handle.read()
        except OSError:
            continue
        for a, b in pattern.findall(data):
            version = (2, int(a), int(b or 0))
            if version > best:
                best, where = version, path
text = f"{best[0]}.{best[1]}"
if best[2]:
    text += f".{best[2]}"
print(text)
print(f"   needs glibc >= {best[0]}.{best[1]} (newest in {where})", file=sys.stderr)
SCAN
step "archive"
tar -C "$BUILD" -czf "$OUTDIR/$NAME.tar.gz" "$NAME"
cp "$BUILD"/channel/linux-64/opengrid-*.conda "$OUTDIR/"
du -sh "$B" "$OUTDIR/$NAME.tar.gz"
step "done"
INSIDE
if [ -n "${OG_NATIVE:-}" ]; then
  if [ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]; then
    echo "OG_NATIVE needs a Linux x86_64 machine (elsewhere, leave it unset: Docker)" >&2
    exit 1
  fi
  OG_VERSION="$VERSION" OG_NATIVE=1 OG_CONDA="${OG_CONDA:-$HOME/.local/share/opengrid/bundle-conda}" \
    OG_WORK="$ROOT" OG_BUILD="$ROOT/build/bundle" OG_OUT="$OUT" \
    bash "$ROOT/build/linux-check/bundle_inside.sh" 2>&1 | tee "$OUT/build.log"
  exit "${PIPESTATUS[0]}"
fi
NAME="opengrid-linux-bundle"
docker rm -f "$NAME" >/dev/null 2>&1 || true
ARGS=(--name "$NAME" --platform "$PLATFORM" -e OG_VERSION="$VERSION"
  -v "opengrid-linux-conda-${PLATFORM//\//-}:/opt/miniforge"
  -v "opengrid-linux-rattler-cache:/root/.cache"
  -v "$ROOT/build/linux-check/source.tar.gz:/source.tar.gz:ro"
  -v "$ROOT/build/linux-check/bundle_inside.sh:/inside.sh:ro"
  -v "$OUT:/out"
  "$IMAGE" bash -c 'bash /inside.sh && echo "OPENGRID BUNDLE: BUILT" || { echo "OPENGRID BUNDLE: FAILED"; exit 1; }')
if [ -n "${OG_DETACH:-}" ]; then
  docker run -d "${ARGS[@]}" >/dev/null
  echo "running in the background: docker logs -f $NAME (output in $OUT)"
else
  docker run --rm "${ARGS[@]}" 2>&1 | tee "$OUT/build.log"
fi
