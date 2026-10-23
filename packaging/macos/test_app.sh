#!/usr/bin/env bash
# Check the macOS DMG as a user gets it: the app copied out of the disk image into another
# folder, its signature, its library links (check_app.py: nothing outside the app and macOS),
# and runs with its own Python and an empty environment (nothing from the shell, conda or the
# build folder): smoke_app.py (imports, 2D layers, exports, 3D Netgen volumes, the window
# offscreen) and examples/09a_cad_import.py (STEP and IGES through OpenCASCADE).
#
#   packaging/macos/test_app.sh [OpenGrid-<version>-macos-<arch>.dmg]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT/pyproject.toml")"
DMG="${1:-$ROOT/packaging/dist/OpenGrid-$VERSION-macos-$(uname -m).dmg}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/opengrid-app-test.XXXXXX")"
LSREGISTER=/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister
cleanup() {
  hdiutil detach "$WORK/mnt" -quiet 2>/dev/null || true
  "$LSREGISTER" -u "$WORK/Applications/OpenGrid.app" >/dev/null 2>&1 || true
  rm -rf "$WORK"
}
trap cleanup EXIT
step() { echo "== $(date +%H:%M:%S) $*"; }

step "$(sw_vers -productName) $(sw_vers -productVersion), $(uname -m): $DMG"
mkdir -p "$WORK/mnt" "$WORK/Applications" "$WORK/home" "$WORK/tmp" "$WORK/run/tests"
hdiutil attach -nobrowse -readonly -mountpoint "$WORK/mnt" "$DMG" >/dev/null
ditto "$WORK/mnt/OpenGrid.app" "$WORK/Applications/OpenGrid.app"
hdiutil detach "$WORK/mnt" -quiet
APP="$WORK/Applications/OpenGrid.app"
PY="$APP/Contents/Resources/env/bin/python"
du -sh "$DMG" "$APP"

step "licences (OpenGrid's and every bundled package's)"
grep -q "MIT License" "$APP/Contents/Resources/LICENSE"
grep -q "^qt6-main " "$APP/Contents/Resources/licenses/README.txt"
test -d "$APP/Contents/Resources/licenses/qt6-main" && test -d "$APP/Contents/Resources/licenses/netgen"
echo "   $(find "$APP/Contents/Resources/licenses" -mindepth 1 -maxdepth 1 -type d | wc -l) packages' licence texts"

step "signature"
codesign --verify --strict "$APP"

step "library links (nothing outside the app and macOS)"
env -i PATH=/usr/bin:/bin "$PY" "$HERE/check_app.py" "$APP"

# A run with nothing from the outside: no environment variables, python -I (no user site, no
# PYTHON* variables, not the current folder).
run() { (cd "$WORK/run" && env -i HOME="$WORK/home" TMPDIR="$WORK/tmp" PATH=/usr/bin:/bin "$PY" -I "$@"); }

step "smoke test (imports, 2D, 3D with Netgen, exports, the window offscreen)"
run "$HERE/smoke_app.py" | tail -8

step "CAD import (OpenCASCADE: STEP, IGES)"
cp -R "$ROOT/examples" "$WORK/run/examples"
rm -rf "$WORK/run/examples/output"
cp -R "$ROOT/tests/data" "$WORK/run/tests/data"
run examples/09a_cad_import.py | tail -3

echo "OPENGRID APP TEST: PASS"
