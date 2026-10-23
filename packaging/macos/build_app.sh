#!/usr/bin/env bash
# Build OpenGrid.app (macOS arm64) and a DMG.
#
# The bundle holds a relocatable conda environment (conda-pack) with the OpenGrid conda package
# (packaging/conda/recipe.yaml) and its conda-forge dependencies (Qt, VTK, Netgen, OpenCASCADE,
# libomp, ...), started by a small launcher: nothing else is needed on the Mac. Its library
# links are relative (relink_macho.py; checked by check_app.py), so it runs wherever it is
# copied. Signed ad hoc (no Developer ID: Gatekeeper asks for confirmation on first launch from
# a download). packaging/macos/test_app.sh checks the DMG as a user would get it.
#
#   packaging/macos/build_app.sh            # builds the conda package first if missing
#   packaging/macos/build_app.sh --rebuild  # always rebuild the conda package
#   OG_RELOCK=1 packaging/macos/build_app.sh   # a fresh solve: rewrites the lock (below)
#
# The environment is the lock packaging/locks/osx-<arch>.txt (exact package builds, installed
# without solving; packaging/lock_env.py), so a rebuild ships the same packages. OG_RELOCK=1
# solves afresh and rewrites that lock: test the app (test_app.sh), then commit the lock. The
# compilers and host packages are pinned in packaging/conda/variants.yaml.
#
# Tools: rattler-build and conda-pack (conda create -n opengrid-pack -c conda-forge
# rattler-build conda-pack), found in $OPENGRID_PACK_ENV (default: that environment).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT/pyproject.toml")"
CONDA_BASE="$(conda info --base)"
TOOLS="${OPENGRID_PACK_ENV:-$CONDA_BASE/envs/opengrid-pack}"
CHANNEL="$ROOT/packaging/output"
DIST="$ROOT/packaging/dist"
WORK="$DIST/work"
APP="$DIST/OpenGrid.app"
ARCH="$(uname -m)"
SUBDIR="osx-$ARCH"
if [[ "$ARCH" == "x86_64" ]]; then
  SUBDIR="osx-64"
fi
LOCK="$ROOT/packaging/locks/$SUBDIR.txt"
# conda's own Python (packaging/lock_env.py uses conda's package records).
CONDA_PYTHON="$CONDA_BASE/bin/python"

if [[ "$(uname)" != "Darwin" ]]; then
  echo "build_app.sh: macOS only" >&2
  exit 1
fi

# 1. The conda package.
if [[ "${1:-}" == "--rebuild" ]] || ! ls "$CHANNEL/$SUBDIR"/opengrid-"$VERSION"-*.conda >/dev/null 2>&1; then
  rm -f "$CHANNEL/$SUBDIR"/opengrid-"$VERSION"-*.conda
  "$TOOLS/bin/rattler-build" build --recipe "$ROOT/packaging/conda/recipe.yaml" \
    -m "$ROOT/packaging/conda/variants.yaml" -c conda-forge --output-dir "$CHANNEL"
fi

# 2. A fresh environment with the package, then packed (relocatable).
rm -rf "$WORK" "$APP"
mkdir -p "$WORK"
# conda reuses a cached package of the same file name: OpenGrid's must be the one just built
# (a rebuilt 0.0.2 has the same name as the 0.0.2 before it).
for pkgs in $(conda info --json | "$TOOLS/bin/python" -c "import json, sys; print(' '.join(json.load(sys.stdin)['pkgs_dirs']))"); do
  rm -rf "$pkgs"/opengrid-"$VERSION"-*
done
if [[ -n "${OG_RELOCK:-}" ]]; then
  conda create -y -q -p "$WORK/env" --override-channels -c "file://$CHANNEL" -c conda-forge \
    "opengrid=$VERSION" python=3.12
  "$CONDA_PYTHON" "$HERE/../lock_env.py" write "$WORK/env" "$LOCK" "a fresh solve on $(date +%F)"
else
  PACKAGE="$(ls -t "$CHANNEL/$SUBDIR"/opengrid-"$VERSION"-*.conda | head -1)"
  "$CONDA_PYTHON" "$HERE/../lock_env.py" explicit "$LOCK" "$PACKAGE" "$WORK/explicit.txt"
  conda create -y -q -p "$WORK/env" --file "$WORK/explicit.txt"
fi
# An explicit list is installed unchecked: every dependency met (OpenGrid's included)?
"$CONDA_PYTHON" "$HERE/../lock_env.py" check "$WORK/env"
# The licence texts of every bundled package (in the package cache, not in a packed env).
"$TOOLS/bin/python" "$HERE/../collect_licenses.py" "$WORK/env" "$WORK/licenses"
"$TOOLS/bin/conda-pack" -p "$WORK/env" -o "$WORK/env.tar.gz" --quiet

# 3. The bundle.
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/env"
cp "$ROOT/LICENSE" "$APP/Contents/Resources/LICENSE"
cp -R "$WORK/licenses" "$APP/Contents/Resources/licenses"
tar -xzf "$WORK/env.tar.gz" -C "$APP/Contents/Resources/env"
# conda-pack / conda-unpack rewrite prefixes inside binaries: re-sign them (ad hoc) before
# anything in the environment runs, and again after conda-unpack.
"$TOOLS/bin/python" "$HERE/sign_macho.py" "$APP/Contents/Resources/env"
# (conda-unpack starts with "#!/usr/bin/env python": run it with the bundled interpreter.)
"$APP/Contents/Resources/env/bin/python" "$APP/Contents/Resources/env/bin/conda-unpack"
# conda-unpack wrote this folder's path into a few libraries' links (libncurses -> libtinfo,
# ...): relative links instead, so the app runs wherever it is copied.
"$TOOLS/bin/python" "$HERE/relink_macho.py" "$APP/Contents/Resources/env"
"$TOOLS/bin/python" "$HERE/sign_macho.py" "$APP/Contents/Resources/env"
# Byte-code caches written by conda-unpack would be rewritten on first run anyway.
find "$APP/Contents/Resources/env" -name "__pycache__" -type d -prune -exec rm -rf {} +

cat >"$APP/Contents/MacOS/OpenGrid" <<'LAUNCHER'
#!/bin/bash
# OpenGrid launcher: the bundled Python environment, nothing from the user's shell.
CONTENTS="$(cd "$(dirname "$0")/.." && pwd)"
ENV="$CONTENTS/Resources/env"
unset PYTHONHOME PYTHONPATH CONDA_PREFIX QT_PLUGIN_PATH
export PATH="$ENV/bin:/usr/bin:/bin"
exec "$ENV/bin/python" -m opengrid "$@"
LAUNCHER
chmod +x "$APP/Contents/MacOS/OpenGrid"

"$APP/Contents/Resources/env/bin/python" "$HERE/make_icon.py" "$WORK/OpenGrid.iconset"
iconutil -c icns "$WORK/OpenGrid.iconset" -o "$APP/Contents/Resources/OpenGrid.icns"

cat >"$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>OpenGrid</string>
  <key>CFBundleDisplayName</key><string>OpenGrid</string>
  <key>CFBundleIdentifier</key><string>org.opengrid.OpenGrid</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>OpenGrid</string>
  <key>CFBundleIconFile</key><string>OpenGrid</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>LSArchitecturePriority</key><array><string>$ARCH</string></array>
  <key>NSHighResolutionCapable</key><true/>
  <key>NSRequiresAquaSystemAppearance</key><false/>
  <key>CFBundleDocumentTypes</key>
  <array>
    <dict>
      <key>CFBundleTypeName</key><string>OpenGrid project</string>
      <key>CFBundleTypeRole</key><string>Editor</string>
      <key>LSHandlerRank</key><string>Owner</string>
      <key>LSItemContentTypes</key><array><string>org.opengrid.project</string></array>
      <key>CFBundleTypeIconFile</key><string>OpenGrid</string>
    </dict>
  </array>
  <key>UTExportedTypeDeclarations</key>
  <array>
    <dict>
      <key>UTTypeIdentifier</key><string>org.opengrid.project</string>
      <key>UTTypeDescription</key><string>OpenGrid project</string>
      <key>UTTypeConformsTo</key><array><string>public.json</string></array>
      <key>UTTypeIconFile</key><string>OpenGrid</string>
      <key>UTTypeTagSpecification</key>
      <dict><key>public.filename-extension</key><array><string>ogp</string></array></dict>
    </dict>
  </array>
</dict>
</plist>
PLIST

# 4. Ad-hoc signature of the bundle itself (its binaries were signed above), and nothing it
# links outside itself and macOS.
codesign --force --sign - "$APP"
codesign --verify "$APP"
"$TOOLS/bin/python" "$HERE/check_app.py" "$APP"

# 5. Disk image.
DMG="$DIST/OpenGrid-$VERSION-macos-$ARCH.dmg"
rm -f "$DMG"
hdiutil create -quiet -volname "OpenGrid $VERSION" -srcfolder "$APP" -ov -format UDZO "$DMG"
du -sh "$APP" "$DMG"
echo "built $APP and $DMG"
