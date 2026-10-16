#!/usr/bin/env bash
# Build OpenGrid.app (macOS arm64) and a DMG.
#
# The bundle holds a relocatable conda environment (conda-pack) with the OpenGrid conda package
# (packaging/conda/recipe.yaml) and its conda-forge dependencies (Qt, VTK, Netgen, OpenCASCADE,
# libomp, ...), started by a small launcher. Signed ad hoc (no Developer ID: Gatekeeper asks for
# confirmation on first launch from a download).
#
#   packaging/macos/build_app.sh            # builds the conda package first if missing
#   packaging/macos/build_app.sh --rebuild  # always rebuild the conda package
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

if [[ "$(uname)" != "Darwin" ]]; then
  echo "build_app.sh: macOS only" >&2
  exit 1
fi

# 1. The conda package.
if [[ "${1:-}" == "--rebuild" ]] || ! ls "$CHANNEL"/*/opengrid-"$VERSION"-*.conda >/dev/null 2>&1; then
  "$TOOLS/bin/rattler-build" build --recipe "$ROOT/packaging/conda/recipe.yaml" \
    -m "$ROOT/packaging/conda/variants.yaml" -c conda-forge --output-dir "$CHANNEL"
fi

# 2. A fresh environment with the package, then packed (relocatable).
rm -rf "$WORK" "$APP"
mkdir -p "$WORK"
conda create -y -q -p "$WORK/env" --override-channels -c "file://$CHANNEL" -c conda-forge \
  "opengrid=$VERSION" python=3.12
"$TOOLS/bin/conda-pack" -p "$WORK/env" -o "$WORK/env.tar.gz" --quiet

# 3. The bundle.
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/env"
tar -xzf "$WORK/env.tar.gz" -C "$APP/Contents/Resources/env"
# conda-pack / conda-unpack rewrite prefixes inside binaries: re-sign them (ad hoc) before
# anything in the environment runs, and again after conda-unpack.
"$TOOLS/bin/python" "$HERE/sign_macho.py" "$APP/Contents/Resources/env"
# (conda-unpack starts with "#!/usr/bin/env python": run it with the bundled interpreter.)
"$APP/Contents/Resources/env/bin/python" "$APP/Contents/Resources/env/bin/conda-unpack"
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

# 4. Ad-hoc signature of the bundle itself (its binaries were signed above).
codesign --force --sign - "$APP"
codesign --verify "$APP"

# 5. Disk image.
DMG="$DIST/OpenGrid-$VERSION-macos-$ARCH.dmg"
rm -f "$DMG"
hdiutil create -quiet -volname "OpenGrid $VERSION" -srcfolder "$APP" -ov -format UDZO "$DMG"
du -sh "$APP" "$DMG"
echo "built $APP and $DMG"
