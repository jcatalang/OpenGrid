#!/usr/bin/env bash
# Install OpenGrid.app, a launcher for the development environment (macOS).
#
# Unlike build_app.sh (a self-contained bundle with a packed copy of the environment), this app
# holds no Python: it starts `python -m opengrid` from the `opengrid` conda environment, so it
# always runs the current code of the editable install (no rebuild after changes). Output goes
# to ~/Library/Logs/OpenGrid.log. .ogp projects open from the Finder (double-click). The
# bundle itself (launcher, Info.plist, icon) changes only with this script and make_icon.py:
# tools/check_launcher.py says whether the installed one is current.
#
#   packaging/macos/install_launcher.sh               # into /Applications
#   packaging/macos/install_launcher.sh ~/Applications
#
# The environment: $OPENGRID_ENV (a prefix), else the conda environment named `opengrid`.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
TARGET="${1:-/Applications}"
ENV_PREFIX="${OPENGRID_ENV:-$(conda info --base)/envs/opengrid}"
PYTHON="$ENV_PREFIX/bin/python"
VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT/pyproject.toml")"
APP="$TARGET/OpenGrid.app"

if [[ "$(uname)" != "Darwin" ]]; then
  echo "install_launcher.sh: macOS only" >&2
  exit 1
fi
if ! "$PYTHON" -c "import opengrid" >/dev/null 2>&1; then
  echo "install_launcher.sh: $PYTHON cannot import opengrid (pip install -e . first)" >&2
  exit 1
fi
if [[ -e "$APP" ]] && ! grep -q "install_launcher" "$APP/Contents/MacOS/OpenGrid" 2>/dev/null; then
  echo "install_launcher.sh: $APP exists and is not this launcher; remove it first" >&2
  exit 1
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"

cat >"$APP/Contents/MacOS/OpenGrid" <<LAUNCHER
#!/bin/bash
# OpenGrid launcher (packaging/macos/install_launcher.sh): the development environment.
export PATH="$ENV_PREFIX/bin:/usr/bin:/bin:/usr/sbin:/sbin"
unset PYTHONHOME PYTHONPATH QT_PLUGIN_PATH
LOG="\$HOME/Library/Logs/OpenGrid.log"
mkdir -p "\$(dirname "\$LOG")"
echo "--- \$(date) OpenGrid started" >>"\$LOG"
exec "$PYTHON" -m opengrid "\$@" >>"\$LOG" 2>&1
LAUNCHER
chmod +x "$APP/Contents/MacOS/OpenGrid"

"$PYTHON" "$HERE/make_icon.py" "$WORK/OpenGrid.iconset"
iconutil -c icns "$WORK/OpenGrid.iconset" -o "$APP/Contents/Resources/OpenGrid.icns"

cat >"$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>OpenGrid</string>
  <key>CFBundleDisplayName</key><string>OpenGrid</string>
  <key>CFBundleIdentifier</key><string>org.opengrid.OpenGrid.dev</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>OpenGrid</string>
  <key>CFBundleIconFile</key><string>OpenGrid</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
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

# What the bundle was made from: checked by tools/check_launcher.py (and shown in the log).
{
  echo "source: $ROOT"
  echo "commit: $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "installer: $(shasum -a 256 "$HERE/install_launcher.sh" | cut -d' ' -f1)"
  echo "icon: $(shasum -a 256 "$HERE/make_icon.py" | cut -d' ' -f1)"
  echo "installed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >"$APP/Contents/Resources/source.txt"

# Refresh the Finder's view of the bundle (icon) and register it with Launch Services.
touch "$APP"
LSREGISTER=/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister
"$LSREGISTER" -f "$APP" >/dev/null 2>&1 || true
# A built self-contained bundle (build_app.sh, packaging/dist) is a frozen copy: keep the
# Finder from choosing it for .ogp files.
"$LSREGISTER" -u "$ROOT/packaging/dist/OpenGrid.app" >/dev/null 2>&1 || true
"$PYTHON" "$HERE/set_default_handler.py" org.opengrid.project org.opengrid.OpenGrid.dev
echo "installed $APP (runs $PYTHON -m opengrid; log: ~/Library/Logs/OpenGrid.log)"
