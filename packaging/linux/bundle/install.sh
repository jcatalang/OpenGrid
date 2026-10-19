#!/bin/bash
# Put OpenGrid in the desktop's application menu and on the PATH (for this user, no root):
#   ~/.local/bin/opengrid, ~/.local/bin/opengrid-python and an "OpenGrid" menu entry that opens
#   .ogp projects. The bundle stays where it is: install from the folder you will keep.
#   ./install.sh            # install
#   ./install.sh --remove   # remove the links and the menu entry (the folder is untouched)
set -euo pipefail
HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
BIN="$HOME/.local/bin"
APPS="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
ICONS="${XDG_DATA_HOME:-$HOME/.local/share}/icons/hicolor/256x256/apps"
MIME="${XDG_DATA_HOME:-$HOME/.local/share}/mime"
if [ "${1:-}" = "--remove" ]; then
  rm -f "$BIN/opengrid" "$BIN/opengrid-python" "$APPS/opengrid.desktop" "$ICONS/opengrid.png" \
    "$MIME/packages/opengrid.xml"
  command -v update-mime-database >/dev/null && update-mime-database "$MIME" || true
  echo "removed the OpenGrid links and menu entry"
  exit 0
fi
"$HERE/opengrid-python" -c "import opengrid; print('OpenGrid', opengrid.__version__, 'ready')"
mkdir -p "$BIN" "$APPS" "$ICONS" "$MIME/packages"
ln -sf "$HERE/opengrid" "$BIN/opengrid"
ln -sf "$HERE/opengrid-python" "$BIN/opengrid-python"
cp "$HERE/opengrid.png" "$ICONS/opengrid.png"
cat > "$MIME/packages/opengrid.xml" <<XML
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-opengrid-project">
    <comment>OpenGrid project</comment>
    <glob pattern="*.ogp"/>
  </mime-type>
</mime-info>
XML
cat > "$APPS/opengrid.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=OpenGrid
Comment=Structured, unstructured and hybrid mesh generator
Exec=$HERE/opengrid %F
Icon=opengrid
Terminal=false
Categories=Science;Engineering;Education;
MimeType=application/x-opengrid-project;
DESKTOP
command -v update-mime-database >/dev/null && update-mime-database "$MIME" || true
command -v update-desktop-database >/dev/null && update-desktop-database "$APPS" || true
echo "installed: $BIN/opengrid, $BIN/opengrid-python and the OpenGrid menu entry"
case ":$PATH:" in *":$BIN:"*) ;; *) echo "note: add $BIN to your PATH to type 'opengrid'";; esac
