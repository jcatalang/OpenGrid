# Sourced by the launchers: prepares the bundled environment once (conda-unpack fixes the paths
# recorded at build time to this folder's location), then activates it (PATH and the packages'
# own settings, e.g. OpenCASCADE's resource paths) without anything from the user's shell.
OG_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OG_ENV="$OG_HERE/env"
if [ -f "$OG_HERE/GLIBC_MIN" ]; then
  og_need="$(cat "$OG_HERE/GLIBC_MIN")"
  og_have="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')"
  if [ -n "$og_have" ] && [ "$(printf '%s\n%s\n' "$og_need" "$og_have" | sort -V | head -1)" != "$og_need" ]; then
    echo "OpenGrid: this bundle needs glibc $og_need or newer; this system has $og_have." >&2
    echo "Build OpenGrid from source instead (docs/INSTALL_LINUX.md in the repository)." >&2
    exit 1
  fi
fi
if [ ! -f "$OG_ENV/.opengrid-location" ]; then
  echo "OpenGrid: first start in $OG_HERE, preparing the environment (once, a few seconds)" >&2
  "$OG_ENV/bin/python" "$OG_ENV/bin/conda-unpack" || {
    echo "OpenGrid: conda-unpack failed in $OG_ENV" >&2
    exit 1
  }
  echo "$OG_HERE" > "$OG_ENV/.opengrid-location"
elif [ "$(cat "$OG_ENV/.opengrid-location")" != "$OG_HERE" ]; then
  echo "OpenGrid: this folder was prepared in $(cat "$OG_ENV/.opengrid-location");" \
    "if something fails after the move, extract the archive again here" >&2
fi
unset PYTHONHOME PYTHONPATH PYTHONSTARTUP QT_PLUGIN_PATH QT_QPA_PLATFORM_PLUGIN_PATH \
  CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL LD_PRELOAD
set +u
source "$OG_ENV/bin/activate" >/dev/null
