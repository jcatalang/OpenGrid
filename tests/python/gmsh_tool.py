"""Locate the Gmsh executable for read-back tests: $OPENGRID_GMSH, then PATH, then the macOS
application bundle."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

CANDIDATES = ("/Applications/Gmsh.app/Contents/MacOS/gmsh",)


def find_gmsh() -> str | None:
    configured = os.environ.get("OPENGRID_GMSH")
    if configured and Path(configured).exists():
        return configured
    found = shutil.which("gmsh")
    if found:
        return found
    for candidate in CANDIDATES:
        if Path(candidate).exists():
            return candidate
    return None


def gmsh_read(path) -> str:
    """Load a mesh with Gmsh in batch mode (no meshing), re-save it, and return the log."""
    gmsh = find_gmsh()
    with tempfile.TemporaryDirectory() as tmp:
        result = subprocess.run(
            [gmsh, str(path), "-0", "-o", str(Path(tmp) / "resaved.msh"), "-format", "msh41"],
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        )
    return result.stdout + result.stderr
