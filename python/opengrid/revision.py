"""The source revision OpenGrid runs from, so a window shows which code it is running."""

import subprocess
from functools import cache
from pathlib import Path


@cache
def source_revision() -> str | None:
    """Short git commit of the source tree this package runs from (an editable install), with
    "+" when the tree has uncommitted changes; None outside a git checkout (an installed
    package) or without git."""
    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists():
        return None
    try:
        # cwd: an application started by the Finder may sit in a folder it cannot read.
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            cwd=root,
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            cwd=root,
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    if status:
        return commit + "+"
    return commit
