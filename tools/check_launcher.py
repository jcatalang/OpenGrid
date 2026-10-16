"""Is the installed OpenGrid.app launcher current?

The launcher (packaging/macos/install_launcher.sh) runs this checkout's code, so code changes
need no reinstall; the bundle itself (launcher script, Info.plist, icon) is current when the
installer and icon scripts it was made from are unchanged. Exit 0: current (or not installed
from this checkout); 1: outdated - rerun packaging/macos/install_launcher.sh.

    python tools/check_launcher.py [/Applications/OpenGrid.app]
"""

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = {
    "installer": ROOT / "packaging" / "macos" / "install_launcher.sh",
    "icon": ROOT / "packaging" / "macos" / "make_icon.py",
}


def stamp(app: Path) -> dict | None:
    path = app / "Contents" / "Resources" / "source.txt"
    if not path.exists():
        return None
    out = {}
    for line in path.read_text().splitlines():
        key, _, value = line.partition(": ")
        out[key] = value
    return out


def problems(app: Path) -> list[str] | None:
    """What makes the installed bundle outdated; None when there is no launcher of this
    checkout to check."""
    found = stamp(app)
    if found is None:
        launcher = app / "Contents" / "MacOS" / "OpenGrid"
        if launcher.exists() and "install_launcher" in launcher.read_text(errors="replace"):
            return ["installed by an older install_launcher.sh (no source.txt)"]
        return None
    if Path(found.get("source", "")).resolve() != ROOT:
        return None
    out = []
    for key, script in SCRIPTS.items():
        digest = hashlib.sha256(script.read_bytes()).hexdigest()
        if found.get(key) != digest:
            out.append(f"{script.relative_to(ROOT)} changed since the bundle was installed")
    return out


def main() -> int:
    app = Path("/Applications/OpenGrid.app")
    if len(sys.argv) > 1:
        app = Path(sys.argv[1])
    issues = problems(app)
    if issues is None:
        print(f"{app}: no launcher of this checkout installed")
        return 0
    if issues:
        for issue in issues:
            print(f"{app}: {issue}")
        print("rerun packaging/macos/install_launcher.sh")
        return 1
    commit = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    print(f"{app}: current (it runs this checkout, now at {commit})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
