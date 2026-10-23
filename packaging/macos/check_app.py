"""Check that an app bundle depends on nothing outside itself and macOS: python check_app.py APP.

Every dependency of every Mach-O file below APP must be a system library (/usr/lib, /System:
the dyld shared cache) or resolve inside APP: @loader_path from the file's folder,
@executable_path from the bundled Python's, @rpath through the file's own rpaths (or the
Python executable's, from its folder), an absolute path only inside APP. Prints what does not
and exits 1. PySide6 modules of Qt components conda-forge's Qt does not ship (WebEngine,
Multimedia, ...: their libQt6*.dylib is missing) cannot be imported anyway: they are counted
apart."""

import os
import sys
from pathlib import Path

from relink_macho import load_commands
from sign_macho import mach_o_files

SYSTEM = ("/usr/lib/", "/System/")


def resolve(dep: str, folder: Path, exe_folder: Path, rpaths, exe_rpaths) -> list[Path]:
    """The candidate files a dependency stands for (an executable's rpaths are relative to
    its own folder)."""

    def expand(text: str, loader: Path) -> Path:
        text = text.replace("@loader_path", str(loader))
        text = text.replace("@executable_path", str(exe_folder))
        return Path(os.path.normpath(text))

    if dep.startswith("@rpath/"):
        rest = dep[len("@rpath/") :]
        found = [expand(rpath, folder) / rest for rpath in rpaths]
        return found + [expand(rpath, exe_folder) / rest for rpath in exe_rpaths]
    return [expand(dep, folder)]


def main() -> int:
    app = Path(sys.argv[1]).resolve()
    env = app / "Contents" / "Resources" / "env"
    exe_folder = env / "bin"
    _own, _deps, exe_rpaths = load_commands(exe_folder / "python3.12")
    problems = []
    absent_qt = set()
    files = mach_o_files(app)
    for path in files:
        own, deps, rpaths = load_commands(path)
        for rpath in rpaths:
            if rpath.startswith("/") and not rpath.startswith(str(app)):
                problems.append(f"{path.relative_to(app)}: rpath {rpath} outside the app")
        for dep in deps:
            if dep == own or dep.startswith(SYSTEM):
                continue
            if "/" not in dep:
                problems.append(f"{path.relative_to(app)}: {dep} (a bare name: not found)")
                continue
            candidates = resolve(dep, path.parent, exe_folder, rpaths, exe_rpaths)
            inside = [c for c in candidates if c.exists() and c.resolve().is_relative_to(app)]
            if inside:
                continue
            if path.parent.name == "PySide6" and Path(dep).name.startswith("libQt6"):
                absent_qt.add(path.name.split(".")[0])
                continue
            problems.append(f"{path.relative_to(app)}: {dep} does not resolve inside the app")
    print(f"{len(files)} Mach-O files; {len(problems)} unresolved dependencies")
    if absent_qt:
        listed = ", ".join(sorted(absent_qt))
        print(f"   (PySide6 modules of Qt components not shipped, not importable: {listed})")
    for line in problems[:40]:
        print("  ", line)
    if problems:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
