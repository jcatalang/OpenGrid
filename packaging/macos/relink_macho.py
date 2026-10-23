"""Make the bundled environment's library links relative: python relink_macho.py ENV.

After conda-unpack, a few conda-forge libraries name their neighbours by absolute path (e.g.
libncurses -> .../env/lib/libtinfo.6.dylib): the path of the folder the app was built in, which
does not exist on another Mac ("import readline" failed in an app copied elsewhere). Every
dependency or rpath below ENV is rewritten relative to the file (@loader_path/...), and a
library's own absolute id to @rpath/<name>, with install_name_tool; the files must be signed
again afterwards (sign_macho.py). Mach-O files are found by their magic number."""

import os
import subprocess
import sys
from pathlib import Path

from sign_macho import mach_o_files


def load_commands(path: Path) -> tuple[str | None, list[str], list[str]]:
    """(the library's own id, its dependencies, its rpaths)."""
    out = subprocess.run(["otool", "-l", str(path)], capture_output=True, text=True).stdout
    lines = [line.strip() for line in out.splitlines()]
    own = None
    deps = []
    rpaths = []
    for k, line in enumerate(lines):
        if line.startswith("cmd "):
            kind = line.split()[1]
            for follow in lines[k + 1 : k + 4]:
                if follow.startswith("name "):
                    value = follow.split()[1]
                    if kind == "LC_ID_DYLIB":
                        own = value
                    elif kind.startswith("LC_LOAD") or kind == "LC_REEXPORT_DYLIB":
                        deps.append(value)
                if follow.startswith("path ") and kind == "LC_RPATH":
                    rpaths.append(follow.split()[1])
    return own, deps, rpaths


def main() -> int:
    env = Path(sys.argv[1]).resolve()
    prefix = str(env) + os.sep
    changed = 0
    for path in mach_o_files(env):
        own, deps, rpaths = load_commands(path)
        args = []
        if own is not None and own.startswith(prefix):
            args += ["-id", f"@rpath/{Path(own).name}"]
        for dep in deps:
            if dep.startswith(prefix):
                relative = os.path.relpath(dep, path.parent)
                args += ["-change", dep, f"@loader_path/{relative}"]
        for rpath in rpaths:
            if rpath.startswith(prefix) or rpath == str(env):
                relative = os.path.relpath(rpath, path.parent)
                args += ["-rpath", rpath, f"@loader_path/{relative}"]
        if args:
            subprocess.run(
                ["install_name_tool", *args, str(path)],
                check=True,
                stderr=subprocess.DEVNULL,
            )
            changed += 1
    print(f"relinked {changed} Mach-O files below {env}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
