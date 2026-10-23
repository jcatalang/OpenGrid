"""The locked runtime environments of the downloadable bundles (packaging/locks/<platform>.txt):
python lock_env.py COMMAND ..., run with conda's own Python (it uses conda's package records
and matching).

  write ENV LOCK [NOTE]     the lock of an environment: every package but OpenGrid's own, in
                            conda's explicit format (URL#md5, in install order)
  explicit LOCK PKG OUT     LOCK plus OpenGrid's freshly built package file PKG: the list that
                            "conda create --file OUT" installs as it is, without solving
  check ENV                 every package's dependencies and constraints are met inside ENV;
                            prints what is not and exits 1

A build from a lock gets the same packages every time (the bundles once changed under us: a
fresh solve picks up whatever conda-forge published since the last build). An explicit list is
not checked by conda, hence `check`: it fails the build when OpenGrid's package asks for
something the lock does not hold (a changed recipe, or host packages that moved: then refresh
the lock on purpose, OG_RELOCK=1 in packaging/linux/build_bundle.sh and
packaging/macos/build_app.sh)."""

import sys
import textwrap
from datetime import date
from pathlib import Path

from conda.core.prefix_data import PrefixData
from conda.models.match_spec import MatchSpec

OWN = "opengrid"


def records(env: Path):
    """The environment's packages in install order (dependencies first)."""
    return list(PrefixData(str(env)).iter_records_sorted())


def write(env: Path, lock: Path, note: str) -> int:
    packages = [r for r in records(env) if r.name != OWN]
    subdirs = sorted({r.subdir for r in packages if r.subdir != "noarch"})
    lines = [
        f"# The runtime environment of OpenGrid's downloadable bundle ({', '.join(subdirs)}):",
        "# every package but OpenGrid's own, exactly as installed (conda's explicit format,",
        "# URL#md5, in install order). Written by packaging/lock_env.py on " + str(date.today()),
    ]
    if note:
        lines += textwrap.wrap(
            f"from {note}.", width=94, initial_indent="# ", subsequent_indent="# "
        )
    lines += [
        "# Refreshed only on purpose: OG_RELOCK=1 packaging/linux/build_bundle.sh (the Linux",
        "# bundle workflow's 'relock' input) or packaging/macos/build_app.sh.",
        "@EXPLICIT",
    ]
    for r in packages:
        if not r.url or not r.md5:
            print(f"{r.name} {r.version} {r.build}: no download URL or md5 recorded")
            return 1
        lines.append(f"{r.url}#{r.md5}")
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.write_text("\n".join(lines) + "\n")
    print(f"{lock}: {len(packages)} packages")
    return 0


def explicit(lock: Path, package: Path, out: Path) -> int:
    text = lock.read_text()
    if "@EXPLICIT" not in text:
        print(f"{lock}: not an explicit package list")
        return 1
    out.write_text(text.rstrip("\n") + "\n" + package.resolve().as_uri() + "\n")
    return 0


def check(env: Path) -> int:
    installed = {r.name: r for r in records(env)}
    problems = []
    for r in installed.values():
        for dep in r.depends:
            spec = MatchSpec(dep)
            if spec.name.startswith("__"):
                continue
            match = installed.get(spec.name)
            if match is None:
                problems.append(f"{r.name} {r.version} needs {dep}: not installed")
            elif not spec.match(match):
                problems.append(f"{r.name} {r.version} needs {dep}: {match.version} {match.build}")
        for constraint in r.constrains:
            spec = MatchSpec(constraint)
            match = installed.get(spec.name)
            if match is not None and not spec.match(match):
                problems.append(
                    f"{r.name} {r.version} allows only {constraint}: {match.version} {match.build}"
                )
    print(f"{env}: {len(installed)} packages; {len(problems)} unmet dependencies")
    for line in problems:
        print("  ", line)
    if problems:
        print("   (the lock no longer fits OpenGrid's package: refresh it with OG_RELOCK=1)")
        return 1
    return 0


def main() -> int:
    command, args = sys.argv[1], [Path(a) for a in sys.argv[2:]]
    if command == "write":
        note = ""
        if len(sys.argv) > 4:
            note = sys.argv[4]
        return write(args[0], args[1], note)
    if command == "explicit":
        return explicit(args[0], args[1], args[2])
    if command == "check":
        return check(args[0])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
