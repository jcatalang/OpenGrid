"""The licence texts of every package in a conda environment, for the downloadable bundles:
python collect_licenses.py ENV OUT.

conda keeps a package's licence files in its extracted folder in the package cache
(info/licenses/, recorded as extracted_package_dir in ENV/conda-meta/*.json), not in the
environment itself, so a packed environment (conda-pack) carries none. They are copied to
OUT/<package>/ and listed in OUT/README.txt with each package's version and licence; packages
without licence files are listed too. OpenGrid's own licence is LICENSE at the bundle's top."""

import json
import shutil
import sys
from pathlib import Path


def main() -> int:
    env, out = Path(sys.argv[1]), Path(sys.argv[2])
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    missing = []
    for record in sorted((env / "conda-meta").glob("*.json")):
        meta = json.loads(record.read_text())
        name = meta["name"]
        licence = meta.get("license") or "(not stated)"
        rows.append((name, meta["version"], licence))
        source = Path(meta.get("extracted_package_dir") or "") / "info" / "licenses"
        if not source.is_dir() or not any(source.iterdir()):
            missing.append(name)
            continue
        shutil.copytree(source, out / name)
    width = max(len(name) for name, _version, _licence in rows)
    lines = [
        "Licences of the third-party packages bundled with OpenGrid (conda-forge builds).",
        "OpenGrid itself is MIT-licensed: ../LICENSE. Each package's licence texts are in the",
        "folder of its name; the sources are on https://conda-forge.org (feedstocks) and the",
        "projects' own sites.",
        "",
    ]
    for name, version, licence in rows:
        lines.append(f"{name.ljust(width)}  {version:<16} {licence}")
    if missing:
        lines += ["", "Without licence files in their package: " + ", ".join(missing)]
    (out / "README.txt").write_text("\n".join(lines) + "\n")
    print(f"licences of {len(rows) - len(missing)} of {len(rows)} packages in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
