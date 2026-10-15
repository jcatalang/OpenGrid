"""Ad-hoc sign every Mach-O file below a directory: python sign_macho.py DIR.

conda-pack and conda-unpack rewrite the environment prefix inside binaries; on Apple silicon a
binary whose signature no longer matches is killed when it starts, so each one is signed again
(codesign --force --sign -). Mach-O files are found by their magic number."""

import subprocess
import sys
from pathlib import Path

MAGIC = {
    bytes.fromhex("cffaedfe"),  # 64-bit, little endian
    bytes.fromhex("cefaedfe"),  # 32-bit, little endian
    bytes.fromhex("cafebabe"),  # universal
}


def mach_o_files(root: Path) -> list[Path]:
    found = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        with path.open("rb") as handle:
            if handle.read(4) in MAGIC:
                found.append(path)
    return found


def main() -> int:
    root = Path(sys.argv[1])
    files = mach_o_files(root)
    for start in range(0, len(files), 200):
        batch = [str(p) for p in files[start : start + 200]]
        subprocess.run(
            ["codesign", "--force", "--sign", "-", *batch],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    print(f"signed {len(files)} Mach-O files below {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
