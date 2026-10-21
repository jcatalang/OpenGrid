"""Enforce the project rule: no inline ternary conditionals (docs/ARCHITECTURE.md, §9).

* Python: any ``ast.IfExp`` node (``a if c else b``) is reported. Exact.
* C++: after removing comments, string/char literals (including raw strings), any remaining
  ``?`` token is reported. Since trigraphs were removed in C++17, ``?`` outside literals and
  comments can only be the conditional operator, so this lexical check is exact in practice.

Vendored third-party code (any path containing a ``third_party`` directory) is skipped.

Usage: python tools/check_no_ternary.py [paths...]   (default: repository source dirs)
Exit status 1 if any violation is found.
"""

import ast
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIRS = ("core", "bindings", "backends", "python", "tests", "examples", "tools")
CPP_SUFFIXES = {".h", ".hpp", ".hh", ".hxx", ".c", ".cc", ".cpp", ".cxx", ".ipp", ".inl"}
SKIP_PARTS = {"third_party", "build", "_skbuild", "__pycache__"}

RAW_STRING = re.compile(r'(?:u8|u|U|L)?R"([^()\\\s]{0,16})\(')


def check_python(path: Path) -> list[str]:
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as error:
        return [f"{path}:{error.lineno}: syntax error: {error.msg}"]
    problems = []
    for node in ast.walk(tree):
        if isinstance(node, ast.IfExp):
            problems.append(
                f"{path}:{node.lineno}: inline conditional expression "
                "('a if c else b'); use an explicit if/else block"
            )
    return problems


def strip_cpp(source: str) -> str:
    """Replace comments and literals by spaces, preserving newlines (so line numbers hold)."""
    out = []
    i = 0
    n = len(source)
    while i < n:
        ch = source[i]
        raw = RAW_STRING.match(source, i)
        if raw is not None:
            terminator = ")" + raw.group(1) + '"'
            end = source.find(terminator, raw.end())
            if end < 0:
                end = n
            else:
                end = end + len(terminator)
            out.append(_blank(source[i:end]))
            i = end
        elif source.startswith("//", i):
            end = source.find("\n", i)
            if end < 0:
                end = n
            out.append(_blank(source[i:end]))
            i = end
        elif source.startswith("/*", i):
            end = source.find("*/", i + 2)
            if end < 0:
                end = n
            else:
                end = end + 2
            out.append(_blank(source[i:end]))
            i = end
        elif ch == '"' or ch == "'":
            end = _end_of_quoted(source, i, ch)
            out.append(_blank(source[i:end]))
            i = end
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _end_of_quoted(source: str, start: int, quote: str) -> int:
    # A quote preceded by a digit is a C++14 digit separator (1'000'000), not a char literal.
    if quote == "'" and start > 0 and source[start - 1].isalnum():
        return start + 1
    i = start + 1
    n = len(source)
    while i < n:
        if source[i] == "\\":
            i += 2
            continue
        if source[i] == quote or source[i] == "\n":
            return i + 1
        i += 1
    return n


def _blank(text: str) -> str:
    return re.sub(r"[^\n]", " ", text)


def check_cpp(path: Path) -> list[str]:
    stripped = strip_cpp(path.read_text(encoding="utf-8"))
    problems = []
    for lineno, line in enumerate(stripped.splitlines(), start=1):
        if "?" in line:
            problems.append(
                f"{path}:{lineno}: conditional operator '?:'; use an explicit if/else block"
            )
    return problems


def iter_files(paths: list[Path]):
    for path in paths:
        if path.is_file():
            yield path
            continue
        for candidate in sorted(path.rglob("*")):
            if not candidate.is_file():
                continue
            if SKIP_PARTS.intersection(candidate.parts):
                continue
            yield candidate


def main(argv: list[str]) -> int:
    if argv:
        roots = [Path(arg) for arg in argv]
    else:
        roots = []
        for name in DEFAULT_DIRS:
            if (REPO_ROOT / name).exists():
                roots.append(REPO_ROOT / name)
    problems = []
    checked = 0
    for path in iter_files(roots):
        if path.suffix == ".py":
            problems.extend(check_python(path))
            checked += 1
        elif path.suffix in CPP_SUFFIXES:
            problems.extend(check_cpp(path))
            checked += 1
    for problem in problems:
        print(problem)
    print(f"check_no_ternary: {checked} files checked, {len(problems)} violation(s)")
    if problems:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
