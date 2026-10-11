"""Point lists from CSV / whitespace-separated text files.

Accepted: comma, semicolon, tab or whitespace delimiters; '#' comments; an optional header
line naming the columns (x, y and optionally z, any order, case-insensitive); 2 or 3 columns
without a header. Problems are reported with file line numbers.
"""

import re
from pathlib import Path

from opengrid.validation import Issue, ValidationError

_SPLIT = re.compile(r"[,;\t ]+")


def read_points_csv(path: str | Path) -> list[tuple[float, float, float]]:
    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as error:
        raise ValidationError([Issue(str(path), f"cannot read file: {error.strerror}")]) from None

    issues: list[Issue] = []
    columns: list[int] | None = None  # indices of x, y, z (z may be -1 = absent)
    points: list[tuple[float, float, float]] = []
    line_numbers: list[int] = []
    for number, line in enumerate(lines, start=1):
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        fields = [field for field in _SPLIT.split(stripped) if field]
        values = []
        numeric = True
        for field in fields:
            try:
                values.append(float(field))
            except ValueError:
                numeric = False
                break
        where = f"{path.name}:{number}"
        if not numeric:
            if columns is None and not points:
                columns = _header_columns(fields, where, issues)
                continue
            issues.append(Issue(where, f"non-numeric value in {stripped!r}"))
            continue
        if columns is None:
            if len(values) not in (2, 3):
                issues.append(Issue(where, f"expected 2 or 3 columns, got {len(values)}"))
                continue
            columns = [0, 1, -1]
            if len(values) == 3:
                columns = [0, 1, 2]
        needed = max(columns) + 1
        if len(values) < needed:
            issues.append(Issue(where, f"expected {needed} columns, got {len(values)}"))
            continue
        z = 0.0
        if columns[2] >= 0:
            z = values[columns[2]]
        points.append((values[columns[0]], values[columns[1]], z))
        line_numbers.append(number)
    if issues:
        raise ValidationError(issues)
    if len(points) < 2:
        raise ValidationError(
            [Issue(path.name, f"at least 2 points are required, got {len(points)}")]
        )
    for index in range(1, len(points)):
        if points[index] == points[index - 1]:
            issues.append(
                Issue(
                    f"{path.name}:{line_numbers[index]}",
                    f"point duplicates the previous one (line {line_numbers[index - 1]})",
                )
            )
    if issues:
        raise ValidationError(issues)
    return points


def _header_columns(fields: list[str], where: str, issues: list[Issue]) -> list[int]:
    lowered = [field.strip().lower() for field in fields]
    result = []
    for axis in ("x", "y"):
        if axis not in lowered:
            issues.append(Issue(where, f"header {fields} has no '{axis}' column"))
            result.append(0)
        else:
            result.append(lowered.index(axis))
    if "z" in lowered:
        result.append(lowered.index("z"))
    else:
        result.append(-1)
    return result
