"""Airfoil coordinate files (.dat) in Selig and Lednicer formats.

* Selig: optional title line, then x y pairs running from the trailing edge over the upper
  surface to the leading edge and back along the lower surface to the trailing edge.
* Lednicer: title line, a line with the point counts of the upper and lower surface, then the
  upper surface LE -> TE and the lower surface LE -> TE, usually separated by blank lines.

The result is always returned in Selig order. Problems are reported with file line numbers.
"""

from dataclasses import dataclass
from pathlib import Path

from opengrid.validation import Issue, ValidationError


@dataclass
class AirfoilData:
    title: str
    format: str  # "selig" or "lednicer"
    points: list[tuple[float, float]]  # Selig order: TE -> upper -> LE -> lower -> TE

    def leading_edge_index(self) -> int:
        """Index of the point with minimum x (first one if several)."""
        best = 0
        for index, (x, _y) in enumerate(self.points):
            if x < self.points[best][0]:
                best = index
        return best


def _parse_numbers(text: str) -> list[float] | None:
    parts = text.replace(",", " ").split()
    values = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            return None
    return values


def read_airfoil_dat(path: str | Path) -> AirfoilData:
    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as error:
        raise ValidationError([Issue(str(path), f"cannot read file: {error.strerror}")]) from None

    issues: list[Issue] = []
    title = path.stem
    rows: list[tuple[int, list[float]]] = []  # (line number, values)
    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values = _parse_numbers(stripped)
        if values is None:
            if not rows and title == path.stem:
                title = stripped
                continue
            issues.append(Issue(f"{path.name}:{number}", f"expected two numbers, got {stripped!r}"))
            continue
        if len(values) != 2:
            issues.append(
                Issue(f"{path.name}:{number}", f"expected two numbers, got {len(values)}")
            )
            continue
        rows.append((number, values))
    if issues:
        raise ValidationError(issues)
    if not rows:
        raise ValidationError([Issue(path.name, "no coordinates found")])

    first_line, first = rows[0]
    is_lednicer = (
        first[0] >= 2.0 and first[1] >= 2.0 and first[0].is_integer() and first[1].is_integer()
    )
    if is_lednicer:
        n_upper = int(first[0])
        n_lower = int(first[1])
        coordinates = rows[1:]
        if len(coordinates) != n_upper + n_lower:
            raise ValidationError(
                [
                    Issue(
                        f"{path.name}:{first_line}",
                        f"header announces {n_upper} + {n_lower} points, "
                        f"file has {len(coordinates)}",
                    )
                ]
            )
        upper = coordinates[:n_upper]
        lower = coordinates[n_upper:]
        ordered = list(reversed(upper))
        if lower and upper and lower[0][1] == upper[0][1]:
            lower = lower[1:]  # shared leading-edge point
        ordered.extend(lower)
        fmt = "lednicer"
    else:
        ordered = rows
        fmt = "selig"

    if len(ordered) < 3:
        raise ValidationError(
            [Issue(path.name, f"at least 3 points are required, got {len(ordered)}")]
        )
    for index in range(1, len(ordered)):
        if ordered[index][1] == ordered[index - 1][1]:
            issues.append(
                Issue(
                    f"{path.name}:{ordered[index][0]}",
                    f"point duplicates the previous one (line {ordered[index - 1][0]})",
                )
            )
    if issues:
        raise ValidationError(issues)
    points = [(values[0], values[1]) for _number, values in ordered]
    return AirfoilData(title=title, format=fmt, points=points)
