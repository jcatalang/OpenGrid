"""Shared helpers for the Python tests."""

import math


def naca0012_selig(n: int = 41) -> list[tuple[float, float]]:
    """NACA 0012 (closed TE) in Selig order, cosine spacing, n points per surface."""
    surface = []
    for i in range(n):
        x = 0.5 * (1.0 - math.cos(math.pi * i / (n - 1)))
        y = 0.6 * (
            0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
        )
        surface.append((x, y))
    upper = list(reversed(surface))
    lower = [(x, -y) for x, y in surface[1:]]
    return upper + lower


def write_selig(path, points, title="NACA 0012") -> None:
    lines = [title]
    for x, y in points:
        lines.append(f"{x:.10f} {y:.10f}")
    path.write_text("\n".join(lines) + "\n")
