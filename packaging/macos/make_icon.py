"""Draw the OpenGrid icon (a stretched O-grid around a circle) into an .iconset folder:
python make_icon.py out.iconset, then iconutil -c icns out.iconset."""

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def draw(path: Path, pixels: int) -> None:
    fig = plt.figure(figsize=(1, 1), dpi=pixels)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(plt.Rectangle((-1, -1), 2, 2, color="#1d3557", zorder=0))
    theta = np.linspace(0, 2 * math.pi, 25)
    radii = 0.22 + 0.68 * (np.linspace(0, 1, 9) ** 1.6)
    width = max(0.4, pixels / 256)
    for r in radii:
        ax.plot(r * np.cos(theta), r * np.sin(theta), color="#a8dadc", lw=width)
    for t in theta[:-1]:
        ax.plot(
            [radii[0] * math.cos(t), radii[-1] * math.cos(t)],
            [radii[0] * math.sin(t), radii[-1] * math.sin(t)],
            color="#a8dadc",
            lw=width,
        )
    ax.add_patch(plt.Circle((0, 0), radii[0], color="#e63946", zorder=3))
    fig.savefig(path, dpi=pixels)
    plt.close(fig)


def main() -> int:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for size in (16, 32, 128, 256, 512):
        draw(out / f"icon_{size}x{size}.png", size)
        draw(out / f"icon_{size}x{size}@2x.png", 2 * size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
