"""Spacing a boundary curve like the inflation layers that run into it (2D sliding layers).

Where a wall's layers slide along a neighbouring boundary curve instead of tapering
(``slide`` in the 2D inflation, core/src/unstructured/inflation.cpp), each layer steps through
one of that curve's nodes, so the curve needs nodes at the layer heights near the wall. This
gives the curve, at each end that meets a wall, explicit nodes at the cumulative layer heights
h, h (1 + g), h (1 + g + g^2), ... over the layers' thickness (a breakpoint there), and a tanh
distribution in between that starts from the next layer height and keeps about the curve's
former spacing. Only the curve's parameters are computed here; its nodes are placed by the C++
distribution as for any curve.
"""

import math

# Where a wall and the boundary next to it meet at an angle (inside the domain) outside these
# limits, the layers do not slide (the C++ inflation's slide_min_angle_deg /
# slide_max_angle_deg): a sharper wedge would squeeze them into thin cells, a wider opening
# has nothing to run into.
SLIDE_MIN_ANGLE = 45.0
SLIDE_MAX_ANGLE = 135.0


def layer_heights(first: float, growth: float, layers: int, total_thickness=None) -> list[float]:
    """Cumulative heights 0, h, h (1 + g), ... of ``layers`` layers (fewer when a total
    thickness is reached first)."""
    heights = [0.0]
    h = first
    for _ in range(layers):
        heights.append(heights[-1] + h)
        h *= growth
        if total_thickness is not None and heights[-1] >= total_thickness:
            break
    return heights


def layered_params(params: dict, length: float, growth: float, begin=None, end=None) -> dict:
    """Connector parameters spaced like layers at its begin and / or end (cumulative heights
    from ``layer_heights``; None: that end keeps no layers). Raises ValueError when the layers
    would take most of the curve."""
    begin_stack = begin or [0.0]
    end_stack = end or [0.0]
    thickness_a = begin_stack[-1]
    thickness_b = end_stack[-1]
    middle = length - thickness_a - thickness_b
    if middle < 0.25 * length:
        raise ValueError(
            f"the layers ({thickness_a + thickness_b:.4g} thick) would fill most of the curve "
            f"({length:.4g} long)"
        )
    spacing = length / max(params["dimension"] - 1, 1)
    intervals = max(2, round(middle / spacing))
    count_a = len(begin_stack) - 1
    count_b = len(end_stack) - 1
    dimension = count_a + intervals + count_b + 1
    # The middle: tanh from the next layer height at each end with layers.
    middle_distribution = {"type": "tanh", "begin_spacing": None, "end_spacing": None}
    if count_a:
        middle_distribution["begin_spacing"] = (thickness_a - begin_stack[-2]) * growth
    if count_b:
        middle_distribution["end_spacing"] = (thickness_b - end_stack[-2]) * growth
    breakpoints = []
    distribution = middle_distribution
    if count_a:
        distribution = {"type": "explicit", "positions": [c / thickness_a for c in begin_stack]}
        breakpoints.append(
            {
                "position": thickness_a / length,
                "index": count_a,
                "distribution": middle_distribution,
            }
        )
    if count_b:
        positions = [1.0 - c / thickness_b for c in reversed(end_stack)]
        breakpoints.append(
            {
                "position": 1.0 - thickness_b / length,
                "index": count_a + intervals,
                "distribution": {"type": "explicit", "positions": positions},
            }
        )
    result = dict(params)
    result.update(
        {
            "dimension": dimension,
            "mode": "arc_length",
            "distribution": distribution,
            "breakpoints": tuple(breakpoints),
        }
    )
    return result


def suggested_first_height(spacing: float) -> float:
    """A first layer height from the wall spacing (a twentieth, as the property panel's)."""
    if not math.isfinite(spacing) or spacing <= 0:
        return 1e-3
    return float(f"{spacing / 20.0:.3g}")
