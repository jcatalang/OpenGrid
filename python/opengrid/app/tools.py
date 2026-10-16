"""Interactive creation tools and snapping (GUI logic without any rendering, unit-testable).

Tools receive world-space clicks (already snapped) and create entities through the Session,
so everything they do is an undoable, journaled command.
"""

import math
from dataclasses import dataclass

from opengrid.api import Session
from opengrid.model.connectors import CURVE_KINDS


@dataclass(frozen=True)
class SnapTarget:
    xyz: tuple[float, float, float]
    label: str  # e.g. "point-1" or "line-1 start"


def snap_targets(session: Session) -> list[SnapTarget]:
    """Point entities and the end points of every curve, in id order (connector ends lie on
    curve ends, so they need no separate targets)."""
    targets = []
    for entity in session.store:
        if entity.kind == "point":
            targets.append(SnapTarget(tuple(entity.params["xyz"]), entity.name))
            continue
        if entity.kind not in CURVE_KINDS:
            continue
        curve = session.geometry(entity.name)
        targets.append(SnapTarget(tuple(curve.eval(0.0)), f"{entity.name} start"))
        targets.append(SnapTarget(tuple(curve.eval(1.0)), f"{entity.name} end"))
    return targets


def snap(xyz, targets: list[SnapTarget], tolerance: float) -> tuple[tuple, SnapTarget | None]:
    """Nearest target within ``tolerance`` (first one on ties), or the input point."""
    best = None
    best_distance = math.inf
    for target in targets:
        d = math.dist(xyz, target.xyz)
        if d <= tolerance and d < best_distance:
            best = target
            best_distance = d
    if best is None:
        return tuple(xyz), None
    return best.xyz, best


class Tool:
    """Base class; ``points`` are the clicks collected so far."""

    name = "select"
    min_points = 0
    max_points = 0  # 0 = unlimited (finish explicitly)

    def __init__(self, session: Session):
        self.session = session
        self.points: list[tuple[float, float, float]] = []

    def prompt(self) -> str:
        return ""

    def click(self, xyz) -> str | None:
        """Add a click; returns the created entity name when the tool completes."""
        self.points.append(tuple(xyz))
        if self.max_points and len(self.points) >= self.max_points:
            return self.finish()
        return None

    def can_finish(self) -> bool:
        return len(self.points) >= self.min_points and self.min_points > 0

    def finish(self) -> str | None:
        if not self.can_finish():
            return None
        points = list(self.points)
        self.points.clear()
        return self.create(points)

    def cancel(self) -> None:
        self.points.clear()

    def create(self, points) -> str:
        raise NotImplementedError

    def preview(self, cursor) -> list[tuple]:
        """Polyline to draw as rubber band given the current cursor position."""
        if not self.points:
            return []
        return [*self.points, tuple(cursor)]


class SelectTool(Tool):
    name = "select"

    def click(self, xyz) -> None:
        return None

    def preview(self, cursor) -> list[tuple]:
        return []


class PointTool(Tool):
    name = "point"
    min_points = 1
    max_points = 1

    def prompt(self) -> str:
        return "Click to place a point"

    def create(self, points) -> str:
        return self.session.create_point(points[0])


class LineTool(Tool):
    name = "line"
    min_points = 2
    max_points = 2

    def prompt(self) -> str:
        if not self.points:
            return "Line: click the start point"
        return "Line: click the end point"

    def create(self, points) -> str:
        return self.session.create_line(points[0], points[1])


class ArcTool(Tool):
    name = "arc"
    min_points = 3
    max_points = 3

    def prompt(self) -> str:
        steps = ("start point", "a point on the arc", "end point")
        return f"Arc: click the {steps[len(self.points)]}"

    def create(self, points) -> str:
        return self.session.create_arc_3pt(points[0], points[1], points[2])


class PolylineTool(Tool):
    name = "polyline"
    min_points = 2

    def prompt(self) -> str:
        return f"Polyline: click points ({len(self.points)} so far); Enter to finish, Esc to cancel"

    def create(self, points) -> str:
        return self.session.create_polyline(points)


class SplineTool(Tool):
    name = "spline"
    min_points = 2

    def prompt(self) -> str:
        return f"Spline: click points ({len(self.points)} so far); Enter to finish, Esc to cancel"

    def create(self, points) -> str:
        return self.session.create_spline(points, method="natural")

    def preview(self, cursor) -> list[tuple]:
        points = [*self.points, tuple(cursor)]
        if len(points) < 3:
            return super().preview(cursor)
        # Draw the spline that would be created (if the points are valid).
        from opengrid.model.entities import build_geometry, normalize_params
        from opengrid.validation import ValidationError

        try:
            params = normalize_params("spline", {"points": points}, "preview")
        except ValidationError:
            return super().preview(cursor)
        curve = build_geometry("spline", params)
        return [tuple(p) for p in curve.tessellate_points(1e-3)]


class NurbsTool(Tool):
    """Click the control points of a cubic (or lower degree) clamped NURBS curve."""

    name = "nurbs"
    min_points = 2

    def prompt(self) -> str:
        return (
            f"NURBS: click control points ({len(self.points)} so far); Enter to finish, Esc to "
            "cancel"
        )

    def create(self, points) -> str:
        degree = min(3, len(points) - 1)
        spans = len(points) - degree
        inner = [k / spans for k in range(spans + 1)]
        knots = [0.0] * degree + inner + [1.0] * degree
        return self.session.create_nurbs(degree, knots, points)

    def preview(self, cursor) -> list[tuple]:
        # The control polygon (the curve is drawn once created).
        return super().preview(cursor)


class FreehandTool(Tool):
    """Press, draw and release: the stroke becomes NURBS curves (Session.create_freehand)."""

    name = "freehand"
    stroke = True
    tolerance = 0.005
    corner_angle = 60.0

    def prompt(self) -> str:
        return "Freehand: press and draw; release to fit NURBS curves (corners kept)"

    def click(self, xyz) -> None:
        return None

    def create_stroke(self, points) -> list[str]:
        if len(points) < 2:
            return []
        return self.session.create_freehand(points, self.tolerance, self.corner_angle)

    def preview(self, cursor) -> list[tuple]:
        return list(self.points)


TOOLS = {
    "select": SelectTool,
    "point": PointTool,
    "line": LineTool,
    "arc": ArcTool,
    "polyline": PolylineTool,
    "spline": SplineTool,
    "nurbs": NurbsTool,
    "freehand": FreehandTool,
}
