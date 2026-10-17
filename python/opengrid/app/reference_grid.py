"""Reference grid: lines of the z = 0 plane at a step that follows the zoom, with the x and y
axes marked and, in a view along z, the coordinates of the major lines written along the
bottom and left edges of the viewport.

Method: the step is the power of ten nearest to a fifth of the visible extent, split into
minor lines (10 per major step, 5 when that would be too dense); the lines cover the visible
part of the plane with a margin, so panning and zooming only move the camera until the next
update, which rebuilds them (cheap: a few hundred line cells). Labels are 2D text actors at
the display positions of the major lines' crossings with the viewport edges.
"""

import math

import numpy as np
import pyvista as pv

MINOR_COLOR = "#3a4249"
MAJOR_COLOR = "#5b6670"
AXIS_X_COLOR = "#8a4a4a"
AXIS_Y_COLOR = "#4a8a5a"
LABEL_COLOR = "#aab2ba"
TARGET_MAJOR_LINES = 5  # major lines across the view, about
MAX_MINOR_LINES = 400  # across the view: beyond this the minor lines are dropped
MARGIN = 1.5  # the grid covers this many visible extents


def grid_step(extent: float) -> tuple[float, int]:
    """(major step, minor divisions per major step) for a visible extent."""
    if not math.isfinite(extent) or extent <= 0.0:
        return 1.0, 10
    raw = extent / TARGET_MAJOR_LINES
    power = 10.0 ** math.floor(math.log10(raw))
    candidates = (1.0, 2.0, 5.0, 10.0)
    step = min(candidates, key=lambda c: abs(c * power - raw)) * power
    divisions = 10
    if step / power == 2.0:
        divisions = 4
    if step / power == 5.0:
        divisions = 5
    return step, divisions


def _format(value: float, step: float) -> str:
    decimals = max(0, -math.floor(math.log10(step)))
    text = f"{value:.{decimals}f}"
    if text in ("-0", "-0.0", "-0.00"):
        text = text[1:]
    return text


class ReferenceGrid:
    def __init__(self, plotter):
        self.plotter = plotter
        self.enabled = False
        self.step = 1.0
        self._actors: list = []
        self._labels: list = []

    def set_enabled(self, enabled: bool) -> None:
        self.enabled = enabled
        if enabled:
            self.update()
        else:
            self.clear()
            self.plotter.render()

    def clear(self) -> None:
        for actor in self._actors:
            self.plotter.remove_actor(actor, render=False)
        self._actors.clear()
        for label in self._labels:
            self.plotter.renderer.RemoveViewProp(label)
        self._labels.clear()

    def visible_extent(self) -> tuple[np.ndarray, float]:
        """(centre of the view on the z = 0 plane, visible extent in world units)."""
        renderer = self.plotter.renderer
        camera = renderer.GetActiveCamera()
        width, height = renderer.GetSize()
        aspect = max(width, 1) / max(height, 1)
        if camera.GetParallelProjection():
            half_height = camera.GetParallelScale()
        else:
            half_height = camera.GetDistance() * math.tan(0.5 * math.radians(camera.GetViewAngle()))
        extent = 2.0 * half_height * max(aspect, 1.0)
        focal = np.asarray(camera.GetFocalPoint(), dtype=float)
        direction = np.asarray(camera.GetDirectionOfProjection(), dtype=float)
        centre = focal.copy()
        if abs(direction[2]) > 1e-9:
            # Where the view axis meets the plane z = 0.
            s = -focal[2] / direction[2]
            if abs(s) <= 10.0 * extent:
                centre = focal + s * direction
        centre[2] = 0.0
        return centre, float(extent)

    def update(self) -> None:
        """Rebuild the lines and labels for the current camera."""
        self.clear()
        if not self.enabled:
            return
        centre, extent = self.visible_extent()
        step, divisions = grid_step(extent)
        self.step = step
        half = MARGIN * extent
        if not self._looking_along_z():
            half *= 2.0  # a tilted plane is foreshortened: cover more of it
        minor = step / divisions
        low = np.floor((centre[:2] - half) / minor) * minor
        high = np.ceil((centre[:2] + half) / minor) * minor
        xs = np.arange(low[0], high[0] + 0.5 * minor, minor)
        ys = np.arange(low[1], high[1] + 0.5 * minor, minor)
        if len(xs) + len(ys) > 2 * MAX_MINOR_LINES:
            minor = step
            xs = np.arange(np.floor(low[0] / step) * step, high[0] + 0.5 * step, step)
            ys = np.arange(np.floor(low[1] / step) * step, high[1] + 0.5 * step, step)
        major_x = np.isclose(np.round(xs / step) * step, xs, atol=1e-9 * step)
        major_y = np.isclose(np.round(ys / step) * step, ys, atol=1e-9 * step)
        axis_x = np.isclose(xs, 0.0, atol=1e-9 * step)
        axis_y = np.isclose(ys, 0.0, atol=1e-9 * step)
        pieces = {
            "minor": (xs[~major_x], ys[~major_y], MINOR_COLOR, 1.0),
            "major": (xs[major_x & ~axis_x], ys[major_y & ~axis_y], MAJOR_COLOR, 1.0),
        }
        for name, (vx, vy, color, width) in pieces.items():
            mesh = _lines(vx, vy, low, high)
            if mesh is not None:
                self._add(mesh, f"reference-grid-{name}", color, width)
        if np.any(axis_x):
            self._add(
                _lines(np.array([0.0]), np.array([]), low, high),
                "reference-grid-y-axis",
                AXIS_Y_COLOR,
                2.0,
            )
        if np.any(axis_y):
            self._add(
                _lines(np.array([]), np.array([0.0]), low, high),
                "reference-grid-x-axis",
                AXIS_X_COLOR,
                2.0,
            )
        if self._looking_along_z():
            self._add_labels(xs[major_x], ys[major_y], step)
        self.plotter.render()

    def _add(self, mesh, name: str, color: str, width: float) -> None:
        actor = self.plotter.add_mesh(
            mesh,
            color=color,
            line_width=width,
            pickable=False,
            reset_camera=False,
            render=False,
            name=None,
        )
        actor.name = name
        self._actors.append(actor)

    def _looking_along_z(self) -> bool:
        direction = self.plotter.renderer.GetActiveCamera().GetDirectionOfProjection()
        return abs(direction[2]) > 0.999

    def _add_labels(self, xs, ys, step: float) -> None:
        import vtk

        renderer = self.plotter.renderer
        width, height = renderer.GetSize()
        if width < 2 or height < 2:
            return
        # World coordinates of the viewport's bottom-left corner on the plane.
        renderer.SetDisplayPoint(0.0, 0.0, 0.5)
        renderer.DisplayToWorld()
        corner = np.asarray(renderer.GetWorldPoint()[:3], dtype=float)

        def display(point):
            renderer.SetWorldPoint(point[0], point[1], 0.0, 1.0)
            renderer.WorldToDisplay()
            return renderer.GetDisplayPoint()[:2]

        for x in xs:
            px, _py = display((x, corner[1]))
            if 0 <= px <= width:
                self._labels.append(_text(_format(float(x), step), px + 3, 4, LABEL_COLOR))
        for y in ys:
            _px, py = display((corner[0], y))
            if 0 <= py <= height:
                self._labels.append(_text(_format(float(y), step), 4, py + 3, LABEL_COLOR))
        for label in self._labels:
            renderer.AddViewProp(label)
        del vtk


def _lines(xs, ys, low, high) -> pv.PolyData | None:
    """Vertical lines at xs and horizontal lines at ys covering [low, high] in the plane."""
    points = []
    cells = []
    for x in xs:
        points.append((x, low[1], 0.0))
        points.append((x, high[1], 0.0))
    for y in ys:
        points.append((low[0], y, 0.0))
        points.append((high[0], y, 0.0))
    if not points:
        return None
    n = len(points) // 2
    cells = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)])
    return pv.PolyData(np.asarray(points, dtype=float), lines=cells.ravel())


def _text(text: str, x: float, y: float, color: str):
    import vtk

    actor = vtk.vtkTextActor()
    actor.SetInput(text)
    actor.SetDisplayPosition(int(x), int(y))
    prop = actor.GetTextProperty()
    prop.SetFontSize(11)
    rgb = pv.Color(color).float_rgb
    prop.SetColor(*rgb)
    prop.SetOpacity(0.9)
    actor.SetPickable(False)
    return actor
