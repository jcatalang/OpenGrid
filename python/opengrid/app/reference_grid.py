"""Reference grid: the world axis planes at a step that follows the zoom.

A view along an axis shows the plane facing it (XY when looking along z, XZ along y, YZ along
x) with the coordinates of the major lines written along the viewport edges. A tilted view
shows the three axis planes, the one facing the view most with its minor lines, with the
coordinates written along the axes (billboard text that faces the camera). The planes pass
through the point with round coordinates nearest the view's focal point, so a model far from
the origin still has a grid; the lines on the true axes are coloured x red, y green, z blue.

Method: the step is 1, 2 or 5 times the power of ten nearest to a fifth of the visible
extent, split into 10, 4 or 5 minor divisions; the lines cover the visible extent with a
margin, so panning and zooming move only the camera until the next update rebuilds them (a
few hundred line cells).
"""

import math

import numpy as np
import pyvista as pv

MINOR_COLOR = "#3a4249"
MAJOR_COLOR = "#5b6670"
AXIS_COLORS = ("#8a4a4a", "#4a8a5a", "#4a6a9a")  # x, y, z
LABEL_COLOR = "#aab2ba"
TARGET_MAJOR_LINES = 5  # major lines across the view, about
MAX_MINOR_LINES = 400  # across the view: beyond this the minor lines are dropped
MARGIN = 1.5  # the grid covers this many visible extents
PLANE_NAMES = {2: "xy", 1: "xz", 0: "yz"}  # by the axis normal to the plane


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
    if text in ("-0", "-0.0", "-0.00", "-0.000"):
        text = text[1:]
    return text


class ReferenceGrid:
    def __init__(self, plotter):
        self.plotter = plotter
        self.enabled = False
        self.step = 1.0
        self.origin = np.zeros(3)  # where the planes meet (round coordinates)
        self._actors: list = []
        self._labels: list = []  # vtkTextActor (edges) or vtkBillboardTextActor3D (axes)

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

    def view(self) -> tuple[np.ndarray, float, np.ndarray]:
        """(focal point, visible extent in world units, viewing direction)."""
        renderer = self.plotter.renderer
        camera = renderer.GetActiveCamera()
        width, height = renderer.GetSize()
        aspect = max(width, 1) / max(height, 1)
        if camera.GetParallelProjection():
            half_height = camera.GetParallelScale()
        else:
            angle = 0.5 * math.radians(camera.GetViewAngle())
            half_height = camera.GetDistance() * math.tan(angle)
        extent = 2.0 * half_height * max(aspect, 1.0)
        focal = np.asarray(camera.GetFocalPoint(), dtype=float)
        direction = np.asarray(camera.GetDirectionOfProjection(), dtype=float)
        return focal, float(extent), direction

    def visible_extent(self) -> tuple[np.ndarray, float]:
        """(the planes' meeting point for the current view, visible extent)."""
        focal, extent, _direction = self.view()
        step, _divisions = grid_step(extent)
        return np.round(focal / step) * step, extent

    def update(self) -> None:
        """Rebuild the lines and labels for the current camera."""
        self.clear()
        if not self.enabled:
            return
        focal, extent, direction = self.view()
        step, divisions = grid_step(extent)
        self.step = step
        origin = np.round(focal / step) * step
        self.origin = origin
        dominant = int(np.argmax(np.abs(direction)))
        aligned = abs(direction[dominant]) > 0.999
        half = MARGIN * extent
        planes = [dominant]
        if not aligned:
            half *= 2.0  # a tilted plane is foreshortened: cover more of it
            planes = [0, 1, 2]
        for normal in planes:
            minor_divisions = 1
            if normal == dominant:
                minor_divisions = divisions
            self._plane(normal, origin, half, step, minor_divisions)
        if aligned:
            self._edge_labels(dominant, origin, step)
        else:
            self._axis_labels(origin, half, step)
        self._corner_text(origin, step)
        self.plotter.render()

    # ------------------------------------------------------------------ lines
    def _plane(self, normal: int, origin, half: float, step: float, divisions: int) -> None:
        u, v = [a for a in range(3) if a != normal]
        minor = step / divisions
        lo = origin[[u, v]] - half
        hi = origin[[u, v]] + half
        us = np.arange(np.floor(lo[0] / minor) * minor, hi[0] + 0.5 * minor, minor)
        vs = np.arange(np.floor(lo[1] / minor) * minor, hi[1] + 0.5 * minor, minor)
        if len(us) + len(vs) > 2 * MAX_MINOR_LINES:
            minor = step
            us = np.arange(np.floor(lo[0] / step) * step, hi[0] + 0.5 * step, step)
            vs = np.arange(np.floor(lo[1] / step) * step, hi[1] + 0.5 * step, step)
        major_u = np.isclose(np.round(us / step) * step, us, atol=1e-9 * step)
        major_v = np.isclose(np.round(vs / step) * step, vs, atol=1e-9 * step)
        on_axis = abs(float(origin[normal])) <= 1e-9 * step  # the plane holds two world axes
        axis_u = np.isclose(us, 0.0, atol=1e-9 * step) & on_axis
        axis_v = np.isclose(vs, 0.0, atol=1e-9 * step) & on_axis
        name = PLANE_NAMES[normal]
        offset = float(origin[normal])
        pieces = (
            ("minor", us[~major_u], vs[~major_v], MINOR_COLOR, 1.0),
            ("major", us[major_u & ~axis_u], vs[major_v & ~axis_v], MAJOR_COLOR, 1.0),
            # The line u = 0 runs along axis v and the line v = 0 along axis u.
            (f"axis-{'xyz'[v]}", us[axis_u], np.array([]), AXIS_COLORS[v], 2.0),
            (f"axis-{'xyz'[u]}", np.array([]), vs[axis_v], AXIS_COLORS[u], 2.0),
        )
        for label, line_us, line_vs, color, width in pieces:
            mesh = _lines(line_us, line_vs, lo, hi, normal, u, v, offset)
            if mesh is not None:
                self._add(mesh, f"reference-grid-{name}-{label}", color, width)

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

    # ----------------------------------------------------------------- labels
    def _edge_labels(self, normal: int, origin, step: float) -> None:
        """Coordinates of the major lines along the bottom (first in-plane axis) and left
        (second) edges of the viewport, in a view along ``normal``."""
        renderer = self.plotter.renderer
        width, height = renderer.GetSize()
        if width < 2 or height < 2:
            return
        u, v = [a for a in range(3) if a != normal]
        corner = _display_to_world(renderer, 0.0, 0.0)
        far = _display_to_world(renderer, float(width), float(height))
        offset = float(origin[normal])
        edges = ((u, v, corner[v], True), (v, u, corner[u], False))
        for axis, other, edge_value, along_x in edges:
            low = min(corner[axis], far[axis])
            high = max(corner[axis], far[axis])
            ticks = np.arange(math.ceil(low / step) * step, high + 0.5 * step, step)
            for tick in ticks:
                point = np.zeros(3)
                point[axis] = tick
                point[other] = edge_value
                point[normal] = offset
                px, py = _world_to_display(renderer, point)
                if along_x and 0 <= px <= width:
                    self._labels.append(_text(_format(float(tick), step), px + 3, 4))
                elif not along_x and 0 <= py <= height:
                    self._labels.append(_text(_format(float(tick), step), 4, py + 3))
        for label in self._labels:
            renderer.AddViewProp(label)

    def _axis_labels(self, origin, half: float, step: float) -> None:
        """Coordinates of the major ticks along the three axes through ``origin`` (tilted
        views), as text facing the camera."""
        renderer = self.plotter.renderer
        for axis in range(3):
            low = origin[axis] - half
            high = origin[axis] + half
            ticks = np.arange(math.ceil(low / step) * step, high + 0.5 * step, step)
            for tick in ticks:
                if axis > 0 and abs(tick - origin[axis]) <= 1e-9 * step:
                    continue  # the meeting point is labelled once (by the x axis)
                point = np.array(origin, dtype=float)
                point[axis] = tick
                self._labels.append(_billboard(_format(float(tick), step), point))
        for label in self._labels:
            renderer.AddViewProp(label)

    def _corner_text(self, origin, step: float) -> None:
        renderer = self.plotter.renderer
        _width, height = renderer.GetSize()
        text = f"grid {step:g}"
        if np.any(np.abs(origin) > 1e-9 * step):
            at = ", ".join(_format(float(c), step) for c in origin)
            text += f"   planes through ({at})"
        label = _text(text, 6, max(height - 18, 4))
        self._labels.append(label)
        renderer.AddViewProp(label)


def _lines(us, vs, lo, hi, normal: int, u: int, v: int, offset: float) -> pv.PolyData | None:
    """Lines at u = us (along v) and v = vs (along u) covering [lo, hi] of the plane whose
    coordinate ``normal`` is ``offset``."""
    points = []
    for value in us:
        for end in (lo[1], hi[1]):
            point = [0.0, 0.0, 0.0]
            point[u] = float(value)
            point[v] = float(end)
            point[normal] = offset
            points.append(point)
    for value in vs:
        for end in (lo[0], hi[0]):
            point = [0.0, 0.0, 0.0]
            point[u] = float(end)
            point[v] = float(value)
            point[normal] = offset
            points.append(point)
    if not points:
        return None
    n = len(points) // 2
    cells = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)])
    return pv.PolyData(np.asarray(points, dtype=float), lines=cells.ravel())


def _display_to_world(renderer, x: float, y: float) -> np.ndarray:
    renderer.SetDisplayPoint(x, y, 0.5)
    renderer.DisplayToWorld()
    wx, wy, wz, w = renderer.GetWorldPoint()
    if w == 0.0:
        w = 1.0
    return np.array([wx / w, wy / w, wz / w])


def _world_to_display(renderer, point) -> tuple[float, float]:
    renderer.SetWorldPoint(float(point[0]), float(point[1]), float(point[2]), 1.0)
    renderer.WorldToDisplay()
    x, y, _z = renderer.GetDisplayPoint()
    return float(x), float(y)


def _text(text: str, x: float, y: float):
    import vtk

    actor = vtk.vtkTextActor()
    actor.SetInput(text)
    actor.SetDisplayPosition(int(x), int(y))
    prop = actor.GetTextProperty()
    prop.SetFontSize(11)
    prop.SetColor(*pv.Color(LABEL_COLOR).float_rgb)
    prop.SetOpacity(0.9)
    actor.SetPickable(False)
    return actor


def _billboard(text: str, point):
    import vtk

    actor = vtk.vtkBillboardTextActor3D()
    actor.SetInput(text)
    actor.SetPosition(float(point[0]), float(point[1]), float(point[2]))
    actor.SetDisplayOffset(4, 4)
    prop = actor.GetTextProperty()
    prop.SetFontSize(10)
    prop.SetColor(*pv.Color(LABEL_COLOR).float_rgb)
    prop.SetOpacity(0.9)
    actor.SetPickable(False)
    return actor
