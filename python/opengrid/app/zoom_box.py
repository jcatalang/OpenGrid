"""Zoom to a window: a left-drag (or a three-finger drag on a trackpad) draws a rectangle that
follows the pointer, and the release zooms the view onto it.

Method: the rectangle is a 2D overlay in display coordinates (vtkActor2D), redrawn on every
mouse move. On release the camera moves so that the rectangle's centre, taken on the focal
plane, becomes the view centre, then zooms by min(view width / box width, view height / box
height): the parallel scale shrinks (parallel projection) or the view angle (perspective), as
VTK's vtkInteractorStyleRubberBandZoom does. Drags shorter than ``min_pixels`` in either
direction are not zooms (a click stays a click).
"""

import numpy as np

BOX_COLOR = (1.0, 0.82, 0.4)


class ZoomBox:
    def __init__(self, renderer, min_pixels: int = 4):
        import vtk

        self.renderer = renderer
        self.min_pixels = min_pixels
        self.start_position = None
        self._points = vtk.vtkPoints()
        self._points.SetNumberOfPoints(4)
        for k in range(4):
            self._points.SetPoint(k, 0.0, 0.0, 0.0)
        outline = vtk.vtkCellArray()
        outline.InsertNextCell(5, [0, 1, 2, 3, 0])
        data = vtk.vtkPolyData()
        data.SetPoints(self._points)
        data.SetLines(outline)
        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(data)
        self.actor = vtk.vtkActor2D()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(*BOX_COLOR)
        self.actor.GetProperty().SetLineWidth(1.5)
        self.actor.VisibilityOff()
        self._data = data
        renderer.AddViewProp(self.actor)

    @property
    def active(self) -> bool:
        return self.start_position is not None

    def start(self, x: float, y: float) -> None:
        self.start_position = (float(x), float(y))

    def update(self, x: float, y: float) -> bool:
        """Redraw the rectangle to the pointer; True once it is large enough to show."""
        if self.start_position is None:
            return False
        x0, y0 = self.start_position
        corners = ((x0, y0), (x, y0), (x, y), (x0, y))
        for k, (px, py) in enumerate(corners):
            self._points.SetPoint(k, px, py, 0.0)
        self._points.Modified()
        self._data.Modified()
        shown = self._large(x, y)
        self.actor.SetVisibility(shown)
        return shown

    def finish(self, x: float, y: float) -> bool:
        """End the drag; zooms (and returns True) when the rectangle is large enough."""
        if self.start_position is None:
            return False
        x0, y0 = self.start_position
        zoomed = self._large(x, y)
        self.cancel()
        if zoomed:
            zoom_to_box(self.renderer, x0, y0, x, y)
        return zoomed

    def cancel(self) -> None:
        self.start_position = None
        self.actor.VisibilityOff()

    def _large(self, x: float, y: float) -> bool:
        x0, y0 = self.start_position
        return abs(x - x0) >= self.min_pixels and abs(y - y0) >= self.min_pixels


def _display_to_world(renderer, x: float, y: float, depth: float) -> np.ndarray:
    renderer.SetDisplayPoint(x, y, depth)
    renderer.DisplayToWorld()
    w = renderer.GetWorldPoint()
    return np.array(w[:3]) / w[3]


def zoom_to_box(renderer, x0: float, y0: float, x1: float, y1: float) -> float:
    """Zoom the renderer's camera onto the display rectangle (x0, y0) - (x1, y1); returns the
    zoom factor."""
    width, height = renderer.GetSize()
    if width <= 0 or height <= 0:
        return 1.0  # no view yet
    camera = renderer.GetActiveCamera()
    focal = np.array(camera.GetFocalPoint())
    renderer.SetWorldPoint(*focal, 1.0)
    renderer.WorldToDisplay()
    depth = renderer.GetDisplayPoint()[2]
    ox, oy = renderer.GetOrigin()
    centre = _display_to_world(renderer, 0.5 * (x0 + x1), 0.5 * (y0 + y1), depth)
    middle = _display_to_world(renderer, ox + 0.5 * width, oy + 0.5 * height, depth)
    shift = centre - middle
    camera.SetFocalPoint(*(focal + shift))
    camera.SetPosition(*(np.array(camera.GetPosition()) + shift))
    factor = min(width / max(abs(x1 - x0), 1.0), height / max(abs(y1 - y0), 1.0))
    camera.Zoom(factor)
    renderer.ResetCameraClippingRange()
    return factor
