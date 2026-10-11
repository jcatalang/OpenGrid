"""Viewport scene: renders model entities with pyvista and maps screen positions to the model.

Picking is done in world space with the core's closest-point projection (not VTK prop
picking), so it is exact for curves, independent of tessellation, and unit-testable without a
render window.
"""

import math

import numpy as np
import pyvista as pv

from opengrid.api import Session

CURVE_COLOR = "#5ad1ff"
POINT_COLOR = "#f2f2f2"
SELECTED_COLOR = "#ffd60a"
PREVIEW_COLOR = "#ff9f43"
SNAP_COLOR = "#2ee59d"
ENDPOINT_COLOR = "#8a96a3"
TESSELLATION_TOLERANCE = 2e-4


class Scene:
    def __init__(self, plotter, session: Session):
        self.plotter = plotter
        self.session = session
        self.selected_ids: set[int] = set()
        self._actors: dict[int, str] = {}  # entity id -> actor name

    # --------------------------------------------------------------- rendering
    def rebuild(self) -> None:
        for actor_name in list(self._actors.values()):
            self.plotter.remove_actor(actor_name, render=False)
        self._actors.clear()
        for entity in self.session.store:
            self._add_entity_actor(entity.id)
        self._update_endpoints()
        self.plotter.render()

    def update_entities(self, ids) -> None:
        for entity_id in ids:
            actor_name = self._actors.pop(entity_id, None)
            if actor_name is not None:
                self.plotter.remove_actor(actor_name, render=False)
            if entity_id in self.session.store:
                self._add_entity_actor(entity_id)
        self._update_endpoints()
        self.plotter.render()

    def set_selection(self, ids) -> None:
        changed = self.selected_ids.symmetric_difference(ids)
        self.selected_ids = set(ids)
        for entity_id in changed:
            if entity_id in self.session.store:
                self._add_entity_actor(entity_id)
        self.plotter.render()

    def show_preview(self, points) -> None:
        if len(points) < 2:
            self.plotter.remove_actor("tool-preview", render=False)
        else:
            self.plotter.add_mesh(
                _polyline(np.asarray(points, dtype=float)),
                color=PREVIEW_COLOR,
                line_width=2,
                name="tool-preview",
                pickable=False,
                reset_camera=False,
            )
        self.plotter.render()

    def show_snap(self, xyz) -> None:
        if xyz is None:
            self.plotter.remove_actor("snap-marker", render=False)
        else:
            self.plotter.add_mesh(
                pv.PolyData(np.asarray([xyz], dtype=float)),
                color=SNAP_COLOR,
                point_size=16,
                render_points_as_spheres=True,
                name="snap-marker",
                pickable=False,
                reset_camera=False,
            )
        self.plotter.render()

    def fit(self) -> None:
        self.plotter.view_xy()
        self.plotter.reset_camera()
        self.plotter.render()

    def _add_entity_actor(self, entity_id: int) -> None:
        entity = self.session.store.get(entity_id)
        actor_name = f"entity-{entity_id}"
        selected = entity_id in self.selected_ids
        if entity.kind == "point":
            color = POINT_COLOR
            if selected:
                color = SELECTED_COLOR
            self.plotter.add_mesh(
                pv.PolyData(np.asarray([entity.params["xyz"]], dtype=float)),
                color=color,
                point_size=12,
                render_points_as_spheres=True,
                name=actor_name,
                pickable=False,
                reset_camera=False,
            )
        else:
            color = CURVE_COLOR
            width = 2.5
            if selected:
                color = SELECTED_COLOR
                width = 4.0
            curve = self.session.geometry(entity.name)
            self.plotter.add_mesh(
                _polyline(curve.tessellate_points(TESSELLATION_TOLERANCE)),
                color=color,
                line_width=width,
                name=actor_name,
                pickable=False,
                reset_camera=False,
            )
        self._actors[entity_id] = actor_name

    def _update_endpoints(self) -> None:
        ends = []
        for entity in self.session.store:
            if entity.kind != "point":
                curve = self.session.geometry(entity.name)
                ends.append(curve.eval(0.0))
                ends.append(curve.eval(1.0))
        if not ends:
            self.plotter.remove_actor("curve-endpoints", render=False)
            return
        self.plotter.add_mesh(
            pv.PolyData(np.asarray(ends, dtype=float)),
            color=ENDPOINT_COLOR,
            point_size=7,
            render_points_as_spheres=True,
            name="curve-endpoints",
            pickable=False,
            reset_camera=False,
        )

    # ------------------------------------------------------- screen <-> model
    def display_to_world(self, x: float, y: float) -> tuple[float, float, float]:
        """Intersect the view ray through display pixel (x, y) with the plane z = 0."""
        renderer = self.plotter.renderer
        near = _display_point_to_world(renderer, x, y, 0.0)
        far = _display_point_to_world(renderer, x, y, 1.0)
        dz = far[2] - near[2]
        if abs(dz) < 1e-300:
            return (near[0], near[1], 0.0)
        s = -near[2] / dz
        return (near[0] + s * (far[0] - near[0]), near[1] + s * (far[1] - near[1]), 0.0)

    def world_per_pixel(self) -> float:
        renderer = self.plotter.renderer
        height = max(renderer.GetSize()[1], 1)
        camera = renderer.GetActiveCamera()
        if camera.GetParallelProjection():
            return 2.0 * camera.GetParallelScale() / height
        distance = camera.GetDistance()
        angle = math.radians(camera.GetViewAngle())
        return 2.0 * distance * math.tan(0.5 * angle) / height


def pick_entity(session: Session, xyz, tolerance: float) -> str | None:
    """Name of the entity closest to ``xyz`` within ``tolerance`` (points win over curves)."""
    best_name = None
    best_score = math.inf
    for entity in session.store:
        if entity.kind == "point":
            # Points sit on top of curve ends: halve their distance so they win the tie.
            score = 0.5 * math.dist(xyz, entity.params["xyz"])
            distance = 2.0 * score
        else:
            _t, _point, distance = session.geometry(entity.name).project(tuple(xyz))
            score = distance
        if distance <= tolerance and score < best_score:
            best_name = entity.name
            best_score = score
    return best_name


def _polyline(points: np.ndarray) -> pv.PolyData:
    n = len(points)
    lines = np.concatenate([[n], np.arange(n)])
    return pv.PolyData(points, lines=lines)


def _display_point_to_world(renderer, x, y, z):
    renderer.SetDisplayPoint(x, y, z)
    renderer.DisplayToWorld()
    wx, wy, wz, w = renderer.GetWorldPoint()
    if w == 0.0:
        return (wx, wy, wz)
    return (wx / w, wy / w, wz / w)
