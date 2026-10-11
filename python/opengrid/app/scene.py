"""Viewport scene: renders model entities with pyvista and maps screen positions to the model.

Picking is done in world space with the core's closest-point projection (not VTK prop
picking), so it is exact for curves, independent of tessellation, and unit-testable without a
render window.
"""

import math

import numpy as np
import pyvista as pv

from opengrid.api import Session
from opengrid.model.connectors import CURVE_KINDS

CURVE_COLOR = "#5ad1ff"
POINT_COLOR = "#f2f2f2"
SELECTED_COLOR = "#ffd60a"
PREVIEW_COLOR = "#ff9f43"
SNAP_COLOR = "#2ee59d"
ENDPOINT_COLOR = "#8a96a3"
CONNECTOR_COLOR = "#ff9f43"
NODE_COLOR = "#e8e8e8"
STRETCHING_CMAP = "coolwarm"
DOMAIN_COLOR = "#9fb3c8"
QUALITY_CMAP = "viridis"
WORST_COLOR = "#ff3b30"
TESSELLATION_TOLERANCE = 2e-4


class Scene:
    def __init__(self, plotter, session: Session):
        self.plotter = plotter
        self.session = session
        self.selected_ids: set[int] = set()
        self._actors: dict[int, list[str]] = {}  # entity id -> actor names
        self.color_by_stretching = True
        self.show_nodes = True
        self.quality_metric: str | None = None  # colour domains by this metric when set

    # --------------------------------------------------------------- rendering
    def rebuild(self) -> None:
        for names in list(self._actors.values()):
            for actor_name in names:
                self.plotter.remove_actor(actor_name, render=False)
        self._actors.clear()
        for entity in self.session.store:
            self._add_entity_actor(entity.id)
        self._update_endpoints()
        self.plotter.render()

    def update_entities(self, ids) -> None:
        # Connectors depend on curves and (through spacing propagation) on each other:
        # redraw all of them together with the changed entities.
        targets = set(ids)
        for entity in self.session.store:
            if entity.kind in ("connector", "domain"):
                targets.add(entity.id)
        for entity_id in sorted(targets):
            self._remove_entity_actor(entity_id)
            if entity_id in self.session.store:
                self._add_entity_actor(entity_id)
        self._update_endpoints()
        self.plotter.render()

    def set_selection(self, ids) -> None:
        changed = self.selected_ids.symmetric_difference(ids)
        self.selected_ids = set(ids)
        for entity_id in changed:
            if entity_id in self.session.store:
                self._remove_entity_actor(entity_id)
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

    def _remove_entity_actor(self, entity_id: int) -> None:
        for actor_name in self._actors.pop(entity_id, []):
            self.plotter.remove_actor(actor_name, render=False)

    def _add_connector_actor(self, entity) -> None:
        nodes = self.session.connector_nodes(entity.name)
        selected = entity.id in self.selected_ids
        line_name = f"entity-{entity.id}"
        mesh = _polyline(nodes.points)
        if selected or not self.color_by_stretching:
            color = CONNECTOR_COLOR
            if selected:
                color = SELECTED_COLOR
            self.plotter.add_mesh(
                mesh,
                color=color,
                line_width=4.0,
                name=line_name,
                pickable=False,
                reset_camera=False,
            )
        else:
            segments = pv.PolyData(nodes.points, lines=_segment_cells(len(nodes.points)))
            segments.cell_data["stretching"] = nodes.interval_ratios()
            threshold = max(self.session.stretching_threshold, 1.0001)
            self.plotter.add_mesh(
                segments,
                scalars="stretching",
                cmap=STRETCHING_CMAP,
                clim=(1.0, threshold),
                line_width=4.0,
                name=line_name,
                pickable=False,
                reset_camera=False,
                show_scalar_bar=False,
            )
        names = [line_name]
        if self.show_nodes:
            node_name = f"entity-{entity.id}-nodes"
            self.plotter.add_mesh(
                pv.PolyData(nodes.points),
                color=NODE_COLOR,
                point_size=5,
                render_points_as_spheres=True,
                name=node_name,
                pickable=False,
                reset_camera=False,
            )
            names.append(node_name)
        self._actors[entity.id] = names

    def _add_domain_actor(self, entity) -> None:
        grid = self.session.domain_grid(entity.name, compute=False)
        mesh = quad_mesh(grid)
        name = f"entity-{entity.id}"
        selected = entity.id in self.selected_ids
        if self.quality_metric is not None:
            values = _core_quality(grid)[self.quality_metric]
            mesh.cell_data[self.quality_metric] = values.reshape(-1)
            self.plotter.add_mesh(
                mesh,
                scalars=self.quality_metric,
                cmap=QUALITY_CMAP,
                show_edges=True,
                edge_color="#30363d",
                name=name,
                pickable=False,
                reset_camera=False,
                scalar_bar_args={"title": self.quality_metric.replace("_", " ")},
            )
        else:
            color = DOMAIN_COLOR
            width = 1.0
            if selected:
                color = SELECTED_COLOR
                width = 1.5
            self.plotter.add_mesh(
                mesh,
                style="wireframe",
                color=color,
                line_width=width,
                name=name,
                pickable=False,
                reset_camera=False,
            )
        self._actors[entity.id] = [name]

    def show_worst_cell(self, grid: np.ndarray | None, j: int = 0, i: int = 0) -> None:
        if grid is None:
            self.plotter.remove_actor("worst-cell", render=False)
            self.plotter.render()
            return
        corners = np.array(
            [grid[j, i], grid[j, i + 1], grid[j + 1, i + 1], grid[j + 1, i], grid[j, i]]
        )
        self.plotter.add_mesh(
            _polyline(corners),
            color=WORST_COLOR,
            line_width=5,
            name="worst-cell",
            pickable=False,
            reset_camera=False,
        )
        center = corners[:4].mean(axis=0)
        size = float(np.max(np.ptp(corners, axis=0)))
        camera = self.plotter.camera
        camera.focal_point = tuple(center)
        camera.position = (center[0], center[1], center[2] + 10.0 * max(size, 1e-12))
        camera.parallel_scale = max(10.0 * size, 1e-9)
        self.plotter.render()

    def _add_entity_actor(self, entity_id: int) -> None:
        entity = self.session.store.get(entity_id)
        actor_name = f"entity-{entity_id}"
        selected = entity_id in self.selected_ids
        if entity.kind == "connector":
            self._add_connector_actor(entity)
            return
        if entity.kind == "domain":
            self._add_domain_actor(entity)
            return
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
        self._actors[entity_id] = [actor_name]

    def _update_endpoints(self) -> None:
        ends = []
        for entity in self.session.store:
            if entity.kind in CURVE_KINDS:
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
    """Name of the entity closest to ``xyz`` within ``tolerance``.

    Ties are broken by priority: points, then connectors (drawn on top of their curves), then
    curves.
    """
    best_name = None
    best_score = math.inf
    for entity in session.store:
        if entity.kind == "point":
            # Points sit on top of curve ends: halve their distance so they win the tie.
            score = 0.5 * math.dist(xyz, entity.params["xyz"])
            distance = 2.0 * score
        elif entity.kind == "connector":
            distance = _distance_to_polyline(
                np.asarray(xyz, dtype=float), session.connector_nodes(entity.name).points
            )
            score = 0.75 * distance
        elif entity.kind == "domain":
            # Inside a domain counts as a hit, at the lowest priority.
            grid = session.domain_grid(entity.name, compute=False)
            if not _inside_loop(np.asarray(xyz, dtype=float), grid):
                continue
            distance = 0.0
            score = tolerance
        elif entity.kind in CURVE_KINDS:
            _t, _point, distance = session.geometry(entity.name).project(tuple(xyz))
            score = distance
        else:
            continue
        if distance <= tolerance and score < best_score:
            best_name = entity.name
            best_score = score
    return best_name


def quad_mesh(grid: np.ndarray) -> pv.PolyData:
    """Quads of a structured grid (nj, ni, 3) as PolyData (cell order: i fastest)."""
    nj, ni, _ = grid.shape
    index = np.arange(ni * nj).reshape(nj, ni)
    quads = np.stack(
        [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
    ).reshape(-1, 4)
    faces = np.hstack([np.full((len(quads), 1), 4), quads]).ravel()
    return pv.PolyData(grid.reshape(-1, 3), faces=faces)


def _core_quality(grid):
    from opengrid import _core

    return _core.quality_2d(grid)


def _distance_to_polyline(point: np.ndarray, vertices: np.ndarray) -> float:
    a = vertices[:-1]
    b = vertices[1:]
    ab = b - a
    lengths = np.einsum("ij,ij->i", ab, ab)
    t = np.zeros(len(a))
    nonzero = lengths > 0.0
    t[nonzero] = np.clip(
        np.einsum("ij,ij->i", point - a[nonzero], ab[nonzero]) / lengths[nonzero], 0.0, 1.0
    )
    closest = a + t[:, None] * ab
    return float(np.min(np.linalg.norm(closest - point, axis=1)))


def _inside_loop(point: np.ndarray, grid: np.ndarray) -> bool:
    """Even-odd test against the boundary loop of a structured grid (xy plane)."""
    loop = np.concatenate([grid[0, :-1], grid[:-1, -1], grid[-1, :0:-1], grid[:0:-1, 0]])
    x = loop[:, 0]
    y = loop[:, 1]
    xn = np.roll(x, -1)
    yn = np.roll(y, -1)
    crosses = (y > point[1]) != (yn > point[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        at = x + (point[1] - y) * (xn - x) / (yn - y)
    return bool(np.count_nonzero(crosses & (point[0] < at)) % 2)


def _segment_cells(n: int) -> np.ndarray:
    """VTK cell array of n - 1 two-point line cells (one per interval)."""
    cells = np.empty((n - 1, 3), dtype=np.int64)
    cells[:, 0] = 2
    cells[:, 1] = np.arange(n - 1)
    cells[:, 2] = np.arange(1, n)
    return cells.ravel()


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
