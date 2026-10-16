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
from opengrid.model.entities import SURFACE_KINDS
from opengrid.validation import ValidationError

CURVE_COLOR = "#5ad1ff"
POINT_COLOR = "#f2f2f2"
SELECTED_COLOR = "#ffd60a"
PREVIEW_COLOR = "#ff9f43"
SNAP_COLOR = "#2ee59d"
HANDLE_COLOR = "#ffb000"
ENDPOINT_COLOR = "#8a96a3"
CONNECTOR_COLOR = "#ff9f43"
NODE_COLOR = "#e8e8e8"
STRETCHING_CMAP = "coolwarm"
DOMAIN_COLOR = "#9fb3c8"
BLOCK_COLOR = "#7fd1b9"
SURFACE_COLOR = "#8a7fd1"
QUALITY_CMAP = "viridis"
WORST_COLOR = "#ff3b30"
TESSELLATION_TOLERANCE = 2e-4


def volume_grid(mesh) -> pv.UnstructuredGrid:
    """A volume mesh (opengrid.model.volumes) as a VTK unstructured grid, cells in the canonical
    order (so cell indices match the quality arrays)."""
    from opengrid.model.volumes import volume_cells

    cells = []
    types = []
    for _key, array, vtk_type in volume_cells(mesh):
        if len(array):
            cells.append(np.hstack([np.full((len(array), 1), array.shape[1]), array]).ravel())
            types.append(np.full(len(array), vtk_type, dtype=np.uint8))
    return pv.UnstructuredGrid(np.concatenate(cells), np.concatenate(types), mesh["points"])


class Scene:
    def __init__(self, plotter, session: Session):
        self.plotter = plotter
        self.session = session
        self.selected_ids: set[int] = set()
        self._actors: dict[int, list[str]] = {}  # entity id -> actor names
        self.color_by_stretching = True
        self.show_nodes = True
        self.quality_metric: str | None = None  # colour domains by this metric when set
        self.patch_colors: dict[int, str] | None = None  # connector id -> colour (patch mode)
        # 3D views of blocks: keep the cells on one side of a plane (origin, normal), and/or
        # shrink every cell towards its centre.
        self.cut_plane: tuple | None = None
        self.shrink = False

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
            if entity.kind in ("connector", "domain", "unstructured_domain", "block", "volume"):
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

    def show_handles(self, points) -> None:
        """Drag handles (control points of the selected entities)."""
        if len(points) == 0:
            self.plotter.remove_actor("drag-handles", render=False)
        else:
            self.plotter.add_mesh(
                pv.PolyData(np.asarray(points, dtype=float)),
                color=HANDLE_COLOR,
                point_size=11,
                render_points_as_spheres=True,
                name="drag-handles",
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
        nodes = self.session.display_connector_nodes(entity.name)
        selected = entity.id in self.selected_ids
        line_name = f"entity-{entity.id}"
        mesh = _polyline(nodes.points)
        if selected or not self.color_by_stretching or self.patch_colors is not None:
            color = CONNECTOR_COLOR
            if self.patch_colors is not None:
                color = self.patch_colors.get(entity.id, CONNECTOR_COLOR)
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
        if entity.kind == "unstructured_domain":
            unstructured = self.session.unstructured_mesh(entity.name)
            mesh = triangle_mesh(
                unstructured["points"], unstructured["triangles"], unstructured["quads"]
            )
        else:
            mesh = quad_mesh(self.session.domain_grid(entity.name, compute=False))
        name = f"entity-{entity.id}"
        selected = entity.id in self.selected_ids
        # In a 3D model the cells are the blocks: domains (their faces) stay wireframes.
        if self.quality_metric is not None and not self.session.is_3d():
            values = self._display_quality(entity)[self.quality_metric]
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
            opacity = 1.0
            if self.patch_colors is not None:
                opacity = 0.25  # let the patch colours of the boundary connectors stand out
            self.plotter.add_mesh(
                mesh,
                style="wireframe",
                color=color,
                line_width=width,
                opacity=opacity,
                name=name,
                pickable=False,
                reset_camera=False,
            )
        self._actors[entity.id] = [name]

    def _add_surface_actor(self, entity) -> None:
        grid = self.session.store.geometry(entity.id).grid(25, 25)
        name = f"entity-{entity.id}"
        color = SURFACE_COLOR
        if entity.id in self.selected_ids:
            color = SELECTED_COLOR
        self.plotter.add_mesh(
            quad_mesh(np.asarray(grid)),
            color=color,
            opacity=0.25,
            name=name,
            pickable=False,
            reset_camera=False,
        )
        self._actors[entity.id] = [name]

    def block_mesh(self, entity) -> pv.DataSet:
        """The block (structured grid) or volume (unstructured grid), with the cut plane
        and shrink applied."""
        if entity.kind == "volume":
            mesh = volume_grid(self.session.volume_mesh(entity.name))
        else:
            grid = np.asarray(self.session.block_grid(entity.name))
            mesh = pv.StructuredGrid(
                np.ascontiguousarray(grid[..., 0].transpose(2, 1, 0)),
                np.ascontiguousarray(grid[..., 1].transpose(2, 1, 0)),
                np.ascontiguousarray(grid[..., 2].transpose(2, 1, 0)),
            )
        if self.quality_metric is not None:
            values = self._display_quality(entity)[self.quality_metric]
            mesh.cell_data[self.quality_metric] = values.reshape(-1)
        if self.cut_plane is not None:
            origin, normal = self.cut_plane
            mesh = mesh.clip(normal=normal, origin=origin, crinkle=True)
        if self.shrink:
            mesh = mesh.shrink(0.8)
        return mesh

    def _add_block_actor(self, entity) -> None:
        try:
            mesh = self.block_mesh(entity)
        except ValidationError:
            return  # an invalid block is reported by the property panel / boundary checks
        name = f"entity-{entity.id}"
        if mesh.n_cells == 0:
            self._actors[entity.id] = []
            return
        cells = self.cut_plane is not None or self.shrink
        if self.quality_metric is not None:
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
        elif cells:
            color = BLOCK_COLOR
            if entity.id in self.selected_ids:
                color = SELECTED_COLOR
            self.plotter.add_mesh(
                mesh,
                color=color,
                show_edges=True,
                edge_color="#30363d",
                name=name,
                pickable=False,
                reset_camera=False,
            )
        else:
            color = BLOCK_COLOR
            if entity.id in self.selected_ids:
                color = SELECTED_COLOR
            self.plotter.add_mesh(
                mesh.extract_surface(algorithm="dataset_surface"),
                style="wireframe",
                color=color,
                line_width=1.0,
                name=name,
                pickable=False,
                reset_camera=False,
            )
        self._actors[entity.id] = [name]

    def _display_quality(self, entity) -> dict:
        """Quality of the displayed mesh (the TFI grid while smoothing is pending); for blocks
        the 3D metrics, with area / area_ratio standing for volume / volume_ratio."""
        from opengrid import _core

        if entity.kind in ("block", "volume"):
            if entity.kind == "volume":
                q = dict(self.session.volume_quality(entity.name))
            else:
                q = dict(self.session.block_quality(entity.name))
            q["area"] = q["volume"]
            q["area_ratio"] = q["volume_ratio"]
            return q
        if entity.kind == "unstructured_domain":
            return self.session.domain_quality(entity.name)
        return _core.quality_2d(self.session.domain_grid(entity.name, compute=False))

    def show_worst_cell(self, corners: np.ndarray | None) -> None:
        """Outline a cell given its corner points (or remove the outline)."""
        if corners is None:
            self.plotter.remove_actor("worst-cell", render=False)
            self.plotter.render()
            return
        corners = np.vstack([corners, corners[:1]])
        self.plotter.add_mesh(
            _polyline(corners),
            color=WORST_COLOR,
            line_width=5,
            name="worst-cell",
            pickable=False,
            reset_camera=False,
        )
        center = corners[:-1].mean(axis=0)
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
        if entity.kind in ("domain", "unstructured_domain"):
            self._add_domain_actor(entity)
            return
        if entity.kind in ("block", "volume"):
            self._add_block_actor(entity)
            return
        if entity.kind in SURFACE_KINDS:
            self._add_surface_actor(entity)
            return
        if entity.kind not in CURVE_KINDS and entity.kind != "point":
            return  # patches and zones are shown through their members
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
        elif entity.kind == "unstructured_domain":
            mesh = session.unstructured_mesh(entity.name)
            # Quads are convex: test them as two triangles each.
            quads = mesh["quads"]
            pieces = np.concatenate([mesh["triangles"], quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
            if not _inside_triangles(np.asarray(xyz, dtype=float), mesh["points"], pieces):
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


def triangle_mesh(points: np.ndarray, triangles: np.ndarray, quads=None) -> pv.PolyData:
    """Triangles, then quads (optional), as PolyData; the cell order matches the quality arrays."""
    faces = [np.hstack([np.full((len(triangles), 1), 3), triangles]).ravel()]
    if quads is not None and len(quads):
        faces.append(np.hstack([np.full((len(quads), 1), 4), quads]).ravel())
    return pv.PolyData(np.asarray(points, dtype=float), faces=np.concatenate(faces))


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


def _inside_triangles(point: np.ndarray, points: np.ndarray, triangles: np.ndarray) -> bool:
    a = points[triangles[:, 0], :2]
    b = points[triangles[:, 1], :2]
    c = points[triangles[:, 2], :2]
    p = point[:2]

    def side(o, v):
        return (v[:, 0] - o[:, 0]) * (p[1] - o[:, 1]) - (v[:, 1] - o[:, 1]) * (p[0] - o[:, 0])

    return bool(np.any((side(a, b) >= 0) & (side(b, c) >= 0) & (side(c, a) >= 0)))


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
