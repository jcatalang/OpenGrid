"""Viewport scene: renders model entities with pyvista and maps screen positions to the model.

Drawing strategy. Curves, points and connectors are many and small: each kind is one merged
actor (plus one for the selected ones, one for the connector nodes and one for the curve end
points), rebuilt from numpy arrays when any of them changes, so a model edit costs a handful
of actor replacements however many entities there are. Domains, blocks, volumes, surfaces and
sketches are few and large: one actor each, replaced only when the entity or something it
depends on changed (the dependents closure of the change) and recoloured in place when
selected. Actors are removed by object, never by name: pyvista's name lookup walks the
renderer's whole prop list.

Picking is done in world space with the core's closest-point projection (not VTK prop
picking), so it is exact for curves, independent of tessellation, and unit-testable without a
render window.
"""

import math

import numpy as np
import pyvista as pv

from opengrid.api import Session
from opengrid.app.reference_grid import ReferenceGrid
from opengrid.model.connectors import CURVE_KINDS
from opengrid.model.entities import SURFACE_KINDS
from opengrid.validation import ValidationError

CURVE_COLOR = "#5ad1ff"
POINT_COLOR = "#f2f2f2"
SELECTED_COLOR = "#ffd60a"
PREVIEW_COLOR = "#ff9f43"
SNAP_COLOR = "#2ee59d"
HANDLE_COLOR = "#ffb000"
HIGH_ORDER_COLOR = "#ff9f1c"  # curved element edges
HIGH_ORDER_NODE_COLOR = "#ffd166"
ENDPOINT_COLOR = "#8a96a3"
CONNECTOR_COLOR = "#ff9f43"
NODE_COLOR = "#e8e8e8"
STRETCHING_CMAP = "coolwarm"
DOMAIN_COLOR = "#9fb3c8"
BLOCK_COLOR = "#7fd1b9"
SURFACE_COLOR = "#8a7fd1"
QUALITY_CMAP = "viridis"
WORST_COLOR = "#ff3b30"
SKETCH_COLOR = "#ff6fd8"
HOVER_COLOR = "#ffffff"
TESSELLATION_TOLERANCE = 2e-4
PICK_PIXELS = 8.0  # ray picking: within this many pixels of a line or point
# A selected curve is drawn bold (in the selection colour, or in its patch colour when the
# curves are coloured by patch: the boundary types stay visible).
SELECTED_CURVE_WIDTH = 9.0
FACE_COLOR = SELECTED_COLOR  # a selected block / volume face
POINTS_PREVIEW_COLOR = "#ff79c6"  # where the points of curves without points would go

# Entity kinds drawn as one actor each; the others (points, curves, connectors) share layers.
ENTITY_KINDS = ("domain", "unstructured_domain", "block", "volume", "sketch", *SURFACE_KINDS)
# Plain and selected (colour, line width) of the per-entity kinds.
STYLES = {
    "domain": ((DOMAIN_COLOR, 1.0), (SELECTED_COLOR, 1.5)),
    "unstructured_domain": ((DOMAIN_COLOR, 1.0), (SELECTED_COLOR, 1.5)),
    "block": ((BLOCK_COLOR, 1.0), (SELECTED_COLOR, 1.0)),
    "volume": ((BLOCK_COLOR, 1.0), (SELECTED_COLOR, 1.0)),
    "sketch": ((SKETCH_COLOR, 2.0), (SELECTED_COLOR, 3.0)),
}
for _kind in SURFACE_KINDS:
    STYLES[_kind] = ((SURFACE_COLOR, 1.0), (SELECTED_COLOR, 1.0))
LAYER_GROUPS = {
    "curves": ("curves", "curves-selected", "curve-endpoints"),
    "points": ("points", "points-selected"),
    "connectors": ("connectors", "connectors-selected", "connector-nodes"),
}


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
        self._entity_actors: dict[int, list] = {}  # entity id -> its actors (per-entity kinds)
        self._layers: dict[str, object] = {}  # merged layer name -> actor
        self._layer_data: dict[str, object] = {}  # merged layer name -> the mesh it draws
        self._overlays: dict[str, object] = {}  # tool preview, snap marker, handles, ...
        self._hidden: list = []  # per-entity actors hidden by a drag preview
        self._previewing: set[int] = set()  # entities left out of the layers during a drag
        self.hidden: set[int] = set()  # entities not drawn (View > Hide selected)
        self.color_by_stretching = True
        self.show_nodes = True
        self.quality_metric: str | None = None  # colour domains by this metric when set
        self.patch_colors: dict[int, str] | None = None  # connector id -> colour (patch mode)
        # 3D views of blocks: keep the cells on one side of a plane (origin, normal), and/or
        # shrink every cell towards its centre.
        self.cut_plane: tuple | None = None
        self.shrink = False
        # The view "fit" starts from: "xy" (2D navigation) or "iso" (3D).
        self.default_view = "xy"
        self.hovered: int | None = None  # entity highlighted under the pointer
        self.grid = ReferenceGrid(plotter)  # the z = 0 plane's grid, when enabled
        # Faces of blocks and volumes: the patch member of each face of their drawn surface
        # (cell data "face" indexes it), the selected ones, the one the last pick hit.
        self._face_members: dict[int, list[dict]] = {}
        self.selected_faces: list[dict] = []
        self.last_face: dict | None = None
        self.patch_faces: list | None = None  # [(face member, colour)] while colouring by patch

    # ---------------------------------------------------------- actor plumbing
    def _add(self, mesh, name: str, pickable: bool = True, **kwargs):
        """Add a mesh without pyvista's name lookup (O(props)); the name is for reading."""
        actor = self.plotter.add_mesh(
            mesh, pickable=pickable, reset_camera=False, render=False, **kwargs
        )
        actor.name = name
        return actor

    def _remove(self, actor) -> None:
        self.plotter.remove_actor(actor, render=False)  # by object: no name lookup

    def _set_layer(self, name: str, mesh, **kwargs) -> None:
        actor = self._layers.pop(name, None)
        self._layer_data.pop(name, None)
        if actor is not None:
            self._remove(actor)
        if mesh is None or mesh.n_points == 0:
            return
        self._layers[name] = self._add(mesh, name, **kwargs)
        self._layer_data[name] = mesh

    def _set_overlay(self, name: str, mesh, **kwargs) -> None:
        actor = self._overlays.pop(name, None)
        if actor is not None:
            self._remove(actor)
        if mesh is None:
            return
        self._overlays[name] = self._add(mesh, name, pickable=False, **kwargs)

    def _remove_entity(self, entity_id: int) -> None:
        for actor in self._entity_actors.pop(entity_id, []):
            self._remove(actor)
        self._face_members.pop(entity_id, None)

    def _refresh_faces(self) -> None:
        """Redraw the selected faces and the patch colours on the redrawn blocks (dropping the
        faces that went)."""
        if self.patch_faces is not None:
            self.show_patch_faces(self.patch_faces)
        if self.selected_faces:
            kept = [m for m in self.selected_faces if self.face_surface(m) is not None]
            self.set_face_selection(kept)

    def layer(self, name: str):
        """The actor of a merged layer ("connectors", "connector-nodes", ...), or None."""
        return self._layers.get(name)

    def layer_data(self, name: str):
        """The mesh a merged layer draws (with its cell data), or None."""
        return self._layer_data.get(name)

    def has_geometry(self) -> bool:
        """Whether anything is drawn (a drag preview counts: its entities are taken out of
        the layers and hidden while it lasts, which must not look like an empty scene)."""
        return bool(self._entity_actors) or bool(self._layers) or bool(self._previewing)

    # --------------------------------------------------------------- rendering
    def rebuild(self) -> None:
        for entity_id in list(self._entity_actors):
            self._remove_entity(entity_id)
        for name in list(self._layers):
            self._remove(self._layers.pop(name))
            self._layer_data.pop(name, None)
        for entity in self.session.store:
            if entity.kind in ENTITY_KINDS:
                self._draw_entity(entity.id)
        self._rebuild_layers(set(LAYER_GROUPS))
        self._refresh_faces()
        self.plotter.render()

    def _closure(self, ids) -> set[int]:
        """The given entities and everything that depends on them, transitively (ids no longer
        in the store are kept: their actors go)."""
        store = self.session.store
        out: set[int] = set()
        stack = list(ids)
        while stack:
            entity_id = stack.pop()
            if entity_id in out:
                continue
            out.add(entity_id)
            if entity_id in store:
                stack.extend(dependent.id for dependent in store.dependents(entity_id))
        return out

    def update_entities(self, ids) -> None:
        """Redraw the changed entities and their dependents. Connectors are redrawn whatever
        changed: they follow their curves, each other (spacing propagation) and the smoothing
        groups that slide them."""
        store = self.session.store
        layers = {"connectors"}
        for entity_id in sorted(self._closure(ids)):
            if entity_id not in store:
                self._remove_entity(entity_id)
                layers.update(("curves", "points"))  # gone: whatever it was
                continue
            kind = store.get(entity_id).kind
            if kind in CURVE_KINDS:
                layers.add("curves")
            elif kind == "point":
                layers.add("points")
            elif kind in ENTITY_KINDS:
                self._draw_entity(entity_id)
        self._rebuild_layers(layers)
        self._refresh_faces()
        self.plotter.render()

    def _rebuild_layers(self, groups) -> None:
        if "curves" in groups:
            self._curve_layers()
        if "points" in groups:
            self._point_layers()
        if "connectors" in groups:
            self._connector_layers()

    def _curve_layers(self) -> None:
        plain = []
        plain_ids = []
        selected = []
        selected_ids = []
        ends = []
        for entity in self.session.store:
            if entity.kind not in CURVE_KINDS or self._skipped(entity.id):
                continue
            curve = self.session.geometry(entity.name)
            points = np.asarray(curve.tessellate_points(TESSELLATION_TOLERANCE), dtype=float)
            if entity.id in self.selected_ids:
                selected.append(points)
                selected_ids.append(entity.id)
            else:
                plain.append(points)
                plain_ids.append(entity.id)
            ends.append(curve.eval(0.0))
            ends.append(curve.eval(1.0))
        self._set_layer("curves", _polylines(plain, plain_ids), color=CURVE_COLOR, line_width=2.5)
        self._set_layer(
            "curves-selected",
            _polylines(selected, selected_ids),
            color=SELECTED_COLOR,
            line_width=SELECTED_CURVE_WIDTH,
        )
        self._set_layer(
            "curve-endpoints",
            _points(ends),
            pickable=False,
            color=ENDPOINT_COLOR,
            point_size=7,
            render_points_as_spheres=True,
        )

    def _point_layers(self) -> None:
        plain = []
        plain_ids = []
        selected = []
        selected_ids = []
        for entity in self.session.store:
            if entity.kind != "point" or self._skipped(entity.id):
                continue
            if entity.id in self.selected_ids:
                selected.append(entity.params["xyz"])
                selected_ids.append(entity.id)
            else:
                plain.append(entity.params["xyz"])
                plain_ids.append(entity.id)
        for name, points, ids, color in (
            ("points", plain, plain_ids, POINT_COLOR),
            ("points-selected", selected, selected_ids, SELECTED_COLOR),
        ):
            self._set_layer(
                name,
                _points(points, ids),
                color=color,
                point_size=12,
                render_points_as_spheres=True,
            )

    def _connector_layers(self) -> None:
        plain = []
        plain_ids = []
        ratios = []
        colors = []
        selected = []
        selected_ids = []
        nodes_all = []
        by_stretching = self.color_by_stretching and self.patch_colors is None
        for entity in self.session.store:
            if entity.kind != "connector" or self._skipped(entity.id):
                continue
            nodes = self.session.display_connector_nodes(entity.name)
            points = np.asarray(nodes.points, dtype=float)
            nodes_all.append(points)
            if entity.id in self.selected_ids:
                selected.append(points)
                selected_ids.append(entity.id)
                continue
            plain.append(points)
            plain_ids.append(entity.id)
            if by_stretching:
                ratios.append(np.asarray(nodes.interval_ratios(), dtype=float))
            else:
                color = CONNECTOR_COLOR
                if self.patch_colors is not None:
                    color = self.patch_colors.get(entity.id, CONNECTOR_COLOR)
                colors.append(pv.Color(color).int_rgb)
        if by_stretching:
            mesh = _segments(plain, plain_ids)
            if mesh is not None:
                mesh.cell_data["stretching"] = np.concatenate(ratios)
            threshold = max(self.session.stretching_threshold, 1.0001)
            self._set_layer(
                "connectors",
                mesh,
                scalars="stretching",
                cmap=STRETCHING_CMAP,
                clim=(1.0, threshold),
                line_width=4.0,
                show_scalar_bar=False,
            )
        else:
            mesh = _polylines(plain, plain_ids)
            if mesh is not None:
                mesh.cell_data["colors"] = np.asarray(colors, dtype=np.uint8)
            self._set_layer("connectors", mesh, scalars="colors", rgb=True, line_width=4.0)
        mesh = _polylines(selected, selected_ids)
        if self.patch_colors is not None and mesh is not None:
            selected_colors = [
                pv.Color(self.patch_colors.get(i, CONNECTOR_COLOR)).int_rgb for i in selected_ids
            ]
            mesh.cell_data["colors"] = np.asarray(selected_colors, dtype=np.uint8)
            self._set_layer(
                "connectors-selected",
                mesh,
                scalars="colors",
                rgb=True,
                line_width=SELECTED_CURVE_WIDTH,
            )
        else:
            self._set_layer(
                "connectors-selected", mesh, color=SELECTED_COLOR, line_width=SELECTED_CURVE_WIDTH
            )
        nodes_mesh = None
        if self.show_nodes:
            nodes_mesh = _points(nodes_all)
        self._set_layer(
            "connector-nodes",
            nodes_mesh,
            pickable=False,
            color=NODE_COLOR,
            point_size=5,
            render_points_as_spheres=True,
        )

    def preview_entities(self, ids) -> None:
        """Drag preview: the given entities (domains, connectors, curves, sketches) drawn as
        one wireframe actor while their own drawing is hidden (cheap to redraw every frame);
        end_preview restores them."""
        ids = {i for i in ids if i in self.session.store}
        if ids != self._previewing:
            # First frame of this drag: leave the entities out of the merged layers and hide
            # their own actors.
            self._previewing = ids
            self._rebuild_layers(set(LAYER_GROUPS))
            for entity_id in sorted(ids):
                for actor in self._entity_actors.get(entity_id, []):
                    actor.SetVisibility(False)
                    self._hidden.append(actor)
        pieces = []
        for entity_id in sorted(ids):
            entity = self.session.store.get(entity_id)
            if entity.kind == "domain":
                pieces.append(quad_mesh(self.session.domain_grid(entity.name, compute=False)))
            elif entity.kind == "connector":
                pieces.append(_polyline(self.session.connector_nodes(entity.name).points))
            elif entity.kind in CURVE_KINDS:
                curve = self.session.geometry(entity.name)
                pieces.append(_polyline(curve.tessellate_points(TESSELLATION_TOLERANCE)))
            elif entity.kind == "sketch":
                pieces.append(_sketch_lines(entity))
        mesh = None
        if pieces:
            mesh = pv.merge(pieces)
        self._set_overlay(
            "drag-preview", mesh, style="wireframe", color=PREVIEW_COLOR, line_width=1.0
        )
        self.plotter.render()

    def end_preview(self) -> None:
        self._set_overlay("drag-preview", None)
        for actor in self._hidden:
            actor.SetVisibility(True)
        self._hidden.clear()
        if self._previewing:
            self._previewing = set()
            self._rebuild_layers(set(LAYER_GROUPS))

    def set_selection(self, ids) -> None:
        changed = self.selected_ids.symmetric_difference(ids)
        self.selected_ids = set(ids)
        store = self.session.store
        layers = set()
        for entity_id in changed:
            if entity_id not in store:
                continue
            kind = store.get(entity_id).kind
            if kind in CURVE_KINDS:
                layers.add("curves")
            elif kind == "point":
                layers.add("points")
            elif kind == "connector":
                layers.add("connectors")
            else:
                self._recolour(entity_id)
        self._rebuild_layers(layers)
        self.plotter.render()

    def _recolour(self, entity_id: int) -> None:
        """Selection colour of a per-entity actor, in place (actors coloured by a quality
        metric keep their scalars)."""
        entity = self.session.store.get(entity_id)
        if entity.kind not in STYLES or self._coloured_by_quality(entity):
            return
        color, width = STYLES[entity.kind][int(entity_id in self.selected_ids)]
        for actor in self._entity_actors.get(entity_id, []):
            actor.prop.color = color
            actor.prop.line_width = width

    def _coloured_by_quality(self, entity) -> bool:
        if self.quality_metric is None:
            return False
        if entity.kind in ("block", "volume"):
            return True
        return entity.kind in ("domain", "unstructured_domain") and not self.session.is_3d()

    def set_hover(self, entity_id: int | None) -> None:
        """Highlight the entity under the pointer: its curve, its nodes' polyline, a domain's
        boundary, a block's outline (none for the selected, or None)."""
        if entity_id is not None and (
            entity_id in self.selected_ids or entity_id not in self.session.store
        ):
            entity_id = None
        if entity_id == self.hovered:
            return
        self.hovered = entity_id
        mesh = None
        if entity_id is not None:
            mesh = self._hover_mesh(self.session.store.get(entity_id))
        self._set_overlay(
            "hover",
            mesh,
            color=HOVER_COLOR,
            line_width=3.0,
            point_size=14,
            render_points_as_spheres=True,
        )
        self.plotter.render()

    def _hover_mesh(self, entity):
        if entity.kind == "point":
            return pv.PolyData(np.asarray([entity.params["xyz"]], dtype=float))
        points = self._entity_points(entity)
        if points is not None:
            return _polyline(points)
        if entity.kind == "domain":
            grid = self.session.domain_grid(entity.name, compute=False)
            loop = np.concatenate([grid[0, :-1], grid[:-1, -1], grid[-1, :0:-1], grid[:0:-1, 0]])
            return _polyline(np.vstack([loop, loop[:1]]))
        if entity.kind == "unstructured_domain":
            mesh = self.session.unstructured_mesh(entity.name)
            boundary = triangle_mesh(mesh["points"], mesh["triangles"], mesh["quads"])
            return boundary.extract_feature_edges(
                boundary_edges=True,
                feature_edges=False,
                manifold_edges=False,
                non_manifold_edges=False,
            )
        actors = self._entity_actors.get(entity.id, [])
        if actors:
            return pv.Box(actors[0].GetBounds()).extract_all_edges()
        return None

    def show_preview(self, points) -> None:
        mesh = None
        if len(points) >= 2:
            mesh = _polyline(np.asarray(points, dtype=float))
        self._set_overlay("tool-preview", mesh, color=PREVIEW_COLOR, line_width=2)
        self.plotter.render()

    def show_snap(self, xyz) -> None:
        mesh = None
        if xyz is not None:
            mesh = pv.PolyData(np.asarray([xyz], dtype=float))
        self._set_overlay(
            "snap-marker", mesh, color=SNAP_COLOR, point_size=16, render_points_as_spheres=True
        )
        self.plotter.render()

    def show_handles(self, points) -> None:
        """Drag handles (control points of the selected entities)."""
        mesh = None
        if len(points):
            mesh = pv.PolyData(np.asarray(points, dtype=float))
        self._set_overlay(
            "drag-handles", mesh, color=HANDLE_COLOR, point_size=11, render_points_as_spheres=True
        )
        self.plotter.render()

    def show_high_order(self, mesh) -> None:
        """Overlay of a high-order mesh (app/high_order_view.py), or none."""
        from opengrid.app.high_order_view import edge_polydata, node_polydata

        edges = None
        nodes = None
        if mesh is not None:
            edges = edge_polydata(mesh)
            nodes = node_polydata(mesh)
        self._set_overlay("high-order-edges", edges, color=HIGH_ORDER_COLOR, line_width=1.5)
        self._set_overlay(
            "high-order-nodes",
            nodes,
            color=HIGH_ORDER_NODE_COLOR,
            point_size=5,
            render_points_as_spheres=True,
        )
        self.plotter.render()

    # ------------------------------------------------------------------ camera
    def fit(self) -> None:
        """The default view (top in 2D navigation, isometric in 3D), centred on everything;
        never the flipped side (only the view keys toggle)."""
        self.view_preset(self.default_view, negative=False)

    def camera_changed(self) -> None:
        """After the camera moved: the reference grid follows, then a render."""
        if self.grid.enabled:
            self.grid.update()  # renders
        else:
            self.plotter.render()

    def reset_camera(self) -> None:
        """Centre the view on everything, keeping the viewing direction."""
        self.plotter.reset_camera()
        self.camera_changed()

    def view_preset(self, name: str, negative: bool | None = None) -> None:
        """Look at the xy plane from +z ("xy", top), the xz plane from -y ("xz", front), the
        yz plane from +x ("yz", side) or from the isometric corner ("iso"), centred on
        everything; ``negative`` looks from the opposite side, and None flips to it when the
        view already looks from the preset's side (press the key again to turn around)."""
        views = {
            "xy": (self.plotter.view_xy, (0.0, 0.0, -1.0)),
            "xz": (self.plotter.view_xz, (0.0, 1.0, 0.0)),
            "yz": (self.plotter.view_yz, (-1.0, 0.0, 0.0)),
            "iso": (self.plotter.view_isometric, None),
        }
        if name not in views:
            raise ValueError(f"unknown view {name!r} (xy, xz, yz or iso)")
        look, direction = views[name]
        if negative is None:
            negative = False
            if direction is not None:
                current = self.plotter.renderer.GetActiveCamera().GetDirectionOfProjection()
                negative = bool(np.allclose(current, direction, atol=1e-9))
        look(negative=negative)
        self.plotter.reset_camera()
        self.camera_changed()

    def view_name(self) -> str | None:
        """ "xy", "-xy", "xz", "-xz", "yz", "-yz" when the view looks along an axis, else None."""
        current = np.asarray(self.plotter.renderer.GetActiveCamera().GetDirectionOfProjection())
        for name, direction in (
            ("xy", (0.0, 0.0, -1.0)),
            ("xz", (0.0, 1.0, 0.0)),
            ("yz", (-1.0, 0.0, 0.0)),
        ):
            if np.allclose(current, direction, atol=1e-9):
                return name
            if np.allclose(current, -np.asarray(direction), atol=1e-9):
                return "-" + name
        return None

    def set_perspective(self, enabled: bool) -> None:
        if enabled:
            self.plotter.disable_parallel_projection()
        else:
            self.plotter.enable_parallel_projection()
        self.camera_changed()

    def zoom_to(self, ids) -> bool:
        """Centre the view on the given entities (keeping the viewing direction); False when
        none of them is drawn."""
        bounds = self.entity_bounds(ids)
        if bounds is None:
            return False
        self.plotter.reset_camera(bounds=bounds)
        self.camera_changed()
        return True

    def entity_bounds(self, ids) -> tuple | None:
        """(xmin, xmax, ymin, ymax, zmin, zmax) of these entities (patches and zones: their
        members), padded so a single point still has an extent."""
        store = self.session.store
        wanted = set()
        for entity_id in ids:
            if entity_id not in store:
                continue
            entity = store.get(entity_id)
            if entity.kind in ("patch", "zone"):
                for item in entity.params["members"]:
                    for key in ("connector", "domain", "block", "volume"):
                        if key in item:
                            wanted.add(item[key])
            else:
                wanted.add(entity_id)
        boxes = []
        for entity_id in sorted(wanted):
            if entity_id not in store:
                continue
            entity = store.get(entity_id)
            for actor in self._entity_actors.get(entity_id, []):
                boxes.append(np.array(actor.GetBounds(), dtype=float))
            points = self._entity_points(entity)
            if points is not None and len(points):
                low = points.min(axis=0)
                high = points.max(axis=0)
                boxes.append(np.array([low[0], high[0], low[1], high[1], low[2], high[2]]))
        if not boxes:
            return None
        stack = np.array(boxes)
        box = [
            stack[:, 0].min(),
            stack[:, 1].max(),
            stack[:, 2].min(),
            stack[:, 3].max(),
            stack[:, 4].min(),
            stack[:, 5].max(),
        ]
        size = float(max(box[1] - box[0], box[3] - box[2], box[5] - box[4]))
        pad = 0.05 * size
        if size == 0.0:
            pad = 0.5
        padded = (
            box[0] - pad,
            box[1] + pad,
            box[2] - pad,
            box[3] + pad,
            box[4] - pad,
            box[5] + pad,
        )
        return tuple(float(v) for v in padded)

    def _entity_points(self, entity) -> np.ndarray | None:
        """The points of an entity drawn in a merged layer (None for the per-entity kinds)."""
        if entity.kind == "point":
            return np.asarray([entity.params["xyz"]], dtype=float)
        if entity.kind in CURVE_KINDS:
            curve = self.session.geometry(entity.name)
            return np.asarray(curve.tessellate_points(TESSELLATION_TOLERANCE), dtype=float)
        if entity.kind == "connector":
            return np.asarray(self.session.display_connector_nodes(entity.name).points)
        return None

    # ----------------------------------------------------------- entity actors
    def _skipped(self, entity_id: int) -> bool:
        """Left out of the drawing: hidden, or taken by a drag preview."""
        return entity_id in self.hidden or entity_id in self._previewing

    def set_hidden(self, ids) -> None:
        """Hide these entities (an empty set shows everything again); redraws."""
        self.hidden = {int(i) for i in ids}
        self.rebuild()

    def _draw_entity(self, entity_id: int) -> None:
        self._remove_entity(entity_id)
        if entity_id in self.hidden:
            return
        entity = self.session.store.get(entity_id)
        if entity.kind in ("domain", "unstructured_domain"):
            self._add_domain_actor(entity)
        elif entity.kind in ("block", "volume"):
            self._add_block_actor(entity)
        elif entity.kind in SURFACE_KINDS:
            self._add_surface_actor(entity)
        elif entity.kind == "sketch":
            self._add_sketch_actor(entity)

    def _add_domain_actor(self, entity) -> None:
        if entity.kind == "unstructured_domain":
            unstructured = self.session.unstructured_mesh(entity.name)
            mesh = triangle_mesh(
                unstructured["points"], unstructured["triangles"], unstructured["quads"]
            )
        else:
            mesh = quad_mesh(self.session.domain_grid(entity.name, compute=False))
        name = f"entity-{entity.id}"
        # In a 3D model the cells are the blocks: domains (their faces) stay wireframes.
        if self._coloured_by_quality(entity):
            values = self._display_quality(entity)[self.quality_metric]
            mesh.cell_data[self.quality_metric] = values.reshape(-1)
            actor = self._add(
                mesh,
                name,
                scalars=self.quality_metric,
                cmap=QUALITY_CMAP,
                show_edges=True,
                edge_color="#30363d",
                scalar_bar_args={"title": self.quality_metric.replace("_", " ")},
            )
        else:
            color, width = STYLES[entity.kind][int(entity.id in self.selected_ids)]
            opacity = 1.0
            if self.patch_colors is not None:
                opacity = 0.25  # let the patch colours of the boundary connectors stand out
            actor = self._add(
                mesh, name, style="wireframe", color=color, line_width=width, opacity=opacity
            )
        self._entity_actors[entity.id] = [actor]

    def _add_surface_actor(self, entity) -> None:
        grid = self.session.store.geometry(entity.id).grid(25, 25)
        color, _width = STYLES[entity.kind][int(entity.id in self.selected_ids)]
        actor = self._add(
            quad_mesh(np.asarray(grid)), f"entity-{entity.id}", color=color, opacity=0.25
        )
        self._entity_actors[entity.id] = [actor]

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
            self._entity_actors[entity.id] = []
            return
        cells = self.cut_plane is not None or self.shrink
        color, width = STYLES[entity.kind][int(entity.id in self.selected_ids)]
        if self.quality_metric is not None:
            actor = self._add(
                mesh,
                name,
                scalars=self.quality_metric,
                cmap=QUALITY_CMAP,
                show_edges=True,
                edge_color="#30363d",
                scalar_bar_args={"title": self.quality_metric.replace("_", " ")},
            )
        elif cells:
            actor = self._add(mesh, name, color=color, show_edges=True, edge_color="#30363d")
        else:
            # Its boundary faces, each cell tagged with the face it is on: a click picks a
            # face (a patch member: a block's imin ... kmax, a volume's ends and side walls).
            surface, members = self._boundary_faces(entity)
            actor = self._add(surface, name, style="wireframe", color=color, line_width=width)
            self._face_members[entity.id] = members
        self._entity_actors[entity.id] = [actor]

    def _boundary_faces(self, entity) -> tuple[pv.PolyData, list[dict]]:
        """The drawn surface of a block or volume with cell data "face" (an index into the
        returned patch members)."""
        from opengrid.model.blocks import face_array
        from opengrid.model.patches import BLOCK_FACES

        if entity.kind == "block":
            grid = np.asarray(self.session.block_grid(entity.name))
            meshes = [quad_mesh(np.ascontiguousarray(face_array(grid, f))) for f in BLOCK_FACES]
            members = [{"block": entity.id, "face": f} for f in BLOCK_FACES]
        else:
            mesh = self.session.volume_mesh(entity.name)
            extruded = entity.params["faces"] is None
            meshes = []
            members = []
            for key, cells in mesh["faces"].items():
                triangles = [c for c in cells if len(c) == 3]
                quads = [c for c in cells if len(c) == 4]
                meshes.append(
                    triangle_mesh(
                        mesh["points"],
                        np.asarray(triangles, dtype=np.int64).reshape(-1, 3),
                        np.asarray(quads, dtype=np.int64).reshape(-1, 4),
                    )
                )
                if isinstance(key, tuple) and key[0] == "side":
                    members.append({"connector": int(key[1])})
                elif isinstance(key, tuple):
                    members.append({"volume": entity.id, "face": "top"})
                elif extruded:
                    members.append({"volume": entity.id, "face": "bottom"})
                else:
                    members.append({"domain": int(key)})
        for index, part in enumerate(meshes):
            part.cell_data["face"] = np.full(part.n_cells, index, dtype=np.int64)
        surface = meshes[0]
        if len(meshes) > 1:
            surface = meshes[0].merge(meshes[1:], merge_points=False)
        return surface, members

    def face_surface(self, member: dict) -> pv.PolyData | None:
        """The cells of one face of a block or volume (as drawn), or None (not drawn)."""
        for owner, members in self._face_members.items():
            actors = self._entity_actors.get(owner, [])
            if member not in members or not actors:
                continue
            surface = actors[0].mapper.dataset
            index = members.index(member)
            cells = np.flatnonzero(np.asarray(surface.cell_data["face"]) == index)
            return surface.extract_cells(cells).extract_surface(algorithm="dataset_surface")
        return None

    def show_patch_faces(self, coloured) -> None:
        """Fill the faces of blocks and volumes with their patch colours ([(member, colour)],
        or None: none)."""
        parts = []
        for member, color in coloured or []:
            part = self.face_surface(member)
            if part is None or not part.n_cells:
                continue
            part.cell_data["colors"] = np.tile(
                np.asarray(pv.Color(color).int_rgb, dtype=np.uint8), (part.n_cells, 1)
            )
            parts.append(part)
        surface = None
        if parts:
            surface = parts[0]
            if len(parts) > 1:
                surface = parts[0].merge(parts[1:])
        self._set_overlay("patch-faces", surface, scalars="colors", rgb=True, opacity=0.6)
        self.patch_faces = coloured

    def set_face_selection(self, members) -> None:
        """Highlight selected faces of blocks and volumes: filled and outlined in bold."""
        self.selected_faces = [dict(m) for m in members]
        parts = [self.face_surface(m) for m in self.selected_faces]
        parts = [p for p in parts if p is not None and p.n_cells]
        if not parts:
            self._set_overlay("selected-faces", None)
            self._set_overlay("selected-faces-outline", None)
            self.plotter.render()
            return
        surface = parts[0]
        if len(parts) > 1:
            surface = parts[0].merge(parts[1:])
        self._set_overlay("selected-faces", surface, color=FACE_COLOR, opacity=0.35)
        outline = surface.extract_feature_edges(
            boundary_edges=True, feature_edges=False, manifold_edges=False, non_manifold_edges=False
        )
        self._set_overlay(
            "selected-faces-outline", outline, color=FACE_COLOR, line_width=SELECTED_CURVE_WIDTH
        )
        self.plotter.render()

    def show_points_preview(self, arrays) -> None:
        """Where curves without points would get them (arrays of points), or nothing."""
        mesh = None
        if arrays:
            mesh = _points(arrays)
        self._set_overlay(
            "points-preview",
            mesh,
            color=POINTS_PREVIEW_COLOR,
            point_size=8,
            render_points_as_spheres=True,
        )
        self.plotter.render()

    def _add_sketch_actor(self, entity) -> None:
        """The sketch's wireframe as drawn (its edges between the vertices) and its vertices."""
        color, width = STYLES["sketch"][int(entity.id in self.selected_ids)]
        lines = self._add(
            _sketch_lines(entity), f"entity-{entity.id}", color=color, line_width=width
        )
        vertices = self._add(
            pv.PolyData(np.asarray(entity.params["vertices"], dtype=float)),
            f"entity-{entity.id}-vertices",
            color=color,
            point_size=9,
            render_points_as_spheres=True,
        )
        self._entity_actors[entity.id] = [lines, vertices]

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
            self._set_overlay("worst-cell", None)
            self.plotter.render()
            return
        corners = np.vstack([corners, corners[:1]])
        self._set_overlay("worst-cell", _polyline(corners), color=WORST_COLOR, line_width=5)
        center = corners[:-1].mean(axis=0)
        size = float(np.max(np.ptp(corners, axis=0)))
        camera = self.plotter.camera
        camera.focal_point = tuple(center)
        camera.position = (center[0], center[1], center[2] + 10.0 * max(size, 1e-12))
        camera.parallel_scale = max(10.0 * size, 1e-9)
        self.camera_changed()

    # ------------------------------------------------------- screen <-> model
    def pick_display(self, x: float, y: float) -> tuple[int | None, tuple | None]:
        """The entity drawn under display pixel (x, y), whatever the viewing direction (VTK
        cell picking on the actors: the closest line or point within PICK_PIXELS along the
        view ray), and the picked point. Blocks, volumes and surfaces are picked by their
        drawn lines or faces."""
        import vtk

        renderer = self.plotter.renderer
        width, height = renderer.GetSize()
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(PICK_PIXELS / max(math.hypot(width, height), 1.0))
        # Lines and points first (the layers and the sketches): a curve seen through a
        # surface, a wing inside its far-field box, is what a click on it means.
        picker.PickFromListOn()
        for layer in self._layers.values():
            picker.AddPickList(layer)
        store = self.session.store
        for entity_id, actors in self._entity_actors.items():
            if entity_id in store and store.get(entity_id).kind == "sketch":
                for actor in actors:
                    picker.AddPickList(actor)
        hit = picker.Pick(x, y, 0.0, renderer)
        if not hit:
            # Then the model's other entities (domains, surfaces, the faces of blocks and
            # volumes) - not the reference grid or the overlays, which would catch the ray.
            picker.InitializePickList()
            for actors in self._entity_actors.values():
                for actor in actors:
                    picker.AddPickList(actor)
            hit = picker.Pick(x, y, 0.0, renderer)
        self.last_face = None
        if not hit:
            return None, None
        actor = picker.GetProp3D()
        cell = picker.GetCellId()
        point = tuple(float(v) for v in picker.GetPickPosition())
        for entity_id, actors in self._entity_actors.items():
            if any(a is actor for a in actors):
                members = self._face_members.get(entity_id)
                faces = actor.mapper.dataset.cell_data.get("face")
                if members is not None and faces is not None and 0 <= cell < len(faces):
                    self.last_face = dict(members[int(faces[cell])])
                return entity_id, point
        for layer in self._layers.values():
            if layer is actor:
                ids = layer.mapper.dataset.cell_data.get("entity")
                if ids is None or cell < 0 or cell >= len(ids):
                    return None, point
                return int(ids[cell]), point
        return None, point

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
    curves and sketch edges.
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
        elif entity.kind == "sketch":
            from opengrid.model.sketch import Lattice

            lattice = Lattice.build(entity.params["vertices"], entity.params["blocks"], "")
            vertices = np.asarray(entity.params["vertices"], dtype=float)
            point = np.asarray(xyz, dtype=float)
            distance = min(
                _distance_to_polyline(point, vertices[[a, b]]) for a, b in sorted(lattice.edges)
            )
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


def _sketch_lines(entity) -> pv.PolyData:
    from opengrid.model.sketch import Lattice

    params = entity.params
    lattice = Lattice.build(params["vertices"], params["blocks"], entity.name)
    edges = np.array(sorted(lattice.edges), dtype=np.int64)
    lines = np.concatenate([np.full((len(edges), 1), 2), edges], axis=1).reshape(-1)
    return pv.PolyData(np.asarray(params["vertices"], dtype=float), lines=lines)


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


def _points(arrays, ids=None) -> pv.PolyData | None:
    """One point cloud from arrays of points (or single points), with the entity id of each
    point as cell data when ``ids`` (one per array) is given; None when empty."""
    rows = [np.asarray(a, dtype=float).reshape(-1, 3) for a in arrays]
    keep = [k for k, r in enumerate(rows) if len(r)]
    if not keep:
        return None
    mesh = pv.PolyData(np.concatenate([rows[k] for k in keep]))
    if ids is not None:
        mesh.cell_data["entity"] = np.repeat(
            np.asarray([ids[k] for k in keep], dtype=np.int64), [len(rows[k]) for k in keep]
        )
    return mesh


def _polylines(arrays, ids=None) -> pv.PolyData | None:
    """One PolyData holding every polyline (one line cell each, tagged with its entity id
    when ``ids`` is given); None when empty."""
    rows = [np.asarray(a, dtype=float).reshape(-1, 3) for a in arrays]
    keep = [k for k, r in enumerate(rows) if len(r) >= 2]
    if not keep:
        return None
    rows = [rows[k] for k in keep]
    counts = np.array([len(a) for a in rows], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    cells = np.empty(int(counts.sum()) + len(rows), dtype=np.int64)
    at = 0
    for start, count in zip(starts, counts, strict=True):
        cells[at] = count
        cells[at + 1 : at + 1 + count] = np.arange(start, start + count)
        at += 1 + count
    mesh = pv.PolyData(np.concatenate(rows), lines=cells)
    if ids is not None:
        mesh.cell_data["entity"] = np.asarray([ids[k] for k in keep], dtype=np.int64)
    return mesh


def _segments(arrays, ids=None) -> pv.PolyData | None:
    """One PolyData with a two-point line cell per interval of every polyline, in order (so a
    per-interval scalar can be attached; the entity id per interval when ``ids`` is given);
    None when empty."""
    rows = [np.asarray(a, dtype=float).reshape(-1, 3) for a in arrays]
    keep = [k for k, r in enumerate(rows) if len(r) >= 2]
    if not keep:
        return None
    rows = [rows[k] for k in keep]
    counts = np.array([len(a) for a in rows], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    pieces = []
    for start, count in zip(starts, counts, strict=True):
        cells = np.empty((count - 1, 3), dtype=np.int64)
        cells[:, 0] = 2
        cells[:, 1] = np.arange(start, start + count - 1)
        cells[:, 2] = cells[:, 1] + 1
        pieces.append(cells.ravel())
    mesh = pv.PolyData(np.concatenate(rows), lines=np.concatenate(pieces))
    if ids is not None:
        mesh.cell_data["entity"] = np.repeat(
            np.asarray([ids[k] for k in keep], dtype=np.int64), counts - 1
        )
    return mesh


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
