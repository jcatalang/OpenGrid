"""Public scripting API (the command layer the GUI uses).

A :class:`Session` owns the model (entity store), the command stack and the journal. Every
mutating method builds a command and runs it through the stack, so everything done here, or
in the GUI (which calls the same methods), is undoable and recorded in the journal as a
replayable Python script::

    import opengrid
    og = opengrid.Session()
    og.create_line((0, 0), (1, 0))
    og.save("geometry.ogp")
"""

import logging
import math
from pathlib import Path

import numpy as np

from opengrid import _core, project
from opengrid.commands import (
    CommandStack,
    CreateEntities,
    CreateEntity,
    DeleteEntities,
    ModifyEntity,
    RenameEntity,
)
from opengrid.commands.base import JOURNAL_HEADER
from opengrid.io import read_airfoil_dat, read_points_csv, write_plot3d
from opengrid.io.cad import CadImport, cad_entities, check_options, read_cad_model
from opengrid.model.blocks import BlockEvaluator
from opengrid.model.connectors import ConnectorEvaluator, ConnectorNodes, order_chain
from opengrid.model.domains import DomainEvaluator, SmoothingResult, assemble_edges
from opengrid.model.store import Entity, Store
from opengrid.model.unstructured import UnstructuredEvaluator, assemble_loops
from opengrid.model.volumes import VolumeEvaluator
from opengrid.validation import Issue, ValidationError

IMPORT_MODES = ("spline", "split", "polyline")
log = logging.getLogger("opengrid")


class Session:
    def __init__(self):
        self.store = Store()
        self.commands = CommandStack(self.store)
        self.connectors = ConnectorEvaluator(self.store)
        self.domains = DomainEvaluator(self.store, self.connectors)
        self.unstructured = UnstructuredEvaluator(self.store, self.connectors)
        self.blocks = BlockEvaluator(self.store, self.domains)
        self.volumes = VolumeEvaluator(self)
        self.commands.add_listener(self._invalidate_derived)
        self.project_path: Path | None = None
        # GUI/report threshold for the stretching-ratio warning (not stored in projects).
        self.stretching_threshold = 1.3
        self.last_warnings: list[str] = []  # warnings of the last topology edit

    # ------------------------------------------------------------- creation
    def create_point(self, xyz, name: str | None = None) -> str:
        """Database point."""
        return self._create("point", {"xyz": xyz}, name)

    def create_line(self, p0, p1, name: str | None = None) -> str:
        """Straight line from p0 to p1 (2 or 3 coordinates each)."""
        return self._create("line", {"p0": p0, "p1": p1}, name)

    def create_arc(
        self,
        center,
        radius: float,
        start_angle: float,
        end_angle: float,
        normal=(0.0, 0.0, 1.0),
        x_axis=(1.0, 0.0, 0.0),
        name: str | None = None,
    ) -> str:
        """Arc of a circle; angles in radians measured from ``x_axis`` about ``normal``."""
        params = {
            "center": center,
            "radius": radius,
            "normal": normal,
            "x_axis": x_axis,
            "start_angle": start_angle,
            "end_angle": end_angle,
        }
        return self._create("arc", params, name)

    def create_arc_3pt(self, start, through, end, name: str | None = None) -> str:
        """Arc from ``start`` through ``through`` to ``end``."""
        params = arc_through_points(start, through, end)
        return self._create("arc", params, name)

    def create_circle(
        self, center, radius: float, normal=(0.0, 0.0, 1.0), name: str | None = None
    ) -> str:
        """Full circle starting on the x axis (or y axis if the normal is along x)."""
        x_axis = (1.0, 0.0, 0.0)
        n = _as3(normal)
        if abs(n[0]) >= 0.9 * math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2):
            x_axis = (0.0, 1.0, 0.0)
        return self.create_arc(center, radius, 0.0, math.tau, normal, x_axis, name)

    def create_polyline(self, points, name: str | None = None) -> str:
        """Polyline through the points (straight segments)."""
        return self._create("polyline", {"points": points}, name)

    def create_spline(
        self,
        points,
        method: str = "natural",
        start_tangent=None,
        end_tangent=None,
        name: str | None = None,
    ) -> str:
        """Interpolating cubic spline; method is 'natural', 'clamped' or 'catmull_rom'."""
        params = {
            "points": points,
            "method": method,
            "start_tangent": start_tangent,
            "end_tangent": end_tangent,
        }
        return self._create("spline", params, name)

    def create_nurbs(
        self, degree: int, knots, control_points, weights=None, name: str | None = None
    ) -> str:
        """NURBS curve: full knot vector (len(control_points) + degree + 1), optional weights."""
        params = {
            "degree": degree,
            "knots": knots,
            "control_points": control_points,
            "weights": weights,
        }
        return self._create("nurbs", params, name)

    # ----------------------------------------------------------- connectors
    # -------------------------------------------------------------- surfaces
    def create_plane(self, origin, u_vector, v_vector, name: str | None = None) -> str:
        """Plane patch x = origin + u u_vector + v v_vector, (u, v) in [0, 1]^2."""
        return self._create(
            "plane", {"origin": origin, "u_vector": u_vector, "v_vector": v_vector}, name
        )

    def create_ruled_surface(self, curve0, curve1, name: str | None = None) -> str:
        """x = (1 - v) curve0(u) + v curve1(u)."""
        return self._create("ruled_surface", {"curve0": curve0, "curve1": curve1}, name)

    def create_extruded_surface(self, curve, direction, name: str | None = None) -> str:
        """x = curve(u) + v direction."""
        return self._create("extruded_surface", {"curve": curve, "direction": direction}, name)

    def create_revolved_surface(
        self,
        curve,
        end_angle: float,
        axis_point=(0.0, 0.0, 0.0),
        axis=(0.0, 0.0, 1.0),
        start_angle: float = 0.0,
        name: str | None = None,
    ) -> str:
        """``curve`` rotated about the axis from start_angle to end_angle (radians)."""
        params = {
            "curve": curve,
            "axis_point": axis_point,
            "axis": axis,
            "start_angle": start_angle,
            "end_angle": end_angle,
        }
        return self._create("revolved_surface", params, name)

    def create_nurbs_surface(
        self,
        degree_u: int,
        degree_v: int,
        knots_u,
        knots_v,
        control_points,
        weights=None,
        name: str | None = None,
    ) -> str:
        """Tensor-product NURBS; ``control_points`` is a list of rows (one per u index)."""
        params = {
            "degree_u": degree_u,
            "degree_v": degree_v,
            "knots_u": knots_u,
            "knots_v": knots_v,
            "control_points": control_points,
            "weights": weights,
        }
        return self._create("nurbs_surface", params, name)

    def project_to_surface(self, name: str, point) -> tuple[float, float, tuple, float]:
        """Closest point on a surface: (u, v, point, distance)."""
        u, v, p, d = self.store.geometry(name).project(tuple(float(c) for c in _as3(point)))
        return u, v, tuple(p), d

    def create_connector(
        self,
        segments,
        dimension: int,
        distribution="uniform",
        mode: str = "arc_length",
        breakpoints=(),
        propagate_spacing: bool = True,
        name: str | None = None,
    ) -> str:
        """Connector on a curve or curve chain.

        ``segments``: a curve name, or a list whose items are curve names or
        ``{"curve": name, "t0": .., "t1": ..}`` / ``(name, t0, t1)`` pieces (t1 < t0 runs a curve
        backwards). ``distribution``: a type name ('uniform', 'geometric', 'tanh', 'sine') or a
        dict with ``type`` and optional ``begin_spacing``, ``end_spacing`` (absolute lengths) and
        ``growth_ratio`` (geometric). ``breakpoints``: list of ``{"position": fraction,
        "index": node, "distribution": ...}``; the breakpoint's distribution applies to the
        sub-connector after it.
        """
        params = {
            "segments": segments,
            "dimension": dimension,
            "mode": mode,
            "distribution": distribution,
            "breakpoints": breakpoints,
            "propagate_spacing": propagate_spacing,
        }
        return self._create("connector", params, name)

    def create_connectors_on_curves(
        self, curves, dimension: int, distribution="uniform"
    ) -> list[str]:
        """One connector per curve (a single undo step)."""
        if isinstance(curves, str):
            curves = [curves]
        items = []
        for curve in curves:
            params = {"segments": curve, "dimension": dimension, "distribution": distribution}
            items.append(("connector", params, None))
        return self.commands.run(CreateEntities(items, label=f"Create {len(items)} connectors"))

    def create_chain_connector(
        self, curves, dimension: int, distribution="uniform", name: str | None = None
    ) -> str:
        """One connector along several curves, ordered and oriented automatically."""
        segments = order_chain(self.store, list(curves))
        return self.create_connector(segments, dimension, distribution, name=name)

    def set_dimension(self, name: str, dimension: int) -> str:
        """Set a connector's number of points."""
        return self.modify(name, dimension=dimension)

    def set_distribution(
        self,
        name: str,
        type: str,
        begin_spacing: float | None = None,
        end_spacing: float | None = None,
        growth_ratio: float | None = None,
    ) -> str:
        """Set a connector's distribution ("uniform", "geometric", "tanh", "sine")."""
        distribution = {
            "type": type,
            "begin_spacing": begin_spacing,
            "end_spacing": end_spacing,
            "growth_ratio": growth_ratio,
        }
        return self.modify(name, distribution=distribution)

    def connector_nodes(self, name: str) -> ConnectorNodes:
        """Nodes and diagnostics of a connector (cached until the model changes)."""
        entity = self.store.get(name)
        if entity.kind != "connector":
            raise ValueError(f"{name!r} is a {entity.kind}, not a connector")
        return self.connectors.nodes(entity.id)

    def connector_report(self, name: str) -> str:
        """One-paragraph summary: points, length, achieved spacings, stretching, warnings."""
        nodes = self.connector_nodes(name)
        lines = [
            f"{name}: {len(nodes.points)} points, length {nodes.length:.10g}",
            f"  spacing begin {nodes.begin_spacing:.6g}, end {nodes.end_spacing:.6g}",
            f"  max stretching ratio {nodes.max_stretching_ratio:.4g}",
        ]
        for end, text in nodes.propagated.items():
            lines.append(f"  {end} spacing propagated: {text}")
        if nodes.max_stretching_ratio > self.stretching_threshold:
            lines.append(f"  WARNING: stretching ratio exceeds {self.stretching_threshold:g}")
        for warning in nodes.warnings:
            lines.append(f"  WARNING: {warning}")
        return "\n".join(lines)

    # -------------------------------------------------------------- domains
    def create_structured_domain(
        self, connectors=None, edges=None, generator=None, name: str | None = None, surface=None
    ) -> str:
        """Structured 2D domain.

        Either ``edges``: four lists (counter-clockwise loop: j = 0 with i increasing, i = imax,
        j = jmax with i decreasing, i = 0) of connector names or ``(name, reversed)`` /
        ``{"connector": name, "reversed": bool}`` items; or ``connectors``: a set of connector
        names assembled automatically (one closed loop, or two closed loops joined by a cut).
        ``generator``: ``{"method": "tfi" | "hermite_tfi", "hermite_walls": (0,),
        "hermite_strength": 1.0, "elliptic": None | {...}}`` (see opengrid.model.domains).
        """
        if edges is None:
            if connectors is None:
                raise ValueError("give either edges or connectors")
            if isinstance(connectors, str):
                connectors = [connectors]
            assembled = assemble_edges(self.store, self.connectors, list(connectors))
            edges = []
            for edge in assembled:
                edges.append(
                    [
                        {
                            "connector": self.store.get(i["connector"]).name,
                            "reversed": i["reversed"],
                        }
                        for i in edge
                    ]
                )
        return self.create_domain(edges, generator, name=name, surface=surface)

    def create_domain(self, edges, generator=None, name: str | None = None, surface=None) -> str:
        """Structured domain from explicit edges (the form recorded in journals). With
        ``surface`` (a surface entity), the grid is the parametric TFI on that surface."""
        params = {"edges": edges, "generator": generator, "surface": surface}
        return self._create("domain", params, name)

    def create_unstructured_domain(
        self, loops=None, generator=None, name: str | None = None, connectors=None, surface=None
    ) -> str:
        """Unstructured (triangle) domain.

        Either ``loops``: the outer loop (counter-clockwise) followed by holes (clockwise), each a
        list of connector names or ``(name, reversed)`` items; or ``connectors``: a set of
        connector names assembled into closed loops automatically (largest loop = outer).
        ``generator``: ``{"algorithm": "delaunay", "min_angle": 25, "growth_rate": 1.2,
        "max_size": None, "sources": [], "max_points": 2000000, "refine": True}``; sources are
        ``{"p0", "p1" (optional; line source), "size", "radius"}``. ``algorithm`` is
        ``"delaunay"`` (constrained Delaunay + refinement) or ``"advancing_front"`` (near
        equilateral triangles that follow the size field closely; ``refine`` and ``min_angle``
        apply to the Delaunay mesher only). Further keys: ``element_type`` ("triangles",
        "quad_dominant", "quads"), ``smoothing`` ("none", "laplacian", "length_weighted",
        "angle_based"), ``smoothing_iterations``, ``edge_swap`` and ``inflation``: boundary
        layers from wall connectors, e.g. ``{"connectors": ["c-wall"], "first_height": 1e-5,
        "growth_rate": 1.2, "layers": 30, "total_thickness": None, "transition": True}``.
        """
        if loops is None:
            if connectors is None:
                raise ValueError("give either loops or connectors")
            if isinstance(connectors, str):
                connectors = [connectors]
            assembled = assemble_loops(self.store, self.connectors, list(connectors), surface)
            loops = []
            for loop in assembled:
                loops.append(
                    [
                        {
                            "connector": self.store.get(i["connector"]).name,
                            "reversed": i["reversed"],
                        }
                        for i in loop
                    ]
                )
        params = {"loops": loops, "generator": generator, "surface": surface}
        return self._create("unstructured_domain", params, name)

    def set_unstructured_generator(self, name: str, **changes) -> str:
        """Change some generator settings of an unstructured domain."""
        generator = dict(self.entity(name).params["generator"])
        generator.update(changes)
        return self.modify(name, generator=generator)

    def unstructured_mesh(self, name: str) -> dict:
        """``{"points" (n, 3), "triangles" (m, 3), "rejected", "warnings"}`` of an unstructured
        domain (boundary points first, in loop order)."""
        entity = self.store.get(name)
        if entity.kind != "unstructured_domain":
            raise ValueError(f"{name!r} is a {entity.kind}, not an unstructured domain")
        return self.unstructured.mesh(entity.id)

    def set_domain_generator(self, name: str, **changes) -> str:
        """Update generator settings, e.g. ``set_domain_generator("d", method="hermite_tfi")`` or
        ``elliptic={"walls": [{"edge": 0, "spacing": 1e-5}]}`` (merged into the current elliptic
        settings; ``elliptic=None`` switches smoothing off). Elliptic keys: ``control``
        ("laplace", "thomas_middlecoff", "ttm"), ``walls``, ``attractions``, ``solver``
        ("line_sor", "point_sor"), ``omega``, ``max_iterations``, ``tolerance`` and
        ``acceleration`` (Anderson depth, default 5; 0: plain relaxation)."""
        current = self.entity(name).params["generator"]
        generator = dict(current)
        for key, value in changes.items():
            if key == "elliptic" and isinstance(value, dict) and current["elliptic"] is not None:
                merged = dict(current["elliptic"])
                merged.update(value)
                generator["elliptic"] = merged
            else:
                generator[key] = value
        return self.modify(name, generator=generator)

    def rotate_domain(self, name: str, steps: int = 1) -> str:
        """Cycle the edges (edge ``steps`` becomes edge 0); orientation is preserved."""
        edges = list(self.entity(name).params["edges"])
        steps = steps % 4
        return self.modify(name, edges=tuple(edges[steps:] + edges[:steps]))

    def _domain_id(self, name: str) -> int:
        entity = self.store.get(name)
        if entity.kind != "domain":
            raise ValueError(f"{name!r} is a {entity.kind}, not a domain")
        return entity.id

    def domain_dimensions(self, name: str) -> tuple[int, int]:
        """(ni, nj) of a structured domain."""
        return self.domains.dimensions(self._domain_id(name))

    def domain_grid(self, name: str, compute: bool = True) -> np.ndarray:
        """Grid points (nj, ni, 3); runs the elliptic smoother if configured and not cached
        (``compute=False`` returns the TFI grid instead while smoothing is pending)."""
        return self.domains.grid(self._domain_id(name), compute)

    def smooth_domain(self, name: str, progress=None) -> SmoothingResult:
        """Run the domain's elliptic smoother (cached); ``progress`` is an optional
        ``opengrid._core.SmoothingProgress`` for monitoring and cancelling from another thread."""
        return self.domains.smooth(self._domain_id(name), progress)

    def domain_status(self, name: str) -> str:
        """One-line summary of a structured or unstructured domain."""
        entity = self.store.get(name)
        if entity.kind == "unstructured_domain":
            mesh = self.unstructured_mesh(name)
            parts = []
            if len(mesh["triangles"]) or not len(mesh["quads"]):
                parts.append(f"{len(mesh['triangles'])} triangles")
            if len(mesh["quads"]):
                parts.append(f"{len(mesh['quads'])} quads")
            text = f"{', '.join(parts)}, {len(mesh['points'])} points"
            if mesh["layers"]:
                text += f", {mesh['layers']} inflation layers"
            if mesh["rejected"]:
                text += f" ({mesh['rejected']} refinement points rejected near the boundary)"
            return text
        return self.domains.status(self._domain_id(name))

    def domain_quality(self, name: str) -> dict:
        """Per-cell quality metrics: arrays of shape (nj - 1, ni - 1) for structured domains,
        (1, m) for unstructured ones (triangles first, then quads)."""
        from opengrid.model.face_quality import face_quality

        entity = self.store.get(name)
        if entity.kind == "unstructured_domain":
            mesh = self.unstructured_mesh(name)
            if np.ptp(mesh["points"][:, 2]) > 0.0 or entity.params.get("surface") is not None:
                return face_quality(mesh["points"], mesh["triangles"], mesh["quads"])
            return _core.quality_cells(
                np.ascontiguousarray(mesh["points"][:, :2]), mesh["triangles"], mesh["quads"]
            )
        grid = self.domain_grid(name)
        if np.ptp(grid[..., 2]) > 0.0 or entity.params.get("surface") is not None:
            # Faces of 3D blocks: 3D metrics, arrays (nj - 1, ni - 1) as for planar domains.
            nj, ni, _ = grid.shape
            index = np.arange(ni * nj).reshape(nj, ni)
            quads = np.stack(
                [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
            ).reshape(-1, 4)
            q = face_quality(grid.reshape(-1, 3), np.zeros((0, 3), dtype=np.int64), quads)
            return {key: value.reshape(nj - 1, ni - 1) for key, value in q.items()}
        return _core.quality_2d(grid)

    def _export_names(self, domains) -> list[str]:
        if domains is None:
            return self.names("domain") + self.names("unstructured_domain")
        if isinstance(domains, str):
            return [domains]
        return list(domains)

    # ------------------------------------------------------ patches and zones
    # --------------------------------------------------------------- volumes
    def create_volume(self, faces, generator=None, name: str | None = None) -> str:
        """Unstructured 3D volume (tetrahedra, pyramids on quad faces) enclosed by the given
        domains (one or more closed shells). ``generator``: ``{"backend": "netgen",
        "max_size": None, "grading": 0.3, "pyramids": "auto" | "always" | "never",
        "inflation": None, "optimize": True}``. ``optimize``: improve cells with a scaled
        Jacobian below 0.3 by moving interior nodes and pyramid apexes and swapping tetrahedra
        (boundary and layer nodes never move). ``inflation``: boundary layers (prisms on
        triangles, hexahedra on quads, fans at convex edges turning more than 60 degrees) from
        wall faces, ``{"faces": [domains], "first_height", "growth_rate": 1.2,
        "layers": 30, "total_thickness": None, "transition": True, "slide": True}``; with
        ``slide`` the layers slide along planar unstructured faces that meet the walls and are
        used by this volume only (those faces are re-meshed around the layers)."""
        return self._create("volume", {"faces": faces, "generator": generator}, name)

    def set_volume_generator(self, name: str, **changes) -> str:
        """Change some generator settings of a volume."""
        generator = dict(self.entity(name).params["generator"])
        generator.update(changes)
        return self.modify(name, generator=generator)

    def volume_mesh(self, name: str) -> dict:
        """``{"points", "tets" (m, 4), "pyramids" (k, 5), "prisms" (k, 6), "hexes" (k, 8),
        "faces" {domain id: boundary cells}, "inflation" (None or {"max_layers", "layers",
        "height", "stops", ...}), "optimization" (None or {"moved", "swaps_23", "swaps_32",
        "worst_before", "worst_after", "below_before", "below_after"}), "warnings", ...}``
        (boundary nodes first, identical to the face meshes)."""
        return self.volumes.mesh(self.store.get(name).id)

    def volume_quality(self, name: str) -> dict:
        """Per-cell quality arrays of a volume mesh (scaled Jacobian, angles, aspect ratio, ...)."""
        from opengrid.model.volumes import volume_quality

        return volume_quality(self.volume_mesh(name))

    def volume_status(self, name: str) -> str:
        """One-line summary of a volume mesh: cell counts, points, shells, layers."""
        mesh = self.volume_mesh(name)
        text = f"{len(mesh['tets'])} tetrahedra"
        for key, label in (("pyramids", "pyramids"), ("prisms", "prisms"), ("hexes", "hexahedra")):
            if len(mesh[key]):
                text += f", {len(mesh[key])} {label}"
        text += f", {len(mesh['points'])} points, {mesh['shells']} shell(s)"
        if mesh["inflation"] is not None:
            text += f"; up to {mesh['inflation']['max_layers']} layers"
        if mesh["warnings"]:
            text += "; " + "; ".join(mesh["warnings"])
        return text

    # ---------------------------------------------------------------- blocks
    def create_block(
        self,
        faces=None,
        extrusion=None,
        name: str | None = None,
        smoothing=None,
        quality_floor=None,
    ) -> str:
        """Structured block from six structured domains (``faces``, any order and orientation)
        or from an ``extrusion`` of a domain: ``{"domain", "mode": "translate" | "rotate" |
        "normal", "vector" (translation or rotation axis), "axis_point", "angle" (radians),
        "distance" (normal), "dimension" (points along the extrusion), "distribution":
        {"type", "begin_spacing", "end_spacing", "growth_ratio"}, "edges": (imin, imax, jmin,
        jmax) conditions of hyperbolic extrusion ("free", "constant_x", "constant_y",
        "constant_z", "periodic"), "dissipation"}``. ``smoothing``: 3D elliptic smoothing
        ``{"control": "laplace" | "thomas_middlecoff", "walls": [{"face": "kmin", "spacing",
        "decay", "relaxation", "max_source"}], "omega", "max_iterations", "tolerance",
        "acceleration" (Anderson depth, default 5; 0: plain relaxation)}``.
        ``quality_floor``: ``{"min_jacobian": 0.05, "max_skewness": 0.95}`` (None: the default);
        cells missing it are reported by the boundary report, and exports refuse the model."""
        params = {
            "faces": faces,
            "extrusion": extrusion,
            "smoothing": smoothing,
            "quality_floor": quality_floor,
        }
        return self._create("block", params, name)

    def set_block_smoothing(self, name: str, smoothing) -> str:
        """Set (or remove with None) the 3D elliptic smoothing of a block."""
        return self.modify(name, smoothing=smoothing)

    def create_extruded_block(
        self,
        domain,
        mode: str = "translate",
        name: str | None = None,
        smoothing=None,
        quality_floor=None,
        **settings,
    ) -> str:
        """Shorthand for create_block(extrusion={"domain": domain, "mode": mode, **settings})."""
        return self.create_block(
            extrusion={"domain": domain, "mode": mode, **settings},
            name=name,
            smoothing=smoothing,
            quality_floor=quality_floor,
        )

    def block_floor(self, name: str) -> dict:
        """Quality-floor check of a block: {"floor", "cells", "below_floor", "worst_jacobian",
        "worst_skewness", "failing" (flat cell indices, i fastest)}."""
        from opengrid.model.quality_floor import block_floor

        return block_floor(self.block_grid(name), self.store.get(name).params["quality_floor"])

    def block_grid(self, name: str) -> np.ndarray:
        """Block grid, shape (nk, nj, ni, 3)."""
        return self.blocks.grid(self.store.get(name).id)

    def block_quality(self, name: str) -> dict:
        """Per-cell quality metrics of a block, arrays of shape (nk - 1, nj - 1, ni - 1)."""
        return _core.quality_3d(np.ascontiguousarray(self.block_grid(name)))

    def block_status(self, name: str) -> str:
        """One-line summary of a block: points, cells, source, smoothing and quality floor."""
        nk, nj, ni, _ = self.block_grid(name).shape
        params = self.store.get(name).params
        source = "6 faces, TFI"
        if params["extrusion"] is not None:
            source = f"{params['extrusion']['mode']} extrusion"
        result = self.blocks.evaluate(self.store.get(name).id)
        text = f"{ni} x {nj} x {nk} points ({(ni - 1) * (nj - 1) * (nk - 1)} cells), {source}"
        if result["smoothing"] is not None:
            state = "converged"
            if not result["smoothing"]["converged"]:
                state = "not converged"
            text += (
                f", elliptic smoothing {state} in {result['smoothing']['iterations']} iterations"
            )
        if result["warnings"]:
            text += "; " + "; ".join(result["warnings"])
        floor = self.block_floor(name)
        if floor["below_floor"]:
            text += (
                f"; {floor['below_floor']} cells miss the quality floor (worst scaled Jacobian "
                f"{floor['worst_jacobian']:.3g}, skewness {floor['worst_skewness']:.3g})"
            )
        return text

    def create_patch(
        self, name: str | None = None, patch_type: str = "wall", members=(), periodic=None
    ) -> str:
        """Named boundary-condition group of connectors; ``patch_type`` in PATCH_TYPES
        (opengrid.model.patches). ``periodic`` is normally set with :meth:`pair_periodic`."""
        params = {"patch_type": patch_type, "members": members, "periodic": periodic}
        return self._create("patch", params, name)

    def pair_periodic(
        self,
        source: str,
        target: str,
        translation=None,
        rotation_angle: float = 0.0,
        rotation_center=(0.0, 0.0, 0.0),
    ) -> str:
        """Pair two periodic patches: target = R(angle)(source - center) + center + translation.

        Without a translation and with no rotation, the translation that maps the centroid of
        the source nodes onto that of the target is computed and stored explicitly.
        """
        from opengrid.model.periodic import automatic_translation

        if translation is None:
            if rotation_angle != 0.0:
                raise ValueError("give the translation explicitly together with a rotation")
            translation = automatic_translation(self, self.entity(source), self.entity(target))
        periodic = {
            "partner": target,
            "rotation_center": rotation_center,
            "rotation_angle": rotation_angle,
            "translation": translation,
        }
        return self.modify(source, periodic=periodic)

    def unpair_periodic(self, source: str) -> str:
        """Remove the periodic pairing of a patch."""
        return self.modify(source, periodic=None)

    def set_patch_type(self, name: str, patch_type: str) -> str:
        """Change a patch's boundary-condition type (wall, inlet, outlet, ...)."""
        return self.modify(name, patch_type=patch_type)

    def add_to_patch(self, name: str, connectors) -> str:
        """Add members (connectors in 2D, domains or block faces in 3D) to a patch."""
        return self._change_members(name, "connector", connectors, add=True)

    def remove_from_patch(self, name: str, connectors) -> str:
        """Remove members from a patch."""
        return self._change_members(name, "connector", connectors, add=False)

    def create_zone(self, name: str | None = None, members=()) -> str:
        """Named group of domains (volume tag)."""
        return self._create("zone", {"members": members}, name)

    def add_to_zone(self, name: str, domains) -> str:
        """Add domains to a zone."""
        return self._change_members(name, "domain", domains, add=True)

    def remove_from_zone(self, name: str, domains) -> str:
        """Remove members from a zone."""
        return self._change_members(name, "domain", domains, add=False)

    def _change_members(self, name: str, key: str, refs, add: bool) -> str:
        if isinstance(refs, (str, int)):
            refs = [refs]
        from opengrid.model.patches import member_label, member_script

        entity = self.entity(name)
        current = [member_script(item, self.store) for item in entity.params["members"]]
        labels = [member_label(item, self.store) for item in entity.params["members"]]
        given = []
        for ref in refs:
            if isinstance(ref, dict) or (isinstance(ref, (list, tuple)) and len(ref) == 2):
                if isinstance(ref, dict):
                    block, face = ref["block"], ref["face"]
                else:
                    block, face = ref
                given.append(
                    (
                        {"block": self.store.get(block).name, "face": face},
                        f"{self.store.get(block).name}.{face}",
                    )
                )
            else:
                given.append((self.store.get(ref).name, self.store.get(ref).name))
        if add:
            members = current + [m for m, label in given if label not in labels]
        else:
            missing = [label for _m, label in given if label not in labels]
            if missing:
                raise ValidationError([Issue(name, f"does not contain {', '.join(missing)}")])
            drop = {label for _m, label in given}
            members = [m for m, label in zip(current, labels, strict=True) if label not in drop]
        return self.modify(name, members=members)

    def is_3d(self) -> bool:
        """True when the model has blocks (exports and boundary checks are then 3D)."""
        return bool(self.store.entities("block")) or bool(self.store.entities("volume"))

    def topology(self):
        """2D: connector uses of domains. 3D models: block faces (boundary / interfaces)."""
        if self.is_3d():
            from opengrid.model.blocks import analyze_blocks

            return analyze_blocks(self.store, self.blocks)
        from opengrid.model.topology import analyze

        return analyze(self.store, self)

    def boundary_report(self):
        """Boundary connectors (2D) or block faces (3D), interfaces (detected automatically),
        boundary entities without a patch and all problems
        (opengrid.model.topology.BoundaryReport)."""
        if self.is_3d():
            from opengrid.model.quality_floor import floor_issue
            from opengrid.model.topology import BoundaryReport
            from opengrid.model.volumes import analyze_3d

            result = analyze_3d(self)
            issues = list(result["issues"])
            # The quality floor of every block and volume.
            for entity in self.store.entities("block") + self.store.entities("volume"):
                try:
                    if entity.kind == "block":
                        summary = self.block_floor(entity.name)
                    else:
                        summary = self.volume_mesh(entity.name)["quality"]
                except (ValidationError, ValueError, RuntimeError) as error:
                    issues.append(Issue(entity.name, f"cannot be meshed: {error}"))
                    continue
                issue = floor_issue(entity.name, summary)
                if issue is not None:
                    issues.append(issue)
            boundary = [e["label"] for e in result["entries"]]
            unassigned = [e["label"] for e in result["entries"] if e["patch"] is None]
            return BoundaryReport(boundary, result["interfaces"], [], unassigned, issues, [])
        from opengrid.model.topology import validate_patches

        return validate_patches(self.store, self.topology(), self)

    def unassigned_boundary_connectors(self) -> list[str]:
        """Boundary entities that belong to no patch."""
        return self.boundary_report().unassigned

    def _require_valid_boundary(self):
        report = self.boundary_report()
        if not report.ok:
            raise ValidationError(report.issues)
        return report

    def export_vtkhdf(self, path, domains=None) -> dict:
        """Write domains to VTKHDF.

        Without patches and zones: one UnstructuredGrid (quads and triangles, shared points
        merged) with ``block_id``, ``zone_id``, cell indices and all quality metrics as cell
        data. With patches or zones (the boundary must then validate): a
        PartitionedDataSetCollection with one cell block per zone and one line block per patch,
        grouped in the assembly as Volume/<zone> and Patches/<patch>.
        """
        from opengrid.io.assembly import assemble
        from opengrid.io.vtkhdf import write_vtkhdf_collection, write_vtkhdf_grid

        if self.is_3d():
            from opengrid.io.assembly3d import assemble_blocks
            from opengrid.model.volumes import analyze_3d

            if not self.store.entities("patch") and not self.store.entities("zone"):
                return write_vtkhdf_grid(path, assemble_blocks(self))
            self._require_valid_boundary()
            return write_vtkhdf_collection(path, assemble_blocks(self, analyze_3d(self)))
        if not self.store.entities("patch") and not self.store.entities("zone"):
            names = None
            if domains is not None:
                names = self._export_names(domains)
            return write_vtkhdf_grid(path, assemble(self, None, domains=names))
        if domains is not None:
            raise ValueError("with patches or zones the whole model is exported (domains=None)")
        report = self._require_valid_boundary()
        return write_vtkhdf_collection(
            path, assemble(self, self.topology(), periodic=report.periodic)
        )

    def export_gmsh(self, path) -> dict:
        """Gmsh MSH 4.1 with physical groups (zones: 2D, patches: 1D); the boundary must
        validate."""
        from opengrid.io.assembly import assemble
        from opengrid.io.gmsh import write_gmsh

        report = self._require_valid_boundary()
        if self.is_3d():
            from opengrid.io.assembly3d import assemble_blocks
            from opengrid.model.volumes import analyze_3d

            return write_gmsh(path, assemble_blocks(self, analyze_3d(self), quality=False))
        return write_gmsh(
            path, assemble(self, self.topology(), quality=False, periodic=report.periodic)
        )

    def export_cgns(self, path, base_name: str = "OpenGrid") -> dict:
        """CGNS: structured zones (domains / blocks) and unstructured zones (unstructured domains /
        volumes) in one base, families, BCs per patch and connectivity; the boundary must
        validate."""
        from opengrid.io.cgns import write_cgns, write_cgns_3d

        report = self._require_valid_boundary()
        if self.is_3d():
            return write_cgns_3d(path, self, base_name)
        return write_cgns(path, self, self.topology(), base_name, periodic=report.periodic)

    def export_plot3d(self, path, domains=None, binary: bool = True, three_d: bool = False) -> None:
        """Write domains as a multi-grid Plot3D file (2D, or 3D with nk = 1); in a model with
        blocks, the blocks (3D)."""
        if self.is_3d() and domains is None:
            if self.store.entities("volume"):
                raise ValidationError(
                    [Issue("Plot3D", "holds structured blocks only; the model has volumes")]
                )
            grids = [self.block_grid(b.name) for b in self.store.entities("block")]
            write_plot3d(path, grids, binary=binary, three_d=True)
            return
        names = self._export_names(domains)
        unstructured = [n for n in names if self.store.get(n).kind == "unstructured_domain"]
        if unstructured:
            raise ValidationError(
                [
                    Issue(
                        n,
                        "Plot3D holds structured blocks only; "
                        "unstructured domains cannot be exported",
                    )
                    for n in unstructured
                ]
            )
        grids = [self.domain_grid(name) for name in names]
        write_plot3d(path, grids, binary=binary, three_d=three_d)

    # -------------------------------------------------------------- editing
    def modify(self, name: str, **params) -> str:
        """Change some parameters of an entity, e.g. ``og.modify("line-1", p1=(2, 0))``."""
        return self.commands.run(ModifyEntity(name, params))

    def rename(self, name: str, new_name: str) -> str:
        """Rename an entity (references follow; journals stay valid)."""
        return self.commands.run(RenameEntity(name, new_name))

    def delete(self, names) -> int:
        """Delete entities (a name or a list); refused while other entities use them."""
        if isinstance(names, str):
            names = [names]
        return self.commands.run(DeleteEntities(list(names)))

    def undo(self) -> None:
        """Undo the last command."""
        self.commands.undo()

    def redo(self) -> None:
        """Redo the last undone command."""
        self.commands.redo()

    # ------------------------------------------------------ topology editing
    def _topology(self, method: str, label: str, *args, **kwargs) -> list[str]:
        from opengrid.commands import TopologyEdit
        from opengrid.model import topology_edit

        function = getattr(topology_edit, method)
        command = TopologyEdit(method, function, args, kwargs, label)
        names = self.commands.run(command)
        self.last_warnings = command.warnings
        for warning in command.warnings:
            log.warning("%s: %s", method, warning)
        return names

    def split_connector(self, name: str, index: int, new_name: str | None = None) -> list[str]:
        """Split a connector at node ``index``: it keeps nodes 0..index, a new connector takes
        the rest (explicit distributions keep every node). Domains, loops and patches using it
        get both. Returns [name, new connector]."""
        kwargs = {}
        if new_name is not None:
            kwargs["new_name"] = new_name
        return self._topology("split_connector", f"Split {name}", name, index, **kwargs)

    def join_connectors(self, first: str, second: str) -> str:
        """Join two connectors sharing an end into ``first`` (all nodes kept); users must take
        them consecutively. Returns the joined connector's name."""
        return self._topology("join_connectors", f"Join {first}", first, second)[0]

    def split_domain(self, name: str, i: int | None = None, j: int | None = None) -> list[str]:
        """Split a structured domain along grid line i (or j). The cut becomes a connector on a
        polyline through the grid line; edge connectors crossing it are split. Returns [name
        (part with the lower index), new domain, cut connector]."""
        kwargs = {}
        if i is not None:
            kwargs["i"] = i
        if j is not None:
            kwargs["j"] = j
        return self._topology("split_domain", f"Split {name}", name, **kwargs)

    def join_domains(self, first: str, second: str) -> str:
        """Join two structured domains sharing a whole edge into ``first``; connectors only
        the shared edge used are deleted. Returns the joined domain's name."""
        return self._topology("join_domains", f"Join {first}", first, second)[0]

    def split_block(self, name: str, axis: str, index: int) -> list[str]:
        """Split a six-face block at ``index`` along ``axis`` ("i", "j" or "k"): its side faces
        are split and a new interior face is created. Returns [name, new block, new face]."""
        return self._topology("split_block", f"Split {name}", name, axis, index)

    def join_blocks(self, first: str, second: str) -> str:
        """Join two six-face blocks sharing a face into ``first`` (side faces joined pairwise,
        the shared face deleted). Returns the joined block's name."""
        return self._topology("join_blocks", f"Join {first}", first, second)[0]

    def ogrid_domain(
        self,
        name: str,
        fraction: float = 0.5,
        radial_points: int | None = None,
        radial_spacing: float | None = None,
        core: bool = True,
    ) -> list[str]:
        """O-grid (butterfly) a structured domain: four domains around a core whose corners are
        centre + fraction * (corner - centre); radial connectors with ``radial_points`` (default
        (ni + nj) / 4) and optional tanh ``radial_spacing`` at the outer boundary. ``core=False``
        leaves the core as a hole. The domain is replaced; returns the new domains."""
        return self._topology(
            "ogrid_domain",
            f"O-grid {name}",
            name,
            fraction=fraction,
            radial_points=radial_points,
            radial_spacing=radial_spacing,
            core=core,
        )

    def ogrid_block(
        self,
        name: str,
        fraction: float = 0.5,
        radial_points: int | None = None,
        radial_spacing: float | None = None,
        core: bool = True,
    ) -> list[str]:
        """O-grid (butterfly) a six-face block: six shell blocks around a core block (or only the
        shell with ``core=False``), all hexahedral. Returns the new blocks."""
        return self._topology(
            "ogrid_block",
            f"O-grid {name}",
            name,
            fraction=fraction,
            radial_points=radial_points,
            radial_spacing=radial_spacing,
            core=core,
        )

    # --------------------------------------------------------------- import
    def import_airfoil_dat(
        self, path, name: str | None = None, mode: str = "spline", method: str = "natural"
    ) -> list[str]:
        """Import a Selig or Lednicer airfoil file.

        mode 'spline': one spline TE -> upper -> LE -> lower -> TE (smooth leading edge);
        mode 'split': '<name>-upper' and '<name>-lower' splines, both running LE -> TE, clamped
        at the leading edge to the tangent of the single-loop spline (G1 at the LE);
        mode 'polyline': the raw points joined by straight segments.
        """
        _check_mode(mode)
        data = read_airfoil_dat(path)
        base = name
        if base is None:
            base = _sanitize(Path(path).stem)
        points = [(x, y, 0.0) for x, y in data.points]
        items: list[tuple[str, dict, str | None]] = []
        if mode == "polyline":
            items.append(("polyline", {"points": points}, base))
        elif mode == "spline":
            items.append(("spline", {"points": points, "method": method}, base))
        else:
            from opengrid.model.entities import build_geometry, normalize_params

            loop = build_geometry(
                "spline", normalize_params("spline", {"points": points, "method": method}, base)
            )
            le = data.leading_edge_index()
            if le == 0 or le == len(points) - 1:
                raise ValidationError(
                    [Issue(str(path), "leading edge is at an end of the point list")]
                )
            tangent = loop.d1(loop.knots[le])
            upper = list(reversed(points[: le + 1]))
            lower = points[le:]
            reverse_tangent = (-tangent[0], -tangent[1], -tangent[2])
            items.append(
                (
                    "spline",
                    {"points": upper, "method": "clamped", "start_tangent": reverse_tangent},
                    f"{base}-upper",
                )
            )
            items.append(
                (
                    "spline",
                    {"points": lower, "method": "clamped", "start_tangent": tangent},
                    f"{base}-lower",
                )
            )
        command = CreateEntities(
            items, label=f"Import {Path(path).name}", comment=f"imported from {Path(path).name}"
        )
        return self.commands.run(command)

    def import_csv(
        self, path, name: str | None = None, mode: str = "spline", method: str = "natural"
    ) -> list[str]:
        """Import a point list as one spline (mode 'spline') or polyline (mode 'polyline')."""
        if mode not in ("spline", "polyline"):
            raise ValueError(f"mode must be 'spline' or 'polyline', got {mode!r}")
        points = read_points_csv(path)
        base = name
        if base is None:
            base = _sanitize(Path(path).stem)
        if mode == "spline":
            item = ("spline", {"points": points, "method": method}, base)
        else:
            item = ("polyline", {"points": points}, base)
        command = CreateEntities(
            [item], label=f"Import {Path(path).name}", comment=f"imported from {Path(path).name}"
        )
        return self.commands.run(command)

    def import_cad(
        self,
        path,
        name: str | None = None,
        spacing: float | None = None,
        max_turn_angle: float = 15.0,
        domains: str = "unstructured",
        domain_generator: dict | None = None,
        tolerance: float = 1e-6,
        scale: float = 1.0,
        heal: bool = True,
        sew: bool = True,
        split_closed: bool = True,
    ) -> CadImport:
        """Import a STEP (.step / .stp) or IGES (.iges / .igs) file (needs the OpenCASCADE
        module, see opengrid.io.cad).

        Creates one ``nurbs`` curve and one connector per CAD edge, and one domain per face
        (``domains``: "unstructured", "auto" = structured for four-sided faces, or "none"),
        on a ``nurbs_surface`` for curved faces. Connector dimensions follow ``spacing``
        (default: 1/40 of the bounding-box diagonal) and ``max_turn_angle`` (degrees of
        tangent turning per segment). ``tolerance`` is the healing / sewing precision (the
        file's own tolerances are used when larger) and ``scale`` multiplies the coordinates
        (OCCT reads lengths in millimetres). ``heal``, ``sew`` and ``split_closed`` (split
        periodic faces and closed edges) switch the healing steps. ``domain_generator`` is
        passed to every unstructured domain. Faces that cannot be meshed are skipped and
        reported in the result (``CadImport.report()``); one undo step undoes the import.
        """
        check_options(spacing, max_turn_angle, domains, tolerance, scale)
        model = read_cad_model(path, tolerance, scale, heal, sew, split_closed)
        base = name
        if base is None:
            base = _sanitize(Path(path).stem)
        items, result = cad_entities(
            model, base, spacing, max_turn_angle, domains, domain_generator, Path(path).name
        )
        comment = (
            f"imported from {Path(path).name} (tolerance {tolerance!r}, scale {scale!r}, "
            f"heal {heal}, sew {sew}, split_closed {split_closed})"
        )
        command = CreateEntities(
            items,
            label=f"Import {Path(path).name}",
            comment=comment,
            optional_kinds=("domain", "unstructured_domain"),
        )
        self.commands.run(command)
        failed = {where: reason for where, reason in command.failed}
        result.domains = [d for d in result.domains if d not in failed]
        result.structured = [d for d in result.structured if d not in failed]
        for where, reason in command.failed:
            number = int(where.rsplit("-", 1)[1])
            result.skipped_faces.append((number, f"invalid domain: {reason}"))
        result.skipped_faces.sort()
        return result

    # -------------------------------------------------------------- queries
    def entity(self, name: str) -> Entity:
        """The stored entity (id, name, kind, canonical params)."""
        return self.store.get(name)

    def names(self, kind: str | None = None) -> list[str]:
        """Entity names in creation order, optionally of one kind ("line", "connector", ...)."""
        return [entity.name for entity in self.store.entities(kind)]

    def geometry(self, name: str):
        """The C++ curve of an entity (None for points)."""
        return self.store.geometry(name)

    def length(self, name: str) -> float:
        """Arc length of a curve (0 for points)."""
        curve = self.geometry(name)
        if curve is None:
            return 0.0
        return curve.length()

    def project_point(self, name: str, point) -> tuple[float, tuple, float]:
        """Closest point on a curve: (t, point, distance)."""
        curve = self.geometry(name)
        if curve is None:
            raise ValueError(f"{name!r} is a point, not a curve")
        return curve.project(_as3(point))

    # --------------------------------------------------------- project files
    def save(self, path) -> None:
        """Write the project file (.ogp, canonical JSON)."""
        project.save(self.store, path)
        self.project_path = Path(path)

    def load(self, path) -> None:
        """Replace the model by a project file. Clears undo history; journal restarts from it."""
        store = project.load(path)
        self._replace_store(store)
        self.project_path = Path(path)
        self.commands.journal.append(f"og.load({str(path)!r})")

    def new(self) -> None:
        """Start an empty model (clears the undo history)."""
        self._replace_store(Store())
        self.project_path = None

    def journal_text(self) -> str:
        """The session journal: a Python script that rebuilds the model (opengrid.replay)."""
        return self.commands.journal.text()

    def save_journal(self, path) -> None:
        """Write the journal script to a file."""
        Path(path).write_text(self.journal_text(), encoding="utf-8")

    # ------------------------------------------------------------- internals
    def _create(self, kind: str, params: dict, name: str | None) -> str:
        return self.commands.run(CreateEntity(kind, params, name))

    def _invalidate_derived(self, _change) -> None:
        self.connectors.invalidate()
        self.domains.invalidate()

    def _replace_store(self, store: Store) -> None:
        listeners = list(self.commands._listeners)
        self.store = store
        self.connectors = ConnectorEvaluator(store)
        self.domains = DomainEvaluator(store, self.connectors)
        self.unstructured = UnstructuredEvaluator(store, self.connectors)
        self.blocks = BlockEvaluator(store, self.domains)
        self.volumes = VolumeEvaluator(self)
        self.commands = CommandStack(store)
        for listener in listeners:
            self.commands.add_listener(listener)
            listener(None)


def replay(source: str) -> Session:
    """Execute a journal script and return the resulting session (bound to ``og``)."""
    import opengrid

    namespace = {"opengrid": opengrid}
    exec(compile(source, "<journal>", "exec"), namespace)
    session = namespace.get("og")
    if not isinstance(session, Session):
        raise ValueError("journal did not create a Session named 'og'")
    return session


def arc_through_points(start, through, end) -> dict:
    """Canonical-ready arc parameters for the circle through three points (start -> end)."""
    p0 = _as3(start)
    p1 = _as3(through)
    p2 = _as3(end)
    a = _sub(p1, p0)
    b = _sub(p2, p0)
    n = _cross(a, b)
    nn = _dot(n, n)
    scale = max(_dot(a, a), _dot(b, b))
    if nn <= 1e-24 * scale * scale or scale == 0.0:
        raise ValidationError([Issue("arc", "the three points are collinear or coincident")])
    # Circumcentre: p0 + (|a|^2 (b x n) + |b|^2 (n x a)) / (2 |n|^2).
    term = _add(_mul(_cross(b, n), _dot(a, a)), _mul(_cross(n, a), _dot(b, b)))
    center = _add(p0, _mul(term, 1.0 / (2.0 * nn)))
    u = _sub(p0, center)
    radius = math.sqrt(_dot(u, u))
    normal = _mul(n, 1.0 / math.sqrt(nn))
    v = _cross(normal, u)
    e = _sub(p2, center)
    sweep = math.atan2(_dot(e, v), _dot(e, u))
    if sweep <= 0.0:
        sweep += math.tau
    return {
        "center": center,
        "radius": radius,
        "normal": normal,
        "x_axis": u,
        "start_angle": 0.0,
        "end_angle": sweep,
    }


def _check_mode(mode: str) -> None:
    if mode not in IMPORT_MODES:
        raise ValueError(f"mode must be one of {', '.join(IMPORT_MODES)}, got {mode!r}")


def _sanitize(stem: str) -> str:
    cleaned = []
    for char in stem:
        if char.isalnum() or char in "_.-":
            cleaned.append(char)
        else:
            cleaned.append("_")
    name = "".join(cleaned).lstrip(".-")
    if not name:
        name = "imported"
    return name


def _as3(v) -> tuple[float, float, float]:
    values = [float(c) for c in v]
    if len(values) == 2:
        values.append(0.0)
    if len(values) != 3:
        raise ValueError(f"expected 2 or 3 coordinates, got {len(values)}")
    return (values[0], values[1], values[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _mul(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


__all__ = ["JOURNAL_HEADER", "Session", "arc_through_points", "replay"]
