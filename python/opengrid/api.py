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
from opengrid.model.connectors import CURVE_KINDS, ConnectorEvaluator, ConnectorNodes, order_chain
from opengrid.model.domains import DomainEvaluator, SmoothingResult, assemble_edges
from opengrid.model.groups import GroupEvaluator, GroupResult
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
        self.groups = GroupEvaluator(self.store, self.connectors, self.domains)
        self.domains.groups = self.groups
        self.unstructured = UnstructuredEvaluator(self.store, self.connectors)
        self.blocks = BlockEvaluator(self.store, self.domains)
        self.blocks.groups = self.groups
        self.groups.blocks = self.blocks
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
    def create_freehand(
        self,
        points,
        tolerance: float = 0.005,
        corner_angle: float | None = 60.0,
        name: str | None = None,
        close_tolerance: float = 0.0,
    ) -> list[str]:
        """NURBS curves fitted to a freehand stroke (pointer samples): split at corners (turns
        of ``corner_angle`` degrees or more; None: never) and fitted by cubic B-splines within
        ``tolerance`` times the stroke length (opengrid.model.fitting). A stroke that ends
        within ``close_tolerance`` (absolute; 0: never) of its start is closed into a loop
        whose curves meet exactly (the GUI passes its pick tolerance), ready for
        ``unstructured_domain_from_curves``. One undo step; the journal records the fitted
        curves. ``name`` is used for a single piece, else ``name-1``, ``name-2``..."""
        from opengrid.model.fitting import fit_stroke

        pieces = fit_stroke(
            [_as3(p) for p in points], tolerance, corner_angle, close_tolerance=close_tolerance
        )
        items = []
        for k, piece in enumerate(pieces):
            params = {key: piece[key] for key in ("degree", "knots", "control_points", "weights")}
            piece_name = name
            if name is not None and len(pieces) > 1:
                piece_name = f"{name}-{k + 1}"
            items.append(("nurbs", params, piece_name))
        return self.commands.run(CreateEntities(items, label="Freehand curve"))

    def domains_from_curves(
        self,
        curves,
        spacing: float | None = None,
        corner_angle: float = 30.0,
        unstructured: bool = False,
    ) -> list[str]:
        """Structured domains from planar curves drawn across each other (a closed loop of
        curves, or a network): the curves are cut where they cross or turn sharply, every
        region with four corners becomes a structured domain with connectors on the curve
        pieces (point counts from ``spacing``, default 1/20 of the median piece length;
        opposite sides made equal so neighbouring domains conform). Other regions are
        reported in ``last_warnings`` (or become unstructured domains with
        ``unstructured``). One undo step; returns the domains created."""
        from opengrid.model.network import build_network, plan_domains

        ids = []
        for curve in curves:
            entity = self.store.get(curve)
            if entity.kind not in CURVE_KINDS:
                raise ValidationError([Issue(entity.name, f"is a {entity.kind}, not a curve")])
            ids.append(entity.id)
        if spacing is None:
            net = build_network(self.store, ids, corner_angle)
            lengths = sorted(edge.length for edge in net.edges)
            if not lengths:
                raise ValidationError([Issue("curves", "no curve pieces")])
            spacing = lengths[len(lengths) // 2] / 20.0
        items, issues, _net = plan_domains(self.store, ids, spacing, corner_angle, unstructured)
        self.last_warnings = [f"{issue.where}: {issue.message}" for issue in issues]
        for warning in self.last_warnings:
            log.warning(warning)
        domains = [name for kind, _params, name in items if kind != "connector"]
        if not domains:
            raise ValidationError(issues or [Issue("curves", "no domain could be formed")])
        self.commands.run(CreateEntities(items, label=f"Domains from {len(ids)} curves"))
        return domains

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

    def unstructured_domain_from_curves(
        self,
        curves,
        dimension: int,
        connectors=(),
        generator=None,
        distribution="uniform",
        name: str | None = None,
    ) -> str:
        """Unstructured (triangle) domain bounded by curves: one connector of ``dimension``
        points on each curve, assembled together with the existing ``connectors`` into closed
        loops (the largest = outer), then the domain, as one undo step (the GUI's
        "Unstructured domain" on selected curves, e.g. a freehand loop). The curve ends must
        meet exactly (``create_freehand`` with ``close_tolerance`` closes a stroke). Returns
        the domain's name; ``generator`` as in ``create_unstructured_domain``."""
        from opengrid.model.unstructured import assemble_polylines

        if isinstance(curves, str):
            curves = [curves]
        if isinstance(connectors, str):
            connectors = [connectors]
        polylines = {}
        labels = {}
        items = []
        counter = self.store.name_counters.get("connector", 0)
        for curve in curves:
            entity = self.store.get(curve)
            if entity.kind not in CURVE_KINDS:
                raise ValidationError([Issue(entity.name, f"is a {entity.kind}, not a curve")])
            while True:
                counter += 1
                connector_name = f"connector-{counter}"
                if connector_name not in self.store:
                    break
            geometry = self.store.geometry(entity.id)
            polylines[connector_name] = geometry.eval_many(np.linspace(0.0, 1.0, 65))
            labels[connector_name] = f"{connector_name} (on {entity.name})"
            params = {"segments": entity.name, "dimension": dimension, "distribution": distribution}
            items.append(("connector", params, connector_name))
        for connector in connectors:
            entity = self.store.get(connector)
            if entity.kind != "connector":
                raise ValidationError([Issue(entity.name, f"is a {entity.kind}, not a connector")])
            polylines[entity.name] = self.connectors.nodes(entity.id).points
            labels[entity.name] = entity.name
        if not polylines:
            raise ValidationError([Issue("unstructured domain", "no curves or connectors given")])
        loops = assemble_polylines(polylines, labels)
        params = {"loops": loops, "generator": generator, "surface": None}
        items.append(("unstructured_domain", params, name))
        label = f"Unstructured domain from {len(curves)} curves"
        return self.commands.run(CreateEntities(items, label=label))[-1]

    def set_unstructured_generator(self, name: str, generator=None, **changes) -> str:
        """Change generator settings of an unstructured domain (``generator``, a dict, and
        ``changes`` merged into the current ones). For a sketch's block the settings go to
        the sketch: the triangulation settings to its ``generator`` (shared by its
        unstructured blocks) and ``inflation`` to walls on the block's edges (first height,
        layers, growth), which is what the property panel does there. One undo step."""
        from opengrid.model.unstructured import normalize_generator

        merged = dict(self.entity(name).params["generator"])
        if generator is not None:
            merged.update(generator)
        merged.update(changes)
        merged = normalize_generator(merged, name, self.store)
        if merged["inflation"] is not None:  # names in the journal, not ids
            inflation = dict(merged["inflation"])
            inflation["connectors"] = tuple(self.store.get(c).name for c in inflation["connectors"])
            merged["inflation"] = inflation
        return self._topology("set_unstructured_generator", f"Generator of {name}", name, merged)[0]

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
        (``compute=False`` returns the TFI grid instead while smoothing is pending). A face of
        a block in a smoothed 3D group is returned where the group moved it."""
        entity_id = self._domain_id(name)
        moved = self.groups.face_grid(entity_id)
        if moved is not None:
            return moved
        return self.domains.grid(entity_id, compute)

    def smooth_domain(self, name: str, progress=None) -> SmoothingResult:
        """Run the domain's elliptic smoother (cached); ``progress`` is an optional
        ``opengrid._core.SmoothingProgress`` for monitoring and cancelling from another thread."""
        return self.domains.smooth(self._domain_id(name), progress)

    # ------------------------------------------------------- smoothing groups
    def create_smoothing_group(
        self,
        domains=(),
        walls=(),
        slide="auto",
        name: str | None = None,
        omega: float = 1.0,
        max_iterations: int = 2000,
        tolerance: float = 1e-6,
        acceleration: int = 20,
        blocks=(),
        joints: str = "fixed",
        pinned=(),
        slide_limit: float = 0.5,
        symmetry=(),
    ) -> str:
        """Smooth planar structured ``domains`` (or six-face ``blocks``, 3D) together as one grid
       : interfaces and block corners move, grid lines continue across
        interfaces, and irregular points (3 or 5+ domains around a node) are relaxed too.
        ``walls``: connectors (or dicts ``{"connector", "spacing" (None: the initial
        first-cell height), "decay", "relaxation", "max_source", "layers", "growth"}``) with
        first-cell height and orthogonality control; they stay fixed. ``layers`` > 0: the
        first layers off the wall are marched along the front normals (heights spacing *
        growth^k, default growth 1.2) and kept, so grid lines stay normal to the wall
        through them; the control then acts on the first cell beyond. ``slide``: boundary
        connectors whose nodes slide along their curves ("auto": every boundary connector
        that is not a wall; () for none). Connectors used by domains outside the group stay
        fixed. Members' own elliptic settings are ignored. The control is Laplace plus the
        wall sources; ``acceleration`` is the Anderson depth (default 20; 0: plain
        relaxation). With ``blocks``: ``walls`` and ``slide`` name block faces (domains);
        sliding faces move on their surface (or plane), boundary edges between faces on
        different surfaces along their connectors; faces shared with blocks or volumes
        outside the group stay fixed. Members' own smoothing settings are ignored.
        ``joints="slide"``: block corners between two sliding connectors on one curve slide
        along it (3D: face corners inside one sliding surface slide on it); ``pinned``:
        connectors (3D: faces) kept fixed even between members, e.g. a wake cut (a wall
        between members is an internal wall: fixed, with wall control on both sides, 2D);
        ``slide_limit``: how far a sliding node may move, in initial intervals; ``symmetry``
        (3D): planar faces that are symmetry planes (nodes move as interior nodes mirrored
        across the plane, without the slide limit)."""
        params = {
            "domains": domains,
            "blocks": blocks,
            "walls": walls,
            "slide": slide,
            "joints": joints,
            "pinned": pinned,
            "symmetry": symmetry,
            "slide_limit": slide_limit,
            "omega": omega,
            "max_iterations": max_iterations,
            "tolerance": tolerance,
            "acceleration": acceleration,
        }
        return self._create("smoothing_group", params, name)

    def smooth_group(self, name: str, progress=None) -> GroupResult:
        """Run (or reuse) a smoothing group; ``progress``: optional
        ``opengrid._core.SmoothingProgress``. Returns grids per member domain id, residuals,
        node counts (regular, irregular, sliding, fixed) and the connector roles."""
        return self.groups.smooth(self._group_id(name), progress)

    def display_connector_nodes(self, name: str) -> ConnectorNodes:
        """Connector nodes as drawn: where a smoothed group moved them (interfaces, sliding
        boundaries), else the connector's own nodes."""
        nodes = self.connector_nodes(name)
        moved = self.groups.moved_connector_points(self.store.get(name).id)
        if moved is None:
            return nodes
        arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(moved, axis=0), axis=1))])
        return ConnectorNodes(moved, arc, list(nodes.warnings), nodes.max_stretching_ratio)

    def group_status(self, name: str) -> str:
        """One-line summary of a smoothing group."""
        return self.groups.status(self._group_id(name))

    def _group_id(self, name: str) -> int:
        entity = self.store.get(name)
        if entity.kind != "smoothing_group":
            raise ValueError(f"{name!r} is a {entity.kind}, not a smoothing group")
        return entity.id

    # ------------------------------------------------------------- sketches
    def create_sketch(
        self,
        vertices,
        blocks,
        dimensions=(),
        spacing: float | None = None,
        curves=(),
        surfaces=(),
        walls=(),
        fixed=(),
        pinned=(),
        symmetry=(),
        distributions=(),
        relax: bool = True,
        tolerance: float = 1e-8,
        max_iterations: int = 4000,
        extrude=None,
        interval_multiple: int = 1,
        connectors=(),
        unstructured=(),
        generator=None,
        name: str | None = None,
    ) -> str:
        """Sketch topology (wireframe of blocks, opengrid.model.sketch): ``vertices``
        (points drawn loosely around the geometry) and ``blocks`` (4 vertex indices each in
        2D, counter-clockwise or not; 8 in 3D, VTK hexahedron order). ``dimensions``: points
        per edge, ``{(a, b): n}`` on any edge of each edge family (or ``spacing`` for the
        families without one); ``curves``: ``{(a, b): curve}`` assigns edges to database
        curves, ``connectors``: ``{(a, b): connector}`` makes 2D edges adopt existing
        connectors (the domain shares their nodes with whatever is on their other side, the
        edge family takes their point count), ``surfaces``: ``{(a, b, c, d): surface}`` 3D
        boundary faces to surfaces;
        ``walls``: ``{(a, b): first-cell height or None}`` (edges, 2D; faces, 3D) or dicts
        with ``"edge"``/``"face"``, ``"spacing"``, ``"decay"``, ``"relaxation"``,
        ``"max_source"``, ``"layers"``, ``"growth"`` (wall layers, as in
        ``create_smoothing_group``; ``{(a, b): {"spacing": h, "layers": n}}`` works too);
        ``fixed``: vertices that must not float; ``distributions``:
        ``{(a, b): distribution}`` from a to b. Generates connectors, domains (3D: also
        blocks) and a smoothing group named ``<sketch>-...`` that the sketch owns; with
        ``relax`` (2D) vertices joining two edges on one curve slide along it until the
        grid stops moving them. ``extrude`` (2D): ``{"vector": (0, 0, L), "points": n,
        "distribution": {...}}`` (spacings as lengths) also makes one extruded block per
        domain: the smoothed 2D grid extruded, the same in every section;
        ``interval_multiple`` m: interval counts rounded up to multiples of m (for
        ``high_order`` meshes of order m). ``unstructured``: indices of the blocks meshed
        as unstructured domains on the same connectors (``generator``: their triangulation
        settings; a wall on their edges becomes inflation layers with the wall's first-cell
        height, layers and growth); they stay out of the smoothing group, and the connectors
        they share with structured blocks are pinned. One undo step; warnings (e.g.
        unassigned edges lying on geometry) are in ``last_warnings``."""
        from opengrid.model.sketch import create_sketch

        kwargs = {
            "vertices": vertices,
            "blocks": blocks,
            "dimensions": dimensions,
            "spacing": spacing,
            "curves": curves,
            "surfaces": surfaces,
            "walls": walls,
            "fixed": fixed,
            "pinned": pinned,
            "symmetry": symmetry,
            "distributions": distributions,
            "relax": relax,
            "tolerance": tolerance,
            "max_iterations": max_iterations,
            "extrude": extrude,
            "interval_multiple": interval_multiple,
            "connectors": connectors,
            "unstructured": unstructured,
            "generator": generator,
            "name": name,
        }
        # The journal records what was given (defaults left out).
        defaults = {
            "relax": True,
            "tolerance": 1e-8,
            "max_iterations": 4000,
            "interval_multiple": 1,
        }
        kwargs = {
            k: v
            for k, v in kwargs.items()
            if k in ("vertices", "blocks")
            or (
                not (v is None or (isinstance(v, (tuple, list, dict)) and not v))
                and defaults.get(k, object()) != v
            )
        }
        return self._sketch_command("create_sketch", create_sketch, "Create sketch", (), kwargs)

    def edit_sketch(self, name: str, **changes) -> str:
        """Change sketch parameters (any of those of create_sketch) and regenerate its
        outputs in place (same names: patches and other references to them stay valid).
        One undo step."""
        from opengrid.model.sketch import edit_sketch

        return self._sketch_command("edit_sketch", edit_sketch, f"Edit {name}", (name,), changes)

    def remove_sketch_blocks(self, name: str, blocks) -> str:
        """Take blocks (indices) out of a 2D sketch and regenerate it: the remaining blocks
        are renumbered (domain names follow), edges left in no block lose their walls,
        pinned state, dimensions and assignments; the vertices stay. The last block cannot
        go (delete the sketch). The GUI deletes a sketch's domain this way. One undo step."""
        from opengrid.model.sketch import remove_blocks

        if isinstance(blocks, int):
            blocks = [blocks]
        return self._sketch_command(
            "remove_sketch_blocks",
            remove_blocks,
            f"Remove blocks from {name}",
            (name, sorted(set(blocks))),
            {},
        )

    def regenerate_sketch(self, name: str) -> str:
        """Generate a sketch's outputs again (e.g. after its geometry changed)."""
        from opengrid.model.sketch import regenerate_sketch

        return self._sketch_command(
            "regenerate_sketch", regenerate_sketch, f"Regenerate {name}", (name,), {}
        )

    def sketch_report(self, name: str) -> str:
        """Vertices (where they went and whether they float), edge families with their
        dimensions, errors and suggested assignments of a sketch."""
        from opengrid.model.sketch import report

        return report(self.store, self._sketch(name))

    def sketch_suggestions(self, name: str, tolerance: float | None = None) -> list[dict]:
        """Unassigned boundary edges (3D: faces) of a sketch lying near a curve (surface):
        ``{"what", "edge" | "face", "curve" | "surface", "message"}``. Nothing is assigned;
        see accept_sketch_suggestions."""
        from opengrid.model.sketch import suggestions

        return suggestions(self.store, self._sketch(name), tolerance)

    def accept_sketch_suggestions(self, name: str, tolerance: float | None = None) -> str:
        """Assign every suggested curve and surface (sketch_suggestions) in one edit."""
        from opengrid.model.sketch import suggestions

        entity = self._sketch(name)
        found = suggestions(self.store, entity, tolerance)
        if not found:
            raise ValidationError([Issue(entity.name, "no suggestions to accept")])
        return self._assign_found(entity, found)

    def assign_sketch_geometry(self, name: str, geometry, tolerance: float | None = None) -> str:
        """Assign the given curves (and 3D: surfaces) to the unassigned boundary edges (faces)
        of a sketch that lie near them - within ``tolerance``, by default the ends within a
        fifth of the edge length and the middle within half of it (the nearest of them; the
        GUI's "assign selected geometry"); given connectors are adopted by the 2D edges
        whose vertices lie at their ends. One edit."""
        from opengrid.model.sketch import suggestions

        entity = self._sketch(name)
        if isinstance(geometry, str):
            geometry = [geometry]
        ids = {self.store.get(g).id for g in geometry}
        found = suggestions(
            self.store, entity, tolerance, only=ids, fraction=0.2, middle_fraction=0.5
        )
        if not found:
            raise ValidationError(
                [Issue(entity.name, "no unassigned boundary edge or face lies near that geometry")]
            )
        return self._assign_found(entity, found)

    def _assign_found(self, entity, found) -> str:
        from opengrid.model.sketch import sketch_script_params

        params = sketch_script_params(entity.params, self.store)
        name = entity.name
        changes = {}
        curves = {item["edge"]: item["curve"] for item in params["curves"]}
        connectors = {item["edge"]: item["connector"] for item in params["connectors"]}
        surfaces = {item["face"]: item["surface"] for item in params["surfaces"]}
        for item in found:
            if "connector" in item:
                connectors[item["edge"]] = item["connector"]
                curves.pop(item["edge"], None)  # the connector carries the geometry
            elif "curve" in item:
                if item["edge"] not in connectors:
                    curves[item["edge"]] = item["curve"]
            else:
                surfaces[item["face"]] = item["surface"]
        if curves != {item["edge"]: item["curve"] for item in params["curves"]}:
            changes["curves"] = curves
        if len(connectors) != len(params["connectors"]):
            changes["connectors"] = connectors
        if len(surfaces) != len(params["surfaces"]):
            changes["surfaces"] = surfaces
        return self.edit_sketch(name, **changes)

    def create_template(self, kind: str, name: str | None = None, **options) -> str:
        """A topology template (opengrid.model.templates): the far-field geometry it needs and
        a sketch with no user-placed vertices. 2D: ``"ogrid"`` (body, radius, blocks,
        radial_points, wall), ``"butterfly"`` (boundary, core, radial_points), ``"cgrid"``
        (body with a sharp trailing edge, radius, downstream, ring, radial_points=(inner,
        outer), wall, wake), ``"hgrid"`` (bottom, right, top, left, blocks=(nx, ny),
        normal_points, wall, wall_sides), ``"body_in_channel"`` (body, bottom, right, top,
        left, margin, radial_points, wall). 3D: ``"pipe"``, ``"cylinder_in_channel"``,
        ``"bump_channel"``, ``"wing_in_box"`` (the 2D template's options plus ``length`` and
        ``span_points``: extruded along z; ``wing_in_box`` also ``taper`` and ``sweep``) and
        ``"sphere_in_box"`` (center, radius, box, cube, points, radial_points, outer_points,
        wall). All take ``spacing`` (default: a twentieth of the median boundary edge),
        ``relax``, ``layers`` / ``growth`` (wall layers on every wall; a C-grid's wake cut
        then takes the wall spacing) and ``interval_multiple`` (interval counts rounded up to
        multiples, for ``high_order``). One undo step; returns the sketch."""
        from opengrid.model.templates import create_template

        kwargs = dict(options)
        kwargs["name"] = name
        return self._sketch_command(
            "create_template", create_template, f"Template {kind}", (kind,), kwargs
        )

    def auto_topology(
        self,
        curves,
        size: float | None = None,
        spacing: float | None = None,
        relax: bool = False,
        name: str | None = None,
    ) -> str:
        """Automatic 2D block topology (opengrid.model.crossfield): a cross field of the region
        bounded by the closed loops of ``curves`` (triangulated with mesh ``size``, default
        1/200 of their total length), its separatrices traced from the singularities and
        corners to the boundary (created as polylines ``separatrix-<k>``) and a sketch whose
        blocks are the four-sided faces they cut (point counts from ``spacing``). One undo
        step; returns the sketch (edit it like any other)."""
        from opengrid.model.crossfield import auto_topology

        kwargs = {"curves": curves, "size": size, "spacing": spacing, "relax": relax, "name": name}
        kwargs = {k: v for k, v in kwargs.items() if v is not None and not (k == "relax" and not v)}
        return self._sketch_command(
            "auto_topology", auto_topology, "Automatic topology", (), kwargs
        )

    def cross_field(self, curves, size: float | None = None):
        """The cross field of a planar region and its separatrices, without creating anything:
        a Layout with ``field`` (points, triangles, complex values z = exp(4 i theta),
        corners, singularities, ``index_sum()`` and ``euler()``) and ``separatrices``."""
        from opengrid.model.crossfield import trace_layout

        if isinstance(curves, str):
            curves = [curves]
        return trace_layout(self.store, [self.store.get(c).id for c in curves], size)

    def _sketch(self, name: str):
        entity = self.store.get(name)
        if entity.kind != "sketch":
            raise ValueError(f"{name!r} is a {entity.kind}, not a sketch")
        return entity

    def _sketch_command(self, method, function, label, args, kwargs) -> str:
        from opengrid.commands import TopologyEdit
        from opengrid.model.sketch import _canonical_input

        args = tuple(_canonical_input(a) for a in args)
        kwargs = {k: _canonical_input(v) for k, v in kwargs.items()}
        command = TopologyEdit(method, function, args, kwargs, label)
        names = self.commands.run(command)
        self.last_warnings = command.warnings
        for warning in command.warnings:
            log.warning("%s: %s", method, warning)
        return names[0]

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
    def create_volume(
        self, faces=None, generator=None, name: str | None = None, extrusion=None
    ) -> str:
        """Unstructured 3D volume (tetrahedra, pyramids on quad faces) enclosed by the given
        domains (one or more closed shells), or - with ``extrusion`` instead of ``faces`` -
        prisms and hexahedra swept from an unstructured domain (see
        ``create_extruded_volume``). ``generator``: ``{"backend": "netgen",
        "max_size": None, "grading": 0.3, "pyramids": "auto" | "always" | "never",
        "inflation": None, "optimize": True}``. ``optimize``: improve cells with a scaled
        Jacobian below 0.3 by moving interior nodes and pyramid apexes and swapping tetrahedra
        (boundary and layer nodes never move). ``inflation``: boundary layers (prisms on
        triangles, hexahedra on quads, fans at convex edges turning more than 60 degrees) from
        wall faces, ``{"faces": [domains], "first_height", "growth_rate": 1.2,
        "layers": 30, "total_thickness": None, "transition": True, "slide": True}``; with
        ``slide`` the layers slide along planar unstructured faces that meet the walls and are
        used by this volume only (those faces are re-meshed around the layers)."""
        params = {"faces": faces, "extrusion": extrusion, "generator": generator}
        return self._create("volume", params, name)

    def create_extruded_volume(
        self,
        domain,
        mode: str = "translate",
        name: str | None = None,
        generator=None,
        **settings,
    ) -> str:
        """Volume of prisms (from the triangles) and hexahedra (from the quads) swept from an
        unstructured ``domain``: ``mode`` "translate" (``vector``), "rotate" (axis ``vector``
        through ``axis_point``, ``angle`` in radians) or "normal" (``distance`` along the plane
        normal of a planar domain), with ``dimension`` points across the layers and a
        connector-style ``distribution`` (as ``create_extruded_block``; the base is mirrored
        when the sweep runs against its normal). Its boundary for patches: the domain (the
        bottom), ``{"volume": name, "face": "top"}`` and the domain's connectors (the side
        walls), so a 2D model's patches carry over to the one-cell-thick 3D mesh some solvers
        want. ``generator`` only sets the quality floor. Structured domains extrude to blocks."""
        return self.create_volume(
            extrusion={"domain": domain, "mode": mode, **settings}, generator=generator, name=name
        )

    def set_volume_generator(self, name: str, **changes) -> str:
        """Change some generator settings of a volume."""
        generator = dict(self.entity(name).params["generator"])
        generator.update(changes)
        return self.modify(name, generator=generator)

    def volume_mesh(self, name: str) -> dict:
        """``{"points", "tets" (m, 4), "pyramids" (k, 5), "prisms" (k, 6), "hexes" (k, 8),
        "faces" {domain id, ("top", domain id) or ("side", connector id): outward boundary
        cells}, "extrusion" (None or {"mirrored", "layers"}), "inflation" (None or
        {"max_layers", "layers", "height", "stops", ...}), "optimization" (None or {"moved",
        "swaps_23", "swaps_32",
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
        text += f", {len(mesh['points'])} points"
        if mesh.get("extrusion") is not None:
            params = self.store.get(name).params["extrusion"]
            source = self.store.get(params["domain"]).name
            text += (
                f"; {params['mode']} extrusion of {source}, {mesh['extrusion']['layers']} layers"
            )
            if mesh["extrusion"]["mirrored"]:
                text += " (base mirrored)"
        else:
            text += f", {mesh['shells']} shell(s)"
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
        from opengrid.model.patches import normalize_patch, refuse_shared_connectors

        params = {"patch_type": patch_type, "members": members, "periodic": periodic}
        refuse_shared_connectors(
            self.store, normalize_patch(params, self.store, "patch")["members"]
        )
        return self._create("patch", params, name)

    def pair_periodic(
        self,
        source: str,
        target: str,
        translation=None,
        rotation_angle: float = 0.0,
        rotation_center=(0.0, 0.0, 0.0),
        rotation_axis=(0.0, 0.0, 1.0),
    ) -> str:
        """Pair two periodic patches: target = R(axis, angle)(source - center) + center +
        translation (the rotation by ``rotation_angle`` radians about ``rotation_axis``
        through ``rotation_center``, right-hand rule).

        Without a translation and with no rotation, the translation that maps the centroid of
        the source nodes onto that of the target is computed and stored explicitly. 2D: the
        source connectors must map node by node onto the target's; 3D (block faces, volume
        faces, the ends and side walls of extruded volumes): the faces must map node by node
        and face by face. The boundary report names what does not match.
        """
        from opengrid.model.periodic import automatic_translation

        if translation is None:
            if rotation_angle != 0.0:
                raise ValueError("give the translation explicitly together with a rotation")
            translation = automatic_translation(self, self.entity(source), self.entity(target))
        periodic = {
            "partner": target,
            "rotation_center": rotation_center,
            "rotation_axis": rotation_axis,
            "rotation_angle": rotation_angle,
            "translation": translation,
        }
        return self.modify(source, periodic=periodic)

    def make_periodic_pair(
        self,
        source,
        target,
        translation=None,
        rotation_angle: float = 0.0,
        rotation_center=(0.0, 0.0, 0.0),
        rotation_axis=(0.0, 0.0, 1.0),
        name: str | None = None,
    ) -> list[str]:
        """Two periodic patches from the two sides (connectors in 2D; domains, block faces or
        extruded-volume ends in 3D) and their pairing, in one step (the GUI's "Pair as
        periodic"): patches ``<name>-a`` and ``<name>-b`` (``periodic-1``, ``periodic-2``, ...
        when no name is given). Without a
        translation and no rotation, the translation maps the centroid of the source side
        onto the target's. Returns [source patch, target patch]; one undo step."""
        if isinstance(source, (str, int, dict)):
            source = [source]
        if isinstance(target, (str, int, dict)):
            target = [target]
        if translation is None:
            if rotation_angle != 0.0:
                raise ValueError("give the translation explicitly together with a rotation")
            delta = self._members_centroid(target) - self._members_centroid(source)
            translation = (float(delta[0]), float(delta[1]), float(delta[2]))
        return self._topology(
            "make_periodic_pair",
            "Periodic pair",
            list(source),
            list(target),
            translation=translation,
            rotation_angle=rotation_angle,
            rotation_center=rotation_center,
            rotation_axis=rotation_axis,
            name=name,
        )

    def _members_centroid(self, members) -> np.ndarray:
        """Centroid of the nodes of boundary members (connectors, domains, block faces,
        extruded-volume ends)."""
        from opengrid.model.blocks import face_array

        pieces = []
        for member in members:
            if isinstance(member, dict) and "block" in member:
                grid = self.block_grid(self.store.get(member["block"]).name)
                pieces.append(face_array(grid, member["face"]).reshape(-1, 3))
            elif isinstance(member, dict) and "volume" in member:
                volume = self.store.get(member["volume"])
                mesh = self.volume_mesh(volume.name)
                key = volume.params["extrusion"]["domain"]
                if member["face"] == "top":
                    key = ("top", key)
                ids = np.unique(np.concatenate([np.asarray(f) for f in mesh["faces"][key]]))
                pieces.append(mesh["points"][ids])
            else:
                ref = member
                if isinstance(member, dict):
                    ref = member.get("connector", member.get("domain"))
                entity = self.store.get(ref)
                if entity.kind == "connector":
                    pieces.append(np.asarray(self.connector_nodes(entity.name).points))
                elif entity.kind == "domain":
                    pieces.append(self.domain_grid(entity.name).reshape(-1, 3))
                elif entity.kind == "unstructured_domain":
                    pieces.append(np.asarray(self.unstructured_mesh(entity.name)["points"]))
                else:
                    raise ValueError(f"{entity.name} is a {entity.kind}, not a boundary")
        if not pieces:
            raise ValueError("no boundary members given")
        return np.concatenate(pieces).mean(axis=0)

    def extrude_model(self, thickness: float, layers: int = 1, ends: str = "periodic") -> list[str]:
        """The one-cell-thick (``layers`` cells) 3D model of a 2D one: every domain extruded
        by ``thickness`` along z (blocks from structured domains, volumes of prisms and
        hexahedra from unstructured ones, named ``<domain>-ext``), the 2D patches carried
        over to the side faces with their periodic pairings, zones extended, and the fronts
        and backs in the patches "front" / "back": a periodic pair (``ends="periodic"``) or
        two symmetry planes (``ends="symmetry"``). One undo step; returns the new names."""
        return self._topology(
            "extrude_model", "Extrude the 2D model", thickness, layers=layers, ends=ends
        )

    def convert_domain(self, name: str) -> str:
        """Turn a structured domain into an unstructured one on the same connectors, or an
        unstructured domain (one loop that assembles into four edges) into a structured one;
        the name stays and zones and 3D patches follow. Refused for a sketch's domains (set
        the sketch's unstructured blocks), members of a smoothing group and faces of blocks or
        volumes. One undo step."""
        return self._topology("convert_domain", f"Convert {name}", name)[0]

    def set_boundary_type(self, members, patch_type: str, name: str | None = None) -> str:
        """Give boundary members (connectors in 2D; domains, ``{"block", "face"}``,
        ``{"volume", "face"}`` and side-wall connectors in 3D) a boundary type: they join
        the patch called ``name`` (created when absent), else the first patch of that type,
        else a new patch named after the type - leaving the patches they were in (a member
        belongs to one patch). The GUI's "Boundary type" menu. One undo step; returns the
        patch name."""
        return self._topology(
            "set_boundary_type", f"Boundary type {patch_type}", members, patch_type, name=name
        )[0]

    def set_dimension(self, name: str, dimension: int, **changes) -> list[str]:
        """A connector's number of points, with the connectors opposite to it in structured
        domains (and so across the faces of blocks) following, so the domains stay valid
        (opengrid.model.dimension); of several connectors on an opposite edge the longest
        takes the change (a warning says so). A sketch's connector sets its edge family in
        the sketch, and its ``distribution`` that edge's distribution there (the sketch owns
        the rest of its parameters). ``changes``: other parameters of the connector
        (distribution, breakpoints, ...) set in the same step. The property panel's
        connector edits go through this; ``modify`` changes one connector only and is
        refused when a domain would become invalid. One undo step; returns the connectors
        changed."""
        return self._topology("set_dimension", f"Dimension of {name}", name, dimension, **changes)

    def clear_boundary_type(self, members) -> list[str]:
        """Take boundary members out of every patch: they are unassigned again (the GUI's
        "Clear type"). Returns the patches that changed; one undo step."""
        return self._topology("clear_boundary_type", "Clear boundary type", members)

    def unpair_periodic(self, source: str) -> str:
        """Remove the periodic pairing of a patch."""
        return self.modify(source, periodic=None)

    def set_patch_type(self, name: str, patch_type: str) -> str:
        """Change a patch's boundary-condition type (wall, inlet, outlet, ...)."""
        return self.modify(name, patch_type=patch_type)

    def add_to_patch(self, name: str, connectors) -> str:
        """Add members (connectors in 2D, domains, block faces or extruded-volume ends
        ``{"volume", "face"}`` in 3D; connectors are also the side walls of extruded volumes)
        to a patch."""
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
        from opengrid.model.patches import member_label, member_script, refuse_shared_connectors

        entity = self.entity(name)
        if add and entity.kind == "patch":
            members = []
            for ref in refs:
                if isinstance(ref, (str, int)) and self.store.get(ref).kind == "connector":
                    members.append({"connector": self.store.get(ref).id})
            refuse_shared_connectors(self.store, members)
        current = [member_script(item, self.store) for item in entity.params["members"]]
        labels = [member_label(item, self.store) for item in entity.params["members"]]
        given = []
        for ref in refs:
            if isinstance(ref, dict) and "volume" in ref:
                volume = self.store.get(ref["volume"]).name
                given.append(({"volume": volume, "face": ref["face"]}, f"{volume}.{ref['face']}"))
            elif isinstance(ref, dict) or (isinstance(ref, (list, tuple)) and len(ref) == 2):
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
            assignment = {e["label"]: e["patch"] for e in result["entries"]}
            ids = {}
            for entry in result["entries"]:
                ids[entry["label"]] = sorted(key[1] for key in entry["keys"])[0]
            return BoundaryReport(
                boundary,
                result["interfaces"],
                [],
                unassigned,
                issues,
                result["periodic"],
                assignment,
                ids,
            )
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

    def _linear_assembly(self):
        """The linear mesh of the whole model as exported: with patches the boundary must
        validate (boundary elements per patch); without, there are no boundary elements."""
        from opengrid.io.assembly import assemble

        periodic = ()
        if self.store.entities("patch"):
            periodic = self._require_valid_boundary().periodic
        if self.is_3d():
            from opengrid.io.assembly3d import assemble_blocks
            from opengrid.model.volumes import analyze_3d

            return assemble_blocks(self, analyze_3d(self), quality=False)
        return assemble(self, self.topology(), quality=False, periodic=periodic)

    def high_order(self, order: int = 2, curved: bool = True, repair: bool = True):
        """Order-``order`` elements (1 .. 6) of every mesh cell: triangles and quadrilaterals
        (2D), tetrahedra, pyramids, prisms and hexahedra (3D), the same cells as the linear
        mesh with equispaced Lagrange nodes, shared nodes merged exactly, and boundary
        elements (lines / triangles / quadrilaterals) per patch (with patches the boundary
        must validate). ``curved``: boundary nodes moved onto the geometry (curves,
        surfaces) and the displacement extended into the neighbouring elements, which are
        checked for validity (ValidationError for folded elements; ``repair``: the elements
        that fold or lose most of their quality are first optimised, see ``curving``
        ["repair"]); False: straight-sided elements (straight edges, planar triangles,
        bilinear quadrilateral faces). Returns an
        ElementMesh (``points``, ``cells`` / ``faces`` per kind in Gmsh node order,
        ``connectivity(kind, fmt)``, ``element_nodes(kind, nodes)``, ``curving`` report)."""
        from opengrid.model.curving import curve_boundary
        from opengrid.model.elements import build_elements

        assembly = self._linear_assembly()
        mesh = build_elements(assembly, order)
        if curved:
            curve_boundary(self, assembly, mesh, repair=repair)
        return mesh

    def export_high_order(
        self,
        path,
        order: int = 2,
        nodes: str = "equispaced",
        base_name: str = "OpenGrid",
        curved: bool = True,
        repair: bool = True,
    ) -> dict:
        """Write ``high_order(order, curved, repair)`` by extension: ``.msh`` (Gmsh 4.1
        high-order element types, physical groups, periodic links), ``.cgns`` (SIDS types up
        to order 4, one unstructured zone per zone, BCs per patch), ``.vtkhdf`` (VTK Lagrange
        cells) or ``.h5``
        (element node coordinates in the node set ``nodes``: "equispaced", or "gll" / "gauss"
        for quadrilaterals and hexahedra)."""
        from opengrid.io.cgns_elements import write_cgns_elements
        from opengrid.io.gmsh import write_gmsh_elements
        from opengrid.io.highorder import write_elements_h5, write_vtkhdf_elements

        suffix = Path(path).suffix.lower()
        writers = {".msh": write_gmsh_elements, ".vtkhdf": write_vtkhdf_elements}
        if suffix not in (".msh", ".cgns", ".vtkhdf", ".h5"):
            raise ValueError(f"unknown high-order format {suffix!r}: .msh, .cgns, .vtkhdf, .h5")
        if nodes != "equispaced" and suffix != ".h5":
            raise ValueError(
                f"{suffix} elements carry equispaced nodes: write the {nodes!r} node set to .h5"
            )
        mesh = self.high_order(order, curved, repair)
        if suffix == ".h5":
            return write_elements_h5(path, mesh, nodes)
        if suffix == ".cgns":
            return write_cgns_elements(path, mesh, base_name)
        return writers[suffix](path, mesh)

    def curved_high_order(
        self, names=None, order: int = 2, nodes: str = "equispaced", snap: bool = True
    ):
        """Curved order-``order`` elements of structured domains (2D) or blocks (3D): p x p
        (x p) grid cells per element (interval counts must be multiples of p), the grid's
        nodes are the equispaced node set, so element edges follow the grid lines and curved
        walls. ``nodes``: "equispaced" (file formats), "gll" (Gauss-Lobatto-Legendre:
        continuous Galerkin / spectral elements, DGSEM with boundary points) or "gauss"
        (Gauss-Legendre: DG solution points); ``snap``: project GLL nodes on boundary sides onto
        their curves (surfaces in 3D). Returns a HighOrderMesh (``xyz`` (elements, (p + 1)^d,
        3) in tensor order, ``scaled_jacobian``, ``merged()``, ``owners``)."""
        from opengrid.model.highorder import high_order

        return high_order(self, names, order, nodes, snap)

    def export_curved_high_order(
        self, path, names=None, order: int = 2, nodes: str = "equispaced", snap: bool = True
    ) -> dict:
        """Write ``curved_high_order(...)``: ``.vtkhdf`` (VTK Lagrange cells, equispaced nodes)
        or ``.h5`` (any node set: element node coordinates, Jacobians, merged numbering)."""
        from opengrid.io.highorder import write_high_order_h5, write_vtkhdf_high_order

        mesh = self.curved_high_order(names, order, nodes, snap)
        if str(path).endswith(".h5"):
            return write_high_order_h5(path, mesh)
        return write_vtkhdf_high_order(path, mesh)

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
    def modify(self, name: str, *, new_name: str | None = None, **params) -> str:
        """Change some parameters of an entity, e.g. ``og.modify("line-1", p1=(2, 0))``, and /
        or rename it (``new_name``) as one undo step; returns the (new) name. A sketch:
        edit_sketch, which regenerates it."""
        if name in self.store and self.store.get(name).kind == "sketch":
            if params:
                name = self.edit_sketch(name, **params)
            if new_name is not None:
                return self.rename(name, new_name)
            return name
        return self.commands.run(ModifyEntity(name, params, new_name))

    def rename(self, name: str, new_name: str) -> str:
        """Rename an entity (references follow; journals stay valid)."""
        return self.commands.run(RenameEntity(name, new_name))

    def delete(self, names, dependents: bool = False) -> int:
        """Delete entities (a name or a list). What depends on them (the connectors on a
        curve, the domains on a connector, the blocks and volumes on a domain, ...) goes too
        with ``dependents=True``; without it such a deletion is refused naming them (see
        :meth:`deletion_dependents`). Patches, zones and smoothing groups never block: they
        lose the members (a group left empty goes, a periodic patch losing its partner is
        unpaired). A sketch takes its outputs along; an output on its own is refused (remove
        its block with :meth:`remove_sketch_blocks`). Returns the number deleted."""
        if isinstance(names, str):
            names = [names]
        return self.commands.run(DeleteEntities(list(names), dependents=dependents))

    def deletion_dependents(self, names) -> list[str]:
        """The entities that deleting ``names`` would take along (``delete(...,
        dependents=True)``), in deletion order; raises for what cannot be deleted at all."""
        from opengrid.commands.geometry import delete_plan

        if isinstance(names, str):
            names = [names]
        _expanded, extra, issues = delete_plan(self.store, list(names))
        if issues:
            raise ValidationError(issues)
        return extra

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

    def move_handle(self, name: str, field: str, index: int, xyz) -> str:
        """Move one defining point of a point or curve (``field``: "xyz", "p0"/"p1",
        "points" or "control_points"; ``index`` -1 for single points), as a drag in the GUI
        does: crossings of the curve with other curves are solved again and the connectors
        meeting there follow, so domains built from a curve network stay closed. One undo
        step."""
        return self._topology("move_handle", f"Move {name}", name, field, index, tuple(xyz))[0]

    def insert_handle(self, name: str, xyz) -> str:
        """Add a defining point to a polyline or spline (on the curve where it is closest to
        ``xyz``, between the neighbouring points) or a control point to a NURBS curve (a knot
        inserted at the closest parameter: the curve does not change), as a double-click on
        a curve does in the GUI. Connectors on the curve keep their places and crossings are
        solved again. One undo step."""
        return self._topology("insert_handle", f"Add a point to {name}", name, tuple(xyz))[0]

    def remove_handle(self, name: str, index: int) -> str:
        """Remove defining point ``index`` of a polyline or spline, or control point ``index``
        of a NURBS curve (knot removal: exact when the curve allows it, otherwise
        ``last_warnings`` says how far the curve moved; the end control points cannot be
        removed), as a double-click on a handle does in the GUI. One undo step."""
        return self._topology("remove_handle", f"Remove a point of {name}", name, int(index))[0]

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
        self.groups = GroupEvaluator(store, self.connectors, self.domains)
        self.domains.groups = self.groups
        self.unstructured = UnstructuredEvaluator(store, self.connectors)
        self.blocks = BlockEvaluator(store, self.domains)
        self.blocks.groups = self.groups
        self.groups.blocks = self.blocks
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
