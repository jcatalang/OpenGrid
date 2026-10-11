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
from opengrid.model.connectors import ConnectorEvaluator, ConnectorNodes, order_chain
from opengrid.model.domains import DomainEvaluator, SmoothingResult, assemble_edges
from opengrid.model.store import Entity, Store
from opengrid.model.unstructured import UnstructuredEvaluator, assemble_loops
from opengrid.validation import Issue, ValidationError

IMPORT_MODES = ("spline", "split", "polyline")


class Session:
    def __init__(self):
        self.store = Store()
        self.commands = CommandStack(self.store)
        self.connectors = ConnectorEvaluator(self.store)
        self.domains = DomainEvaluator(self.store, self.connectors)
        self.unstructured = UnstructuredEvaluator(self.store, self.connectors)
        self.commands.add_listener(self._invalidate_derived)
        self.project_path: Path | None = None
        # GUI/report threshold for the stretching-ratio warning (not stored in projects).
        self.stretching_threshold = 1.3

    # ------------------------------------------------------------- creation
    def create_point(self, xyz, name: str | None = None) -> str:
        return self._create("point", {"xyz": xyz}, name)

    def create_line(self, p0, p1, name: str | None = None) -> str:
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
        params = {
            "degree": degree,
            "knots": knots,
            "control_points": control_points,
            "weights": weights,
        }
        return self._create("nurbs", params, name)

    # ----------------------------------------------------------- connectors
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
        return self.modify(name, dimension=dimension)

    def set_distribution(
        self,
        name: str,
        type: str,
        begin_spacing: float | None = None,
        end_spacing: float | None = None,
        growth_ratio: float | None = None,
    ) -> str:
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
        self, connectors=None, edges=None, generator=None, name: str | None = None
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
        return self.create_domain(edges, generator, name=name)

    def create_domain(self, edges, generator=None, name: str | None = None) -> str:
        """Structured domain from explicit edges (the form recorded in journals)."""
        params = {"edges": edges, "generator": generator}
        return self._create("domain", params, name)

    def create_unstructured_domain(
        self, loops=None, generator=None, name: str | None = None, connectors=None
    ) -> str:
        """Unstructured (triangle) domain.

        Either ``loops``: the outer loop (counter-clockwise) followed by holes (clockwise), each a
        list of connector names or ``(name, reversed)`` items; or ``connectors``: a set of
        connector names assembled into closed loops automatically (largest loop = outer).
        ``generator``: ``{"min_angle": 25, "growth_rate": 1.2, "max_size": None, "sources": [],
        "max_points": 2000000, "refine": True}``; sources are ``{"p0", "p1" (optional; line
        source), "size", "radius"}``.
        """
        if loops is None:
            if connectors is None:
                raise ValueError("give either loops or connectors")
            if isinstance(connectors, str):
                connectors = [connectors]
            assembled = assemble_loops(self.store, self.connectors, list(connectors))
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
        return self._create("unstructured_domain", {"loops": loops, "generator": generator}, name)

    def set_unstructured_generator(self, name: str, **changes) -> str:
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
        settings; ``elliptic=None`` switches smoothing off)."""
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
        entity = self.store.get(name)
        if entity.kind == "unstructured_domain":
            mesh = self.unstructured_mesh(name)
            text = f"{len(mesh['triangles'])} triangles, {len(mesh['points'])} points"
            if mesh["rejected"]:
                text += f" ({mesh['rejected']} refinement points rejected near the boundary)"
            return text
        return self.domains.status(self._domain_id(name))

    def domain_quality(self, name: str) -> dict:
        """Per-cell quality metrics: arrays of shape (nj - 1, ni - 1) for structured domains,
        (1, m) for unstructured ones."""
        entity = self.store.get(name)
        if entity.kind == "unstructured_domain":
            mesh = self.unstructured_mesh(name)
            return _core.quality_triangles(
                np.ascontiguousarray(mesh["points"][:, :2]), mesh["triangles"]
            )
        return _core.quality_2d(self.domain_grid(name))

    def _export_names(self, domains) -> list[str]:
        if domains is None:
            return self.names("domain") + self.names("unstructured_domain")
        if isinstance(domains, str):
            return [domains]
        return list(domains)

    # ------------------------------------------------------ patches and zones
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
        return self.modify(source, periodic=None)

    def set_patch_type(self, name: str, patch_type: str) -> str:
        return self.modify(name, patch_type=patch_type)

    def add_to_patch(self, name: str, connectors) -> str:
        return self._change_members(name, "connector", connectors, add=True)

    def remove_from_patch(self, name: str, connectors) -> str:
        return self._change_members(name, "connector", connectors, add=False)

    def create_zone(self, name: str | None = None, members=()) -> str:
        """Named group of domains (volume tag)."""
        return self._create("zone", {"members": members}, name)

    def add_to_zone(self, name: str, domains) -> str:
        return self._change_members(name, "domain", domains, add=True)

    def remove_from_zone(self, name: str, domains) -> str:
        return self._change_members(name, "domain", domains, add=False)

    def _change_members(self, name: str, key: str, refs, add: bool) -> str:
        if isinstance(refs, (str, int)):
            refs = [refs]
        current = [self.store.get(item[key]).name for item in self.entity(name).params["members"]]
        names = [self.store.get(ref).name for ref in refs]
        if add:
            members = current + [n for n in names if n not in current]
        else:
            missing = [n for n in names if n not in current]
            if missing:
                raise ValidationError([Issue(name, f"does not contain {', '.join(missing)}")])
            members = [n for n in current if n not in names]
        return self.modify(name, members=members)

    def topology(self):
        from opengrid.model.topology import analyze

        return analyze(self.store, self)

    def boundary_report(self):
        """Boundary connectors, interfaces (detected automatically), unused connectors, boundary
        connectors without a patch and all problems (opengrid.model.topology.BoundaryReport)."""
        from opengrid.model.topology import validate_patches

        return validate_patches(self.store, self.topology(), self)

    def unassigned_boundary_connectors(self) -> list[str]:
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
        return write_gmsh(
            path, assemble(self, self.topology(), quality=False, periodic=report.periodic)
        )

    def export_cgns(self, path, base_name: str = "OpenGrid") -> dict:
        """Structured CGNS: zones, families, BCs per patch and 1-to-1 connectivity; the boundary
        must validate."""
        from opengrid.io.cgns import write_cgns

        report = self._require_valid_boundary()
        return write_cgns(path, self, self.topology(), base_name, periodic=report.periodic)

    def export_plot3d(self, path, domains=None, binary: bool = True, three_d: bool = False) -> None:
        """Write domains as a multi-grid Plot3D file (2D, or 3D with nk = 1)."""
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
        return self.commands.run(RenameEntity(name, new_name))

    def delete(self, names) -> int:
        if isinstance(names, str):
            names = [names]
        return self.commands.run(DeleteEntities(list(names)))

    def undo(self) -> None:
        self.commands.undo()

    def redo(self) -> None:
        self.commands.redo()

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

    # -------------------------------------------------------------- queries
    def entity(self, name: str) -> Entity:
        return self.store.get(name)

    def names(self, kind: str | None = None) -> list[str]:
        return [entity.name for entity in self.store.entities(kind)]

    def geometry(self, name: str):
        """The C++ curve of an entity (None for points)."""
        return self.store.geometry(name)

    def length(self, name: str) -> float:
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
        project.save(self.store, path)
        self.project_path = Path(path)

    def load(self, path) -> None:
        """Replace the model by a project file. Clears undo history; journal restarts from it."""
        store = project.load(path)
        self._replace_store(store)
        self.project_path = Path(path)
        self.commands.journal.append(f"og.load({str(path)!r})")

    def new(self) -> None:
        self._replace_store(Store())
        self.project_path = None

    def journal_text(self) -> str:
        return self.commands.journal.text()

    def save_journal(self, path) -> None:
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
