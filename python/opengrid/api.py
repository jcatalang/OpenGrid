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

from opengrid import project
from opengrid.commands import (
    CommandStack,
    CreateEntities,
    CreateEntity,
    DeleteEntities,
    ModifyEntity,
    RenameEntity,
)
from opengrid.commands.base import JOURNAL_HEADER
from opengrid.io import read_airfoil_dat, read_points_csv
from opengrid.model.store import Entity, Store
from opengrid.validation import Issue, ValidationError

IMPORT_MODES = ("spline", "split", "polyline")


class Session:
    def __init__(self):
        self.store = Store()
        self.commands = CommandStack(self.store)
        self.project_path: Path | None = None

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

    def _replace_store(self, store: Store) -> None:
        listeners = list(self.commands._listeners)
        self.store = store
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
