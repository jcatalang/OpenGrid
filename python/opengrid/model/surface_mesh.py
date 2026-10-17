"""Surface meshes of database surfaces, so a closed set of surfaces (a box, a cylinder, a
wing) becomes a watertight shell of domains for a volume.

Method
- Each surface's boundary is its four parametric sides: a plane's four edge lines; a ruled
  surface's two curves and the straight lines joining their ends; an extruded surface's curve,
  its translate and the two extrusion lines; a surface of revolution's generating curve at the
  start and end angles (rotated copies) and the circular arcs swept by the curve's ends. Sides
  of zero length (a ruled cap closing at the leading edge, a disc's centre) are dropped.
- Every side is matched against the existing curves by sampling nine points (either way
  round), and a full circle by existing arcs tiling it; a side with no match becomes a new
  curve, so neighbouring surfaces share their edge curves.
- A connector lying on a side's curve alone is reused (the neighbours then share it, which is
  what makes the shell watertight); otherwise one is made with the points of ``spacing`` (or
  ``points``).
- Flat surfaces (planes, flat caps, discs) are meshed as planar domains assembled from the
  sides' connectors; curved ones as domains on the surface, the loop in parametric order
  (reversed when the frame reports it clockwise).
Not meshable here: a full turn of a curved surface of revolution (split it in halves, as
create_cylinder does), a surface of revolution closing at its axis on both ends (a sphere's
poles), NURBS surfaces.
"""

import math
from dataclasses import dataclass

import numpy as np

from opengrid.model.connectors import CURVE_KINDS, ConnectorEvaluator
from opengrid.model.entities import normalize_params
from opengrid.validation import Issue, ValidationError

SAMPLES = 9
FULL_TURN = 2.0 * math.pi


@dataclass(frozen=True)
class Side:
    """One boundary side of a surface: ``spec`` is ("entity", id), ("line", p0, p1),
    ("arc", center, radius, normal, x_axis, a0, a1) or ("copy", id, transform) where the
    transform maps points of the source curve; ``forward``: traversed along the spec's own
    direction in the surface's loop."""

    spec: tuple
    forward: bool


def _unit(v):
    v = np.asarray(v, dtype=float)
    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError("zero vector")
    return v / norm


def _rotation(axis_point, axis, angle: float):
    """A point map: rotation by ``angle`` about the axis through ``axis_point``."""
    p0 = np.asarray(axis_point, dtype=float)
    w = _unit(axis)
    c = math.cos(angle)
    s = math.sin(angle)

    def rotate(point):
        r = np.asarray(point, dtype=float) - p0
        out = r * c + np.cross(w, r) * s + w * float(np.dot(w, r)) * (1.0 - c)
        return tuple(float(x) for x in out + p0)

    def rotate_vector(vector):
        r = np.asarray(vector, dtype=float)
        out = r * c + np.cross(w, r) * s + w * float(np.dot(w, r)) * (1.0 - c)
        return tuple(float(x) for x in out)

    return rotate, rotate_vector


def _translation(direction):
    d = np.asarray(direction, dtype=float)

    def move(point):
        return tuple(float(x) for x in np.asarray(point, dtype=float) + d)

    def keep(vector):
        return tuple(float(x) for x in vector)

    return move, keep


def transformed_curve_params(entity, move, move_vector) -> dict:
    """The raw parameters of ``entity`` (a curve) with its defining points moved by
    ``move`` and its direction vectors by ``move_vector``."""
    p = entity.params
    if entity.kind == "line":
        return {"p0": move(p["p0"]), "p1": move(p["p1"])}
    if entity.kind == "arc":
        return {
            "center": move(p["center"]),
            "radius": p["radius"],
            "normal": move_vector(p["normal"]),
            "x_axis": move_vector(p["x_axis"]),
            "start_angle": p["start_angle"],
            "end_angle": p["end_angle"],
        }
    if entity.kind == "polyline":
        return {"points": [move(q) for q in p["points"]]}
    if entity.kind == "spline":
        out = {"points": [move(q) for q in p["points"]], "method": p["method"]}
        for key in ("start_tangent", "end_tangent"):
            out[key] = None
            if p.get(key) is not None:
                out[key] = move_vector(p[key])
        return out
    if entity.kind == "nurbs":
        return {
            "degree": p["degree"],
            "knots": p["knots"],
            "control_points": [move(q) for q in p["control_points"]],
            "weights": p["weights"],
        }
    raise ValidationError([Issue(entity.name, f"a {entity.kind} cannot be copied as a curve")])


def _arc_point(center, radius, normal, x_axis, angle):
    c = np.asarray(center, dtype=float)
    x = _unit(x_axis)
    y = np.cross(_unit(normal), x)
    return c + radius * (math.cos(angle) * x + math.sin(angle) * y)


def samples(store, spec, count: int = SAMPLES) -> np.ndarray:
    """``count`` points along a side spec, in its own direction."""
    kind = spec[0]
    ts = np.linspace(0.0, 1.0, count)
    if kind == "entity":
        curve = store.geometry(store.get(spec[1]).name)
        return np.asarray([curve.eval(float(t)) for t in ts], dtype=float)
    if kind == "line":
        p0 = np.asarray(spec[1], dtype=float)
        p1 = np.asarray(spec[2], dtype=float)
        return np.asarray([p0 + t * (p1 - p0) for t in ts])
    if kind == "arc":
        _k, center, radius, normal, x_axis, a0, a1 = spec
        return np.asarray(
            [_arc_point(center, radius, normal, x_axis, a0 + t * (a1 - a0)) for t in ts]
        )
    if kind == "copy":
        move = spec[2][0]
        return np.asarray([move(q) for q in samples(store, ("entity", spec[1]), count)])
    raise ValueError(f"unknown side spec {kind!r}")


def is_flat(surface) -> bool:
    grid = np.asarray(
        [surface.eval(u, v) for u in np.linspace(0, 1, 6) for v in np.linspace(0, 1, 6)],
        dtype=float,
    )
    size = float(np.max(np.ptp(grid, axis=0))) or 1.0
    centred = grid - grid.mean(axis=0)
    return float(np.linalg.svd(centred, compute_uv=False)[-1]) <= 1e-7 * size * math.sqrt(len(grid))


def surface_sides(store, entity) -> tuple[list[Side], bool, bool]:
    """(sides in loop order, flat, full turn) of a surface. Degenerate sides are left out."""
    p = entity.params
    sides: list[Side] = []
    full_turn = False
    if entity.kind == "plane":
        o = np.asarray(p["origin"], dtype=float)
        u = np.asarray(p["u_vector"], dtype=float)
        v = np.asarray(p["v_vector"], dtype=float)
        corners = [o, o + u, o + u + v, o + v]
        for k in range(4):
            a = tuple(float(x) for x in corners[k])
            b = tuple(float(x) for x in corners[(k + 1) % 4])
            sides.append(Side(("line", a, b), True))
    elif entity.kind == "ruled_surface":
        g0 = store.geometry(store.get(p["curve0"]).name)
        g1 = store.geometry(store.get(p["curve1"]).name)
        sides.append(Side(("entity", p["curve0"]), True))
        sides.append(Side(("line", tuple(g0.eval(1.0)), tuple(g1.eval(1.0))), True))
        sides.append(Side(("entity", p["curve1"]), False))
        sides.append(Side(("line", tuple(g0.eval(0.0)), tuple(g1.eval(0.0))), False))
    elif entity.kind == "extruded_surface":
        g = store.geometry(store.get(p["curve"]).name)
        d = np.asarray(p["direction"], dtype=float)
        end0 = np.asarray(g.eval(0.0), dtype=float)
        end1 = np.asarray(g.eval(1.0), dtype=float)
        sides.append(Side(("entity", p["curve"]), True))
        sides.append(Side(("line", tuple(end1), tuple(end1 + d)), True))
        sides.append(Side(("copy", p["curve"], _translation(d)), False))
        sides.append(Side(("line", tuple(end0), tuple(end0 + d)), False))
    elif entity.kind == "revolved_surface":
        g = store.geometry(store.get(p["curve"]).name)
        axis_point = np.asarray(p["axis_point"], dtype=float)
        w = _unit(p["axis"])
        a0 = float(p["start_angle"])
        a1 = float(p["end_angle"])
        full_turn = abs(abs(a1 - a0) - FULL_TURN) < 1e-9
        ends = [np.asarray(g.eval(0.0), dtype=float), np.asarray(g.eval(1.0), dtype=float)]
        arcs = []
        for end in ends:
            centre = axis_point + float(np.dot(end - axis_point, w)) * w
            radius = float(np.linalg.norm(end - centre))
            arcs.append((centre, radius, end))
        scale = max(1.0, float(np.max(np.abs(np.asarray([g.eval(t) for t in (0, 0.5, 1)])))))
        for k, (centre, radius, end) in enumerate(arcs):
            if radius <= 1e-9 * scale:
                arcs[k] = None
                continue
            x_axis = tuple(float(x) for x in (end - centre) / radius)
            spec = ("arc", tuple(float(x) for x in centre), radius, tuple(w), x_axis, a0, a1)
            arcs[k] = spec
        if arcs[0] is not None:
            sides.append(Side(arcs[0], True))
        if not full_turn:
            sides.append(Side(("copy", p["curve"], _rotation(axis_point, w, a1)), True))
        if arcs[1] is not None:
            sides.append(Side(arcs[1], False))
        if not full_turn:
            if a0 == 0.0:
                sides.append(Side(("entity", p["curve"]), False))
            else:
                sides.append(Side(("copy", p["curve"], _rotation(axis_point, w, a0)), False))
    else:
        raise ValidationError(
            [Issue(entity.name, f"a {entity.kind} cannot be meshed here (not yet supported)")]
        )
    flat = is_flat(store.geometry(entity.name))
    kept = []
    for side in sides:
        pts = samples(store, side.spec, 3)
        if float(np.linalg.norm(pts[-1] - pts[0])) + float(np.linalg.norm(pts[1] - pts[0])) > 0:
            kept.append(side)
    return kept, flat, full_turn


def _scale(store) -> float:
    scale = 1.0
    for entity in store.entities("point"):
        scale = max(scale, float(np.max(np.abs(entity.params["xyz"]))))
    for entity in store:
        if entity.kind in CURVE_KINDS:
            pts = samples(store, ("entity", entity.id), 3)
            scale = max(scale, float(np.max(np.abs(pts))))
    return scale


def find_curves(store, spec, tolerance: float) -> list[tuple[str, bool]]:
    """Existing curves covering the side: [(name, reversed relative to the spec)], usually
    one; several arcs for a full circle they tile; [] when none."""
    wanted = samples(store, spec)
    for entity in store:
        if entity.kind not in CURVE_KINDS:
            continue
        have = samples(store, ("entity", entity.id))
        if float(np.max(np.linalg.norm(have - wanted, axis=1))) <= tolerance:
            return [(entity.name, False)]
        if float(np.max(np.linalg.norm(have[::-1] - wanted, axis=1))) <= tolerance:
            return [(entity.name, True)]
    if spec[0] == "arc" and abs(abs(spec[6] - spec[5]) - FULL_TURN) < 1e-9:
        return _tiling_arcs(store, spec, tolerance)
    return []


def _tiling_arcs(store, spec, tolerance: float) -> list[tuple[str, bool]]:
    """Existing arcs about the same centre, radius and axis whose spans tile the full
    circle (the two halves of a rim)."""
    _k, center, radius, normal, x_axis, _a0, _a1 = spec
    centre = np.asarray(center, dtype=float)
    w = _unit(normal)
    x = _unit(x_axis)
    y = np.cross(w, x)
    pieces = []
    for entity in store.entities("arc"):
        p = entity.params
        if float(np.linalg.norm(np.asarray(p["center"]) - centre)) > tolerance:
            continue
        if abs(p["radius"] - radius) > tolerance:
            continue
        if abs(abs(float(np.dot(_unit(p["normal"]), w))) - 1.0) > 1e-9:
            continue
        pts = samples(store, ("entity", entity.id), 5)
        angles = []
        for q in pts:
            r = q - centre
            angles.append(math.atan2(float(np.dot(r, y)), float(np.dot(r, x))))
        unwrapped = np.unwrap(angles)
        pieces.append((float(min(unwrapped)), float(max(unwrapped)), entity.name))
    if not pieces:
        return []
    pieces.sort()
    covered = 0.0
    start = pieces[0][0]
    cursor = start
    names = []
    for lo, hi, name in pieces:
        if lo > cursor + 1e-6:
            return []
        covered += max(0.0, hi - max(lo, cursor))
        cursor = max(cursor, hi)
        names.append((name, False))
    if covered < FULL_TURN - 1e-6:
        return []
    return names


def make_curve(store, spec, name: str) -> str:
    """A new curve entity for a side spec."""
    kind = spec[0]
    if kind == "line":
        raw = {"p0": spec[1], "p1": spec[2]}
        kind = "line"
    elif kind == "arc":
        _k, center, radius, normal, x_axis, a0, a1 = spec
        raw = {
            "center": center,
            "radius": radius,
            "normal": normal,
            "x_axis": x_axis,
            "start_angle": a0,
            "end_angle": a1,
        }
        kind = "arc"
    elif kind == "copy":
        source = store.get(spec[1])
        move, move_vector = spec[2]
        raw = transformed_curve_params(source, move, move_vector)
        kind = source.kind
    else:
        raise ValueError(f"a {kind} side needs no new curve")
    unique = name
    number = 1
    while unique in store:
        number += 1
        unique = f"{name}-{number}"
    return store.add(kind, normalize_params(kind, raw, unique, store), unique).name


def connector_on(store, connectors: ConnectorEvaluator, curve: str, spacing, points) -> str:
    """The connector lying on ``curve`` alone, else a new one with ``points`` or the points
    of ``spacing`` along it."""
    curve_id = store.get(curve).id
    for entity in store.entities("connector"):
        segments = entity.params["segments"]
        if len(segments) == 1 and segments[0]["curve"] == curve_id:
            return entity.name
    if points is None:
        length = float(connectors_length(store, curve))
        points = max(2, round(length / spacing) + 1)
    name = f"c-{curve}"
    number = 1
    while name in store:
        number += 1
        name = f"c-{curve}-{number}"
    raw = {"segments": curve, "dimension": int(points), "distribution": "uniform"}
    return store.add("connector", normalize_params("connector", raw, name, store), name).name


def connectors_length(store, curve: str) -> float:
    pts = samples(store, ("entity", store.get(curve).id), 65)
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))


def mesh_surfaces(store, names, spacing=None, points=None, generator=None):
    """Unstructured domains ``<surface>-mesh`` on the given surfaces, their boundary curves
    and connectors found or made (shared between neighbours). Returns (domain names,
    warnings). Curved surfaces first, so the flat ones (caps) find the rims they made."""
    from opengrid.model.unstructured import assemble_loops

    if isinstance(names, str):
        names = [names]
    if spacing is None and points is None:
        raise ValidationError([Issue("mesh_surfaces", "give a spacing or a number of points")])
    if spacing is not None and (not isinstance(spacing, (int, float)) or spacing <= 0):
        raise ValidationError([Issue("spacing", f"must be a positive length, got {spacing!r}")])
    if points is not None and (
        isinstance(points, bool) or not isinstance(points, int) or points < 2
    ):
        raise ValidationError([Issue("points", f"must be a whole number >= 2, got {points!r}")])
    surfaces = []
    for name in names:
        entity = store.get(name)
        if entity.kind not in ("plane", "ruled_surface", "extruded_surface", "revolved_surface"):
            raise ValidationError(
                [
                    Issue(
                        name,
                        f"is a {entity.kind}: mesh_surfaces takes planes, ruled, extruded "
                        "and revolved surfaces",
                    )
                ]
            )
        sides, flat, full_turn = surface_sides(store, entity)
        if not sides:
            reason = "has no boundary sides of any length"
            if entity.kind == "revolved_surface":
                reason = (
                    "closes on its axis at both ends (a sphere): its poles leave no boundary "
                    "to mesh from; a sphere surface mesh is not available yet"
                )
            raise ValidationError([Issue(name, reason)])
        if full_turn and not flat:
            raise ValidationError(
                [
                    Issue(
                        name,
                        "is a full turn of a curved surface of revolution: split it in two "
                        "halves (create_cylinder does), one domain cannot wrap around",
                    )
                ]
            )
        surfaces.append((entity, sides, flat))
    surfaces.sort(key=lambda item: item[2])  # curved first
    tolerance = 1e-6 * _scale(store)
    connectors = ConnectorEvaluator(store)
    settings = {}
    if generator is not None:
        settings.update(generator)
    if spacing is not None and "max_size" not in settings:
        settings["max_size"] = float(spacing)
    made = []
    warnings = []
    for entity, sides, flat in surfaces:
        loop = []  # (connector name, reversed) in loop order
        for k, side in enumerate(sides):
            if side.spec[0] == "entity":
                curves = [(store.get(side.spec[1]).name, False)]
            else:
                curves = find_curves(store, side.spec, tolerance)
                if not curves:
                    curves = [(make_curve(store, side.spec, f"{entity.name}-e{k}"), False)]
            if not side.forward:
                curves = [(c, not r) for c, r in reversed(curves)]
            for curve, reversed_ in curves:
                connector = connector_on(store, connectors, curve, spacing, points)
                connectors.invalidate()
                loop.append((connector, reversed_))
        domain_name = f"{entity.name}-mesh"
        number = 1
        while domain_name in store:
            number += 1
            domain_name = f"{entity.name}-mesh-{number}"
        if flat:
            assembled = assemble_loops(store, connectors, [c for c, _r in loop], None)
            loops = [
                [(store.get(item["connector"]).name, item["reversed"]) for item in one]
                for one in assembled
            ]
            raw = {"loops": loops, "generator": settings, "surface": None}
            params = normalize_params("unstructured_domain", raw, domain_name, store)
        else:
            raw = {"loops": [loop], "generator": settings, "surface": entity.name}
            try:
                params = normalize_params("unstructured_domain", raw, domain_name, store)
            except ValidationError as error:
                if not any("counter-clockwise" in issue.message for issue in error.issues):
                    raise
                raw["loops"] = [[(c, not r) for c, r in reversed(loop)]]
                params = normalize_params("unstructured_domain", raw, domain_name, store)
        made.append(store.add("unstructured_domain", params, domain_name).name)
    return made, warnings
