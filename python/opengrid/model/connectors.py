"""Connectors: parameter schema, validation, curve chains and node evaluation.

A connector is a chain of curve pieces with a number of points ("dimension") and a 1D
distribution. Besides the stretching functions, ``{"type": "explicit", "positions": (...)}``
gives the node positions directly (fractions of the length, 0 to 1, one per point); splitting
and joining connectors produce it so that the new connectors keep the old nodes. Breakpoints
split it into sub-connectors at given positions and node indices; the top-level distribution
applies to the first sub-connector and each breakpoint's distribution to the sub-connector
after it. Spacings are absolute lengths.

Spacing propagation (``propagate_spacing``): an end of a tanh/sine (sub-)connector whose
spacing is unspecified adopts the explicit spacing that other connectors prescribe at the same
node (the smallest one if several). Propagation is one level deep (explicit -> unspecified),
so the result never depends on evaluation order.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.validation import Issue, ValidationError

DISTRIBUTION_TYPES = ("uniform", "geometric", "tanh", "sine")  # the stretching functions
EXPLICIT = "explicit"  # given positions (from split / join)
SPACING_MODES = ("arc_length", "parametric")
CURVE_KINDS = ("line", "arc", "polyline", "spline", "nurbs")
SPACING_TYPES = ("tanh", "sine")  # types whose free ends can adopt propagated spacings
CHAIN_TOLERANCE = 1e-9  # relative to the largest curve length (and at least absolute 1e-9)
DEFAULT_DISTRIBUTION = {
    "type": "uniform",
    "begin_spacing": None,
    "end_spacing": None,
    "growth_ratio": None,
}


def _number(value, field_name, issues, where, positive=False, optional=False):
    if value is None and optional:
        return None
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        value = value.item()
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        issues.append(Issue(f"{where}: {field_name}", f"expected a finite number, got {value!r}"))
        return None
    if positive and value <= 0:
        issues.append(Issue(f"{where}: {field_name}", f"must be positive, got {value!r}"))
        return None
    return float(value)


def normalize_distribution(raw, where: str, field_name: str = "distribution") -> dict | None:
    """Canonical distribution dict; raises ValidationError."""
    issues: list[Issue] = []
    if isinstance(raw, str):
        raw = {"type": raw}
    if not isinstance(raw, dict):
        raise ValidationError(
            [Issue(f"{where}: {field_name}", f"expected a dict or type name, got {raw!r}")]
        )
    allowed = set(DEFAULT_DISTRIBUTION)
    kind = raw.get("type", "uniform")
    if kind == EXPLICIT:
        allowed.add("positions")
    unknown = sorted(set(raw) - allowed)
    for key in unknown:
        issues.append(Issue(f"{where}: {field_name}.{key}", "unexpected key"))
    if kind not in DISTRIBUTION_TYPES and kind != EXPLICIT:
        issues.append(
            Issue(
                f"{where}: {field_name}.type",
                f"must be one of {', '.join((*DISTRIBUTION_TYPES, EXPLICIT))}, got {kind!r}",
            )
        )
    result = {"type": kind}
    for key in ("begin_spacing", "end_spacing", "growth_ratio"):
        result[key] = _number(
            raw.get(key), f"{field_name}.{key}", issues, where, positive=True, optional=True
        )
    if issues:
        raise ValidationError(issues)
    given = [
        key for key in ("begin_spacing", "end_spacing", "growth_ratio") if result[key] is not None
    ]
    if kind == "uniform" and given:
        issues.append(Issue(f"{where}: {field_name}", f"uniform takes no {', '.join(given)}"))
    if kind == "geometric" and len(given) != 1:
        issues.append(
            Issue(
                f"{where}: {field_name}",
                "geometric needs exactly one of begin_spacing, end_spacing, growth_ratio",
            )
        )
    if kind in SPACING_TYPES and result["growth_ratio"] is not None:
        issues.append(Issue(f"{where}: {field_name}.growth_ratio", f"not used by {kind}"))
    if kind == EXPLICIT:
        if given:
            issues.append(Issue(f"{where}: {field_name}", f"explicit takes no {', '.join(given)}"))
        positions = _positions(raw.get("positions"), f"{where}: {field_name}.positions", issues)
        result["positions"] = positions
    if issues:
        raise ValidationError(issues)
    return result


def _positions(raw, label: str, issues: list[Issue]) -> tuple | None:
    """Explicit node positions: at least two finite, strictly increasing numbers from 0 to 1."""
    if not isinstance(raw, (list, tuple, np.ndarray)) or len(raw) < 2:
        issues.append(Issue(label, "expected a list of at least 2 positions from 0 to 1"))
        return None
    values = []
    for item in raw:
        if hasattr(item, "item") and not isinstance(item, (int, float)):
            item = item.item()
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
            issues.append(Issue(label, f"expected numbers, got {item!r}"))
            return None
        values.append(float(item))
    if values[0] != 0.0 or values[-1] != 1.0:
        issues.append(Issue(label, "must start at 0 and end at 1"))
        return None
    for k in range(1, len(values)):
        if not values[k] > values[k - 1]:
            issues.append(Issue(label, f"must increase strictly (position {k})"))
            return None
    return tuple(values)


def _segment(raw, index: int, store, where: str, issues: list[Issue]) -> dict | None:
    label = f"{where}: segments[{index}]"
    t0 = 0.0
    t1 = 1.0
    if isinstance(raw, dict):
        curve = raw.get("curve")
        t0 = raw.get("t0", 0.0)
        t1 = raw.get("t1", 1.0)
    elif isinstance(raw, (list, tuple)) and len(raw) == 3:
        curve, t0, t1 = raw
    else:
        curve = raw
    if isinstance(curve, bool) or not isinstance(curve, (str, int)):
        issues.append(Issue(label, f"expected a curve name or id, got {curve!r}"))
        return None
    if curve not in store:
        issues.append(Issue(label, f"no entity {curve!r}"))
        return None
    entity = store.get(curve)
    if entity.kind not in CURVE_KINDS:
        issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a curve"))
        return None
    t0 = _number(t0, "t0", issues, label)
    t1 = _number(t1, "t1", issues, label)
    if t0 is None or t1 is None:
        return None
    if not (0.0 <= t0 <= 1.0 and 0.0 <= t1 <= 1.0) or t0 == t1:
        issues.append(Issue(label, f"t0, t1 must lie in [0, 1] and differ (got {t0}, {t1})"))
        return None
    return {"curve": entity.id, "t0": t0, "t1": t1}


def normalize_connector(raw: dict, store, where: str) -> dict:
    """Validate connector parameters against the store; returns canonical parameters."""
    issues: list[Issue] = []
    segments_raw = raw["segments"]
    if isinstance(segments_raw, (str, int, dict)):
        segments_raw = [segments_raw]
    segments = []
    for index, item in enumerate(segments_raw):
        segment = _segment(item, index, store, where, issues)
        segments.append(segment)
    if not segments:
        issues.append(Issue(f"{where}: segments", "at least one curve is required"))
    dimension = raw["dimension"]
    if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 2:
        issues.append(Issue(f"{where}: dimension", f"must be an integer >= 2, got {dimension!r}"))
    mode = raw.get("mode", "arc_length")
    if mode not in SPACING_MODES:
        issues.append(
            Issue(f"{where}: mode", f"must be one of {', '.join(SPACING_MODES)}, got {mode!r}")
        )
    propagate = raw.get("propagate_spacing", True)
    if not isinstance(propagate, bool):
        issues.append(
            Issue(f"{where}: propagate_spacing", f"must be true or false, got {propagate!r}")
        )
    try:
        distribution = normalize_distribution(raw.get("distribution", "uniform"), where)
    except ValidationError as error:
        issues.extend(error.issues)
        distribution = None
    breakpoints = []
    for index, item in enumerate(raw.get("breakpoints", ())):
        label = f"{where}: breakpoints[{index}]"
        if not isinstance(item, dict) or "position" not in item or "index" not in item:
            issues.append(
                Issue(label, "needs 'position' and 'index' (and optionally 'distribution')")
            )
            continue
        position = _number(item["position"], "position", issues, label)
        node = item["index"]
        if isinstance(node, bool) or not isinstance(node, int):
            issues.append(Issue(label, f"index must be an integer, got {node!r}"))
            continue
        try:
            dist = normalize_distribution(item.get("distribution", "uniform"), label)
        except ValidationError as error:
            issues.extend(error.issues)
            continue
        if position is not None:
            breakpoints.append({"position": position, "index": node, "distribution": dist})
    if issues:
        raise ValidationError(issues)

    previous_position = 0.0
    previous_index = 0
    for index, item in enumerate(breakpoints):
        label = f"{where}: breakpoints[{index}]"
        if not previous_position < item["position"] < 1.0:
            issues.append(Issue(label, f"position {item['position']} must increase within (0, 1)"))
        if not previous_index < item["index"] < dimension - 1:
            issues.append(
                Issue(label, f"index {item['index']} must increase within (0, {dimension - 1})")
            )
        previous_position = item["position"]
        previous_index = item["index"]
    if not issues:
        indices = [0] + [b["index"] for b in breakpoints] + [dimension - 1]
        dists = [distribution] + [b["distribution"] for b in breakpoints]
        for k, dist in enumerate(dists):
            explicit = dist.get("positions")
            needed = indices[k + 1] - indices[k] + 1
            if dist["type"] == EXPLICIT and explicit is not None and len(explicit) != needed:
                label = "distribution"
                if k > 0:
                    label = f"breakpoints[{k - 1}].distribution"
                issues.append(
                    Issue(
                        f"{where}: {label}.positions",
                        f"has {len(explicit)} positions for {needed} points",
                    )
                )
    params = {
        "segments": tuple(segments),
        "dimension": dimension,
        "mode": mode,
        "distribution": distribution,
        "breakpoints": tuple(breakpoints),
        "propagate_spacing": propagate,
    }
    if not issues:
        try:
            build_chain(params, store)
        except ValueError as error:
            issues.append(Issue(where, str(error)))
    if issues:
        raise ValidationError(issues)
    return params


def build_chain(params: dict, store) -> "_core.CurveChain":
    pieces = []
    scale = 1.0
    for segment in params["segments"]:
        curve = store.geometry(segment["curve"])
        scale = max(scale, curve.length())
        pieces.append((curve, segment["t0"], segment["t1"]))
    return _core.CurveChain(pieces, CHAIN_TOLERANCE * scale)


def node_specs(subs: list[dict]) -> list[tuple]:
    """The C++ node distribution's input for the sub-connectors (sub_connectors)."""
    specs = []
    for sub in subs:
        dist = sub["distribution"]
        specs.append(
            (
                sub["start"],
                sub["end"],
                sub["n_points"],
                dist["type"],
                dist["begin_spacing"],
                dist["end_spacing"],
                dist["growth_ratio"],
                dist.get("positions"),
            )
        )
    return specs


def preview_points(params: dict, store) -> np.ndarray:
    """Where a connector with these canonical parameters puts its points (the C++ node
    distribution, without spacing propagation): previews of parameters not applied yet."""
    chain = build_chain(params, store)
    points, _arc, _warnings, _ratio = _core.connector_nodes(
        chain, params["mode"], node_specs(sub_connectors(params))
    )
    return np.asarray(points)


def script_params(params: dict, store) -> dict:
    """Connector parameters with curve ids replaced by current names (for the journal)."""
    result = dict(params)
    segments = []
    for segment in params["segments"]:
        named = dict(segment)
        named["curve"] = store.get(segment["curve"]).name
        segments.append(named)
    result["segments"] = tuple(segments)
    return result


def order_chain(store, names: list[str]) -> list[dict]:
    """Order and orient curves end to end, starting from a free end if there is one.

    Raises ValidationError if the curves do not form a single open or closed chain.
    """
    if not names:
        raise ValidationError([Issue("chain", "no curves given")])
    curves = {}
    for name in names:
        entity = store.get(name)
        if entity.kind not in CURVE_KINDS:
            raise ValidationError([Issue(name, f"is a {entity.kind}, not a curve")])
        geometry = store.geometry(name)
        curves[name] = (
            np.asarray(geometry.eval(0.0)),
            np.asarray(geometry.eval(1.0)),
            geometry.length(),
        )
    scale = max(1.0, max(length for _a, _b, length in curves.values()))
    tolerance = CHAIN_TOLERANCE * scale

    def touches(point, name):
        start, end, _length = curves[name]
        return (
            np.linalg.norm(point - start) <= tolerance or np.linalg.norm(point - end) <= tolerance
        )

    # Start from a curve with a free end (open chain); otherwise from the first name.
    first = names[0]
    first_reversed = False
    for name in names:
        start, end, _length = curves[name]
        start_free = not any(touches(start, other) for other in names if other != name)
        end_free = not any(touches(end, other) for other in names if other != name)
        if start_free or end_free:
            first = name
            first_reversed = end_free and not start_free
            break
    if first_reversed:
        ordered = [{"curve": first, "t0": 1.0, "t1": 0.0}]
        tail = curves[first][0]
    else:
        ordered = [{"curve": first, "t0": 0.0, "t1": 1.0}]
        tail = curves[first][1]
    remaining = [name for name in names if name != first]
    while remaining:
        for name in remaining:
            start, end, _length = curves[name]
            if np.linalg.norm(tail - start) <= tolerance:
                ordered.append({"curve": name, "t0": 0.0, "t1": 1.0})
                tail = end
                break
            if np.linalg.norm(tail - end) <= tolerance:
                ordered.append({"curve": name, "t0": 1.0, "t1": 0.0})
                tail = start
                break
        else:
            raise ValidationError(
                [
                    Issue(
                        "chain",
                        f"curves {', '.join(remaining)} do not continue the chain "
                        f"after {ordered[-1]['curve']}",
                    )
                ]
            )
        remaining.remove(ordered[-1]["curve"])
    return ordered


@dataclass
class ConnectorNodes:
    """Evaluated connector: nodes and diagnostics."""

    points: np.ndarray  # (n, 3)
    arc_positions: np.ndarray  # (n,)
    warnings: list[str] = field(default_factory=list)
    max_stretching_ratio: float = 1.0
    propagated: dict[str, str] = field(default_factory=dict)  # "begin"/"end" -> source text

    @property
    def length(self) -> float:
        return float(self.arc_positions[-1])

    @property
    def begin_spacing(self) -> float:
        return float(self.arc_positions[1] - self.arc_positions[0])

    @property
    def end_spacing(self) -> float:
        return float(self.arc_positions[-1] - self.arc_positions[-2])

    def interval_ratios(self) -> np.ndarray:
        """Per interval: the largest ratio to its neighbouring intervals (>= 1)."""
        spacing = np.diff(self.arc_positions)
        ratios = np.ones_like(spacing)
        if len(spacing) > 1:
            pair = np.maximum(spacing[1:] / spacing[:-1], spacing[:-1] / spacing[1:])
            ratios[1:] = np.maximum(ratios[1:], pair)
            ratios[:-1] = np.maximum(ratios[:-1], pair)
        return ratios


def sub_connectors(params: dict) -> list[dict]:
    """Sub-connector specs (start, end, n_points, distribution) from canonical parameters."""
    positions = [0.0] + [b["position"] for b in params["breakpoints"]] + [1.0]
    indices = [0] + [b["index"] for b in params["breakpoints"]] + [params["dimension"] - 1]
    distributions = [params["distribution"]] + [b["distribution"] for b in params["breakpoints"]]
    subs = []
    for k in range(len(distributions)):
        subs.append(
            {
                "start": positions[k],
                "end": positions[k + 1],
                "n_points": indices[k + 1] - indices[k] + 1,
                "distribution": dict(distributions[k]),
            }
        )
    return subs


class ConnectorEvaluator:
    """Evaluates connectors of a store (with spacing propagation); results are cached until
    ``invalidate`` is called (the Session does that after every model change)."""

    def __init__(self, store):
        self.store = store
        self._chains: dict[int, object] = {}
        self._nodes: dict[int, ConnectorNodes] = {}
        self._canonical: np.ndarray | None = None  # (m, 3) distinct connector ends
        self._ends: tuple | None = None  # every connector end, see _end_index

    def invalidate(self) -> None:
        self._chains.clear()
        self._nodes.clear()
        self._canonical = None
        self._ends = None

    def _end_index(self) -> tuple:
        """Every connector end, built once per invalidation: points (2 n, 3) and, per end, the
        connector id, "begin" / "end", the explicit spacing there (nan: none) and the
        connector name. Lookups are one vectorised distance each."""
        if self._ends is None:
            points = []
            ids = []
            sides = []
            spacings = []
            names = []
            for other in self.store.entities("connector"):
                chain = self.chain(other.id)
                subs = sub_connectors(other.params)
                for side, point, spacing in (
                    ("begin", chain.start(), subs[0]["distribution"]["begin_spacing"]),
                    ("end", chain.end(), subs[-1]["distribution"]["end_spacing"]),
                ):
                    points.append(point)
                    ids.append(other.id)
                    sides.append(side)
                    if spacing is None:
                        spacings.append(np.nan)
                    else:
                        spacings.append(float(spacing))
                    names.append(other.name)
            self._ends = (
                np.asarray(points, dtype=float).reshape(-1, 3),
                np.asarray(ids, dtype=np.int64),
                sides,
                np.asarray(spacings, dtype=float),
                names,
            )
        return self._ends

    def canonical_point(self, point) -> np.ndarray:
        """Shared coordinate of a connector end: coincident ends of different connectors are
        represented by one bitwise-identical point (the first in id order), so that domains
        sharing a node get identical coordinates."""
        if self._canonical is None:
            kept = np.zeros((0, 3))
            for end in self._end_index()[0]:
                tolerance = CHAIN_TOLERANCE * max(1.0, float(np.max(np.abs(end))))
                if len(kept) == 0 or np.min(np.linalg.norm(kept - end, axis=1)) > tolerance:
                    kept = np.vstack([kept, end[None, :]])
            self._canonical = kept
        point = np.asarray(point, dtype=float)
        if len(self._canonical) == 0:
            return point
        tolerance = CHAIN_TOLERANCE * max(1.0, float(np.max(np.abs(point))))
        distances = np.linalg.norm(self._canonical - point, axis=1)
        nearest = int(np.argmin(distances))
        if distances[nearest] <= tolerance:
            return self._canonical[nearest]
        return point

    def chain(self, entity_id: int):
        if entity_id not in self._chains:
            self._chains[entity_id] = build_chain(self.store.get(entity_id).params, self.store)
        return self._chains[entity_id]

    def nodes(self, entity_id: int) -> ConnectorNodes:
        if entity_id not in self._nodes:
            self._nodes[entity_id] = self._evaluate(entity_id)
        return self._nodes[entity_id]

    def _explicit_spacings_at(
        self, point: np.ndarray, exclude: int, tolerance: float
    ) -> list[tuple[float, str]]:
        """Explicit spacings of the other connectors' ends at ``point``: (spacing, "name
        begin|end"), sorted."""
        points, ids, sides, spacings, names = self._end_index()
        if len(points) == 0:
            return []
        near = (np.linalg.norm(points - point, axis=1) <= tolerance) & (ids != exclude)
        near &= np.isfinite(spacings)
        return sorted((float(spacings[k]), f"{names[k]} {sides[k]}") for k in np.flatnonzero(near))

    def _evaluate(self, entity_id: int) -> ConnectorNodes:
        entity = self.store.get(entity_id)
        params = entity.params
        chain = self.chain(entity_id)
        subs = sub_connectors(params)
        propagated = {}
        if params["propagate_spacing"]:
            tolerance = CHAIN_TOLERANCE * max(1.0, chain.length())
            ends = (
                ("begin", subs[0], "begin_spacing", chain.start()),
                ("end", subs[-1], "end_spacing", chain.end()),
            )
            for label, sub, key, point in ends:
                dist = sub["distribution"]
                if dist["type"] in SPACING_TYPES and dist[key] is None:
                    candidates = self._explicit_spacings_at(np.asarray(point), entity_id, tolerance)
                    if candidates:
                        dist[key] = candidates[0][0]
                        propagated[label] = f"{candidates[0][0]:.6g} from {candidates[0][1]}"
        specs = node_specs(subs)
        points, arc_positions, warnings, ratio = _core.connector_nodes(chain, params["mode"], specs)
        points[0] = self.canonical_point(points[0])
        points[-1] = self.canonical_point(points[-1])
        points.setflags(write=False)
        return ConnectorNodes(points, arc_positions, list(warnings), ratio, propagated)
