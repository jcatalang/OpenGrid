"""Drag handles: the points of entities that can be moved with the mouse (GUI logic without
rendering, unit-testable).

Handles: point entities; the end points of lines; the defining points of polylines and
splines; the control points of NURBS curves; connector break points (moved along the
connector). A drag previews the change inside a store transaction that is rolled back after
drawing, so nothing is recorded until the release, which runs one Modify command (one undo
step, one journal line with the final position: Session.move_handle, which also keeps the
connectors of a curve network on the moved curve's crossings).
"""

from dataclasses import dataclass

import numpy as np

from opengrid.model.connectors import CURVE_KINDS


@dataclass(frozen=True)
class Handle:
    entity: str
    field: str  # parameter name
    index: int  # item in a list parameter (-1 for a single point)
    xyz: tuple[float, float, float]


def entity_handles(session, name: str) -> list[Handle]:
    entity = session.store.get(name)
    p = entity.params
    if entity.kind == "point":
        return [Handle(name, "xyz", -1, tuple(p["xyz"]))]
    if entity.kind == "line":
        return [Handle(name, "p0", -1, tuple(p["p0"])), Handle(name, "p1", -1, tuple(p["p1"]))]
    if entity.kind in ("polyline", "spline"):
        return [Handle(name, "points", k, tuple(q)) for k, q in enumerate(p["points"])]
    if entity.kind == "nurbs":
        return [
            Handle(name, "control_points", k, tuple(q)) for k, q in enumerate(p["control_points"])
        ]
    if entity.kind == "connector":
        chain = session.connectors.chain(entity.id)
        out = []
        for k, item in enumerate(p["breakpoints"]):
            point = tuple(chain.point_at_length(item["position"] * chain.length()))
            out.append(Handle(name, "breakpoints", k, point))
        return out
    return []


def handles_for(session, names) -> list[Handle]:
    out = []
    for name in names:
        entity = session.store.get(name)
        if entity.kind in CURVE_KINDS or entity.kind in ("point", "connector"):
            out.extend(entity_handles(session, name))
    return out


def pick_handle(handles: list[Handle], xyz, tolerance: float) -> Handle | None:
    best = None
    best_distance = tolerance
    for handle in handles:
        d = float(np.linalg.norm(np.asarray(handle.xyz) - np.asarray(xyz, dtype=float)))
        if d <= best_distance:
            best = handle
            best_distance = d
    return best


def moved_params(session, handle: Handle, xyz) -> dict:
    """Parameter changes that put ``handle`` at ``xyz`` (break points: the closest point of
    the connector, kept between their neighbours)."""
    entity = session.store.get(handle.entity)
    point = tuple(float(c) for c in xyz)
    if handle.field == "breakpoints":
        nodes = session.connector_nodes(handle.entity)
        positions = nodes.arc_positions
        k = int(np.argmin(np.linalg.norm(nodes.points - np.asarray(point), axis=1)))
        fraction = float(positions[k] / positions[-1])
        items = [dict(b) for b in entity.params["breakpoints"]]
        lo = 0.0
        hi = 1.0
        if handle.index > 0:
            lo = items[handle.index - 1]["position"]
        if handle.index + 1 < len(items):
            hi = items[handle.index + 1]["position"]
        if not lo < fraction < hi:
            return {}
        items[handle.index]["position"] = fraction
        items[handle.index]["index"] = k
        return {"breakpoints": tuple(items)}
    if handle.index < 0:
        return {handle.field: point}
    values = list(entity.params[handle.field])
    values[handle.index] = point
    return {handle.field: tuple(values)}


def dependents_closure(store, name: str) -> list[int]:
    """Ids of the entity and everything that depends on it, transitively."""
    seen = []
    pending = [store.get(name).id]
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.append(current)
        pending.extend(d.id for d in store.dependents(current))
    return sorted(seen)


class Drag:
    """One drag of a handle: preview(xyz) while moving, finish(xyz) on release."""

    def __init__(self, session, handle: Handle):
        self.session = session
        self.handle = handle
        self.affected = dependents_closure(session.store, handle.entity)
        self.changes: dict = {}

    def preview(self, xyz, draw) -> bool:
        """Apply the move temporarily, call ``draw(affected ids)``, undo it. Returns False when
        the move is invalid (nothing drawn)."""
        from opengrid.model.entities import normalize_params
        from opengrid.model.topology_edit import move_handle
        from opengrid.validation import ValidationError

        store = self.session.store
        changes = moved_params(self.session, self.handle, xyz)
        if not changes:
            return False
        entity = store.get(self.handle.entity)
        merged = dict(entity.params)
        merged.update(changes)
        store.begin()
        try:
            try:
                if self.handle.field == "breakpoints":
                    params = normalize_params(entity.kind, merged, entity.name, store)
                    store.update(entity.id, params=params)
                else:
                    move_handle(
                        store, entity.name, self.handle.field, self.handle.index, xyz, check=False
                    )
            except (ValidationError, ValueError):
                return False
            self.session._invalidate_derived(None)
            try:
                draw(self.affected)
            except (ValidationError, ValueError):
                return False
            self.changes = changes
            return True
        finally:
            store.rollback()
            self.session._invalidate_derived(None)

    def finish(self, xyz) -> str | None:
        """Commit the final position as one command (Session.move_handle, or a Modify of the
        break points); returns the entity name."""
        if self.handle.field == "breakpoints":
            changes = moved_params(self.session, self.handle, xyz)
            if not changes:
                changes = self.changes
            if not changes:
                return None
            return self.session.modify(self.handle.entity, **changes)
        return self.session.move_handle(
            self.handle.entity, self.handle.field, self.handle.index, xyz
        )
