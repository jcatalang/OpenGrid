"""Curves as the window shows them: one curve, its shape and its points.

The model keeps two entities for a meshed curve: the database curve (a line, arc, polyline,
spline or NURBS: the shape, what the handles move) and the connector on it (the points: their
number and distribution, what domains and patches are built on). The window shows them as one
"curve" when they are one-to-one: the connector has a single segment covering the whole curve
and no other connector uses that curve. That pair is shown, selected, edited and deleted as
one item, named after the connector.

A connector on a part of a curve or along several curves is a curve of its own (its shape is
those pieces); a database curve that no connector uses is a curve without points yet; a
database curve used by such connectors is not listed (its connectors are).
"""

from opengrid.model.connectors import CURVE_KINDS
from opengrid.model.store import name_prefix


class CurveIndex:
    """The curve pairing of a store, at the time it is built (cheap: one pass over the
    connectors)."""

    def __init__(self, store):
        self.store = store
        self._generated: set[int] | None = None
        users: dict[int, list[int]] = {}
        for connector in store.entities("connector"):
            for segment in connector.params["segments"]:
                curve_users = users.setdefault(segment["curve"], [])
                if connector.id not in curve_users:
                    curve_users.append(connector.id)
        self.users = users  # database curve id -> ids of the connectors on it
        self.geometry_of: dict[int, int] = {}  # paired connector id -> its curve id
        self.connector_of: dict[int, int] = {}  # paired curve id -> its connector id
        for connector in store.entities("connector"):
            segments = connector.params["segments"]
            if len(segments) != 1 or not _whole(segments[0]):
                continue
            curve = segments[0]["curve"]
            if users[curve] == [connector.id]:
                self.geometry_of[connector.id] = curve
                self.connector_of[curve] = connector.id

    def representative(self, entity_id: int) -> int:
        """The id that stands for the entity's curve (a paired curve: its connector)."""
        return self.connector_of.get(entity_id, entity_id)

    def members(self, entity_id: int) -> list[int]:
        """The entities of the curve an id stands for: [curve, connector] for a pair."""
        curve = self.geometry_of.get(entity_id)
        if curve is None:
            connector = self.connector_of.get(entity_id)
            if connector is None:
                return [entity_id]
            return [entity_id, connector]
        return [curve, entity_id]

    def generated(self) -> set[int]:
        """Ids of the sketches' outputs (their names are made up from the sketch's)."""
        if self._generated is None:
            from opengrid.model.sketch import output_owner

            self._generated = set(output_owner(self.store))
        return self._generated

    def is_listed(self, entity) -> bool:
        """Whether the entity has a row of its own among the curves (a database curve used
        by connectors has none: its connectors stand for it)."""
        if entity.kind in CURVE_KINDS:
            return entity.id not in self.users
        return True

    def describe(self, entity) -> str:
        """A short type text: 'line, 21 pts' ('line rectangle-1-bottom, 21 pts' when the
        shape has a name of its own), 'arc, no points', '21 pts, part of a line', '41 pts,
        along 3 curves'."""
        store = self.store
        if entity.kind in CURVE_KINDS:
            if entity.id in self.users:
                return f"{entity.kind}, in {len(self.users[entity.id])} curve(s)"
            return f"{entity.kind}, no points"
        points = f"{entity.params['dimension']} pts"
        curve = self.geometry_of.get(entity.id)
        if curve is not None:
            shape = store.get(curve)
            if is_automatic_name(shape) or shape.id in self.generated():
                return f"{shape.kind}, {points}"
            return f"{shape.kind} {shape.name}, {points}"  # e.g. a shape's side
        segments = entity.params["segments"]
        if len(segments) == 1:
            curve = segments[0]["curve"]
            kind = store.get(curve).kind
            if _whole(segments[0]):  # the whole shape, which other points share
                return f"{kind}, {points} (1 of {len(self.users[curve])} on it)"
            return f"{points}, part of {_article(kind)} {kind}"
        return f"{points}, along {len(segments)} curves"


def is_automatic_name(entity) -> bool:
    """Whether the entity has an automatic name ("line-3"), not one given to it."""
    prefix = f"{name_prefix(entity.kind)}-"
    suffix = entity.name[len(prefix) :]
    return entity.name.startswith(prefix) and suffix.isdigit()


def _whole(segment) -> bool:
    ends = sorted((float(segment["t0"]), float(segment["t1"])))
    return ends == [0.0, 1.0]


def _article(word: str) -> str:
    if word[:1] in "aeiou":
        return "an"
    return "a"
