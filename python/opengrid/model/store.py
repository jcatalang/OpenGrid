"""Entity store: the authoritative model of a session, with transactional change tracking.

Every mutation must happen inside a transaction (``begin`` ... ``commit``/``rollback``); the
commit returns a :class:`Change` holding the before/after state of every touched entity, which
is all the command stack needs for exact undo and redo.
"""

import re
from collections.abc import Iterator
from dataclasses import dataclass, field

from opengrid.model.entities import SURFACE_CURVES, SURFACE_KINDS, build_geometry, build_surface
from opengrid.validation import Issue, ValidationError

NAME_PATTERN = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.\-]*$")


@dataclass(frozen=True)
class Entity:
    """Immutable snapshot of one entity. ``params`` is canonical and must not be mutated."""

    id: int
    name: str
    kind: str
    params: dict


def references(entity: Entity) -> set[int]:
    """Ids of the entities this entity refers to."""
    if entity.kind == "connector":
        return {segment["curve"] for segment in entity.params["segments"]}
    if entity.kind in SURFACE_CURVES:
        return {entity.params[key] for key in SURFACE_CURVES[entity.kind]}
    if entity.kind == "domain":
        refs = {item["connector"] for edge in entity.params["edges"] for item in edge}
        if entity.params.get("surface") is not None:
            refs.add(entity.params["surface"])
        return refs
    if entity.kind == "block":
        from opengrid.model.blocks import block_references

        return block_references(entity.params)
    if entity.kind == "unstructured_domain":
        return {item["connector"] for loop in entity.params["loops"] for item in loop}
    if entity.kind == "patch":
        refs = {
            item.get("connector", item.get("domain", item.get("block")))
            for item in entity.params["members"]
        }
        if entity.params.get("periodic") is not None:
            refs.add(entity.params["periodic"]["partner"])
        return refs
    if entity.kind == "zone":
        return {item["domain"] for item in entity.params["members"]}
    return set()


@dataclass
class Change:
    """Before/after state of the entities touched by one transaction (None = absent)."""

    before: dict[int, Entity | None] = field(default_factory=dict)
    after: dict[int, Entity | None] = field(default_factory=dict)

    def is_empty(self) -> bool:
        return not self.before

    def ids(self) -> list[int]:
        return sorted(self.before)


class Store:
    def __init__(self):
        self._entities: dict[int, Entity] = {}
        self._names: dict[str, int] = {}
        self._geometry: dict[int, object] = {}
        self._transaction: dict[int, Entity | None] | None = None
        self.next_id = 1
        self.name_counters: dict[str, int] = {}

    # ------------------------------------------------------------------ queries
    def __len__(self) -> int:
        return len(self._entities)

    def __iter__(self) -> Iterator[Entity]:
        for entity_id in sorted(self._entities):
            yield self._entities[entity_id]

    def __contains__(self, key: str | int) -> bool:
        if isinstance(key, str):
            return key in self._names
        return key in self._entities

    def get(self, key: str | int) -> Entity:
        """Entity by name or id; KeyError with a readable message if absent."""
        if isinstance(key, str):
            if key not in self._names:
                raise KeyError(f"no entity named {key!r}")
            return self._entities[self._names[key]]
        if key not in self._entities:
            raise KeyError(f"no entity with id {key}")
        return self._entities[key]

    def entities(self, kind: str | None = None) -> list[Entity]:
        result = []
        for entity in self:
            if kind is None or entity.kind == kind:
                result.append(entity)
        return result

    def geometry(self, key: str | int):
        """C++ curve or surface of an entity (cached), or None for points."""
        entity = self.get(key)
        if entity.id not in self._geometry:
            if entity.kind in SURFACE_KINDS:
                self._geometry[entity.id] = build_surface(entity.kind, entity.params, self.geometry)
            else:
                self._geometry[entity.id] = build_geometry(entity.kind, entity.params)
        return self._geometry[entity.id]

    def dependents(self, key: str | int) -> list[Entity]:
        """Entities that reference the given one (connectors -> curves, domains -> connectors)."""
        entity_id = self.get(key).id
        return [entity for entity in self if entity_id in references(entity)]

    def snapshot(self) -> dict[int, Entity]:
        return dict(self._entities)

    # ------------------------------------------------------------- transactions
    @property
    def in_transaction(self) -> bool:
        return self._transaction is not None

    def begin(self) -> None:
        if self._transaction is not None:
            raise RuntimeError("transaction already open")
        self._transaction = {}

    def commit(self) -> Change:
        transaction = self._require_transaction()
        change = Change()
        for entity_id in sorted(transaction):
            before = transaction[entity_id]
            after = self._entities.get(entity_id)
            if before != after:
                change.before[entity_id] = before
                change.after[entity_id] = after
        self._transaction = None
        return change

    def rollback(self) -> None:
        transaction = self._require_transaction()
        self._transaction = None
        self.apply_states(transaction)

    def apply_states(self, states: dict[int, Entity | None]) -> None:
        """Force the given entities into the given states (used by undo/redo/rollback)."""
        for entity_id in states:
            self._detach(entity_id)
        for entity_id in sorted(states):
            entity = states[entity_id]
            if entity is not None:
                self._attach(entity)

    # ---------------------------------------------------------------- mutations
    def add(self, kind: str, params: dict, name: str | None = None) -> Entity:
        """Add an entity with canonical params; an automatic name is chosen if none given."""
        self._require_transaction()
        if name is None:
            name = self.unique_name(kind)
        else:
            self._check_name(name)
            self._bump_counter(kind, name)
        entity = Entity(self.next_id, name, kind, params)
        self.next_id += 1
        self._touch(entity.id)
        self._attach(entity)
        return entity

    def update(self, key: str | int, params: dict | None = None, name: str | None = None) -> Entity:
        self._require_transaction()
        old = self.get(key)
        new_name = old.name
        if name is not None and name != old.name:
            self._check_name(name)
            new_name = name
        new_params = old.params
        if params is not None:
            new_params = params
        self._touch(old.id)
        self._detach(old.id)
        entity = Entity(old.id, new_name, old.kind, new_params)
        self._attach(entity)
        return entity

    def remove(self, key: str | int) -> None:
        self._require_transaction()
        entity = self.get(key)
        self._touch(entity.id)
        self._detach(entity.id)

    def unique_name(self, kind: str) -> str:
        """Next free automatic name ``<kind>-<n>``; counters only ever increase."""
        counter = self.name_counters.get(kind, 0)
        while True:
            counter += 1
            name = f"{kind}-{counter}"
            if name not in self._names:
                break
        self.name_counters[kind] = counter
        return name

    # ------------------------------------------------------------------ helpers
    def _bump_counter(self, kind: str, name: str) -> None:
        """An explicit name of the automatic form '<kind>-<n>' advances that counter, so a
        journal replay (which always passes explicit names) ends with identical counters."""
        prefix = f"{kind}-"
        suffix = name[len(prefix) :]
        if name.startswith(prefix) and suffix.isdigit() and suffix.isascii():
            number = int(suffix)
            if number > self.name_counters.get(kind, 0):
                self.name_counters[kind] = number

    def _check_name(self, name: str) -> None:
        if not isinstance(name, str) or NAME_PATTERN.match(name) is None:
            raise ValidationError(
                [
                    Issue(
                        repr(name),
                        "invalid name (letters, digits, '_', '.', '-'; "
                        "must not start with '.' or '-')",
                    )
                ]
            )
        if name in self._names:
            raise ValidationError([Issue(name, "name already in use")])

    def _require_transaction(self) -> dict[int, Entity | None]:
        if self._transaction is None:
            raise RuntimeError("model mutations must run inside a transaction (use a Command)")
        return self._transaction

    def _touch(self, entity_id: int) -> None:
        if self._transaction is not None and entity_id not in self._transaction:
            self._transaction[entity_id] = self._entities.get(entity_id)

    def _attach(self, entity: Entity) -> None:
        self._entities[entity.id] = entity
        self._names[entity.name] = entity.id

    def _detach(self, entity_id: int) -> None:
        entity = self._entities.pop(entity_id, None)
        if entity is not None:
            del self._names[entity.name]
        self._geometry.pop(entity_id, None)
        # Surfaces hold the C++ curves they were built from: rebuild them too.
        for other in self._entities.values():
            if other.kind in SURFACE_CURVES and entity_id in references(other):
                self._geometry.pop(other.id, None)
