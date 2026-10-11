"""Commands creating, modifying, renaming and deleting entities (geometry and connectors).

Each command renders its journal line during ``execute`` (entity references are written by
their names at that moment, so a journal stays valid across later renames).
"""

from opengrid.commands.base import Command
from opengrid.model.connectors import script_params
from opengrid.model.entities import normalize_params
from opengrid.model.store import Store
from opengrid.validation import Issue, ValidationError


def format_value(value) -> str:
    """Python literal for a canonical parameter value (repr of floats round-trips exactly)."""
    return repr(value)


def format_kwargs(params: dict) -> str:
    parts = []
    for key, value in params.items():
        parts.append(f"{key}={format_value(value)}")
    return ", ".join(parts)


def _for_script(kind: str, params: dict, store: Store) -> dict:
    if kind == "connector":
        return script_params(params, store)
    if kind in ("domain", "unstructured_domain"):
        key = "edges"
        if kind == "unstructured_domain":
            key = "loops"
        result = dict(params)
        edges = []
        for edge in params[key]:
            items = []
            for item in edge:
                items.append(
                    {"connector": store.get(item["connector"]).name, "reversed": item["reversed"]}
                )
            edges.append(tuple(items))
        result[key] = tuple(edges)
        return result
    if kind in ("patch", "zone"):
        key = "connector"
        if kind == "zone":
            key = "domain"
        result = dict(params)
        result["members"] = tuple(store.get(item[key]).name for item in params["members"])
        if kind == "patch" and params.get("periodic") is not None:
            periodic = dict(params["periodic"])
            periodic["partner"] = store.get(periodic["partner"]).name
            result["periodic"] = periodic
        return result
    return params


def create_call(kind: str, params: dict, name: str, store: Store) -> str:
    return f"og.create_{kind}({format_kwargs(_for_script(kind, params, store))}, name={name!r})"


def _check_dependents_still_valid(store: Store, entity_id: int) -> None:
    """After an entity changed, everything that depends on it (transitively) must still be
    valid: connectors must form continuous chains, domains closed loops with matching counts."""
    issues = []
    pending = [entity_id]
    seen = set()
    while pending:
        current = pending.pop(0)
        for dependent in store.dependents(current):
            if dependent.id in seen:
                continue
            seen.add(dependent.id)
            try:
                normalize_params(dependent.kind, dict(dependent.params), dependent.name, store)
            except ValidationError as error:
                details = "; ".join(issue.message for issue in error.issues)
                issues.append(Issue(dependent.name, f"would become invalid: {details}"))
                continue
            pending.append(dependent.id)
    if issues:
        raise ValidationError(issues)


class CreateEntity(Command):
    def __init__(self, kind: str, params: dict, name: str | None = None):
        self.kind = kind
        self.raw_params = dict(params)
        self.requested_name = name
        self.label = f"Create {kind}"
        self.name: str | None = None
        self._script = ""

    def execute(self, store: Store) -> str:
        where = self.requested_name
        if where is None:
            where = f"new {self.kind}"
        params = normalize_params(self.kind, self.raw_params, where, store)
        entity = store.add(self.kind, params, self.requested_name)
        self.name = entity.name
        self._script = create_call(self.kind, params, entity.name, store)
        return entity.name

    def script(self) -> str:
        return self._script


class CreateEntities(Command):
    """Create several entities as one undo step (e.g. a file import).

    The journal records explicit ``create_*`` calls, so replay does not depend on the source
    file still existing or being unchanged.
    """

    def __init__(self, items: list[tuple[str, dict, str | None]], label: str, comment: str = ""):
        self.items = items
        self.label = label
        self.comment = comment
        self._script = ""

    def execute(self, store: Store) -> list[str]:
        lines = []
        if self.comment:
            lines.append(f"# {self.comment}")
        names = []
        for kind, raw_params, name in self.items:
            where = name
            if where is None:
                where = f"new {kind}"
            params = normalize_params(kind, raw_params, where, store)
            entity = store.add(kind, params, name)
            names.append(entity.name)
            lines.append(create_call(kind, params, entity.name, store))
        self._script = "\n".join(lines)
        return names

    def script(self) -> str:
        return self._script


class ModifyEntity(Command):
    def __init__(self, name: str, changes: dict):
        self.target = name
        self.changes = dict(changes)
        self.label = f"Modify {name}"
        self._script = ""

    def execute(self, store: Store) -> str:
        if not self.changes:
            raise ValidationError([Issue(self.target, "no parameters to modify")])
        entity = store.get(self.target)
        merged = dict(entity.params)
        merged.update(self.changes)
        params = normalize_params(entity.kind, merged, entity.name, store)
        store.update(entity.id, params=params)
        _check_dependents_still_valid(store, entity.id)
        recorded = {}
        shown = _for_script(entity.kind, params, store)
        for key in self.changes:
            recorded[key] = shown[key]
        self._script = f"og.modify({self.target!r}, {format_kwargs(recorded)})"
        return entity.name

    def script(self) -> str:
        return self._script


class RenameEntity(Command):
    def __init__(self, name: str, new_name: str):
        self.target = name
        self.new_name = new_name
        self.label = f"Rename {name}"

    def execute(self, store: Store) -> str:
        store.update(self.target, name=self.new_name)
        return self.new_name

    def script(self) -> str:
        return f"og.rename({self.target!r}, {self.new_name!r})"


class DeleteEntities(Command):
    def __init__(self, names: list[str]):
        self.names = list(names)
        self.label = f"Delete {len(self.names)} entities"
        if len(self.names) == 1:
            self.label = f"Delete {self.names[0]}"

    def execute(self, store: Store) -> int:
        issues = []
        seen = set()
        for name in self.names:
            if name not in store:
                issues.append(Issue(repr(name), "no such entity"))
            elif name in seen:
                issues.append(Issue(repr(name), "listed twice"))
            seen.add(name)
        if issues:
            raise ValidationError(issues)
        for name in self.names:
            for dependent in store.dependents(name):
                if dependent.name not in seen:
                    issues.append(Issue(name, f"is used by {dependent.kind} {dependent.name!r}"))
        if issues:
            raise ValidationError(issues)
        for name in self.names:
            store.remove(name)
        return len(self.names)

    def script(self) -> str:
        return f"og.delete({self.names!r})"
