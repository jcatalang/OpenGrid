"""Commands creating, modifying, renaming and deleting database entities."""

from opengrid.commands.base import Command
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


def create_call(kind: str, params: dict, name: str) -> str:
    return f"og.create_{kind}({format_kwargs(params)}, name={name!r})"


class CreateEntity(Command):
    def __init__(self, kind: str, params: dict, name: str | None = None):
        self.kind = kind
        self.raw_params = dict(params)
        self.requested_name = name
        self.label = f"Create {kind}"
        self.name: str | None = None
        self.params: dict | None = None

    def execute(self, store: Store) -> str:
        where = self.requested_name
        if where is None:
            where = f"new {self.kind}"
        self.params = normalize_params(self.kind, self.raw_params, where)
        entity = store.add(self.kind, self.params, self.requested_name)
        self.name = entity.name
        return entity.name

    def script(self) -> str:
        return create_call(self.kind, self.params, self.name)


class CreateEntities(Command):
    """Create several entities as one undo step (e.g. a file import).

    The journal records explicit ``create_*`` calls, so replay does not depend on the source
    file still existing or being unchanged.
    """

    def __init__(self, items: list[tuple[str, dict, str | None]], label: str, comment: str = ""):
        self.items = items
        self.label = label
        self.comment = comment
        self.created: list[tuple[str, dict, str]] = []

    def execute(self, store: Store) -> list[str]:
        self.created = []
        for kind, raw_params, name in self.items:
            where = name
            if where is None:
                where = f"new {kind}"
            params = normalize_params(kind, raw_params, where)
            entity = store.add(kind, params, name)
            self.created.append((kind, params, entity.name))
        names = []
        for _kind, _params, name in self.created:
            names.append(name)
        return names

    def script(self) -> str:
        lines = []
        if self.comment:
            lines.append(f"# {self.comment}")
        for kind, params, name in self.created:
            lines.append(create_call(kind, params, name))
        return "\n".join(lines)


class ModifyEntity(Command):
    def __init__(self, name: str, changes: dict):
        self.target = name
        self.changes = dict(changes)
        self.label = f"Modify {name}"
        self.recorded: dict = {}

    def execute(self, store: Store) -> str:
        if not self.changes:
            raise ValidationError([Issue(self.target, "no parameters to modify")])
        entity = store.get(self.target)
        merged = dict(entity.params)
        merged.update(self.changes)
        params = normalize_params(entity.kind, merged, entity.name)
        store.update(entity.id, params=params)
        self.recorded = {}
        for key in self.changes:
            self.recorded[key] = params[key]
        return entity.name

    def script(self) -> str:
        return f"og.modify({self.target!r}, {format_kwargs(self.recorded)})"


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
        for name in self.names:
            if name not in store:
                issues.append(Issue(repr(name), "no such entity"))
        if issues:
            raise ValidationError(issues)
        for name in self.names:
            store.remove(name)
        return len(self.names)

    def script(self) -> str:
        return f"og.delete({self.names!r})"
