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
    from opengrid.model.entities import SURFACE_CURVES

    if kind in SURFACE_CURVES:
        result = dict(params)
        for key in SURFACE_CURVES[kind]:
            result[key] = store.get(params[key]).name
        return result
    if kind == "volume":
        from opengrid.model.volumes import volume_script_params

        return volume_script_params(params, store)
    if kind == "smoothing_group":
        from opengrid.model.groups import group_script_params

        return group_script_params(params, store)
    if kind == "block":
        from opengrid.model.blocks import block_script_params

        return block_script_params(params, store)
    if kind == "sketch":
        from opengrid.model.sketch import sketch_script_params

        return sketch_script_params(params, store)
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
        if kind in ("domain", "unstructured_domain") and params.get("surface") is not None:
            result["surface"] = store.get(params["surface"]).name
        if kind == "unstructured_domain" and params["generator"]["inflation"] is not None:
            generator = dict(params["generator"])
            inflation = dict(generator["inflation"])
            inflation["connectors"] = tuple(store.get(c).name for c in inflation["connectors"])
            generator["inflation"] = inflation
            result["generator"] = generator
        return result
    if kind in ("patch", "zone"):
        key = "connector"
        if kind == "zone":
            key = "domain"
        result = dict(params)
        if kind == "patch":
            from opengrid.model.patches import member_script

            result["members"] = tuple(member_script(item, store) for item in params["members"])
        else:
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

    def __init__(
        self,
        items: list[tuple[str, dict, str | None]],
        label: str,
        comment: str = "",
        optional_kinds: tuple[str, ...] = (),
    ):
        self.items = items
        self.label = label
        self.comment = comment
        # Items of these kinds that fail validation are skipped and listed in ``failed``
        # (name, reason) instead of failing the whole command (e.g. one unmeshable CAD face).
        self.optional_kinds = optional_kinds
        self.failed: list[tuple[str, str]] = []
        self._script = ""

    def execute(self, store: Store) -> list[str]:
        lines = []
        if self.comment:
            lines.append(f"# {self.comment}")
        names = []
        self.failed = []
        for kind, raw_params, name in self.items:
            where = name
            if where is None:
                where = f"new {kind}"
            if kind in self.optional_kinds:
                try:
                    params = normalize_params(kind, raw_params, where, store)
                except ValidationError as error:
                    reason = "; ".join(issue.message for issue in error.issues)
                    self.failed.append((where, reason))
                    continue
            else:
                params = normalize_params(kind, raw_params, where, store)
            entity = store.add(kind, params, name)
            names.append(entity.name)
            lines.append(create_call(kind, params, entity.name, store))
        self._script = "\n".join(lines)
        return names

    def script(self) -> str:
        return self._script


class ModifyEntity(Command):
    """Change parameters and / or rename, as one undo step (the journal gets a modify line
    and a rename line)."""

    def __init__(self, name: str, changes: dict, new_name: str | None = None):
        self.target = name
        self.changes = dict(changes)
        self.new_name = new_name
        self.label = f"Modify {name}"
        if new_name is not None and not changes:
            self.label = f"Rename {name}"
        self._script = ""

    def execute(self, store: Store) -> str:
        if not self.changes and self.new_name is None:
            raise ValidationError([Issue(self.target, "no parameters to modify")])
        entity = store.get(self.target)
        lines = []
        if self.changes:
            merged = dict(entity.params)
            merged.update(self.changes)
            params = normalize_params(entity.kind, merged, entity.name, store)
            store.update(entity.id, params=params)
            from opengrid.model.sketch import regenerate_adopters

            regenerate_adopters(store, entity.id)  # sketches follow an adopted connector
            _check_dependents_still_valid(store, entity.id)
            recorded = {}
            shown = _for_script(entity.kind, params, store)
            for key in self.changes:
                recorded[key] = shown[key]
            lines.append(f"og.modify({self.target!r}, {format_kwargs(recorded)})")
        name = entity.name
        if self.new_name is not None and self.new_name != entity.name:
            store.update(entity.id, name=self.new_name)
            lines.append(f"og.rename({self.target!r}, {self.new_name!r})")
            name = self.new_name
        self._script = "\n".join(lines)
        return name

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


from opengrid.model.membership import MEMBERSHIP_KINDS, drop_members  # noqa: E402


def delete_plan(store: Store, names) -> tuple[list[str], list[str], list[Issue]]:
    """What deleting ``names`` means: (those entities plus the outputs of the sketches among
    them, the dependents that would go with them - transitively: the connectors on a curve,
    the domains on a connector, the blocks and volumes on a domain - and the refusals: a
    sketch's output named on its own, an entity a sketch is assigned to). Patches, zones and
    smoothing groups are not dependents: they lose the members."""
    from opengrid.model.sketch import output_owner

    issues: list[Issue] = []
    seen: set[str] = set()
    expanded: list[str] = []
    for name in names:
        if name not in store:
            issues.append(Issue(repr(name), "no such entity"))
        elif name in seen:
            issues.append(Issue(repr(name), "listed twice"))
        else:
            expanded.append(name)
            seen.add(name)
    if issues:
        return [], [], issues
    for name in list(expanded):
        entity = store.get(name)
        if entity.kind == "sketch":
            for _role, output in entity.params["outputs"]:
                if output in store and store.get(output).name not in seen:
                    expanded.append(store.get(output).name)
                    seen.add(store.get(output).name)
    owners = output_owner(store)
    for name in expanded:
        owner = owners.get(store.get(name).id)
        if owner is not None and owner[0].name not in seen:
            issues.append(
                Issue(
                    name,
                    f"is generated by sketch {owner[0].name!r}: remove its block from the "
                    "sketch (remove_sketch_blocks; in the GUI, delete the domain) or delete "
                    "the sketch",
                )
            )
    extra: list[str] = []
    todo = list(expanded)
    while todo:
        name = todo.pop(0)
        for dependent in store.dependents(name):
            if dependent.name in seen or dependent.kind in MEMBERSHIP_KINDS:
                continue
            if dependent.kind == "sketch":
                issues.append(
                    Issue(name, f"is assigned to sketch {dependent.name!r}: unassign it first")
                )
                continue
            owner = owners.get(dependent.id)
            if owner is not None and owner[0].name not in seen:
                issues.append(
                    Issue(
                        name,
                        f"is used by {dependent.kind} {dependent.name!r} of sketch "
                        f"{owner[0].name!r}: delete the sketch first",
                    )
                )
                continue
            seen.add(dependent.name)
            extra.append(dependent.name)
            todo.append(dependent.name)
    return expanded, extra, issues


class DeleteEntities(Command):
    """Delete entities. With ``dependents`` what depends on them goes too (delete_plan);
    without, such a deletion is refused naming the dependents. Patches, zones and smoothing
    groups lose the deleted members (a group left without domains is deleted; a periodic
    patch whose partner goes is unpaired). A sketch takes its outputs along."""

    def __init__(self, names: list[str], dependents: bool = False):
        self.names = list(names)
        self.dependents = dependents
        self.label = f"Delete {len(self.names)} entities"
        if len(self.names) == 1:
            self.label = f"Delete {self.names[0]}"

    def execute(self, store: Store) -> int:
        expanded, extra, issues = delete_plan(store, self.names)
        if issues:
            raise ValidationError(issues)
        if extra and not self.dependents:
            going = set(extra)
            for name in expanded:
                for dependent in store.dependents(name):
                    if dependent.name in going:
                        issues.append(
                            Issue(
                                name,
                                f"is used by {dependent.kind} {dependent.name!r} (delete it "
                                "too, or delete with dependents=True)",
                            )
                        )
            raise ValidationError(issues)
        deleting = expanded + extra
        deleted_ids = {store.get(name).id for name in deleting}
        for group_id in drop_members(store, deleted_ids):
            deleting.append(store.get(group_id).name)
        for name in deleting:
            store.remove(name)
        return len(deleting)

    def script(self) -> str:
        if self.dependents:
            return f"og.delete({self.names!r}, dependents=True)"
        return f"og.delete({self.names!r})"
