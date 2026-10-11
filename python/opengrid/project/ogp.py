"""Read and write OpenGrid project files (.ogp, versioned JSON; see docs/ARCHITECTURE.md §6).

The writer is deterministic: sorted keys, one entity record per line, ``repr``-exact floats.
Hence save -> load -> save is byte-identical.
"""

import json
from pathlib import Path

from opengrid.model.entities import DATABASE_KINDS, normalize_params
from opengrid.model.store import Entity, Store
from opengrid.validation import Issue, ValidationError

FORMAT = "opengrid-project"
VERSION = 1

# Sections reserved for later phases; written empty so the schema is stable.
_FUTURE_SECTIONS = ("blocks", "patches", "zones")
_SECTIONS = {"database": DATABASE_KINDS, "connectors": ("connector",), "domains": ("domain",)}
_RECORD_KEYS = ("id", "name", "type")


def _record(entity: Entity) -> dict:
    record = {"id": entity.id, "name": entity.name, "type": entity.kind}
    for key, value in entity.params.items():
        record[key] = value
    return record


def _dump(value) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def to_json(store: Store) -> str:
    header = {
        "format": FORMAT,
        "version": VERSION,
        "units": "m",
        "tolerance": 1e-9,
        "next_id": store.next_id,
        "name_counters": dict(sorted(store.name_counters.items())),
        "settings": {},
    }
    sections = {}
    for section, kinds in _SECTIONS.items():
        sections[section] = [_record(entity) for entity in store if entity.kind in kinds]
    for name in _FUTURE_SECTIONS:
        sections[name] = []
    lines = ["{"]
    items = []
    for key in sorted(list(header) + list(sections)):
        if key in sections:
            records = sections[key]
            if not records:
                items.append(f' "{key}": []')
                continue
            body = ",\n".join(f"  {_dump(record)}" for record in records)
            items.append(f' "{key}": [\n{body}\n ]')
        else:
            items.append(f' "{key}": {_dump(header[key])}')
    lines.append(",\n".join(items))
    lines.append("}")
    return "\n".join(lines) + "\n"


def from_json(text: str, source: str = "<project>") -> Store:
    """Parse project JSON into a new Store. Raises ValidationError listing every problem."""
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        issue = Issue(f"{source}:{error.lineno}", f"invalid JSON: {error.msg}")
        raise ValidationError([issue]) from None
    if not isinstance(data, dict):
        raise ValidationError([Issue(source, "top level must be a JSON object")])
    if data.get("format") != FORMAT:
        raise ValidationError(
            [Issue(source, f"not an OpenGrid project (format={data.get('format')!r})")]
        )
    version = data.get("version")
    if version != VERSION:
        raise ValidationError(
            [Issue(source, f"unsupported project version {version!r} (expected {VERSION})")]
        )

    issues: list[Issue] = []
    seen_ids: set[int] = set()
    seen_names: set[str] = set()
    records: dict[str, list[tuple[str, dict]]] = {}
    for section, kinds in _SECTIONS.items():
        records[section] = []
        for index, record in enumerate(data.get(section, [])):
            where = f"{source}: {section}[{index}]"
            if not isinstance(record, dict) or any(key not in record for key in _RECORD_KEYS):
                issues.append(Issue(where, "record needs 'id', 'name' and 'type'"))
                continue
            entity_id = record["id"]
            name = record["name"]
            kind = record["type"]
            if isinstance(entity_id, bool) or not isinstance(entity_id, int) or entity_id < 1:
                issues.append(Issue(where, f"invalid id {entity_id!r}"))
                continue
            if entity_id in seen_ids:
                issues.append(Issue(where, f"duplicate id {entity_id}"))
                continue
            if name in seen_names:
                issues.append(Issue(where, f"duplicate name {name!r}"))
                continue
            if kind not in kinds:
                issues.append(Issue(where, f"unknown type {kind!r} in section {section!r}"))
                continue
            seen_ids.add(entity_id)
            seen_names.add(name)
            records[section].append((where, record))

    for section in _FUTURE_SECTIONS:
        if data.get(section):
            issues.append(Issue(f"{source}: {section}", "not supported in this version yet"))
    next_id = data.get("next_id", 1)
    if seen_ids and next_id <= max(seen_ids):
        issues.append(
            Issue(f"{source}: next_id", f"{next_id} must exceed the largest id {max(seen_ids)}")
        )
    if issues:
        raise ValidationError(issues)

    # Database entities first; connectors are validated against them.
    store = Store()
    for section in _SECTIONS:
        entities = []
        for _where, record in records[section]:
            raw = {key: value for key, value in record.items() if key not in _RECORD_KEYS}
            try:
                params = normalize_params(record["type"], raw, f"{source}: {record['name']}", store)
            except ValidationError as error:
                issues.extend(error.issues)
                continue
            entities.append(Entity(record["id"], record["name"], record["type"], params))
        if issues:
            raise ValidationError(issues)
        store.begin()
        for entity in sorted(entities, key=lambda e: e.id):
            store.next_id = entity.id
            store.add(entity.kind, entity.params, entity.name)
        store.commit()
    store.next_id = next_id
    store.name_counters = dict(data.get("name_counters", {}))
    return store


def save(store: Store, path: str | Path) -> None:
    Path(path).write_text(to_json(store), encoding="utf-8")


def load(path: str | Path) -> Store:
    path = Path(path)
    return from_json(path.read_text(encoding="utf-8"), source=path.name)
