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
_FUTURE_SECTIONS = ()
_SECTIONS = {
    "database": DATABASE_KINDS,
    "connectors": ("connector",),
    "domains": ("domain", "unstructured_domain"),
    "blocks": ("block", "volume"),
    "groups": ("smoothing_group",),
    "patches": ("patch",),
    "zones": ("zone",),
}
_RECORD_KEYS = ("id", "name", "type")
# Kinds also accepted (when reading) in an earlier section: smoothing groups of domains were
# written to "domains" before groups of blocks needed their own section after "blocks".
_LEGACY = {"domains": ("smoothing_group",)}


def _record(entity: Entity) -> dict:
    record = {"id": entity.id, "name": entity.name, "type": entity.kind}
    for key, value in entity.params.items():
        if key in record:
            # Would silently overwrite the record header (a Phase 4 bug stored the patch type
            # under "type"); parameters must never use the reserved keys.
            raise ValueError(f"{entity.kind} parameter {key!r} clashes with a project record key")
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
            if kind not in kinds and kind not in _LEGACY.get(section, ()):
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

    # Database entities first; connectors are validated against them. Periodic pairings refer
    # to other patches of the same section, so patches are added without them first and the
    # pairings applied in a second pass.
    # Within a section entities are added one by one in id order, so references inside the
    # section (surfaces -> curves) resolve.
    store = Store()
    deferred = []
    for section in _SECTIONS:
        store.begin()
        for _where, record in sorted(records[section], key=lambda item: item[1]["id"]):
            raw = {key: value for key, value in record.items() if key not in _RECORD_KEYS}
            if record["type"] == "patch" and raw.get("periodic") is not None:
                deferred.append((record, raw))
                raw = dict(raw)
                raw["periodic"] = None
            try:
                params = normalize_params(record["type"], raw, f"{source}: {record['name']}", store)
            except ValidationError as error:
                issues.extend(error.issues)
                continue
            store.next_id = record["id"]
            store.add(record["type"], params, record["name"])
        store.commit()
        if issues:
            raise ValidationError(issues)
    if deferred:
        store.begin()
        for record, raw in deferred:
            try:
                params = normalize_params("patch", raw, f"{source}: {record['name']}", store)
            except ValidationError as error:
                issues.extend(error.issues)
                continue
            store.update(record["id"], params=params)
        store.commit()
        if issues:
            raise ValidationError(issues)
    store.next_id = next_id
    store.name_counters = dict(data.get("name_counters", {}))
    return store


def save(store: Store, path: str | Path) -> None:
    Path(path).write_text(to_json(store), encoding="utf-8")


def load(path: str | Path) -> Store:
    path = Path(path)
    return from_json(path.read_text(encoding="utf-8"), source=path.name)
