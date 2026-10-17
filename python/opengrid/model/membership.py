"""Membership entities - patches, zones and smoothing groups - and what happens to them when
entities they list are deleted: they lose the members rather than block the deletion (a group
left without domains or blocks goes too; a periodic patch whose partner goes is unpaired).
Used by the delete command and by sketch regeneration (outputs no longer generated)."""

from opengrid.model.store import references

MEMBERSHIP_KINDS = ("patch", "zone", "smoothing_group")


def without_members(entity, deleted: set[int]) -> dict | None:
    """The params of a patch, zone or smoothing group without the deleted ids; None when a
    group would keep no domain or block."""
    params = dict(entity.params)
    if entity.kind == "patch":
        params["members"] = tuple(
            m for m in entity.params["members"] if not any(v in deleted for v in m.values())
        )
        pairing = entity.params.get("periodic")
        if pairing is not None and pairing["partner"] in deleted:
            params["periodic"] = None
        return params
    if entity.kind == "zone":
        params["members"] = tuple(m for m in entity.params["members"] if m["domain"] not in deleted)
        return params
    from opengrid.model.groups import _wall_ref

    params["domains"] = tuple(d for d in entity.params["domains"] if d not in deleted)
    params["blocks"] = tuple(b for b in entity.params["blocks"] if b not in deleted)
    if not params["domains"] and not params["blocks"]:
        return None
    params["walls"] = tuple(w for w in entity.params["walls"] if _wall_ref(w) not in deleted)
    if entity.params["slide"] != "auto":
        params["slide"] = tuple(c for c in entity.params["slide"] if c not in deleted)
    for key in ("pinned", "symmetry"):
        if key in entity.params:
            params[key] = tuple(c for c in entity.params[key] if c not in deleted)
    return params


def drop_members(store, deleted: set[int]) -> list[int]:
    """Take the deleted ids out of every patch, zone and smoothing group that lists them
    (inside the open transaction). Returns the ids of the groups left empty, for the caller
    to delete as well."""
    emptied = []
    for entity in list(store):
        if entity.kind not in MEMBERSHIP_KINDS or entity.id in deleted:
            continue
        if not (references(entity) & deleted):
            continue
        params = without_members(entity, deleted)
        if params is None:
            emptied.append(entity.id)
        else:
            store.update(entity.id, params=params)
    return emptied
