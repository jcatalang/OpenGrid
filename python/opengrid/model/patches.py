"""Patches (named boundary-condition groups of connectors) and zones (named groups of domains).

Patch types and their CGNS boundary-condition types:
    wall -> BCWall, inlet -> BCInflow, outlet -> BCOutflow, symmetry -> BCSymmetryPlane,
    farfield -> BCFarfield, periodic / interface / user -> BCTypeUserDefined.
Periodic patches are paired (see opengrid.model.periodic); paired periodic boundaries are
exported as periodic connectivity, not as boundary conditions.
"""

import math

from opengrid.validation import Issue, ValidationError

PATCH_TYPES = ("wall", "inlet", "outlet", "symmetry", "farfield", "periodic", "interface", "user")
CGNS_BC_TYPES = {
    "wall": "BCWall",
    "inlet": "BCInflow",
    "outlet": "BCOutflow",
    "symmetry": "BCSymmetryPlane",
    "farfield": "BCFarfield",
    "periodic": "BCTypeUserDefined",
    "interface": "BCTypeUserDefined",
    "user": "BCTypeUserDefined",
}


def _members(raw, key: str, kind, store, where: str, issues: list[Issue]) -> tuple:
    """Member references; ``kind`` is an entity kind or a tuple of accepted kinds."""
    kinds = kind
    if isinstance(kinds, str):
        kinds = (kinds,)
    if isinstance(raw, (str, int, dict)):
        raw = [raw]
    result = []
    seen = set()
    for index, item in enumerate(raw):
        label = f"{where}: members[{index}]"
        ref = item
        if isinstance(item, dict):
            ref = item.get(key)
        if isinstance(ref, bool) or not isinstance(ref, (str, int)):
            issues.append(Issue(label, f"expected a {kind} name or id, got {ref!r}"))
            continue
        if ref not in store:
            issues.append(Issue(label, f"no entity {ref!r}"))
            continue
        entity = store.get(ref)
        if entity.kind not in kinds:
            issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a {kinds[0]}"))
            continue
        if entity.id in seen:
            issues.append(Issue(label, f"{entity.name!r} listed twice"))
            continue
        seen.add(entity.id)
        result.append({key: entity.id})
    return tuple(result)


def _vec3(value, field, where, issues):
    if isinstance(value, (str, bytes)) or not hasattr(value, "__len__") or len(value) not in (2, 3):
        issues.append(
            Issue(f"{where}: periodic.{field}", f"expected 2 or 3 numbers, got {value!r}")
        )
        return None
    numbers = [float(v) for v in value]
    if len(numbers) == 2:
        numbers.append(0.0)
    return tuple(numbers)


def _periodic(raw, store, where: str, issues: list[Issue]) -> dict | None:
    if raw is None:
        return None
    if not isinstance(raw, dict) or "partner" not in raw:
        issues.append(Issue(f"{where}: periodic", "needs 'partner' (and a transform)"))
        return None
    for key in sorted(set(raw) - {"partner", "rotation_center", "rotation_angle", "translation"}):
        issues.append(Issue(f"{where}: periodic.{key}", "unexpected key"))
    partner = raw["partner"]
    if isinstance(partner, bool) or not isinstance(partner, (str, int)) or partner not in store:
        issues.append(Issue(f"{where}: periodic.partner", f"no patch {partner!r}"))
        return None
    entity = store.get(partner)
    if entity.kind != "patch":
        issues.append(
            Issue(f"{where}: periodic.partner", f"{entity.name!r} is a {entity.kind}, not a patch")
        )
        return None
    angle = raw.get("rotation_angle", 0.0)
    if isinstance(angle, bool) or not isinstance(angle, (int, float)) or not math.isfinite(angle):
        issues.append(
            Issue(f"{where}: periodic.rotation_angle", f"expected a number, got {angle!r}")
        )
        return None
    return {
        "partner": entity.id,
        "rotation_center": _vec3(
            raw.get("rotation_center", (0.0, 0.0, 0.0)), "rotation_center", where, issues
        ),
        "rotation_angle": float(angle),
        "translation": _vec3(raw.get("translation", (0.0, 0.0, 0.0)), "translation", where, issues),
    }


BLOCK_FACES = ("imin", "imax", "jmin", "jmax", "kmin", "kmax")


def _patch_members(raw, store, where: str, issues: list[Issue]) -> tuple:
    """Patch members: connectors (2D boundaries), structured domains (faces of 3D blocks) or
    block faces ({"block": name or id, "face": "imin" ...}, or a (block, face) pair)."""
    if isinstance(raw, (str, int, dict)):
        raw = [raw]
    result = []
    seen = set()
    for index, item in enumerate(raw):
        label = f"{where}: members[{index}]"
        block_face = None
        ref = item
        if isinstance(item, dict) and "block" in item:
            block_face = (item.get("block"), item.get("face"))
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            block_face = tuple(item)
        elif isinstance(item, dict):
            ref = item.get("connector", item.get("domain"))
        if block_face is not None:
            block, face = block_face
            if isinstance(block, bool) or not isinstance(block, (str, int)) or block not in store:
                issues.append(Issue(label, f"no block {block!r}"))
                continue
            entity = store.get(block)
            if entity.kind != "block":
                issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a block"))
                continue
            if face not in BLOCK_FACES:
                issues.append(Issue(label, f"face must be one of {', '.join(BLOCK_FACES)}"))
                continue
            member = {"block": entity.id, "face": face}
            key = (entity.id, face)
        else:
            if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
                issues.append(Issue(label, f"no entity {ref!r}"))
                continue
            entity = store.get(ref)
            if entity.kind == "connector":
                member = {"connector": entity.id}
            elif entity.kind in ("domain", "unstructured_domain"):
                member = {"domain": entity.id}
            else:
                issues.append(
                    Issue(
                        label,
                        f"{entity.name!r} is a {entity.kind}; patch members are connectors, "
                        "structured domains or block faces",
                    )
                )
                continue
            key = entity.id
        if key in seen:
            issues.append(Issue(label, "listed twice"))
            continue
        seen.add(key)
        result.append(member)
    return tuple(result)


def member_label(item: dict, store) -> str:
    if "block" in item:
        return f"{store.get(item['block']).name}.{item['face']}"
    ref = item.get("connector", item.get("domain"))
    return store.get(ref).name


def member_script(item: dict, store):
    if "block" in item:
        return {"block": store.get(item["block"]).name, "face": item["face"]}
    return member_label(item, store)


def normalize_patch(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"patch_type", "members", "periodic"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a patch"))
    kind = raw.get("patch_type", "wall")
    if kind not in PATCH_TYPES:
        issues.append(
            Issue(f"{where}: patch_type", f"must be one of {', '.join(PATCH_TYPES)}, got {kind!r}")
        )
    members = _patch_members(raw.get("members", ()), store, where, issues)
    periodic = _periodic(raw.get("periodic"), store, where, issues)
    if periodic is not None and kind != "periodic":
        issues.append(
            Issue(f"{where}: periodic", f"only periodic patches can be paired (type is {kind!r})")
        )
    if issues:
        raise ValidationError(issues)
    return {"patch_type": kind, "members": members, "periodic": periodic}


def normalize_zone(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"members"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a zone"))
    members = _members(
        raw.get("members", ()),
        "domain",
        ("domain", "unstructured_domain", "block", "volume"),
        store,
        where,
        issues,
    )
    if issues:
        raise ValidationError(issues)
    return {"members": members}


def member_names(params: dict, key: str, store) -> tuple:
    return tuple(member_label(item, store) for item in params["members"])
