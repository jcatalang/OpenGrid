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


def _members(raw, key: str, kind: str, store, where: str, issues: list[Issue]) -> tuple:
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
        if entity.kind != kind:
            issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a {kind}"))
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


def normalize_patch(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"patch_type", "members", "periodic"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a patch"))
    kind = raw.get("patch_type", "wall")
    if kind not in PATCH_TYPES:
        issues.append(
            Issue(f"{where}: patch_type", f"must be one of {', '.join(PATCH_TYPES)}, got {kind!r}")
        )
    members = _members(raw.get("members", ()), "connector", "connector", store, where, issues)
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
    members = _members(raw.get("members", ()), "domain", "domain", store, where, issues)
    if issues:
        raise ValidationError(issues)
    return {"members": members}


def member_names(params: dict, key: str, store) -> tuple:
    return tuple(store.get(item[key]).name for item in params["members"])
