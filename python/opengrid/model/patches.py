"""Patches (named boundary-condition groups of connectors) and zones (named groups of domains).

Patch types and their CGNS boundary-condition types:
    wall -> BCWall, inlet -> BCInflow, outlet -> BCOutflow, symmetry -> BCSymmetryPlane,
    farfield -> BCFarfield, periodic / interface / user -> BCTypeUserDefined.
Periodic patches are paired (see opengrid.model.periodic: a rotation about an axis through a
centre, then a translation; 3D faces are matched by opengrid.model.periodic_faces); paired
periodic boundaries are exported as periodic connectivity, not as boundary conditions.
"""

import math

from opengrid.model.periodic import DEFAULT_AXIS
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
    known = {"partner", "rotation_center", "rotation_axis", "rotation_angle", "translation"}
    for key in sorted(set(raw) - known):
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
    axis = _vec3(raw.get("rotation_axis", DEFAULT_AXIS), "rotation_axis", where, issues)
    if axis is not None:
        length = math.sqrt(sum(v * v for v in axis))
        if length == 0.0 or not math.isfinite(length):
            issues.append(Issue(f"{where}: periodic.rotation_axis", "must not be zero"))
        elif abs(length - 1.0) > 1e-12:
            axis = tuple(v / length for v in axis)  # a unit vector (kept once it is one)
    return {
        "partner": entity.id,
        "rotation_center": _vec3(
            raw.get("rotation_center", (0.0, 0.0, 0.0)), "rotation_center", where, issues
        ),
        "rotation_axis": axis,
        "rotation_angle": float(angle),
        "translation": _vec3(raw.get("translation", (0.0, 0.0, 0.0)), "translation", where, issues),
    }


BLOCK_FACES = ("imin", "imax", "jmin", "jmax", "kmin", "kmax")
VOLUME_FACES = ("bottom", "top")  # of an extruded volume (its sides are its connectors)


def _patch_members(raw, store, where: str, issues: list[Issue]) -> tuple:
    """Patch members: connectors (2D boundaries, and the side walls of extruded volumes),
    domains (faces of 3D blocks and volumes, the bottom of an extruded volume), block faces
    ({"block": name or id, "face": "imin" ...}, or a (block, face) pair) and the ends of
    extruded volumes ({"volume": name or id, "face": "top" | "bottom"})."""
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
        elif isinstance(item, dict) and "volume" in item:
            volume = item.get("volume")
            face = item.get("face")
            if (
                isinstance(volume, bool)
                or not isinstance(volume, (str, int))
                or volume not in store
                or store.get(volume).kind != "volume"
            ):
                issues.append(Issue(label, f"no volume {volume!r}"))
                continue
            entity = store.get(volume)
            if entity.params["faces"] is not None:
                issues.append(
                    Issue(label, f"{entity.name!r} is not extruded: patch its face domains")
                )
                continue
            if face not in VOLUME_FACES:
                issues.append(Issue(label, f"face must be one of {', '.join(VOLUME_FACES)}"))
                continue
            key = (entity.id, face)
            if key in seen:
                issues.append(Issue(label, "listed twice"))
                continue
            seen.add(key)
            result.append({"volume": entity.id, "face": face})
            continue
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
    if "volume" in item:
        return f"{store.get(item['volume']).name}.{item['face']}"
    ref = item.get("connector", item.get("domain"))
    return store.get(ref).name


def member_script(item: dict, store):
    if "block" in item:
        return {"block": store.get(item["block"]).name, "face": item["face"]}
    if "volume" in item:
        return {"volume": store.get(item["volume"]).name, "face": item["face"]}
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


def set_boundary_type(store, members, patch_type: str, name: str | None = None) -> list[str]:
    """Put boundary members (connector names, domains, ``{"block", "face"}``, ``{"volume",
    "face"}``) into one patch of ``patch_type``: the patch called ``name`` (created when
    absent), else the first existing patch of that type (one named after the type first),
    else a new patch named after the type. The members leave whatever patches they were in
    (a member belongs to one patch). Returns [patch name]."""
    issues: list[Issue] = []
    if patch_type not in PATCH_TYPES:
        issues.append(
            Issue("patch_type", f"must be one of {', '.join(PATCH_TYPES)}, got {patch_type!r}")
        )
    if isinstance(members, (str, int, dict)):
        members = [members]
    wanted = _patch_members(members, store, "set_boundary_type", issues)
    if issues:
        raise ValidationError(issues)
    if not wanted:
        raise ValidationError([Issue("set_boundary_type", "no boundary members given")])
    refuse_shared_connectors(store, wanted)
    target = None
    if name is not None and name in store:
        target = store.get(name)
        if target.kind != "patch":
            raise ValidationError([Issue(name, f"is a {target.kind}, not a patch")])
        if target.params["patch_type"] != patch_type:
            raise ValidationError(
                [Issue(name, f"is a {target.params['patch_type']} patch, not {patch_type}")]
            )
    elif name is None:
        candidates = [p for p in store.entities("patch") if p.params["patch_type"] == patch_type]
        named = [p for p in candidates if p.name == patch_type]
        if named:
            target = named[0]
        elif candidates:
            target = candidates[0]
    for patch in store.entities("patch"):
        if target is not None and patch.id == target.id:
            continue
        kept = tuple(m for m in patch.params["members"] if m not in wanted)
        if len(kept) != len(patch.params["members"]):
            params = dict(patch.params)
            params["members"] = kept
            store.update(patch.id, params=params)
    if target is None:
        new_name = name
        if new_name is None:
            new_name = patch_type
            number = 1
            while new_name in store:
                number += 1
                new_name = f"{patch_type}-{number}"
        params = normalize_patch(
            {"patch_type": patch_type, "members": wanted, "periodic": None}, store, new_name
        )
        target = store.add("patch", params, new_name)
        return [target.name]
    existing = list(target.params["members"])
    for member in wanted:
        if member not in existing:
            existing.append(member)
    params = dict(target.params)
    params["members"] = tuple(existing)
    store.update(target.id, params=params)
    return [target.name]


def partner_of(patch, store) -> str | None:
    """The name of a periodic patch's partner, paired from either side; None when unpaired."""
    pairing = patch.params.get("periodic")
    if pairing is not None:
        return store.get(pairing["partner"]).name
    for other in store.entities("patch"):
        other_pairing = other.params.get("periodic")
        if other_pairing is not None and other_pairing["partner"] == patch.id:
            return other.name
    return None


def shared_connectors(store) -> dict[int, list[str]]:
    """{connector id: [domain names]} for the connectors used by two domains of a 2D model:
    interfaces, found automatically, that take no boundary type. Empty for 3D models (there a
    connector is the side wall of an extruded volume)."""
    if store.entities("block") or store.entities("volume"):
        return {}
    uses: dict[int, list[str]] = {}
    for domain in store:
        if domain.kind == "domain":
            loops = domain.params["edges"]
        elif domain.kind == "unstructured_domain":
            loops = domain.params["loops"]
        else:
            continue
        for loop in loops:
            for item in loop:
                uses.setdefault(item["connector"], []).append(domain.name)
    return {c: names for c, names in uses.items() if len(names) >= 2}


def refuse_shared_connectors(store, members) -> None:
    """A connector shared by two domains is an interface, not a boundary: refuse to put it in
    a patch, saying so."""
    shared = shared_connectors(store)
    issues = []
    for member in members:
        if "connector" in member and member["connector"] in shared:
            domains = " and ".join(shared[member["connector"]])
            issues.append(
                Issue(
                    store.get(member["connector"]).name,
                    f"is shared by {domains}: it is an interface between them, found "
                    "automatically, and takes no boundary type",
                )
            )
    if issues:
        raise ValidationError(issues)


def clear_boundary_type(store, members) -> list[str]:
    """Take boundary members out of every patch (they are then unassigned). Returns the names
    of the patches that changed."""
    issues: list[Issue] = []
    if isinstance(members, (str, int, dict)):
        members = [members]
    wanted = _patch_members(members, store, "clear_boundary_type", issues)
    if issues:
        raise ValidationError(issues)
    changed = []
    for patch in store.entities("patch"):
        kept = tuple(m for m in patch.params["members"] if m not in wanted)
        if len(kept) != len(patch.params["members"]):
            params = dict(patch.params)
            params["members"] = kept
            store.update(patch.id, params=params)
            changed.append(patch.name)
    return changed


def pair_names(store, name: str | None = None) -> tuple[str, str]:
    """The names of a new periodic pair: ``<name>-a`` / ``<name>-b`` (``periodic-1a`` /
    ``periodic-1b``, ``periodic-2a`` / ``-2b``, ... by default), both unused."""
    number = 0
    while True:
        number += 1
        base = name
        if base is None:
            base = f"periodic-{number}"
        elif number > 1:
            base = f"{name}-{number}"
        joint = "-"
        if base[-1].isdigit():
            joint = ""  # periodic-1a / periodic-1b, span-2a / span-2b
        candidates = (f"{base}{joint}a", f"{base}{joint}b")
        if not any(c in store for c in candidates):
            return candidates


def make_periodic_pair(
    store,
    source_members,
    target_members,
    translation=(0.0, 0.0, 0.0),
    rotation_angle: float = 0.0,
    rotation_center=(0.0, 0.0, 0.0),
    rotation_axis=(0.0, 0.0, 1.0),
    name: str | None = None,
) -> list[str]:
    """Two periodic patches (``<name>-a`` / ``<name>-b``, ``periodic-1`` and up by default, see
    ``pair_names``) from the two sides' members, paired source -> target with the transform.
    The members leave the patches they were in. Returns [source patch, target patch]."""
    names = []
    for members, candidate in zip(
        (source_members, target_members), pair_names(store, name), strict=True
    ):
        names.append(set_boundary_type(store, members, "periodic", candidate)[0])
    source = store.get(names[0])
    params = dict(source.params)
    params["periodic"] = {
        "partner": store.get(names[1]).id,
        "rotation_center": rotation_center,
        "rotation_axis": rotation_axis,
        "rotation_angle": rotation_angle,
        "translation": translation,
    }
    store.update(source.id, params=normalize_patch(params, store, source.name))
    return names
