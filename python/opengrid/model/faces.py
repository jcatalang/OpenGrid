"""The faces of blocks and volumes and their identifiers (what boundary types are given to in
3D, as curves are in 2D).

A block has six faces, ``<block>.imin`` ... ``<block>.kmax``. An extruded volume has
``<volume>.bottom`` (its base domain), ``<volume>.top`` and one side wall per curve of its base
domain, ``<volume>.side-1`` ... ``side-n`` in the order of the base's loops. A volume bounded by
domains has those domains as faces (their own names). Each face is also a patch member
(``{"block": id, "face": "imin"}``, ``{"volume": id, "face": "top"}``, ``{"connector": id}`` for
a side wall, ``{"domain": id}``): ``face_member`` and ``face_identifier`` convert between the
two, and patches accept identifiers wherever they take members.
"""

from dataclasses import dataclass

from opengrid.model.patches import BLOCK_FACES


@dataclass(frozen=True)
class Face:
    identifier: str  # "block-1.imin", "volume-1.side-2", "dom-3"
    short: str  # "imin", "side-2", "dom-3" (under its block or volume)
    member: dict  # the patch member, with ids
    owner: int  # the block or volume


def faces_of(store, entity) -> list[Face]:
    """The faces of a block or volume (none for other entities)."""
    if entity.kind == "block":
        return [
            Face(f"{entity.name}.{face}", face, {"block": entity.id, "face": face}, entity.id)
            for face in BLOCK_FACES
        ]
    if entity.kind != "volume":
        return []
    if entity.params["faces"] is not None:
        faces = []
        for domain in entity.params["faces"]:
            name = store.get(domain).name
            faces.append(Face(name, name, {"domain": domain}, entity.id))
        return faces
    faces = [
        Face(f"{entity.name}.bottom", "bottom", {"volume": entity.id, "face": "bottom"}, entity.id),
        Face(f"{entity.name}.top", "top", {"volume": entity.id, "face": "top"}, entity.id),
    ]
    for k, connector in enumerate(side_curves(store, entity), start=1):
        faces.append(
            Face(f"{entity.name}.side-{k}", f"side-{k}", {"connector": connector}, entity.id)
        )
    return faces


def side_curves(store, volume) -> list[int]:
    """The curves (connector ids) an extruded volume's side walls rise from, in loop order."""
    base = store.get(volume.params["extrusion"]["domain"])
    curves = []
    for loop in base.params.get("loops", ()):
        for item in loop:
            if item["connector"] not in curves:
                curves.append(item["connector"])
    return curves


def all_faces(store) -> list[Face]:
    faces = []
    for entity in store:
        faces.extend(faces_of(store, entity))
    return faces


def face_identifier(store, member: dict) -> str:
    """The identifier of a patch member that is a face (or the name of the entity it names:
    a 2D curve, a domain that bounds no volume)."""
    if "block" in member:
        return f"{store.get(member['block']).name}.{member['face']}"
    if "volume" in member:
        return f"{store.get(member['volume']).name}.{member['face']}"
    for face in all_faces(store):
        if face.member == member:
            return face.identifier
    if "domain" in member:
        for entity in store.entities("volume"):
            extrusion = entity.params.get("extrusion")
            if extrusion is not None and extrusion["domain"] == member["domain"]:
                return f"{entity.name}.bottom"
    return store.get(member.get("connector", member.get("domain"))).name


def face_member(store, identifier: str) -> dict | None:
    """The patch member a face identifier stands for, or None (not a face identifier)."""
    owner, _dot, short = identifier.rpartition(".")
    if not owner or owner not in store:
        return None
    for face in faces_of(store, store.get(owner)):
        if face.short == short:
            return dict(face.member)
    return None
