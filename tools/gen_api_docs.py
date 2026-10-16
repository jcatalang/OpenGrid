"""Generate docs/API.md, the reference of the scripting API (opengrid.Session), from the
docstrings. Every public method must be listed in exactly one group below
(tests/python/test_api_docs.py checks that and that docs/API.md is up to date).

    python tools/gen_api_docs.py          # rewrite docs/API.md
    python tools/gen_api_docs.py --check  # exit 1 if docs/API.md is stale
"""

import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "API.md"

GROUPS = [
    (
        "Geometry (database entities)",
        [
            "create_point",
            "create_line",
            "create_arc",
            "create_arc_3pt",
            "create_circle",
            "create_polyline",
            "create_spline",
            "create_nurbs",
            "create_plane",
            "create_ruled_surface",
            "create_extruded_surface",
            "create_revolved_surface",
            "create_nurbs_surface",
            "geometry",
            "length",
            "project_point",
            "project_to_surface",
        ],
    ),
    (
        "Connectors",
        [
            "create_connector",
            "create_connectors_on_curves",
            "create_chain_connector",
            "set_dimension",
            "set_distribution",
            "connector_nodes",
            "connector_report",
        ],
    ),
    (
        "Domains (2D, and faces in 3D)",
        [
            "create_structured_domain",
            "create_domain",
            "create_unstructured_domain",
            "set_domain_generator",
            "set_unstructured_generator",
            "rotate_domain",
            "smooth_domain",
            "domain_dimensions",
            "domain_grid",
            "unstructured_mesh",
            "domain_quality",
            "domain_status",
        ],
    ),
    (
        "Freehand curves, curve networks, live editing",
        [
            "create_freehand",
            "domains_from_curves",
            "move_handle",
        ],
    ),
    (
        "Multi-block smoothing groups",
        [
            "create_smoothing_group",
            "smooth_group",
            "group_status",
            "display_connector_nodes",
        ],
    ),
    (
        "Sketch topology (floating wireframe)",
        [
            "create_sketch",
            "edit_sketch",
            "regenerate_sketch",
            "sketch_report",
            "sketch_suggestions",
            "accept_sketch_suggestions",
            "assign_sketch_geometry",
        ],
    ),
    (
        "Blocks (structured 3D)",
        [
            "create_block",
            "create_extruded_block",
            "set_block_smoothing",
            "block_grid",
            "block_quality",
            "block_floor",
            "block_status",
        ],
    ),
    (
        "Volumes (unstructured / hybrid 3D)",
        [
            "create_volume",
            "set_volume_generator",
            "volume_mesh",
            "volume_quality",
            "volume_status",
        ],
    ),
    (
        "Topology editing",
        [
            "split_connector",
            "join_connectors",
            "split_domain",
            "join_domains",
            "split_block",
            "join_blocks",
            "ogrid_domain",
            "ogrid_block",
        ],
    ),
    (
        "Patches, zones and topology",
        [
            "create_patch",
            "set_patch_type",
            "add_to_patch",
            "remove_from_patch",
            "pair_periodic",
            "unpair_periodic",
            "create_zone",
            "add_to_zone",
            "remove_from_zone",
            "is_3d",
            "topology",
            "boundary_report",
            "unassigned_boundary_connectors",
        ],
    ),
    ("Import", ["import_airfoil_dat", "import_csv", "import_cad"]),
    ("Export", ["export_vtkhdf", "export_cgns", "export_gmsh", "export_plot3d"]),
    (
        "Editing, undo and the model",
        ["modify", "rename", "delete", "undo", "redo", "entity", "names"],
    ),
    ("Projects and journals", ["new", "save", "load", "journal_text", "save_journal"]),
]

HEADER = """# OpenGrid scripting API reference

Generated from the docstrings by `tools/gen_api_docs.py`; do not edit by hand. For concepts and
worked examples see `docs/SCRIPTING.md`; the golden cases in `examples/golden/` are complete
scripts.

```python
import opengrid

og = opengrid.Session()  # every method below is called on a Session
```

Module functions: `opengrid.replay(source)` runs a journal script and returns its session;
`opengrid.ValidationError` (with `.issues`, each an `opengrid.Issue(where, message)`) is raised
for invalid input; nothing is changed when it is raised.
"""


def public_methods(session_class) -> list[str]:
    return sorted(
        name
        for name, member in inspect.getmembers(session_class)
        if not name.startswith("_") and callable(member)
    )


def render() -> str:
    from opengrid.api import Session

    lines = [HEADER]
    for title, names in GROUPS:
        lines.append(f"## {title}\n")
        for name in names:
            method = getattr(Session, name)
            signature = (
                str(inspect.signature(method)).replace("(self, ", "(").replace("(self)", "()")
            )
            doc = inspect.getdoc(method) or ""
            lines.append(f"### `{name}{signature}`\n")
            lines.append(doc + "\n")
    return "\n".join(lines).rstrip() + "\n"


def check_groups() -> list[str]:
    from opengrid.api import Session

    listed = [name for _title, names in GROUPS for name in names]
    problems = []
    for name in sorted(set(public_methods(Session)) - set(listed)):
        problems.append(f"{name}: not listed in tools/gen_api_docs.py GROUPS")
    for name in sorted(set(listed) - set(public_methods(Session))):
        problems.append(f"{name}: listed but not a public Session method")
    for name in sorted({n for n in listed if listed.count(n) > 1}):
        problems.append(f"{name}: listed twice")
    return problems


def main() -> int:
    problems = check_groups()
    if problems:
        print("\n".join(problems))
        return 1
    text = render()
    if "--check" in sys.argv:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != text:
            print("docs/API.md is out of date: run python tools/gen_api_docs.py")
            return 1
        return 0
    OUTPUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
