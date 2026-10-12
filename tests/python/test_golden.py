"""The golden validation cases of SPEC.md, reproduced from scripts alone (examples/golden): each
script builds its model through the public API and checks its invariants."""

import importlib.util
import sys
from pathlib import Path

import pytest

GOLDEN = Path(__file__).resolve().parents[2] / "examples" / "golden"
CASES = sorted(GOLDEN.glob("g[0-9][0-9]_*.py"))
NEEDS_NETGEN = {"g09_sphere_in_box", "g10_wing_stub"}


def test_all_ten_cases_exist():
    assert len(CASES) == 10


@pytest.mark.parametrize("path", CASES, ids=[p.stem for p in CASES])
def test_golden_case(path):
    if path.stem in NEEDS_NETGEN:
        pytest.importorskip("netgen.meshing")
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    spec.loader.exec_module(module)
    og = module.build()
    result = module.check(og)
    assert isinstance(result, dict) and result
    # The journal of the script replays to the same model.
    import opengrid
    from opengrid.project import to_json

    assert to_json(opengrid.replay(og.journal_text()).store) == to_json(og.store)
