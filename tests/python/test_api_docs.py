"""The scripting API reference (docs/API.md) covers every public Session method and is current."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "gen_api_docs", ROOT / "tools" / "gen_api_docs.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_method_is_documented_and_listed():
    import inspect

    from opengrid.api import Session

    tool = load_tool()
    assert tool.check_groups() == []
    for name in tool.public_methods(Session):
        assert inspect.getdoc(getattr(Session, name)), f"{name} has no docstring"


def test_api_reference_is_up_to_date():
    tool = load_tool()
    assert (ROOT / "docs" / "API.md").read_text(encoding="utf-8") == tool.render()
