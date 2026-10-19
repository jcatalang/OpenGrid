"""The complete cases of examples/12a_airfoil_cases.py and examples/12b_bodies_in_a_box.py build
as the examples say (they use the newest API: shapes, boundary types, periodic pairs, the
one-cell extrusion, surface meshes), so the recipes the README points to stay current."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_airfoil_cases():
    module = _load("12a_airfoil_cases")
    far = module.far_field_case()
    report = far.boundary_report()
    assert report.ok and len(report.boundary) == 6
    mesh = far.unstructured_mesh("flow")
    assert mesh["layers"] == 40 and len(mesh["quads"]) > 5000
    cascade = module.cascade_case()
    report = cascade.boundary_report()
    assert report.ok and len(report.periodic) == 1
    assert cascade.entity("pitch-a").params["periodic"]["translation"] == pytest.approx(
        (0.0, 1.0, 0.0)
    )
    names = cascade.extrude_model(0.05)
    assert names == ["passage-flow-ext", "front", "back"]
    report = cascade.boundary_report()
    assert report.ok and len(report.periodic) == 2


def test_wing_in_a_box():
    pytest.importorskip("netgen.meshing")
    module = _load("12b_bodies_in_a_box")
    og = module.wing_case()
    assert og.boundary_report().ok
    mesh = og.volume_mesh("fluid")
    assert len(mesh["prisms"]) > 10000 and len(mesh["tets"]) > 10000
    assert np.asarray(og.volume_quality("fluid")["scaled_jacobian"]).min() > 0.05
