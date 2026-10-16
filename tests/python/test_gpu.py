"""Phase 11d: GPU compute (OpenCL) - detection, modes, and the experimental GPU sweeps."""

import numpy as np
import pytest

import opengrid
from opengrid.app.main_window import compute_summary


def test_gpu_info_and_modes():
    info = opengrid.gpu_info()
    assert set(info) == {
        "built",
        "available",
        "platform",
        "device",
        "version",
        "compute_units",
        "float64",
        "reason",
    }
    if not info["available"]:
        assert info["reason"]
    assert opengrid.gpu_mode() == "off"  # the default
    with pytest.raises(ValueError, match="expected 'off', 'auto' or 'on'"):
        opengrid.set_gpu("fast")
    assert "CPU threads" in compute_summary()


def test_gpu_group_smoothing_matches_the_cpu():
    if not opengrid.gpu_info()["available"]:
        pytest.skip("no GPU")
    from test_templates import channel

    og = opengrid.Session()
    og.create_template(
        "bump_channel",
        **channel(og, bump=True),
        blocks=(3, 1),
        wall=5e-3,
        length=1.0,
        three_d_blocks=True,
    )
    group = og.names("smoothing_group")[0]
    cpu = og.smooth_group(group)
    try:
        opengrid.set_gpu("on")
        og.groups._results.clear()
        gpu = og.smooth_group(group)
    finally:
        opengrid.set_gpu("off")
    assert gpu.converged
    assert any("on the GPU" in w for w in gpu.warnings)
    worst = max(float(np.abs(gpu.grids[b] - cpu.grids[b]).max()) for b in cpu.grids)
    assert worst < 1e-3  # the channel is 3 x 1 x 1; the GPU stops at a relative 1e-5
