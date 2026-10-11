import opengrid
from opengrid import _core


def test_core_version_matches_package():
    assert _core.version() == opengrid.__version__
    assert opengrid.core_version() == opengrid.__version__


def test_add():
    assert _core.add(2, 3) == 5
    assert _core.add(a=-4, b=4) == 0


def test_add_is_64_bit():
    big = 2**40
    assert _core.add(big, big) == 2**41
