import math

import pytest
from helpers import naca0012_selig, write_selig

import opengrid
from opengrid.io import read_airfoil_dat, read_points_csv


def test_selig(tmp_path):
    path = tmp_path / "naca0012.dat"
    points = naca0012_selig(41)
    write_selig(path, points)
    data = read_airfoil_dat(path)
    assert data.format == "selig"
    assert data.title == "NACA 0012"
    assert len(data.points) == 81
    assert data.leading_edge_index() == 40
    assert data.points[40] == pytest.approx((0.0, 0.0))


def test_lednicer_matches_selig(tmp_path):
    surface = naca0012_selig(21)[20:]  # LE -> TE lower surface
    upper = [(x, -y) for x, y in surface]
    lines = ["NACA 0012 Lednicer", f"{len(upper)}. {len(surface)}.", ""]
    lines += [f"{x:.10f} {y:.10f}" for x, y in upper]
    lines.append("")
    lines += [f"{x:.10f} {y:.10f}" for x, y in surface]
    path = tmp_path / "led.dat"
    path.write_text("\n".join(lines) + "\n")
    data = read_airfoil_dat(path)
    assert data.format == "lednicer"
    assert len(data.points) == len(upper) + len(surface) - 1
    assert data.points[0] == pytest.approx(upper[-1])
    assert data.points[-1] == pytest.approx(surface[-1])


def test_airfoil_errors_report_line_numbers(tmp_path):
    path = tmp_path / "bad.dat"
    path.write_text("Title\n1.0 0.0\n0.5 abc\n0.0 0.0\n0.5 -0.1 7\n")
    with pytest.raises(opengrid.ValidationError) as error:
        read_airfoil_dat(path)
    text = str(error.value)
    assert "bad.dat:3" in text
    assert "bad.dat:5" in text


def test_airfoil_duplicate_point(tmp_path):
    path = tmp_path / "dup.dat"
    path.write_text("1 0\n0.5 0.1\n0.5 0.1\n0 0\n0.5 -0.1\n1 0\n")
    with pytest.raises(opengrid.ValidationError, match=r"dup.dat:3.*line 2"):
        read_airfoil_dat(path)


def test_lednicer_count_mismatch(tmp_path):
    path = tmp_path / "led.dat"
    path.write_text("Title\n3. 3.\n0 0\n0.5 0.1\n1 0\n0 0\n0.5 -0.1\n")
    with pytest.raises(opengrid.ValidationError, match=r"led.dat:2.*3 \+ 3"):
        read_airfoil_dat(path)


def test_import_modes(tmp_path):
    path = tmp_path / "naca0012.dat"
    write_selig(path, naca0012_selig(41))
    og = opengrid.Session()
    assert og.import_airfoil_dat(path) == ["naca0012"]
    names = og.import_airfoil_dat(path, name="split", mode="split")
    assert names == ["split-upper", "split-lower"]
    og.import_airfoil_dat(path, name="raw", mode="polyline")
    upper = og.geometry("split-upper")
    lower = og.geometry("split-lower")
    # Both run LE -> TE and meet at the LE with opposite (G1) tangents.
    assert upper.eval(0.0) == pytest.approx([0.0, 0.0, 0.0], abs=1e-15)
    assert lower.eval(0.0) == pytest.approx([0.0, 0.0, 0.0], abs=1e-15)
    tu = upper.d1(0.0)
    tl = lower.d1(0.0)
    cos = (tu[0] * tl[0] + tu[1] * tl[1]) / (math.hypot(*tu[:2]) * math.hypot(*tl[:2]))
    assert cos == pytest.approx(-1.0, abs=1e-12)
    # The single spline length is close to the NACA 0012 perimeter (~2.0393 chords).
    assert og.length("naca0012") == pytest.approx(2.0393, abs=2e-3)
    og.undo()
    assert "raw" not in og.store


def test_csv_variants(tmp_path):
    path = tmp_path / "a.csv"
    path.write_text("# comment\nY;X\n1;0\n2;1\n0;2\n")
    assert read_points_csv(path) == [(0.0, 1.0, 0.0), (1.0, 2.0, 0.0), (2.0, 0.0, 0.0)]
    path.write_text("0 0 1\n1\t1\t2  # trailing comment\n2,0,3\n")
    assert read_points_csv(path) == [(0.0, 0.0, 1.0), (1.0, 1.0, 2.0), (2.0, 0.0, 3.0)]


def test_csv_errors_report_line_numbers(tmp_path):
    path = tmp_path / "b.csv"
    path.write_text("x,y\n0,0\n1,oops\n2\n")
    with pytest.raises(opengrid.ValidationError) as error:
        read_points_csv(path)
    text = str(error.value)
    assert "b.csv:3" in text
    assert "b.csv:4" in text
    path.write_text("a,b\n0,0\n")
    with pytest.raises(opengrid.ValidationError, match="no 'x' column"):
        read_points_csv(path)


def test_import_csv_creates_spline(tmp_path):
    path = tmp_path / "pts.csv"
    path.write_text("0,0\n1,1\n2,0\n")
    og = opengrid.Session()
    assert og.import_csv(path, name="hump") == ["hump"]
    assert og.entity("hump").kind == "spline"
