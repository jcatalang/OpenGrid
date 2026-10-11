import importlib.util
from pathlib import Path

_TOOL = Path(__file__).resolve().parents[2] / "tools" / "check_no_ternary.py"
_spec = importlib.util.spec_from_file_location("check_no_ternary", _TOOL)
check_no_ternary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_no_ternary)


def test_python_ifexp_detected(tmp_path):
    bad = tmp_path / "bad.py"
    bad.write_text("x = 1\ny = 2 if x else 3\n")
    problems = check_no_ternary.check_python(bad)
    assert len(problems) == 1
    assert ":2:" in problems[0]


def test_python_if_block_ok(tmp_path):
    good = tmp_path / "good.py"
    good.write_text("x = 1\nif x:\n    y = 2\nelse:\n    y = 3\ns = 'a ? b : c'\n")
    assert check_no_ternary.check_python(good) == []


def test_cpp_conditional_detected(tmp_path):
    bad = tmp_path / "bad.cpp"
    bad.write_text("int f(int a) {\n  return a > 0 ? a : -a;\n}\n")
    problems = check_no_ternary.check_cpp(bad)
    assert len(problems) == 1
    assert ":2:" in problems[0]


def test_cpp_literals_and_comments_ignored(tmp_path):
    good = tmp_path / "good.cpp"
    good.write_text(
        "// is it? yes\n"
        "/* multi\n line ? */\n"
        'const char* s = "what?";\n'
        "char c = '?';\n"
        'const char* r = R"x(raw ? ")x";\n'
        "long big = 1'000'000;\n"
        "int g() { return 1; }\n"
    )
    assert check_no_ternary.check_cpp(good) == []


def test_repository_is_clean():
    assert check_no_ternary.main([]) == 0
