"""The locked runtime environments of the downloadable bundles (packaging/locks): well formed,
without OpenGrid's own package, holding what the recipe asks for at run time, and agreeing
with the compilers and host packages pinned in packaging/conda/variants.yaml (OpenGrid's
package depends on the OCCT / CGNS / C++ runtime it was compiled with: a pin moved without a
new lock would fail the bundle builds, packaging/lock_env.py check, so it fails here first)."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
LOCKS = {"linux-64": "linux", "osx-arm64": "osx"}
# Run exports of the host packages: how many leading version parts the runtime must share
# with the build (occt >=8.0.1,<8.0.2; cgns >=4.5.2,<4.6), and at least the pinned version.
RUN_EXPORTS = {"occt": 3, "cgns": 2}
# The C++ runtime of each platform's compiler (libstdcxx >=16, libcxx >=23).
CXX_RUNTIME = {"linux": "libstdcxx", "osx": "libcxx"}
LINE = re.compile(
    r"^https://conda\.anaconda\.org/conda-forge/(?P<subdir>[\w-]+)/"
    r"(?P<file>[^/#]+)\.(conda|tar\.bz2)#[0-9a-f]{32}$"
)


def read_lock(subdir):
    """{name: version} of a lock."""
    text = (ROOT / "packaging" / "locks" / f"{subdir}.txt").read_text()
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    assert lines[0] == "@EXPLICIT"
    packages = {}
    for line in lines[1:]:
        match = LINE.match(line)
        assert match, line
        assert match["subdir"] in (subdir, "noarch"), line
        name, version, _build = match["file"].rsplit("-", 2)
        assert name not in packages, f"{name} twice"
        packages[name] = version
    return packages


def variant_pins():
    """{key: {platform or 'all': value}} of variants.yaml (its plain lists and if / then)."""
    pins = {}
    key = None
    platform = "all"
    for line in (ROOT / "packaging" / "conda" / "variants.yaml").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if not line.startswith(" "):
            key = line.strip().rstrip(":")
            pins[key] = {}
            continue
        item = line.strip()
        if item.startswith("- if:"):
            platform = item.split(":", 1)[1].strip()
            continue
        if item.startswith("then:"):
            pins[key][platform] = item.split(":", 1)[1].strip().strip('"')
            platform = "all"
            continue
        pins[key]["all"] = item.removeprefix("- ").strip('"')
    return pins


def run_requirements(platform):
    """The names of the recipe's run requirements on a platform ('linux' / 'osx')."""
    names = []
    inside = False
    condition = None
    for line in (ROOT / "packaging" / "conda" / "recipe.yaml").read_text().splitlines():
        if line.startswith("  run:"):
            inside = True
            continue
        if not inside or line.lstrip().startswith("#") or not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        if indent <= 2:
            break
        item = line.strip()
        if item.startswith("- if:"):
            condition = item.split(":", 1)[1].strip()
            continue
        if item.startswith("then:"):
            if condition == platform:
                names.append(item.split(":", 1)[1].strip())
            continue
        names.append(item.removeprefix("- ").split()[0])
    return names


def parts(version):
    return tuple(int(p) for p in re.findall(r"\d+", version))


@pytest.mark.parametrize("subdir", sorted(LOCKS))
def test_lock_holds_the_run_requirements(subdir):
    packages = read_lock(subdir)
    assert "opengrid" not in packages
    assert packages["python"].startswith("3.12.")
    missing = [name for name in run_requirements(LOCKS[subdir]) if name not in packages]
    assert not missing, f"{subdir}: {missing} not in the lock"


@pytest.mark.parametrize("subdir", sorted(LOCKS))
def test_lock_fits_the_pinned_host_packages(subdir):
    platform = LOCKS[subdir]
    packages = read_lock(subdir)
    pins = variant_pins()
    for name, shared in RUN_EXPORTS.items():
        pin = pins[name]["all"]
        assert parts(packages[name])[:shared] == parts(pin)[:shared], (name, packages[name])
        assert parts(packages[name]) >= parts(pin), (name, packages[name])
    compiler = pins["cxx_compiler_version"][platform]
    runtime = CXX_RUNTIME[platform]
    assert parts(packages[runtime])[0] >= parts(compiler)[0], (runtime, packages[runtime])
    if platform == "osx":
        pin = pins["llvm-openmp"]["all"]
        assert parts(packages["llvm-openmp"]) >= parts(pin)
