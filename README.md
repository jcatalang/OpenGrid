# OpenGrid

Interactive block-structured / unstructured / hybrid mesh generator for CFD (2D and 3D). C++20 core, Python bindings (nanobind),
PySide6 + VTK GUI, and a Python scripting API.

Status: all phases (0–9c) done. What it does:

- **Geometry:** lines, arcs, polylines, splines, NURBS curves and surfaces; airfoil `.dat`, CSV and
  STEP / IGES import (OpenCASCADE, optional) with healing and automatic connectors and domains.
- **Connectors:** uniform, geometric, tanh and sine distributions with end spacings, breakpoints,
  spacing propagation; split and join.
- **2D / surface meshing:** structured domains (TFI, Hermite TFI, elliptic smoothing with wall
  control), unstructured domains (constrained Delaunay, advancing front, quads), boundary layers.
- **3D:** structured blocks (six faces, extrusion, hyperbolic extrusion, elliptic smoothing) and
  unstructured / hybrid volumes (Netgen tetrahedra, prism and hexahedron layers with fans,
  pyramids, own optimiser), with a quality floor every cell must meet.
- **Multi-block topology:** split / join domains and blocks, O-grid (butterfly) domains and
  blocks for all-hexahedral meshes.
- **Export:** CGNS (structured, unstructured, mixed), VTKHDF, Gmsh 4.1, Plot3D; patches and
  zones, boundary checks before export.
- **Everything is scriptable:** the GUI and scripts share one API; every action is an undoable
  command recorded in a replayable Python journal.
- **Fast and deterministic:** parallel (OpenMP) hot loops whose results do not depend on the
  thread count; 1M-cell 2D meshes in 2–5 s, a 500 x 200 elliptic grid with wall control in
  about 1 s (`examples/09c_performance.py`; `opengrid.set_threads(n)`).

See `docs/SCRIPTING.md` (guide), `docs/API.md` (reference), `examples/` (one demo per phase)
and `examples/golden/` (the validation cases).

```python
import opengrid

og = opengrid.Session()
og.create_circle((0, 0), 0.5, name="cylinder")
og.create_circle((0, 0), 10.0, name="far")
og.create_line((0.5, 0), (10, 0), name="cut")
og.create_connector("cylinder", 129, name="c-cyl")
og.create_connector("far", 129, name="c-far")
og.create_connector(
    "cut", 81, {"type": "tanh", "begin_spacing": 1e-4, "end_spacing": 1.0}, name="c-cut"
)
og.create_domain([[("c-cyl", True)], ["c-cut"], ["c-far"], [("c-cut", True)]], name="ogrid")
og.create_patch("wall", "wall", ["c-cyl"])
og.create_patch("farfield", "farfield", ["c-far"])
og.export_cgns("cylinder.cgns")
print(og.journal_text())  # replayable script of everything done
```

## Install
- **macOS (Apple silicon):** `OpenGrid.app` / the DMG built by `packaging/macos/build_app.sh`
  (self-contained; signed ad hoc, so the first launch of a downloaded copy needs right-click →
  Open).
- **Linux / macOS, conda package:** build it with
  `rattler-build build --recipe packaging/conda/recipe.yaml -m packaging/conda/variants.yaml
  -c conda-forge --output-dir packaging/output`, then
  `conda create -n og -c file://$PWD/packaging/output -c conda-forge opengrid` and run
  `opengrid` (GUI) or `python -c "import opengrid"`.

## Setup (conda-forge, development)
```bash
conda env create -f environment.yml
conda activate opengrid
pip install --no-build-isolation -e .
python -m opengrid
```

### macOS compiler note
Apple Clang ≥ 15 is required. If `clang++ --version` reports an older version because
`xcode-select` points to an old Xcode, either switch globally:
```bash
sudo xcode-select -s /Library/Developer/CommandLineTools
```
or only inside the environment:
```bash
conda env config vars set DEVELOPER_DIR=/Library/Developer/CommandLineTools -n opengrid
conda deactivate && conda activate opengrid
```

## Setup (uv, alternative)
Needs CMake ≥ 3.28, Ninja, a C++20 compiler, and (from Phase 3/4) HDF5 and CGNS from your
system package manager (Homebrew: `brew install cmake ninja googletest hdf5 cgns`;
Ubuntu: `apt install cmake ninja-build libgtest-dev libhdf5-dev libcgns-dev`).
```bash
uv venv --python 3.12 && source .venv/bin/activate
uv pip install scikit-build-core nanobind
uv pip install --no-build-isolation -e ".[gui,test]"
```

## C++ development
```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
```

## Tests and lint
```bash
pytest -q
ruff check . && ruff format --check .
python tools/check_no_ternary.py
```

## Documentation
- `docs/SCRIPTING.md`: scripting guide; `docs/API.md`: API reference (generated:
  `python tools/gen_api_docs.py`)
- `docs/SPEC.md`: requirements and the phase prompts
- `docs/ARCHITECTURE.md`: design decisions and licensing
- `docs/ROADMAP.md`: phases and acceptance criteria
- `docs/PROGRESS.md`: what has been done
- `docs/OPEN_ISSUES.md`: limitations and caveats
