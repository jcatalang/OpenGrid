# OpenGrid

Interactive block-structured / unstructured / hybrid mesh generator for CFD (2D and 3D). C++20 core, Python bindings (nanobind),
PySide6 + VTK GUI, and a Python scripting API.

Status: Phase 0 (scaffold). See `docs/ROADMAP.md`.

## Setup (conda-forge, primary)
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
- `docs/SPEC.md`: requirements and the phase prompts
- `docs/ARCHITECTURE.md`: design decisions and licensing
- `docs/ROADMAP.md`: phases and acceptance criteria
- `docs/PROGRESS.md`: what has been done
