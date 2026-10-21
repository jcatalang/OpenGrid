#!/usr/bin/env bash
# Every style check of the repository: ruff (lint and format), clang-format (the style in
# tools/clang-format.yaml) and the no-inline-conditional rule. Run from anywhere, in the
# opengrid environment:  tools/lint.sh
set -euo pipefail
cd "$(dirname "$0")/.."
ruff check .
ruff format --check .
clang-format --style=file:tools/clang-format.yaml --dry-run --Werror \
  $(git ls-files '*.cpp' '*.hpp' '*.h')
python tools/check_no_ternary.py
echo "lint: all clean"
