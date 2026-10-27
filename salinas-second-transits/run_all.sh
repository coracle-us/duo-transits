#!/usr/bin/env bash
# One-shot run (Mac, Linux or Git Bash on Windows): makes a venv, installs, runs every step.
# Usage: ./run_all.sh            (all candidates of type single + multiplanet)
#        ./run_all.sh --workers 2
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d .venv ]; then
  python3 -m venv .venv
fi
PY=.venv/bin/python
[ -x "$PY" ] || PY=.venv/Scripts/python.exe   # Windows venv layout
if [ ! -f .venv/.installed ]; then
  "$PY" -m pip install -q --upgrade pip
  "$PY" -m pip install -q -r requirements.txt
  touch .venv/.installed
fi
"$PY" tests/test_synthetic.py
"$PY" st2_pipeline.py all --candidates candidates.csv --work work --out results "$@"
echo "Done. Open results/shortlist.md"
