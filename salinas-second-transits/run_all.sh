#!/usr/bin/env bash
# One-shot run on a Mac: makes a venv, installs, runs every step.
# Usage: ./run_all.sh            (all candidates of type single + multiplanet)
#        ./run_all.sh --workers 2
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d .venv ]; then
  python3 -m venv .venv
  .venv/bin/pip install -q --upgrade pip
  .venv/bin/pip install -q -r requirements.txt
fi
.venv/bin/python tests/test_synthetic.py
.venv/bin/python st2_pipeline.py all --candidates candidates.csv --work work --out results "$@"
echo "Done. Open results/shortlist.md"
