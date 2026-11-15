#!/usr/bin/env bash
# Build the arXiv source package for AAS82207 from the accepted text (tag aas82207R2-submitted).
# Usage: build_arxiv.sh /path/to/aastex701.cls OUTDIR
# aastex701.cls: AASTeX 7.0.1 (2025/05/09), e.g. from https://mirrors.ctan.org/macros/latex/contrib/aastex.zip
# Test compile: cd OUTDIR/build && tectonic -k --keep-logs rnaas.tex
set -euo pipefail
cls=$1; out=$2
repo=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
grep -q 'ProvidesClass{aastex701}\[2025/05/09 Version 7.0.1' "$cls" || { echo "not AASTeX 7.0.1: $cls" >&2; exit 1; }
mkdir -p "$out/src" "$out/build"
git -C "$repo" show aas82207R2-submitted:singles/rnaas/rnaas.tex > "$out/src/rnaas.tex"
cp "$cls" "$out/src/aastex701.cls"
cp "$out/src/"* "$out/build/"
tar -C "$out/src" -czf "$out/McBee_AAS82207_arXiv.tar.gz" rnaas.tex aastex701.cls
tar -tzf "$out/McBee_AAS82207_arXiv.tar.gz"
