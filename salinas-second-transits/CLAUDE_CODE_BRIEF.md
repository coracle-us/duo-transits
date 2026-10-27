# Brief for Claude Code on Jon's Mac: Salinas second-transit search

Paste this into a Claude Code session opened in a clone of coracle-us/duo-transits.

## Goal
Find a second TESS transit (Sector 61 or later) for the new single-transit candidates in
Salinas et al., arXiv 2610.24882. Each confirmed hit becomes a new duotransit whose period
aliases Jon resolves with his usual method, and a set of them could make a Research Note
of the AAS (RNAAS).

## Starting point
The pipeline is in `salinas-second-transits/` (unzip it into the repo root). Read its
README.md first; it explains every step, the files and the pass criteria. It was built in
a cloud session that could not reach MAST or ExoFOP, so it has only been tested on
synthetic light curves (`tests/test_synthetic.py`). Expect the first real run to hit
problems in the MAST download, TIC query, ExoFOP column handling or the CNN step.

## Steps
1. Put `salinas-second-transits/` in the repo as its own folder, add a .gitignore for
   `.venv/`, `work/` caches of raw light curves and other large downloads, and commit.
2. Run `./run_all.sh` (or the steps one by one as the README shows). Use `--workers 1`
   if MAST throttles.
3. When a step errors on real data, fix the code, rerun from that step, and commit each
   fix with a clear message. Keep the synthetic test passing.
4. Sanity check: the report should flag TIC 302377820 (transit in S68 matching
   candidate b) and TIC 404736082 (more transits in S70), which the paper itself
   mentions. If either is missing, find out why before trusting the rest.
5. Optional: run the CNN step with helemysm/singletransiters (needs its own Python 3.10
   env, numpy<2, torch). That repo has no license file: run and cite it, but don't copy
   its code into duo-transits.
6. Commit `results/shortlist.csv`, `results/shortlist.md`, the per-hit plots,
   `work/events_all.csv` and `work/aliases_summary.csv` (not raw light curves), then push.

## Report back
Finish with a short summary Jon can paste into his astrophysics channel: how many
candidates had Sector 61+ data, how many passing hits, whether the two sanity-check
systems came out, any hits with a `supported` alias (a third transit), what you had to
fix, and the commit hash. Don't vet or rank hits beyond what the pipeline reports;
that happens in the channel.

## Limits
- Read-only use of MAST, ExoFOP and GitHub data; don't submit anything to ExoFOP or
  anywhere else.
- Don't push to any repo other than coracle-us/duo-transits, and don't force-push.
- Install only packages in requirements.txt (plus the CNN repo's own requirements if
  running step 5).
