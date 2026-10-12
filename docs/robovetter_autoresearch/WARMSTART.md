# Warm start: autoresearch on the Kepler DR25 Robovetter

Copy this file into the new repository as `CLAUDE.md`. It is the brief for a fresh Claude Code session. Written
2026-10-12 by the session that ran the `duo-transits` project (github.com/coracle-us/duo-transits).

## 1. Who and why

- **Who.** Jon McBee: independent researcher, hobby astronomer, ORCID 0009-0000-7524-4470, Windows 11, an
  NVIDIA RTX A4500 (20 GB) and 32 CPU cores. He has two Research Notes of the AAS under review (AAS82207,
  AAS82225) from `duo-transits`. He wants blunt go/no-go calls, and he presses submit/send buttons himself.
- **Goal.** Run a Karpathy-style **autoresearch** loop (github.com/karpathy/autoresearch) on a problem where an
  agent can make many small, *verifiable* improvements overnight. Inspiration:
  github.com/CrocSwap/integer-mult-bounds, which has an exact scalar metric, cheap trustworthy checking, and
  stacking wins.
- **Chosen problem** (scoping: `duo-transits/docs/autoresearch_scoping.txt`): beat the **Kepler DR25 Robovetter**
  at separating real planets from false alarms in the hard regime of **long periods and low signal-to-noise**.
  - This regime drives the uncertainty in Earth-analogue occurrence rates.
  - The labels are *exact* because they are simulated: we know which signals were injected planets.
  - The baseline is frozen and public.
  - Experiments take seconds, not hours.
- **Not chosen, and why.** Headline ExoMiner/Astronet metrics (ROC-AUC ~0.998) are saturated and limited by label
  noise; a loop would learn the noise.
- **Phase 2, later:** a TESS single-event benchmark (see the scoping file, option 2).

## 2. Data (verified 2026-10-12)

Source: **github.com/nasa/kepler-robovetter** (last push 2019). It contains:

- `Kepler-RoboVetter-DR25.cpp`, the reference Robovetter (C++11). Compile it: `g++ -std=c++11 -O2 -o robovet
  Kepler-RoboVetter-DR25.cpp`. The repo includes Windows instructions.
- One `kplr_dr25_<set>_robovetter_input.tar.gz` and one `..._output.txt` per data set. The input files hold **every
  metric the Robovetter uses**, so the problem is a tabular classifier, not a light-curve one. Sets:

| Set | What it is | Role |
|---|---|---|
| `obs` | real Kepler Q1–Q17 DR25 TCEs | real data, no ground-truth labels |
| `inj1` | on-target injected planets | **positives** |
| `inj2` | off-target injected planets (centroid offset) | astrophysical false positives |
| `inj3` | injected eclipsing binaries | astrophysical false positives |
| `inv` | inverted light curves | false alarms (systematics) |
| `scr1`, `scr2`, `scr3` | scrambled-cadence light curves | false alarms |
| `sup` | supplemental TCEs | read the docs before using |

The column definitions and the meaning of every test are in **KSCI-19114** (Coughlin 2017,
https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19114-001.pdf) and the DR25 catalogue paper (Thompson et al.
2018, "Planetary Candidates Observed by Kepler. VIII", DOI 10.3847/1538-4365/aab4f9, arXiv:1710.06758). Both were
checked on 2026-10-12: the PDF URL resolves and the arXiv metadata matches. Confirm the journal volume and page in
ADS before citing. Read KSCI-19114 before designing features.

Each `*_output.txt` gives the Robovetter's own disposition (PC/FP), flags and score: this is the **baseline**.
Reproduce at least one output file exactly with the compiled code as a sanity check (README: `./robovet INPUT
OUTPUT 10000`).

## 3. Harness design (build this first; it is the whole game)

The agent may edit **only `train.py`**. Everything else is frozen and hash-checked at the start of every run.

- **`prepare.py`** (frozen; run once):
  - Download, unpack and parse all sets into one table (`data/tces.parquet`) with `set`, `kic`, `tce_id`,
    `period`, `mes`, the Robovetter metrics, and the Robovetter baseline columns (disposition, score).
  - **Split by star (KIC), never by TCE,** using a fixed hash of the KIC ID: `train` 60%, `val` 20%, `locked` 20%.
  - **Hold out `scr3` completely,** to test generalisation to an unseen false-alarm realisation. `obs` is never used
    for training.
  - Write `data/split_manifest.json` with SHA-256 hashes of every file.
- **`evaluate.py`** (frozen): given a model's scores on `val` (or `locked`, human-run only):
  - **Regime:** 100 ≤ P ≤ 500 d and MES ≤ 10. Also report P 200–500 d and the full range.
  - **Primary metric: completeness at matched effectiveness.**
    - Choose the score threshold so that **effectiveness** (fraction of `inv`+`scr1`+`scr2` TCEs in the regime
      dispositioned FP) equals the Robovetter's effectiveness on the same rows.
    - Then report **completeness**, the fraction of `inj1` dispositioned PC.
    - Higher is better. This puts the comparison at the Robovetter's own operating point.
  - Secondary metrics, reported but never optimised:
    - completeness at effectiveness 0.99 and 0.995;
    - the rejection rate of astrophysical false positives (`inj2`, `inj3`);
    - an estimated reliability on `obs` with the Thompson et al. 2018 formula (from the obs PC count and the
      effectiveness; check the exact equation in the paper);
    - the PC fraction on `obs`.
  - **Uncertainty:** a bootstrap over stars (KICs), 1,000 resamples, gives σ for the primary metric.
  - Output one line of JSON plus a row appended to `results.tsv` (commit hash, metric, σ, wall time, short
    description).
- **`train.py`** (the only file the agent edits): reads the train split and writes scores for `val`. It may use any
  features from the input metrics, any model, and any training procedure. Fixed **5-minute** wall-clock budget;
  record seeds. The GPU is optional (gradient boosting on CPU is the obvious start).
- **`program.md`** (Jon iterates on this; the agent reads it): the autoresearch instructions. Edit `train.py`, run
  `evaluate.py` on `val`, and **keep a change only if the primary metric improves by > 2σ (bootstrap)**; otherwise
  revert. Log every attempt, including failures. Never read or score `locked`. Never edit the frozen files.

## 4. Known traps (address them in the harness, not after the fact)

- **Set-identity leakage.** Inverted and scrambled light curves differ systematically from real ones, so a model can
  learn "is this simulated?" instead of "is this a planet?".
  - Run **adversarial validation** early: train a classifier to tell `obs` from `inv`/`scr` using only the
    features.
  - Features that make that task easy are suspect. Consider dropping them, or report results with and without
    them.
  - The `scr3` hold-out and the obs-based reliability are the checks.
- **The Robovetter was hand-tuned on these same sets.** Its numbers are in-sample, and so are ours on `val`. That is
  fair, but say so. Headroom may be smaller than it looks.
- **Regime size.** The 100–500 d, MES ≤ 10 regime may contain only hundreds of `inj1` TCEs. Check the counts first,
  and make sure the bootstrap σ isn't larger than any plausible gain. Widen the regime if needed, and document it.
- **Multiple comparisons.** Overnight runs try about 100 changes, so a few will pass 2σ by luck. The **locked**
  split, scored once by a human at the end, is the only number that goes in a paper.
- **Don't rediscover published work.** Check prior ML on DR25 Robovetter metrics before claiming novelty (e.g.
  Armstrong et al. 2021 on Kepler validation; ExoMiner trained on Kepler). Verify every citation in ADS/arXiv.

## 5. Running the loop

- Option A (closest to Karpathy, cheaper): a separate headless Claude Code session that follows `program.md` for N
  iterations. Restrict its permissions so it can **edit only `train.py`** and run only `python train.py` /
  `python evaluate.py`. Do NOT disable all permissions, despite the autoresearch README.
- Option B: a multi-agent Workflow inside an interactive session. It costs more tokens; Jon must explicitly opt in.
- Start with a **short pilot** (about 10 iterations) to confirm the harness, timing and logging before an overnight
  run.

## 6. Working rules carried over from duo-transits

- **Verification over volume.** Every claimed number traces to a file and a script in the repo. Before writing
  anything up, get an independent check: a blind re-run, or a separate reviewer agent.
- Prefer "no significant gain" (publishable as a null) over a shaky win.
- Keep a `numbers.md` provenance file from day one.
- **Windows notes:**
  - Git Bash heredocs mangle backslashes, so write Python and LaTeX with the file tools.
  - Use the venv's Python explicitly.
  - Large data goes in `data/` and is gitignored.
- **Publishing.** RNAAS is the likely venue. It has a 1,350-word limit, one table *or* one figure, and AASTeX v7. Add
  an AI disclosure statement at the end (see journals.aas.org/author-llm-guidelines). Jon is the author; Claude is
  not.
- Never put private correspondence in the public repo.

## 7. First-session checklist

1. Create the repo layout: `prepare.py`, `evaluate.py`, `train.py`, `program.md`, `numbers.md`, `data/`
   (gitignored), and `.venv` with numpy, pandas, scikit-learn, lightgbm or xgboost, pyarrow.
2. Download and parse the data. Print the per-set TCE counts and the counts in the regime.
3. Reproduce one Robovetter output file with the compiled C++ code.
4. Compute the **baseline** (Robovetter) primary metric on `val` with bootstrap σ.
5. Run adversarial validation and report the leakage risk.
6. Write `train.py` v0 (a plain gradient-boosted model on all metrics) and score it.
7. Stop and report to Jon: the counts, the baseline, v0, σ, the leakage findings, and a go/no-go for the overnight
   loop.
