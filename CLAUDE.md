# Duotransit Period Recovery — project briefing

Hobby research project (Jon McBee + Claude). Goal: resolve or narrow the orbital periods of
TESS **duotransit** exoplanet candidates using TESS data taken *after* the discovery paper,
and publish the results as a **Research Note of the AAS (RNAAS)**.

Read this file first. It is the source of truth for scope, method and rules.

---

## 1. The science in one paragraph

TESS observes a patch of sky for ~27 days, then returns roughly two years later. A planet
with a long period may transit once in each visit: two dips separated by a large gap Δt.
The true period must be **P_n = Δt / n** for some integer n ≥ 1. Most n are excluded
because TESS was observing when that period would have produced another transit and saw
none; the survivors are the **aliases**. One more transit collapses the list to a single
period. Resolved periods turn candidates into followable targets for RV mass measurement
and characterization. Most of these candidates are warm Jupiters (P ≳ 20 d), which keep
the orbital fingerprints of how giant planets migrate, which is why they matter.

## 2. The source catalog

**Hawthorn, Gill, Bayliss et al., "TESS Duotransit Candidates from the Southern Ecliptic
Hemisphere", MNRAS 528, 1841 (2024), arXiv:2310.17268 (v2).** ADS: 2024MNRAS.528.1841H,
DOI 10.1093/mnras/stad3783.

- 85 candidates (25 previously known, 60 new), from SPOC FFI light curves.
- Data used: Cycle 1 (Sectors 1–13) and Cycle 3 (Sectors 27–39; one Tc2 is in S42) only. Cycle 5 was **not**
  used; the authors explicitly name "search Cycle 5 for third transits" as a way to cut
  aliases. TESS has revisited the south since then (Cycle 5, Cycle 7, possibly later):
  **that newer data is our opportunity.**
- Appendix tables: A1 host properties, A2 transit events (Tc1, Tc2 in TBJD, depths,
  durations), A3 aliases (N_alias, most probable period P_marg; "shown in full online").
  In the LaTeX these are `table_1.tex`, `table_4.tex`, `table_5.tex`.
- N_alias behaves as a *count* of allowed aliases (despite the caption); for 71/85 rows it
  equals n_max+1 where P_min = dT/n_max, which looks like an off-by-one in their counting.
- Notable long-period candidates: TIC 294097549 (~200 d), TIC 25194908 (~338 d),
  TIC 349091983 (~549 d).

## 3. Data provenance — IMPORTANT

`data/catalog/candidates.csv` (all 85 rows, Tables A1–A3) is parsed from the arXiv LaTeX
source by `scripts/parse_catalog.py`. **Use it, not `data/transcribed/`.** The web
transcription agreed with the LaTeX everywhere it had data; differences (distance strings
in A1 rows 33/36/38/41/53, A2 rows 77–85 missing) are in `data/catalog/transcription_diff.csv`.

## 4. Plan

1. **Get data** — `python scripts/fetch_data.py` (paper source, NASA Exoplanet Archive
   `pscomppars`, ExoFOP TOI/CTOI tables, all TESS light curves for the 85 TICs via
   lightkurve: SPOC 2-min/20-s, TESS-SPOC FFI, QLP, TGLC). Resumable; logs to
   `data/raw/fetch_log.txt`. Windows machine; the script is cross-platform.
2. **Triage: which candidates are still open?** For each TIC, check confirmed planets
   (pscomppars by `tic_id`), TOI dispositions (CP/KP/FP/FA), CTOI notes, and literature,
   especially NGTS/Warwick follow-up papers (e.g. NGTS-30 b, NGTS-36 b/TOI-2449 b,
   NGTS-38 b, and the Ulmer-Moll et al. 2026 review arXiv:2604.09254). Record status +
   source per candidate in `data/catalog/status.csv`. Already-solved candidates are
   excluded from claims (they can serve as validation cases; see step 4).
3. **Recompute aliases ourselves** from Tc1, Tc2 and the actual time coverage of *all*
   TESS data (old + new), rather than trusting A3. For each n: predict transit times
   T = Tc1 + k·P_n over every observed baseline; mark the alias *excluded* if a predicted
   full transit falls in good-quality data with no dip at the expected depth/duration,
   *supported* if a dip appears, *unconstrained* if it falls in gaps. Partial transits
   (ingress/egress only) count, with care.
4. **Validate the pipeline** before trusting it: run it on candidates that already have
   published periods; it must recover their true alias and exclude the rest. Also confirm
   our Cycle 1/3-only alias lists match the paper's N_alias.
5. **Search new data for third transits**: per open candidate, per surviving alias, fit a
   transit model at the predicted time (depth and duration fixed or tightly bounded from
   Tc1/Tc2). Report Δχ² / BIC for transit vs. no-transit.
6. **Vet any detection**: depth and duration consistent with the first two; odd/even and
   secondary-eclipse checks (eclipsing binaries); not coincident with momentum dumps,
   sector edges or scattered light; check SAP vs PDCSAP; consider TGLC/QLP as an
   independent reduction; nearby-star contamination (TIC/Gaia neighbours).
7. **Results table**: TIC | status | aliases before → after | resolved P (± σ) or
   surviving alias list | data used | notes.
8. **Write the RNAAS note** (see §6), with one table.

## 5. Rules of the road

- **Verification over volume.** A claimed period must survive step 6 and be independently
  re-derived (fresh script or a separate agent that hasn't seen the first analysis) before
  it goes in the note. Prefer "alias narrowed to {a, b}" over a shaky "resolved".
- Time system: TESS BTJD/TBJD = BJD_TDB − 2457000. Keep everything in BTJD internally.
- Never edit `data/raw/`; derived products go in `data/catalog/` or `results/`.
- Every number in the note must trace to a file and a script in this repo.
- Check for duplication of others' work before claiming novelty: if a period is already
  published, cite it, don't claim it.
- Null results are publishable (RNAAS accepts them); don't over-interpret marginal dips.

## 6. Publication target: Research Notes of the AAS

- Moderated, not peer reviewed (Lead Editor Chris Lintott); fast publication; indexed in
  ADS, citable. No author publication charges (as of Oct 2026).
- ≤ 1,500 words, abstract required, **one figure or one table, not both**.
- AASTeX (v7+) preferred; Overleaf/Word accepted.
- Author: Jon McBee (independent; consider listing Coracle as affiliation). Claude is
  not an author: add an acknowledgment disclosing AI assistance, and **check current AAS
  journals policy on AI use before submitting.**
- Natural follow-ons: report solved candidates to ExoFOP/TFOP; if results are strong, grow
  into a peer-reviewed AJ/MNRAS paper.

## 7. Useful tools/references

- `lightkurve` (data access, light curves), `astroquery` (archives), `batman` or
  `exoplanet`/`juliet` (transit models), `tess-point` (sector coverage predictions).
- MonoTools (H. Osborn) — the duo/monotransit alias-probability code used in the field.
- Hawthorn et al. 2024 for method details (their alias definitions and vetting flags).

## 8. Status (updated 2026-10-11)

Done (steps 1–2 + catalog + coverage):
- `.venv` with requirements; `scripts/fetch_data.py` fixed and run: paper source, pscomppars,
  TOI/CTOI tables, 1,354 light curves for all 85 TICs (`data/raw/`, gitignored, ~750 MB).
- `data/catalog/candidates.csv` from LaTeX (see §3).
- Triage → `data/catalog/status.csv` (`scripts/build_status.py`, sources per row):
  49 open, 26 resolved, 7 unclear, 3 false-positive. Resolved tiers: 9 peer-reviewed,
  1 catalog paper (549 d), 5 Hawthorn PhD thesis (Warwick 2024, NGTS-solved), 8 ExoFOP/TFOP
  community, 3 TOI periods from QLP faint-star search.
- Coverage beyond the paper's data → `data/catalog/coverage.csv` (`scripts/coverage_summary.py`):
  78/85 have new sectors (up to S107); 7 have none.

Lessons that change the plan:
- **The novelty gap is narrower than assumed.** ExoFOP community users (2025–26) already used
  S87–S103 to pin periods for ≥8 candidates, and Hawthorn's thesis solved 5 more with NGTS.
  Always check per-target ExoFOP uploads/notes (`scripts/fetch_exofop.py`), not just the bulk
  TOI table, which can be stale (e.g. still lists max periods for solved TOIs).
- Those community periods are our **validation set** (step 4); independent confirmation is
  reportable, but cite the upload and don't claim discovery.
- 104 of 276 new target-sectors are QLP-only (TESS-SPOC stops at S72, TGLC has no recent
  sectors) → independent-reduction vetting needs our own TESSCut photometry.
- Sectors 97 and 98 are ~55 d long (genuine). MAST has a mislabelled "Sector 1751" product.
- Open conflict: TIC 13072758 thesis 23.794 d (n=31) vs ExoFOP 49.17 d (n=15).

Later the same day (steps 3–6 done, project wound down):
- Alias pipeline built and validated: `results/validation_notes.md`.
- Vetting of the 9 single-period leads: `results/vetting/vetting_log.csv`.
  - 3 are off-target, likely nearby eclipsing binaries (a 4th, TIC 265466589, was a
    reference-point artefact in diff_image.py, corrected; see results/vetting/vetting_log.csv).
  - 1 is source-ambiguous (TOI-6559).
  - 3 are on-target periods confirmed by an independent BLS: 207783865, 38138512, 153838604.
- **Verdict: dead end as a research direction.** The catalogue is picked over and the yield
  is too small. Results were emailed to the catalogue authors instead of an RNAAS:
  `results/outreach/`.
- Jon cannot get an ExoFOP account (institutional affiliation required). Outputs go via RNAAS
  or to the groups that own the candidates.
- Side quest: flares on the 2-min hosts, `sidequests/flares/` (1 flare; a sentence at most).

**Next problem set:** see `docs/problem_set_survey.md`. Current recommendation: stale
single-transit CTOIs plus difference-image localisation. A pilot is planned but not started.
