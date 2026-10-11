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
Hemisphere", MNRAS (2024), arXiv:2310.17268 (v2).**

- 85 candidates (25 previously known, 60 new), from SPOC FFI light curves.
- Data used: Cycle 1 (Sectors 1–13) and Cycle 3 (Sectors 27–38) only. Cycle 5 was **not**
  used; the authors explicitly name "search Cycle 5 for third transits" as a way to cut
  aliases. TESS has revisited the south since then (Cycle 5, Cycle 7, possibly later):
  **that newer data is our opportunity.**
- Appendix tables: A1 host properties, A2 transit events (Tc1, Tc2 in TBJD, depths,
  durations), A3 aliases (N_alias, most probable period P_marg; "shown in full online").
- Notable long-period candidates: TIC 294097549 (~200 d), TIC 25194908 (~338 d),
  TIC 349091983 (~549 d).

## 3. Data provenance — IMPORTANT

`data/transcribed/` was produced in a cloud session by a web-reading model, not from the
LaTeX source. Treat it as **unverified**:

- `hosts_tableA1.csv`: all 85 rows. Row 41 distance garbled.
- `transits_tableA2_partial.csv`: rows 1–76 only (77–85 missing).
- Table A3 (aliases) not obtained.

**First task:** run `scripts/fetch_data.py`, which downloads the arXiv LaTeX source into
`data/raw/paper_src/`. Parse Tables A1–A3 from the `.tex`, write
`data/catalog/candidates.csv`, and diff it against the transcribed files. Any disagreement
is resolved in favor of the LaTeX source. Do not use the transcribed files for analysis.

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

## 8. Status at handoff (2026-10-11)

- Topic chosen and literature gap checked (active area; newer TESS sectors unexploited
  for this catalog as far as we found — re-check during triage).
- Fetch script written, syntax-checked, **not yet run**.
- Nothing analyzed yet.
