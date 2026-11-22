# Gaia DR4 Sweep (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Before Gaia DR4 (2026-12-02), freeze the host list, have a cross-match that is tested on DR3
and switches to DR4 by one argument, and have a pre-registration ready for Jon to approve.

**Architecture:** A new folder `gaia_dr4/`. Pure decision functions in `core.py` (tested offline);
archive access in `xmatch.py` (schema is a parameter: `gaiadr3` now, `gaiadr4` on the day); the frozen
inputs are CSVs built by `build_hosts.py` from tables already in this repo.

**Tech Stack:** Python 3.13 (uv venv), astroquery (Gaia TAP, MAST TIC), pandas, pytest, PyYAML.

**Spec:** `docs/STRATEGY_living_alias_catalogue.md` (§5.4, §6 phase 1, §6a).

## Global Constraints

- Read-only toward the outside world: archive queries only. Nothing is posted or sent.
- Times in BTJD (BJD_TDB − 2457000).
- Gaia `source_id` stays an exact int64 everywhere (a float column corrupts 19-digit IDs).
- The Gaia ESA Archive is offline from 2026-11-25 12:00 CET until 2026-12-02 12:00 CET. The DR3 run
  must finish before then.
- Thresholds live in `prereg_dr4.yaml`, not in code defaults that differ from it.
- No raw downloads committed; tables, code and logs only.

## Review Focus

1. A TIC with no Gaia DR2 ID in the TIC: expect a positional fallback, flagged, not a dropped star.
2. A DR2 ID that maps to two DR3 sources: expect the nearest, with the count recorded.
3. A star with several NSS solutions: expect every solution kept, not the first.
4. A Gaia period so imprecise that any Δt matches: expect "uninformative", never "match".
5. A catalogue row with a missing TIC or epoch: expect it listed in a dropped file with the reason.

---

### Task 1: Decision functions

**Files:** Create `gaia_dr4/core.py`; Test `gaia_dr4/tests/test_core.py`.

**Produces:**
- `rv_variable(p, gof, n) -> bool` (Katz et al. 2023 cut: p ≤ 0.01, gof > 4, n ≥ 10; missing → False)
- `phase_match(dt, dt_err, p, p_err, nsig=3.0, max_window=0.5) -> dict` with keys `n`, `frac`, `dphi`,
  `window`, `verdict` in {`match`, `no_match`, `uninformative`}
- `period_match(p_hyp, p, p_err, nsig=3.0, rel_floor=1e-3, harmonics=(0.5, 1.0, 2.0)) -> dict` with
  keys `verdict` in {`match`, `no_match`}, `harmonic`, `ratio`
- `binarity_flags(row, ruwe_max=1.4) -> dict` with `ruwe_flag`, `nss_flag`, `rv_flag`, `any_flag`

- [ ] Tests: the cases in `test_core.py` (exact values there).
- [ ] Implement; `pytest gaia_dr4/tests/test_core.py` passes.
- [ ] Commit.

### Task 2: Frozen candidate and host lists

**Files:** Create `gaia_dr4/build_hosts.py`; outputs `gaia_dr4/candidates_frozen.csv`,
`gaia_dr4/hosts_frozen.csv`, `gaia_dr4/candidates_dropped.csv`; Test `gaia_dr4/tests/test_build_hosts.py`.

**Produces:** `candidates_frozen.csv` columns `tic, catalogue, cand_id, t1_btjd, t2_btjd, period_d`;
`hosts_frozen.csv` columns `tic, catalogues, n_cand`.

Inputs: `data/catalog/candidates.csv` (hawthorn), `singles/data/ctoi_singles.csv` (ctoi),
`salinas-second-transits/candidates.csv` (salinas), `k2singles/data/k2_singles.csv` (k2),
`toi_singles_k2/targets.csv` (toi_single).

- [ ] Tests: row counts per catalogue equal the source tables minus dropped rows; no duplicate
  `(catalogue, cand_id)`; every `tic` is a positive integer.
- [ ] Implement, run, commit the three CSVs.

### Task 3: Cross-match

**Files:** Create `gaia_dr4/xmatch.py`; outputs `gaia_dr4/results_dr3/{ids,source,nss,report}.csv`;
Test `gaia_dr4/tests/test_xmatch_network.py` (marker `network`).

**Consumes:** `hosts_frozen.csv`, `core.py`.
**Produces:** `python gaia_dr4/xmatch.py --schema gaiadr3 --out gaia_dr4/results_dr3`.

- [ ] TIC → DR2 ID from the MAST TIC; DR2 → DR3 through `dr2_neighbourhood`; positional fallback.
- [ ] `gaia_source` flags and every `nss_two_body_orbit` row per source.
- [ ] `report.csv`: per candidate, flags and the `phase_match` / `period_match` verdict per NSS solution.
- [ ] Network test: TIC 147035057 → 2947647054662200192; TIC 159614322 → 1669046625055549312 with
  RUWE 1.51 ± 0.01; TIC 231969683 and 230130556 each have an NSS orbit.
- [ ] Compare with `singles/data/gaia_binarity_flags.csv` and report differences.
- [ ] Commit.

### Task 4: Pre-registration draft

**Files:** Create `gaia_dr4/prereg_dr4.yaml`, `gaia_dr4/README.md`.

- [ ] Claim types, frozen inputs with SHA-256, DR4 tables to query, thresholds, outcome definitions,
  what is human review. Status `DRAFT - not approved`.
- [ ] Commit; open a PR.
