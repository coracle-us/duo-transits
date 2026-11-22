# Gaia sweep of long-period TESS candidate hosts

Phase 1 of `docs/STRATEGY_living_alias_catalogue.md`: have everything ready before Gaia DR4 opens on
2026-12-02 12:00 CET. Plan: `docs/superpowers/plans/2026-11-22-gaia-dr4-sweep.md`.

| File | What |
|---|---|
| `build_hosts.py` | Builds the frozen lists from tables already in the repo. |
| `candidates_frozen.csv` | One row per candidate, plus the periods this project derived. |
| `hosts_frozen.csv` | One row per star. |
| `core.py` | The decision functions (flags, phase test, period test). No network. |
| `xmatch.py` | TIC to Gaia source, `gaia_source` flags, every NSS orbit, and the report. |
| `prereg_dr4.yaml` | Pre-registration. **Draft; not approved.** Thresholds are read from it. |
| `results_dr3/` | The DR3 baseline run. |
| `t16_overlap.py`, `data_t16/` | T16 Planet Hunt singles: source table, overlap check, ID disagreements. |
| `tests/` | `pytest gaia_dr4/tests` (offline); `pytest -m network gaia_dr4/tests` checks the DR3 run. |

## Run

```bash
uv venv --python 3.13 .venv && uv pip install --python .venv/bin/python astroquery pandas pytest pyyaml
.venv/bin/python gaia_dr4/build_hosts.py
.venv/bin/python gaia_dr4/xmatch.py --schema gaiadr3 --out gaia_dr4/results_dr3
```

On release day, after the open items in `prereg_dr4.yaml` are settled:
`xmatch.py --schema gaiadr4 --out gaia_dr4/results_dr4` (copy `ids.csv` across first if DR3 source IDs
carry over).

## The two tests

- **Phase test** (duotransits): is the gap between the two transits a whole number of Gaia orbital
  periods? `window` is the fraction of all phases the test would accept, so it is also the chance an
  unrelated orbit passes. Above 0.5 the test is called uninformative.
- **Period test** (a stated period): does the Gaia period agree within 3σ, or 0.1% if larger, directly
  or at half or twice the period?

## DR3 baseline (run 2026-11-22, `results_dr3/`)

1368 candidate rows on 1301 stars (rerun the same day after the T16 singles were added). 1300 stars have a
DR3 source: 1297 through the TIC's DR2 ID and `dr2_neighbourhood`, 3 by position (TIC 3181884, 197736350,
456862427; check these by hand). TIC 252803603 (Tmag 6.9) has no Gaia source within 5.4 arcsec. 35 DR2 IDs
map to more than one DR3 source; the nearest is used.

| | Stars |
|---|---|
| Any flag | 277 |
| RUWE > 1.4 | 202 |
| `non_single_star` > 0 | 78 |
| RV-variable | 141 |
| With an NSS orbit | 55 (59 solutions) |

Flagged stars per catalogue (first catalogue listed for the star): ctoi 104 of 354, hawthorn 21 of 85,
k2 43 of 163, salinas 34 of 162, t16 48 of 395, toi_single 27 of 142.

Nothing here is new. Every orbit verdict below repeats what the repo already records, which is the point
of the baseline:

| TIC | Catalogue | Gaia orbit | Verdict | Already known (source) |
|---|---|---|---|---|
| 32179255 | hawthorn | SB1 72.13 d | gap = 10 periods; window 0.09 | `data/catalog/status.csv`: possible binary |
| 100776118 | hawthorn | AstroSpectroSB1 251.3 d | gap = 3 periods; window 0.06 | `status.csv`: possible binary |
| 70561926 | hawthorn | SB1 21.535 d | gap = 35 periods; window 0.24 | `status.csv`: false positive (KELT, 21.4926 d) |
| 77437543 | hawthorn | SB1 4.44 d | gap = 163 periods; window 0.16 | **A false match.** The real period is 60.33 d (TOI-2490 b, Henderson et al. 2024). |
| 157119927 | hawthorn | SB1 21.865 d | no match, by a hair (frac −0.0138 against a limit of 0.0130) | `status.csv`: false positive, same Gaia orbit, "~3 sigma" |
| 207783865 | hawthorn | SB1 7.545 d | no match | `status.csv`: not an alias |
| 230130556, 231969683 | own_ctoi | SB1 22.245 d, 21.969 d | period match | `singles/scripts/gaia_flags.py` |
| 51238317 | ctoi | SB1 14.2495 d | matches the listed 14.261 d | — |
| 264333501 | own_ctoi | SB1 15.770 d | matches at half our 31.5515 d | classified EB-on-target (Gaia) |

Two lessons for the pre-registration:
- TIC 77437543 shows the phase test passing by chance at a window of 0.16. The draft sends any match with
  a window above 0.1 to human review, which catches it.
- TIC 157119927 fails by less than the assumed timing error could explain. `dt_err_d` is one assumed
  value for every candidate; per-candidate timing errors would be better.

Compared with `singles/data/gaia_binarity_flags.csv` (359 stars, matched by position there): RUWE agrees
to 0.01 for all but 2 stars, one star resolves to a different source, and 2 more stars are flagged here
(TIC 147855066, TIC 130370378).

## T16 Planet Hunt singles (added 2026-11-22)

Roth et al. (arXiv 2604.18579, accepted in ApJS) publish 411 single-transit candidates from Cycle 1 FFIs
of stars to T = 16, with no period constraints: "even the transit times may be inaccurate". The table is
`data_t16/singleTransit_candidate_table.tab`, from doi:10.7910/DVN/CWEUGW (CC0; file ID 13605256, sha256
d357f60a…0a87a1). Depth is in magnitudes, duration in days, Tc in BJD − 2457000.

`t16_overlap.py` checks each star against ExoFOP's TOI and CTOI tables and pscomppars as downloaded on
2026-11-22, the Prša et al. 2022 TESS EB catalogue, and our lists (`data_t16/t16_overlap.csv`):

| Check | Stars of 411 |
|---|---|
| Has a TOI | 10 (7 PC, 1 APC, 1 CP, 1 FP) |
| Confirmed-planet host | 2 |
| Has a CTOI | 13 (hawthorn 5, salinas 3, six others 1 each); 4 carry a period |
| TESS EB (Prša et al. 2022) | 0 (the T16 authors had already removed those) |
| Already in our lists | 8 (5 are Hawthorn duotransits) |
| Transit time within 0.1 d of another T16 single on the same sector, camera and CCD | 72 |
| Deeper than 0.05 mag | 42 |

**Kept: 401.** Dropped, with the reason in `candidates_dropped.csv`: 8 with a TOI and 2 confirmed-planet
hosts. The 9 kept stars with someone else's CTOI stay in, flagged in the overlap table.

Things to know before using these:
- **72 kept stars share their transit time with another T16 single on the same CCD.** That pattern points
  to a systematic, not a transit. They are flagged (`n_coincident`), not dropped.
- **For 11 stars the Gaia DR3 ID that T16 gives differs from the one we get from the TIC**
  (`data_t16/t16_gaia_id_disagreements.csv`). The two are neighbouring sources. T16 made its light curve
  for its own source, so on these stars the TIC number may point at the neighbour. Not resolved here.
- The hosts are faint (median G 14.4; 336 fainter than G 13), so Gaia RVs exist for about half and the
  RUWE and NSS flags carry most of the weight.
- In DR3, 48 of the 395 T16-only stars carry a flag and 5 have an SB1 orbit (TIC 99839636, 317096831,
  332483493, 340562993, 1695417; periods 7.9 to 70 d). With one transit each, no period test is possible.
- No search of the literature beyond those tables was made for these stars.

## Limits

- The Gaia ESA Archive is offline from 2026-11-25 12:00 CET until the release.
- The TSTPC and Planet Hunters lists are not in the host list.
- Transit-time errors are not carried per candidate; one assumed value (`dt_err_d`) is used.
- Single-transit candidates with no stated period get flags only. No period test is possible.
- A flag (RUWE, NSS, RV variability) says the star may have a companion. It does not say the transit
  signal comes from that companion.
