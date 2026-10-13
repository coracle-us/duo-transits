# EPIC 202126877 RNAAS draft: provenance of every number

**Draft:** `eb_depth_evolution/rnaas/rnaas.tex`. NOT submitted.
**Figure:** `depth_history_summary.png`, copied from `eb_depth_evolution/figures/` (`scripts/25_summary_figure.py`).
**Word count:** abstract 137, body 734 (`rnaas/wordcount.py`, 2026-10-12, after fixing the `\%` comment-regex bug).

## Abstract and depths

| Number | Source |
|---|---|
| P = 21.775 d; linear TESS ephemeris P = 21.774799 d, T0 = BTJD 2488.58926 | `REPORT.txt` §3; `scripts/13_oc_basic.py` |
| 2004–08 primary < 4.5% (3σ worst case), typically < 1.6%, for offsets −0.6…+0.3 d | `scripts/19_ground_combined_limits.py`, `12_ground_eclipse_scan.py`; `REPORT.txt` §2 |
| K2 2014 primary 6.18% (EVEREST/K2SFF); secondary < 0.03% (3σ) | `results_depth_events.csv` (`scripts/20_depth_table.py`); `08_k2_events.py` |
| TESS primary 20.58 ± 0.12% (2021), 27.40 ± 0.16% (2023), 35.64% (2026) | `REPORT.txt` §2 table; `results_depth_events.csv` |
| TESS secondary 4.36 ± 0.14% (2021), 12.52 ± 0.10% (2023); S/P 0.212 → 0.457 | same |
| Depth = median flux within ±15 min of mid-eclipse, own TESScut photometry, r = 1.5 px | `REPORT.txt` §2; `scripts/04_tesscut_photometry.py`, `06_measure_tess_events.py` |
| Neighbour G 14.9 at 21″, ~2.7% of the TESS flux | `scripts/03_gaia_neighbours.py`; `REPORT.txt` §2 |
| 2026 segment BTJD 4055.8–4063.0 | `REPORT.txt` (Sector "1751", BTJD 4055.76–4063.00) |
| DASCH: no excess faint points versus 30 comparison stars | `scripts/22_dasch.py`–`24_dasch_sensitivity.py`; `REPORT.txt` §2 |
| ZTF and ATLAS saturated; ASAS-SN V sparse | `REPORT.txt` §2 and Addenda 1–4 |

## Inclination model

| Number | Source |
|---|---|
| Joint fit to 11 TESS eclipses, reduced χ² 1.43 | `scripts/16_physical_fit.py`, `ebmodel.py`; `REPORT.txt` §4 |
| R1/a = 0.0297, R2/a = 0.0391, e = 0.35, ω ≈ 138° | same |
| i: 2014.4 85.69–85.84; 2021.8 87.24–87.39; 2023.9 87.88–88.10; 2026.0 88.72–89.08 (text: 85.7, 87.3, 88.0, 88.9) | `REPORT.txt` §4; `data/inclination_trend.csv` |
| di/dt 0.19–0.22, 0.31–0.34, 0.39–0.45 deg/yr (text: 0.21, 0.33, 0.42); average 0.25 | `scripts/18_inclination_trend.py` |
| Thresholds: primary i > 84.47°, secondary i > 86.55°, totality i > 89.25° | `REPORT.txt` §4 |
| Onset: primary 2008–2011, secondary 2018–19; totality 2027–29; a turnaround predicts 4–7% in 2004–08 (excluded) | `REPORT.txt` §4 |
| R1 = 1.6, R2 = 2.1 R☉; M_tot ≈ 4.3 M☉ | `REPORT.txt` §4 (sanity check) |

## Timing

| Number | Source |
|---|---|
| O−C rms 10.3 min; χ² 483 for 4 dof | `REPORT.txt` §3; `13_oc_basic.py`, `14_timing_check.py` |
| S−P interval −26.2 ± 2.5 min in one orbit (2021); −104 min from 2021 to 2023 | `scripts/21_timing_models.py` |
| Outer period ~100–600 d (order of magnitude) | `REPORT.txt` §4, third body |
| Gaia RUWE 1.007, no NSS solution, no co-moving companion | `scripts/17_gaia_astrometry.py`, `03_gaia_neighbours.py` |

## Predictions

| Number | Source |
|---|---|
| 2026-10-25 22:39, 11-16 17:12, 12-08 11:46, 12-30 06:21 UTC (geocentric, light-travel corrected) | `scripts/28_predictions.py` → `predictions_2026_27.csv` |

`REPORT.txt` §3 lists the same epochs about 3–8 min later; those were BJD_TDB read as UTC (the same bug as in the K2 note).

## References

All 14 journal references were verified against ADS on 2026-10-12: `scripts/27_verify_refs.py` → `rnaas/refs_check.txt`.

UAT keywords were checked in the official UAT.json:
- Eclipsing binary stars 444;
- Trinary stars 1714 (not "Triple stars");
- Orbital evolution 1178.

## Still pending before any submission

1. The blind independent re-derivation (`eb_depth_evolution/independent/`) must agree.
2. A Fable fact-check of the draft.
3. Optional: an ASAS-SN g-band retry when Sky Patrol V2 is back.
4. The "McBee 2026, submitted" citation refers to AAS82225. Update it when that note is accepted.

## Status decision (2026-10-12, Jon)

**Finish, then HOLD.** Bring the draft to submission-ready (blind check plus Fable fact-check), but do not submit
until AAS82225 (the K2 note, which already mentions this star's deepening) has an editor decision. Then cite it
properly and submit. Reason: avoid the appearance of salami-slicing with an overlapping note weeks later.

## Blind re-derivation (2026-10-12, `eb_depth_evolution/independent/`; agent without access to our analysis)

| Quantity | Ours (REPORT.txt) | Blind (independent/REPORT.txt, results.csv) |
|---|---|---|
| i 2014.36 | 85.69–85.84° | 85.77 ± 0.04 (stat) ± 0.19 (syst) |
| i 2021.8 | 87.24–87.39 | 87.234 ± 0.004 ± 0.21 |
| i 2023.9 | 87.88–88.10 | 87.860 ± 0.004 ± 0.29 |
| i 2026.04 | 88.72–89.08 | 88.549 ± 0.011 ± 0.32 |
| di/dt 2014→21 | 0.19–0.22 | 0.197 ± 0.023 |
| di/dt 2021→23 | 0.31–0.34 | 0.305 ± 0.003 |
| di/dt 2023→26 | 0.39–0.45 | 0.316 ± 0.005 (0.30–0.37) |
| Primary onset threshold | i > 84.47 → 2008–2011 | i > 84.41 (84.33–84.56) → 2007–2010 |
| Totality threshold | 89.25 → 2027–29 | 89.59 (89.20–89.95) → ~2029.5 (2027.5–2031) |
| S−P interval change | −26.2 min/orbit (2021) | −27 (2021), −52 (2023); ~4 min/orbit between seasons |
| K2 primary early | 102 min | 80–100 ± 29 min |
| Model depths (r = 2 px) | median ±15 min (r = 1.5): 20.6 / 27.4 / 35.6 | 20.1 / 27.0 / 34.5; reduction systematics 1.6–2.2% |

**Changes made to the draft:**
- "Accelerating 0.21 → 0.33 → 0.42" became "0.20 (2014–21) to 0.31 deg/yr (2021–26); a constant rate is excluded,
  continued increase not established".
- Onset: 2007–2011.
- Totality: i ≈ 89.2–89.9°, ~2028–2031, 40–45% deep.
- The S−P interval sentence now has both seasons.
- Abstract updated to match.
