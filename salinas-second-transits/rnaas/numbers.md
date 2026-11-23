# RNAAS draft (`rnaas.tex`): where every number comes from, and what's left before submitting

**Status: draft, not submitted.** The CLAUDE.md rule is that every number in the note must trace to a
file and a script. Paths below are relative to `salinas-second-transits/`.

## Sample and search

| Quantity in the note | Value | File | Script |
|---|---|---|---|
| Single-transit candidates and multi-planet systems searched | 108 and 4 (8 events); 116 events on 112 stars | `candidates.csv` (`type` single / multiplanet) | `parse_salinas_tables.py` |
| Sectors in the paper | 27–60 | arXiv 2610.24882v1, §6 | — |
| Stars with light curves from S61 on; last sector | 101; S107 | `work/coverage.csv` (`n_new_lc`, `lc_sectors`) | `st2_pipeline.py coverage` |
| QLP clock offset, S100+ | up to 423 s (≈7 min) | `results/round2/SUMMARY_2026-10-27.md`, addendum "all S100+ times"; `results/round2/qlp_s94_timing.csv` | `round2/s6_qlp_times.py` |
| Search thresholds | SNR ≥ 7; 0.5–2× depth and duration; 1.5× the strongest brightening | `README.md`, step 2 | `st2_pipeline.py search` |
| Sector-wide dip windows masked | 4, in S94–S97 | `st2core.COMMON_MODE_WINDOWS`; `results/round2/common_mode_scan.csv` | `round2/s9_common_mode.py` |
| Passing events; candidates with one | 191; 53 | `work/events_all.csv` (`passes`); `results/shortlist.csv` (53 rows) | `st2_pipeline.py search`, `report` |
| Candidates with 2+ distinct events; periodic | 36; 8 | `work/common_period.csv`; `results/shortlist.csv` (`cp_flag`) | `st2_pipeline.py aliases` |
| TOI-7717.01 (317132251) and TOI-7743.01 (454039549) | 25.21 d, 10.97 d | `results/shortlist.csv` (`toi`, `toi_period`, `cp_period_d`) | `st2_pipeline.py toi` |
| 141423461: no clean event, 2× duration | `n_clean_events` 0, `dur_ratio` 2.0 | `results/shortlist.csv` | same |
| 198009543: 11.71 d against a 42.8 d circular minimum; 9.7 h | `cp_period_d`, `pmin_circ_d`, `dur_ref_h` | `results/shortlist.csv` | same |
| Four stars: not a TOI or CTOI, no ExoFOP follow-up, no literature | | `results/round3/NOVELTY_2026-11-22.md`, `novelty_2026-11-22.csv` | `round3/s13_novelty.py` |

## Vetting (round 3 and the blind run)

| Quantity | Value | File | Script |
|---|---|---|---|
| P, T0, depth, T14, N_tr (table) | see table | `results/round3/oddeven_summary.csv` (`P_d`, `T0_btjd`, `depth_ppm`, `t14_h`, `n_transits_timed`) | `round3/s11_oddeven.py` |
| Reduced χ² of the timing fits; no point > 3σ | 0.30–3.23; `n_oc_over_3sigma` = 0 | same (`timing_chi2_red`) | same |
| Shorter aliases excluded | ≥ 15σ (lowest: 107113048 at P/5, 15.3σ) | `results/round3/independent/out/aliases_*.txt`; `RESULTS.md` table | `independent/scripts/s07_aliases.py` |
| Odd/even | ≤ 1.3σ (round 3); 2.2σ blind for 107113048 (2 + 2 transits) | `oddeven_summary.csv` (`odd_even_sigma`); `independent/COMPARISON.md` | `s11_oddeven.py`; `independent/scripts/s09_checks.py` |
| Phase-0.5 3σ limits | 44 / 479 / 705 / 50 ppm | `oddeven_summary.csv` (`sec_phase05_upper3s_ppm`) | `s11_oddeven.py` |
| TESS EB catalogue and VSX within 30″ | none | `oddeven_summary.csv` (`tess_eb`, `vsx`) | same |
| Difference-image offsets | 2.2–3.7″ (round 3); 1.9–3.9″ for 230133616 (round 2 §5); within ~2″ / ~1″ (blind) | `results/round3/PIXELS_2026-11-22.md`; `results/round2/SUMMARY_2026-10-27.md` §5; `independent/RESULTS.md` | `round3/s12_pixels.py`; `independent/scripts/s13_diffimg.py` |
| Blind re-derivation within 1.6σ | | `results/round3/independent/COMPARISON.md` | `independent/scripts/` |

## Per star

| Quantity | Value | File | Script |
|---|---|---|---|
| 241550976: R★, companion size | 0.96 R☉; ≈1.2 R_J | `independent/out/radius_241550976.txt`; `oddeven_summary.csv` (`rp_rj` 1.19) | `independent/scripts/s14_radius.py`; `s11_oddeven.py` |
| 107113048: R★, size | 1.0 R☉; ≈0.9 R_J | `independent/out/radius_107113048.txt`; `rp_rj` 0.91 | same |
| 107113048: EB neighbour 97″, eclipse 0.68 h from E29, 7–18%, O−C −1.1 min | | `independent/out/eb_neighbour_107113048.txt` | `independent/scripts/s16_eb_neighbour_107113048.py` |
| TIC 107116357 not in Gaia DR3 EBs, TESS EB catalogue, VSX | | `results/round3/NOVELTY_2026-11-22.md` | `round3/s14_eb_neighbour_catalogs.py` |
| 107113048: E29 alone rules out 74.55 and 149.10 d | | `independent/RESULTS.md`, 107113048 follow-up table | blind analysis |
| 107113048: next transit | 2026-12-25 21:44 UTC ± 11 min | `results/round3/predictions_2026-11-22.csv` | `round3/s15_predictions.py` |
| 275647095: 1.42× circular duration; R★ 1.14 ± 0.06 R☉; e ≳ 0.34 | | `independent/out/radius_275647095.txt`; `independent/RESULTS.md` | `independent/scripts/s14_radius.py` |
| 275647095: G = 14.3 star at 19.6″, 42σ from the dip | | `independent/out/diffimg_275647095.txt` | `independent/scripts/s13_diffimg.py` |
| 230133616: ecliptic latitude | +82.7° | `work/round3_novelty/exofop_target_230133616.json` (not committed); `NOVELTY_2026-11-22.md` | `round3/s13_novelty.py` |
| 230133616: 37 transits in S16–S86 | | `independent/out/ephem_230133616.txt` | `independent/scripts/s06_ephemeris.py` |
| 230133616: R★ 2.45 R☉, RUWE 3.79, Δχ² ≈ 49, Rp 1.7–1.9 R_J | | `independent/out/radius_230133616.txt`, `gaia_230133616.txt`, `rhofixed_230133616.txt`; `oddeven_summary.csv` | `s14_radius.py`, `s12_gaia.py`, `s15_rho_fixed.py` |
| 230133616: 3100 K M dwarf excluded at > 2σ on 218 of 228 orbits | | `independent/out/secondary_230133616.txt`; `independent/RESULTS.md` | `independent/scripts/s17_secondary_search.py` |
| Three southern stars, T = 12.2–12.4, depths 0.8–1.6% | | `work/coverage.csv` (`tmag`, `dec`); `oddeven_summary.csv` | — |

## Before submitting

- [ ] **The Salinas team has replied, or had a fair chance to.** They are following up their own singles
  (`outreach/2026-11-22-salinas-team-email-draft.md`, not sent).
- [ ] **New sectors.** Fold in the blind session's check of sectors after S107, and the 2026-12-25 transit
  of 107113048 if it was observed.
- [ ] **References.** Check every entry against ADS and arXiv (not done yet). Salinas et al. 2026 is
  arXiv-only so far; check for a journal reference.
- [ ] **AAS AI policy.** Re-check the current AAS journals policy; the disclosure text is copied from
  AAS82207.
- [ ] **Word count.** `python ../singles/rnaas/aas_wordcount.py rnaas.tex` gives an abstract of 143
  (limit 150) and 1,014 counted words (limit 1,350), as of this draft.
- [ ] **Compile.** Compile on Overleaf and read the PDF. `rnaas_draft_reading_copy.pdf` is a check build with
  the older aastex631 class (tectonic), because aastex701.cls only ships with the Overleaf template. It
  compiles with no errors or undefined citations.
