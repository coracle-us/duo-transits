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
| 275647095: dip at phase 0.872, 694 ± 147 ppm, 8 of 9 orbits, no background or centroid change; six events in S111–S114 | | `independent/out/secondary_275647095.txt`, `phasedip_275647095.txt`; `independent/RESULTS.md` | `independent/scripts/s17_secondary_search.py`, `s20_phase_dip.py` |
| 107113048: E42 outside TESS coverage; S114 tests 74.55 and 29.82 d | | `independent/out/predictions.txt`; `independent/RESULTS.md`, "Predictions" | `independent/scripts/s19_predictions.py` (tess-point) |
| No TESS data after S102 for any of the four (MAST, 2026-11-23) | | `independent/RESULTS.md`, "Predictions" | same |
| 230133616: ecliptic latitude | +82.7° | `work/round3_novelty/exofop_target_230133616.json` (not committed); `NOVELTY_2026-11-22.md` | `round3/s13_novelty.py` |
| 230133616: 37 transits in S16–S86 | | `independent/out/ephem_230133616.txt` | `independent/scripts/s06_ephemeris.py` |
| 230133616: R★ 2.45 R☉, RUWE 3.79, Δχ² ≈ 49, Rp 1.7–1.9 R_J | | `independent/out/radius_230133616.txt`, `gaia_230133616.txt`, `rhofixed_230133616.txt`; `oddeven_summary.csv` | `s14_radius.py`, `s12_gaia.py`, `s15_rho_fixed.py` |
| 230133616: 3100 K M dwarf excluded at > 2σ on 218 of 228 orbits | | `independent/out/secondary_230133616.txt`; `independent/RESULTS.md` | `independent/scripts/s17_secondary_search.py` |
| Three southern stars, T = 12.2–12.4, depths 0.8–1.6% | | `work/coverage.csv` (`tmag`, `dec`); `oddeven_summary.csv` | — |

## Before submitting

- [ ] **The Salinas team has replied, or had a fair chance to.** They are following up their own singles
  (`outreach/2026-11-22-salinas-team-email-draft.md`, not sent).
- [x] **New sectors (2026-11-23).** There is no TESS data after S102, and none of the four is on silicon in
  S103–S110.
- [ ] **S111–S114.** Fold these in when public:
  - the 275647095 phase-0.872 dip (six events);
  - 107113048 E43 and the 74.55 / 29.82 d tests in S114;
  - a ground-based 107113048 E42 (2026-12-25), if anyone observes it.
- [x] **References (2026-11-24).** `checks/verify_refs.py` → `checks/refs_check.txt`: each of the 11 entries matches
  exactly one ADS record (authors, year, journal, volume, page), and the two ASCL entries match their ascl ids.
  Salinas et al. 2026 is still arXiv-only: arXiv:2610.24882 has no journal reference or DOI (last updated
  2026-10-17). Re-check just before submitting.
- [x] **AAS AI policy (2026-11-24).** Read https://journals.aas.org/author-llm-guidelines (posted 9 Sep 2026,
  no later update shown). It asks for a disclosure statement at the end of the article naming the tool and
  what it was used for, with responsibility resting on the author; AI cannot be an author; every reference
  must be checked to exist (ADS recommended); AI-made figures need a caption note. The draft's "AI disclosure
  statement" names Claude (Anthropic) and its uses, takes responsibility, and the note has no figure, so it
  complies. Nothing specific to Research Notes. Optional: the page suggests telling the editor at submission
  which sections the LLM helped with.
- [x] **Word count (2026-11-24).** The RNAAS page (https://journals.aas.org/research-note-preparation-guidelines/,
  read 2026-11-24) counts title, headers, captions and references towards 1,500 words, with 150 reserved for
  the abstract; the Editorial Office's count is final.
  - **Abstract trimmed by 3 words** ("Two of these are" → "Two are"; "TIC 241550976 and TIC 107113048" →
    "TIC 241550976 and 107113048"). `aas_wordcount.py` counts `\citet{Salinas2026}` as one word, but it prints
    as "Salinas et al. (2026)", so the printed abstract was 152. Now: 146 by `aas_wordcount.py`, 150 as
    printed (PDF text) and 150 in `abstract_plaintext.txt` by `wc -w` (limit 150).
  - Text after the abstract, without references: 1,109 by `aas_wordcount.py` (the 1,350 AAS enforced on
    AAS82207R1).
  - Total by `checks/total_words.py` (same word rule): 1,486 with the table body, 1,430 without it (the page
    says table contents don't count), 1,491 with the author line. From the PDF text: 1,529 raw, 1,488 scaled
    by the AAS82207R1 calibration (1706/1753). Limit 1,500: it fits, but with little room; any added
    sentence (for example after the Salinas reply or S111) needs a matching cut.
- [x] **Compile (2026-11-24).** tectonic 0.17.0 with `aastex701.cls` v7.0.1 (2025/05/09, from
  `singles/rnaas/arxiv/McBee_AAS82207_arXiv.tar.gz`, as for AAS84040/41/43): 3 pages, the table on page 2, no
  undefined references or citations, no overfull boxes. `rnaas.pdf` is that build. It replaces
  `rnaas_draft_reading_copy.pdf` (aastex631) as the copy to read.
- [x] **Keywords and data DOIs (2026-11-24).** UAT 489 "Exoplanet detection methods" and 1709 "Transit
  photometry" match their names in the current UAT.json (github.com/astrothesaurus/UAT), with no change notes.
  All six `\dataset` DOIs resolve (doi.org handle API): SPOC t9-nmc8-f686, TESS-SPOC t9-wpz1-8s54, QLP
  t9-r086-e880 and TESS FFI 0cp4-2j79 at MAST; ExoFOP5 and NEA13 at IPAC.
- [x] **Plain-text abstract for the AAS form:** `abstract_plaintext.txt` (markup removed; hyphens for en
  dashes, "sigma" for σ), 150 words by `wc -w`. Affiliation "Independent researcher, USA", as in the other
  notes.
