# Be-star RNAAS draft: provenance of every number

Draft: `rnaas/ms.tex`, v3 of 2026-11-24. **Submitted 2026-11-24 as AAS84043** (tag `aas84043-submitted`; see the end of this file). v3 = v2 (revised after
`REVIEW_mythos.md` and `REVIEW_mythos_v2.md`) plus an AI note in the figure caption and a trim to
fit the 1,500-word total (see "v3" at the end of this file).

Word count (v3):
- abstract 146 words (limit 150), from `singles/rnaas/aas_wordcount.py`;
- text after the abstract, excluding references, 1,132 words (the 1,350 limit AAS enforced on
  AAS82207R1), from the same script;
- total with title, abstract, headers, caption and references: 1,473 by the same word rule, and
  1,455 from the compiled PDF's text scaled by the AAS82207R1 calibration (1706/1753). Limit 1,500
  per https://journals.aas.org/research-note-preparation-guidelines/ (read 2026-11-24).

All paths below are under `be-star-pulsations/`. Scripts run with `BE_RUN=full`. Abbreviations:
"stage-4 test" is `stats4.run_test`, and "the c5 log" is the printed output of
`src/c5_review_checks.py`.

## Sample and method

| Number | Source |
|---|---|
| 3,300 candidates; Gaia OBA selection; W1 and W1−W2 changing together; mostly unconfirmed | Jia et al. 2026 (arXiv:2611.02450); `REVIEW_mythos.md` M1 |
| 3,237 non-EB stars; 2,541 pass | `results/full/sample.csv` (`src/s1_sample.py`) |
| Products: 2,020 / 219 / 198 / 104 stars | `results/full/sample.csv` |
| 30-min bins, 5-d knots, 0.25–24 c/d | `src/metrics.py` (`prepare`, `BANDS`, `sector_metrics`) |
| Metric definitions | `src/metrics.py` |
| Flux columns: PDCSAP (SPOC), SAP_FLUX (QLP), CORR_FLUX (ELEANOR) | `src/s1_sample.py` `FLUX_COL` |
| Spline s = 0.04 | `src/wise.py` (`fit_spline`) |
| Gap ≤ 250 d; overshoot ≤ 0.1 mag; W1_quiet | `src/wise.py` (`MAX_GAP`, `MAX_OVERSHOOT`, `DiskCurve`) |
| TESS and WISE about 3 months apart | Observing geometry: TESS near opposition, WISE at 90° elongation (`REVIEW_mythos.md` S5) |
| CROWDSAP ≥ 0.9 / contamination ratio ≤ 0.1 | `src/s3_metrics.py` (`CROWD_MIN`, `CONTRATIO_MAX`) |
| 1,084 stars, 4,108 sectors, 755 QLP | `results/full/test_results.csv`, main rows |
| 10,000 draws per null; bootstrap over stars | `src/stats4.py`, `src/s4_test.py` |
| Injection tests | `tests/test_synthetic.py` |

## Results

### Main test (`results/full/test_results.csv`, main, k = 0)

| Number | Value |
|---|---|
| B1 | 0.80 [0.70, 0.91] |
| red noise | 0.74 |
| flicker | 0.84 [0.73, 0.96] |
| A1 | 0.40 |
| white noise | 0.33 |
| ρ | 0.09–0.20 |
| Noise floor (A1) | 0.54 rms |

The abstract's "about 8%" is B1: e^0.080 ≈ 1.083. The "about 3%" is the 5–24 c/d β of 0.31–0.37 in
`freq_beta.csv`.

### Robustness

| Number | Source |
|---|---|
| No crowding cut: 2,449 stars, B1 0.59 | `test_results.csv`, variant `no crowding cut` |
| Drop mixed-pipeline sectors; drop W1 < 8 | `test_results.csv` variants |
| Drop the top 100 leverage stars: B1 1.03 | `results/checks/review_checks.csv`, M5 rows |
| Linear interpolation | `results/checks/linear_interp.csv` |
| Sector fixed effects: every β changes by < 0.03 (largest +0.025, B1 and red noise) | `results/checks/review_extras.csv`, "two-way star+sector FE" |
| Gaia ESP-ELS beStar: 256 stars, B1 0.58 [0.29, 0.87], p_swap 5e-4 | `review_checks.csv`, M1 rows |

Gaia labels come from `classlabel_espels` in the cached `gaia_dr3_params.parquet`. Of the main
stars, 265 are labelled beStar and 256 have at least 2 usable sectors.

### Frequency bins

| Number | Source |
|---|---|
| Absolute: 0.71–0.92 below 1 c/d; 0.31–0.37 above 5 c/d | `results/full/freq_beta.csv` (`src/s6_freq_beta.py`), main variant |
| Relative, recentred bins: 0.74–1.17 below 0.9 f1; 0.43 [0.34, 0.53] at 0.891–1.122 f1; 0.28–0.49 above | `review_checks.csv`, "S1 relative bin main" rows |

### A1

| Number | Source |
|---|---|
| Stable f1: 510 stars, A1 0.23 [0.10, 0.36] | `review_checks.csv`, M2 rows |
| Stable f1 with f1 > 0.8 c/d: A1 0.08, p_perm 0.37, p_swap 0.29 | `review_checks.csv`, M2 rows |

### White noise

| Number | Source |
|---|---|
| Noisiest tercile: white −0.002 (p_perm 0.98, p_swap 0.87); B1 0.63 | `review_checks.csv`, M4 rows |
| Within-star correlation 0.47 | `results/checks/review_extras.csv` |
| Excess coefficient −0.033 [−0.123, 0.058] with ln B1 in the fit | `review_extras.csv` |
| SPOC formal flux error 0.09, p 0.56 | `results/full/control_tests.csv` |
| Ratios 0.41–0.51; SPOC p_swap ≥ 0.13 | `control_tests.csv` |

### Dilution

| Number | Source |
|---|---|
| QLP adds the TIC-expected flux to the difference flux | Huang et al. 2020, RNAAS 4, 204, §II item 7 and §III (arXiv HTML read 2026-11-24) |
| QLP SAP_FLUX median ≈ 1 per sector | Checked in the cached files |
| k = 0.21 (0.11–0.30), 125 stars | `results/checks/dilution_k.csv` |
| SPOC with k = 0.213: A1 0.369, B1 0.881, flicker 0.936; p_swap ≤ 0.0042 | `review_checks.csv`, M3 rows |
| QLP: 0.52, 0.86, 0.90 | `test_results.csv`, QLP variant |
| 4 shell candidates; changes < 0.01 | `results/checks/dilution_k.csv`, `shell_candidates.csv` |

### Scatter and the discussion

| Number | Source |
|---|---|
| Nearly half (519 of 1,084 = 48%) change by < 0.1 mag | c5 log; `results/checks/review_signs.csv` |
| Of the 177 stars with a range > 0.4 mag, 71% (0.712) have a positive B1 slope | `review_signs.csv` |
| 48% of stars have a negative A1 slope | `src/s7b_bexrb_slopes.py` output; `results/full/LITERATURE_CHECK.md` |
| RX J0440.9+4431 = TIC 410336237, slope −0.86 | `results/full/bexrb_slopes.csv` |

## Figure

The figure is `figures/rnaas_fig1.pdf`, copied to `rnaas/` and unchanged in v2. It is made by
`src/c4_rnaas_fig.py`.

## References and DOIs checked (2026-11-23/24)

**arXiv API**
- Jia et al. 2026 (5 authors);
- Niwano et al. 2024;
- Labadie-Bartz et al. 2025;
- Huang et al. 2020.

**Crossref**
- Niwano, MNRAS 534, 2168;
- Labadie-Bartz et al. 2022, AJ 163, 226;
- Mainzer et al. 2014, ApJ 792, 30;
- Wright et al. 2010, AJ 140, 1868;
- Huat et al. 2009, A&A 506, 95;
- Semaan et al. 2018, A&A 613, A70;
- Creevey et al. 2023, A&A 674, A26 (Apsis I).

**DataCite**
- GSFC-ELEANOR-LITE: 10.17909/j2yt-t417;
- NEOWISE-R single-exposure source table: 10.26131/IRSA144. That is the table queried in
  `src/s2_wise.py` (`neowiser_p1bs_psd`).

**UAT keyword IDs** (checked 2026-11-24 against `UAT.json` in github.com/astrothesaurus/UAT,
release 6.0.0 of 2025-12-17): Be stars 142, Circumstellar disks 235, Stellar pulsations 1625,
Time domain astronomy 2109. All four match their names and none has a change note.

## v3 (2026-11-24): AI note in the caption, trimmed to the 1,500-word total

- The caption now ends "Plotting code was written with Claude (Anthropic).", following the AAS LLM
  guidelines (https://journals.aas.org/author-llm-guidelines, posted 9 Sep 2026). matplotlib is
  added to `\software`.
- The AAS page counts title, headers, captions and references towards 1,500 words. v2 came to
  1,626 on that rule, so the text was tightened. A number-by-number diff against v2 shows these
  removals only (their sources above still hold):
  - the product counts (QLP 2,020, TESS-SPOC 219, SPOC 120-s 198, GSFC-ELEANOR-LITE 104);
  - "755 of the stars from QLP" in Method (still given in the dilution paragraph);
  - the mixed-pipeline and W1 < 8 robustness checks;
  - "2,449 stars" for the no-crowding-cut run (B1 0.59 kept);
  - the SPOC formal flux error check (beta 0.09, p 0.56);
  - "15--24 c/d" in the white-noise paragraph (it is defined in Method).
- Other edits are wording only: shorter caption, swap-null sentence, white-noise, dilution,
  scatter and Discussion paragraphs, the AI statement's list of uses, and the Lightkurve reference
  without its title (ASCL id kept).
- Test build with `tectonic` and `aastex701.cls`: 4 pages, no undefined references, figure on one page.

## Submitted 2026-11-24 as AAS84043

- Submitted through the AAS Peer Review System (aas.msubmit.net) from Jon's signed-in Chrome, driven by Claude on
  Jon's instruction; the ethics statement was ticked on Jon's explicit yes. The system answered
  "Manuscript Approved (AAS84043)".
- **Submitted text:** `rnaas/ms.tex` at main commit `104032a` (tag `aas84043-submitted`; the file was last changed in
  `35ba2d0`). Files: `ms.tex` (Manuscript file), `rnaas_fig1.pdf` (Figure file), `aastex701.cls` v7.0.1 (Other LaTeX
  file, the copy in `singles/rnaas/arxiv/McBee_AAS82207_arXiv.tar.gz`). `ms.tex` and the figure are byte-identical to
  main at that commit.
- **Merged PDF:** built by AAS, 4 pages; checked page by page against our tectonic build (title, author with ORCID
  icon, abstract, keywords, Sections 1-3, Figure 1 with the Claude caption line, acknowledgments with the five DOIs,
  facilities, software with matplotlib, AI statement, 13 references). Approved.
- **Form fields:**
  - corridor: Stars and Stellar Physics;
  - author: Mr. Jon McBee, Independent researcher, United States (from Jon's AAS profile);
  - title: from `ms.tex`; abstract: `rnaas/abstract_plaintext.txt` (form count 148 of 150);
  - ethics: "Yes, I agree";
  - AI usage: Yes, with the note's AI disclosure paragraph verbatim (no model version given);
  - MAST DOIs: 10.17909/t9-nmc8-f686, 10.17909/t9-wpz1-8s54, 10.17909/t9-r086-e880, 10.17909/j2yt-t417;
  - IPAC: Yes (NEOWISE-R from IRSA; the form takes no DOI, and 10.26131/IRSA144 is cited in the note);
  - digital assets: none; funding: No;
  - UAT: Be stars (142), Circumstellar disks (235), Stellar pulsations (1625), Time domain astronomy (2109).

## Shortened for the AAS word count (2026-11-25), not uploaded

The 2026-11-25 audit (`rnaas/audit-2026-11-25.md`, branch `rnaas-audit-2026-11-25`) found the submitted text at 1366
words by the gate (references included), AAS range 1356-1396, over AAS's 1350. It passed the quality check and is
Under Review. Branch `be-aas84043-length` cuts 94 words of wording: no result, number or citation changed; the
abstract is identical; a numeric-token diff against tag `aas84043-submitted` shows nothing lost (two "at least two"
became "$\geq2$"). Gate (d385d67): abstract 146, main text 1272, AAS range 1263-1300. `ms.pdf` is the tectonic
build (4 pages; page 1 checked by eye). The pre-existing "Float too large for page by 28.4pt" warning for the figure
is unchanged (the figure was not touched). To reach AAS it needs the manuscript back in Author Approval.
