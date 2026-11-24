# Boyajian's star RNAAS note: build, counts and provenance

## Source and build (2026-11-24)

- **Text:** the Claude Doc "RNAAS draft: Boyajian's star companion period"
  (https://claude.ai/artifact/9mW8wn3XzkgUs2u1NcQTye), rev 20, converted to AASTeX by hand. The doc's
  "Notes for Jon" section is left out, and it has been deleted from the doc too.
- **Added in the LaTeX version** (not in the doc): ORCID and email (as in the three earlier notes), the UAT
  keywords, the MAST data DOIs in the acknowledgments, the facilities line and the AI disclosure statement (the
  EB note's wording without its blind re-analysis clause).
- **Affiliation:** "Independent researcher, USA" (Jon, 2026-11-24), as in the three earlier notes (AAS82207,
  AAS84033 and the EB note). The doc says the same.
- **Plain-text abstract for the AAS form:** `abstract_plaintext.txt`, from the LaTeX abstract with the markup
  removed (~ for $\sim$, M_J, hyphens for en dashes); 143 words by `wc -w`.
- **Figure:** `boyajian_alias_figure.png`, a copy of `../boyajian_alias_figure.png` (byte-identical; made by
  `../make_figure_cloud.py`). It is the same image as the doc's Figure 1.
- **Compile:** tectonic with `aastex701.cls` v7.0.1 (from `singles/rnaas/arxiv/McBee_AAS82207_arXiv.tar.gz`):
  3 pages, no undefined references or overfull boxes. `rnaas.pdf` is that build.

## Word counts

- `singles/rnaas/aas_wordcount.py` (the count AAS's quality check applied to AAS82207R1): abstract **144**
  (limit 150), text after the abstract without references **1086** (limit 1350).
- Plain count of the whole PDF (`pdftotext rnaas.pdf - | wc -w`, so title, author block, headers, captions,
  acknowledgments and references; the figure is a raster image, so its labels are not counted): **1458**
  (limit 1500). The doc ran about 1350 by the same kind of count; the additions above take about 100 words.

## Checks done for submission

- **Figure 1 vs text, window 4:** the bottom-right panel (Window 4, P = 1345.65 d, n = -2, BKJD 1205.9)
  shows a dip of about 0.4% against the expected 1.1% transit (red), so the core holds data dimmed by less
  than 0.8%: excluded, as the text says. The top panel's "allowed" row shows the four surviving ranges near
  780 d (776.2-776.6, 790.5-790.7, 792.1-792.9, 793.3-793.4 d), as in the text.
- **References:** all nine checked in ADS (2026-11-24). Zhu & Liu = 2026arXiv261034806Z, "The anomalously
  large transiting planet around Boyajian's star: evidence for a circumplanetary ring" (Wei Zhu, Quanyi Liu);
  cited in Sections 1, 3 and 4. Madurga-Favieres et al. = 2026MNRAS.550g1273M.
- **MAST DOIs:** 10.17909/T9488N resolves to "Kepler LC, Q0-Q17", 10.17909/t9-nmc8-f686 to "TESS Light
  Curves - All Sectors".

## Where the numbers come from

| Numbers | File |
|---|---|
| Kepler survivors 776.17-776.60, 790.50-790.66, 792.08-792.92, 793.25-793.40 d; 781.6 and 800.9 d excluded | `../kepler_alias_survivors.csv`, `../kepler_alias_scan.py`, `../kepler_alias_scan_local.log` |
| Window 4 n = -2 at BKJD 1205.9; windows 1, 2, 3, 5 recovered | `../window4_check.py` → `../window4_check_2026-11-24.log`, `../window_scan_1025_1375.csv` (rerun 2026-11-24, same data and rule as `kepler_alias_scan.py`). All 19 core cadences are dimmed by less than 0.8% at P = 1345.640, 1345.648 and 1345.656 d; median core deficit 0.37%, deepest within +-1 d **0.44%**. The note said 0.45% (the cloud run's `../README.md`, no log) until 2026-11-24, when Jon had it changed to 0.44% in the abstract, Section 3, `abstract_plaintext.txt` and the doc. No surviving range overlaps 1345.0-1346.3 d. Windows 1, 2, 3, 5 survive as 1029.62-1031.54, 1034.51-1035.46, 1300.48-1303.76 and 1366.16-1368.21 d. Kepler alone also leaves 33 other short ranges in 1025-1375 d, each with its transits in gaps |
| TESS S15 fit: BTJD 1730.154 +- 0.003, depth 1.16%, T14 20.3 h | `../README.md` (cloud run) |
| No TESS sector in window 3; next S120-S121 | `../RESULTS_LOCAL.md`, `../tess_sectors.log` |
| AAVSO: 22,234 points; 19 visits on six nights from three observers; 43% coverage; dchi2 -3.4 to +13.3; window 3's 2023 March return +19.4 on 15% | `../RESULTS_AAVSO.md`, `../ground_check_aavso_run.log` |
| ASAS-SN curve ends 2025 June | `../RESULTS_ASASSN.md` |

## Submitted 2026-11-24 as AAS84040

- Submitted through the AAS Peer Review System (aas.msubmit.net) from Jon's signed-in Chrome, driven by Claude with
  Jon's approval at each decision; the system answered "Manuscript Approved (AAS84040)".
- **Submitted text:** `rnaas.tex` at commit `7ef61ff` (tag `aas84040-submitted`; the file was last changed in
  `a74c9a0`). Files: `rnaas.tex` (Manuscript file), `boyajian_alias_figure.png` (Figure file), `aastex701.cls` v7.0.1
  (Other LaTeX file). All three are byte-identical to main at that commit (the class file to the copy in
  `singles/rnaas/arxiv/McBee_AAS82207_arXiv.tar.gz`). No bst file (the note uses `thebibliography`) and no ORCID icon
  file; AAS's build showed the icon anyway.
- **Merged PDF:** built by AAS, 3 pages; all three pages checked against our build (title, affiliation, abstract with
  0.44%, Sections 1-4, Figure 1, acknowledgments with both DOIs, AI statement, nine references). Approved on Jon's yes.
- **Form fields:**
  - corridor: The Solar System, Exoplanets, and Astrobiology (Jon's choice; the K2 notes used Stars and Stellar Physics);
  - author: Mr. Jon McBee, Independent researcher, filled from Jon's AAS profile;
  - title: the LaTeX title with `$\sim$` written as `~`;
  - abstract: `abstract_plaintext.txt`;
  - ethics: "Yes, I agree" (ticked on Jon's yes);
  - AI usage: Yes, described as "Anthropic's Claude (Claude Code and claude.ai, with Claude Opus 5.5) wrote and ran
    the analysis code under my supervision and drafted and edited the note. I reviewed the text, checked every number
    against the saved outputs and every reference against ADS, and take full responsibility for the content. The
    figure was made by the analysis code from the data; no AI-generated figures.";
  - MAST DOIs: 10.17909/T9488N, 10.17909/t9-nmc8-f686; IPAC: No; digital assets: none; funding: No;
  - UAT: Exoplanet detection methods (489), Transit photometry (1709).
