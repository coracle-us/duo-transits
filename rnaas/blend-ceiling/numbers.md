# Blend-ceiling RNAAS note: build, counts and provenance

## Source and build (2026-11-25)
- **Text:** `DRAFT_RNAAS.md` at 800412e, converted to AASTeX by hand as `rnaas.tex`. The introduction was reworded
  to match the two arXiv abstracts (below).
- **Compile:** tectonic with `aastex701.cls` v7.0.1 (from `singles/rnaas/arxiv/McBee_AAS82207_arXiv.tar.gz`; not
  committed here). Result: 3 pages; no undefined references, no overfull or underfull boxes. The only warnings are
  lineno.sty's invalid UTF-8 byte, from the class package. `rnaas.pdf` is that build.
- **Gate:** `scripts/rnaas_gate.py` from rnaas-word-gate d385d67:
  `rnaas_gate PASS: rnaas.tex: abstract 135/150; main text with references 1190/1300; AAS range 1182-1216 (limit 1350); page 1: rnaas_page1.png, 28 words before ABSTRACT`
- **Page 1** (`rnaas_page1.png`), checked by eye: draft header, title, author, affiliation, abstract, keywords,
  then section 1. No body text above the title. Pages 2 and 3 were also checked: figure 1 and its caption on page 3,
  and nine references.

## References (checked on ADS, 2026-11-25)
All nine resolve on ADS with the authors, year, journal, volume and page given:
- 2023A&A...674A...1G
- 2021ApJS..254...39G
- 2020RNAAS...4..204H
- 2016SPIE.9913E..3EJ
- 2018ascl.soft12013L
- 2026arXiv261042009P (single author, Pagés Navarrete, Juan F.)
- 2011JMLR...12.2825P
- 2022ApJS..258...16P
- 2026arXiv261024882S (six authors: Salinas, Brahm, Olmschenk, Soares-Furtado, Pichara, Schlieder)

Claims about the two arXiv papers:
- **Pagés Navarrete:**
  - From the abstract: ultra-short-period planets; AUC 0.97 vs catalogued EBs and 0.50 vs TFOP false positives;
    at least 42% blends of the 145 examined; bound 1 − f/2 ≈ 0.79.
  - The P < 1.5 d test-set cut (Sect. 3.2) was seen only through a summary of the HTML, so the note no longer
    states it.
- **Salinas et al.:**
  - From the abstract: the CNN predicts depth, duration and mid-time with calibrated uncertainties; it is trained
    on synthetic light curves; it is evaluated on 120 confirmed exoplanets.
  - "P > 10 d" and the CP/KP TOI selection (Sect. 5.1) were seen only through a summary of the HTML.
  - That the 169 new candidates came from their earlier Transformer was also seen only through the summary.

## Every number in the note
| Number | Source |
|---|---|
| 336 / 215 scored, 216 FP/FA, median 8 windows | `results/auc.txt`, `data/features.csv` |
| 410 EBs, 1146 TCE-like, 50 sanity windows | `results/auc.txt` |
| 1182 quiet windows | `data/features_quiet.csv` (make_quiet.py) |
| f, ceilings, Wilson intervals, rule G rates (70%, 60%) | `results/blend_fraction.txt` |
| Model B AUCs, sanity medians, depth bins, KS, Tmag | `results/auc.txt` |
| Model A AUCs (brackets), sanity median 0.91 | `results/amendmentA/auc.txt`; A no-transit AUC 0.786 from `results/amendmentA/object_scores.csv` (also plotted in fig1) |
| Box SNR medians in RULES amendment B | the amendment-A run (summarised in RULES.md) |
