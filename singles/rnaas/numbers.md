# RNAAS note (`rnaas.tex`): where every number comes from, and the pre-submission checklist

The CLAUDE.md rule is that every number in the note must trace to a file and a script.

## Sample and counts

| Quantity in the note | Value | File | Script |
|---|---|---|---|
| Open single/duo CTOIs (not active groups) | 450 | `singles/data/ctoi_singles.csv` + `..._dropped.csv` | `singles/scripts/build_ctoi_singles.py` |
| Removed: missing params / confirmed hosts / TOI / TESS EB | 46 / 8 / 4 / 24 | `singles/data/ctoi_singles_dropped.csv` | same |
| CTOIs / stars kept | 368 / 359 | `singles/data/ctoi_singles.csv` | same |
| Stars skipped (community period on ExoFOP) | 32 | `singles/results/full/skipped.csv` | `singles/scripts/run_remaining.py` |
| Stars classified | 306 | `singles/results/final_classification_2026-10-11_rnaas.csv` (frozen copy for the note) | `singles/scripts/classify_final.py` |
| Event not recovered / single | 34 / 193 | same | same |
| With repeats (event recovered) | 79 | same | same |
| Repeats: Gaia binary / firmly off-target / inconsistent / hand-vetted | 21 / 18 / 23 / 17 | same | same |
| Hand-vetted: survive / tentative / possibly off-target / rejected | 2 / 3 / 5 / 7 | `singles/results/vetting/leads_full_run.md`, MANUAL notes in `classify_final.py` | `scripts/vet_candidate.py` |
| Median excluded log-P fraction (singles) | 0.52 | `final_classification.csv` (`frac_logP_excluded`) | `classify_final.py` |
| Validation: 26/26 kept, 19 unique, 0 wrong | | `singles/results/validation/` (see `singles/README.md`) | `singles/scripts/` validation run |
| Own photometry recovers a known eclipse to 2% | 3.21 ± 0.07 vs 3.28 ppt (S33, TIC 147035057) | `singles/writeups/tic147035057/build_pdf.py` text; `singles/README.md` | `singles/scripts/tesscut_photometry.py` |
| Gaia SB1 matches | 21.978/21.969, 22.259/22.245, 18.36/18.355 d | `singles/data/gaia_binarity_flags.csv`; batch-2 notes in `singles/README.md` | `singles/scripts/gaia_flags.py` |
| Partial-transit sweep: 2-product dips vs bumps | 48 vs 37 (27 vs 24 excl. 236770234) | `singles/results/partial_sweep/summary.csv` | `singles/scripts/partial_sweep.py` |

## Individual systems

| Quantity | Value | File | Script |
|---|---|---|---|
| 147035057: P, T0 | 37.31497 ± 0.00003 d | `singles/results/vetting/147035057_ephemeris.json` | vetting scripts (see README) |
| 147035057: ΔBIC, k, R₂ limit | ≈10; 0.34 (0.09–0.75); > 1.7 R_J | `singles/writeups/tic147035057/shape_fit.json` | `shape_fit.py` |
| 147035057: RV amplitude, M₂, RUWE | 50 km/s; ≥ 0.64 M☉; 1.12 | `singles/writeups/tic147035057/gaia_tension.json` | `gaia_tension.py` |
| 159614322: P, T0 | 89.9651 ± 0.0006 / 179.9302 ± 0.0012; 2458703.8277 ± 0.0041 | `singles/writeups/tic159614322/ephemeris.json` | `ephemeris.py` |
| 159614322: S21 O−C, inflation | +13 to +26 min; ×2.1 | same | same |
| 159614322: next transit | 2026-12-28 06:42 UTC ± 45 min | `singles/writeups/tic159614322/predictions.csv` | `ephemeris.py` |
| 159614322: depth, duration, Rp | ≈3.2 ppt, ≈9.7 h, ≈0.9 R_J | `analysis.json`, `independent/results.json` | `analysis.py`, independent scripts |
| 229055790: aliases, Rp | 38, ≥ 26.7 d, ≈3.8 R⊕ | `singles/results/full/229055790.json`, `leads_full_run.md` | `pilot.py` |
| Table event sectors | | per-target JSONs in `singles/results/<run>/` | sector lookup from `data/raw/tess_orbit_times.csv` |

## Before submitting

- [x] **References** (verified 2026-10-12). All 12 entries were checked for first and co-authors,
      journal, volume and page against arXiv API metadata, ASCL, the NASA TESS BibTeX file
      (Huang 2020, Caldwell 2020) and the A&A/arXiv listing (Heidari 2024). All match `rnaas.tex`;
      no changes were needed.
- [x] **ORCID** added: 0009-0000-7524-4470.
- [x] **Affiliation:** "Independent researcher, USA" (confirmed by Jon 2026-10-12).
- [x] **Repository** pushed to `github.com/coracle-us/duo-transits` (public) on 2026-10-12, with
      correspondence and family documents kept private. Optional: a Zenodo DOI for the release.
- [x] **AI policy** read by Jon (2026-10-12). A dedicated AI disclosure statement was added
      at the end of the note, following the AAS examples.
- [x] **Compiled** on Overleaf, using the official AAS template (`aastex701.cls`), with no errors.
- [x] **Submitted** by Jon on 2026-10-12 via Overleaf's 'Submit to AAS' (not waiting for Gaia
      DR4). Next: AAS emails instructions to finish the metadata on its site (choose RNAAS;
      paste the AI disclosure). Gaia DR4 (2 Dec) can go in a follow-up note.

## After submission

- **Received** by AAS on 12 Oct 2026 as manuscript **AAS82207** (Research Notes of the AAS,
  corridor "The Solar System, Exoplanets, and Astrobiology").
- **Submission site fields.** The abstract was entered as plain text and is the same as the
  manuscript. The UAT concepts are those in the note. Funding: none.
- **MAST dataset DOIs entered:** 10.17909/t9-nmc8-f686 (SPOC LCs), 10.17909/t9-wpz1-8s54
  (TESS-SPOC), 10.17909/t9-r086-e880 (QLP) and 10.17909/0cp4-2j79 (FFIs). All four were verified
  on DataCite on 2026-10-12.
- **IPAC: yes.** We cite ExoFOP 10.26134/ExoFOP5 and NASA Exoplanet Archive PSCompPars
  10.26133/NEA13.
- **Revision-ready.** The `\dataset` tags for all six DOIs are in the acknowledgments of
  `rnaas.tex`, as AAS asks for revised manuscripts. They are not in the submitted version. Use
  this file if the editor asks for changes.

