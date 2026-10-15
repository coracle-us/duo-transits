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


- **Sample arithmetic** (Jon's reviewer, 2026-10-12). The 359 stars break down as:
  - 32 skipped because a community period already exists;
  - 22 outside the cuts in `run_remaining.py` (Tmag ≥ 7 and depth 200–50,000 ppm): 7 brighter than
    T = 7 and 15 with depths outside that range, mostly placeholder values near 100%;
  - 305 analysed.

  The 306 classified also include TIC 53847727, a pilot target later removed from the sample as a
  Prša+2022 EB. The fix is a sentence in `rnaas.tex` (revision-ready); the submitted version
  lacks it.

**Known minor issue (found 2026-10-12, K2 fact-check).** `singles/writeups/tic159614322/ephemeris.py`
converts BJD_TDB to UTC without removing the barycentric light-travel time. For TIC 159614322 this
is 1.2 min: the 2026-12-28 transit is 06:40 geocentric UTC, not 06:42. That is negligible against
the ±45 min uncertainty. Fix it in any revision: quote BJD_TDB, or apply the correction.

## AAS82207 revision, after the Fable adversarial review (2026-10-13)

Review: `singles/review_fable_revision/FINDINGS.txt`. New or changed numbers in `rnaas.tex`:

- **Event test used pipeline light curves** (SPOC 2-min for 285 stars): `singles/revision/audit_34.csv`
  (`pipeline_product`), from `singles/revision/audit_34.py`.
- **18 of the 31 failures from one uploader batch** (ExoFOP user `amarino`, UAN), including all
  six half-day slips: `audit_34.csv`, `uploader` column, after dropping the three truncated-TIC
  duplicates (1461133, 3181884, 4689632).
- **Three near-threshold cases** (TIC 38005437, 355095316, 425585592) that pass in another
  pipeline product: `audit_34.csv`, `depth_snr_per_product`.
- **Any-depth search, 14 stars with vetted extra eclipses** (9 binaries, 4 no-single-period, 1
  single = TIC 250492521, SNR 10.0, three reductions); periods for 4 stars including TIC 188835620
  (22.217 d) and TIC 83384234 (4.33 d): `singles/revision_search/manual_vetting.csv`,
  `per_star_summary.csv`, `period_check.py`.
- **234/272 re-detected; 21 of the 38 misses below SNR 8, 6 momentum dumps**:
  `singles/revision_search/recovery_check.csv`, `REPORT.txt` §3.
- **TIC 147035057 primary depth 5.3–5.6%**: 5.3% minimum in S6 and S87 from
  `singles/revision_search/verify_147035057.py`; 5.6% PDCSAP/QLP point minimum. **Neighbour flux:**
  target G = 8.83; Gaia DR3 sources within 60″ supply 3.5% of the target flux in total, the
  brightest G = 13.50 at 16.5″ and G = 13.53 at 24.4″ (~1.3% each):
  `singles/revision_search/neighbours_147035057.py` (Gaia DR3 cone search, run 2026-10-13).
- **Photometry validation: ten usable stars, T = 8.5–14.1**; the 11th (T = 15.5) is 96%
  contaminated and excluded from the 1.07 median; worst depth ratio 0.69 ("about 30% shallower"):
  `singles/revision/REPORT.txt`.
- **Period range 5–1500 d** for the repeat search: `singles/scripts` alias grid (as for the 52%
  log-period exclusion in the abstract).
- Sample rule wording: "794 open CTOIs with no period or tagged single, duo or mono", from
  `singles/scripts/build_ctoi_singles.py` and `singles/revision/count_owner_exclusion.py`.
- Removed: the HD 88986 b sentence and the Heidari et al. 2024 reference; "offsets on known
  planets are 1–3″"; the backgrounds 1–9% clause (still in the letter). Body 1339 words, abstract 146.

## Resubmitted 2026-10-13 as AAS82207R1

- Received by AAS 13-Oct-2026 10:37, manuscript AAS82207R1.
- Submitted text = commit 99afd3f (`singles/rnaas/rnaas.tex`; git tag `aas82207R1-submitted`), uploaded as a
  single self-contained file `McBee_AAS82207_R1.tex` with aastex701.cls, aasjournalv7.bst and orcid-ID.png.
  The AAS-built merged PDF was checked word for word against the Overleaf build: identical, 4 pages.
- The plain-text abstract in the submission form (`private/aas82207_R1_abstract_plaintext.txt`, 149 words by
  the form's space-split count) differs from the manuscript only typographically ("37.315-d", "P=89.965") and
  in "or none" for "or no dip".
- Response to the editor pasted in "Reply to Referee": `private/aas82207_response_to_editor.txt`.

## AAS82207R1 bounced on length (2026-10-13), shortened

AAS quality check: "over the 1350 word limit by 356 words" (= 1706). AAS counts everything after the abstract
except the references: body + table + acknowledgments + facilities/software + AI statement. The old counter counted
the body only (1345). New counter: `singles/rnaas/aas_wordcount.py`, which gives 1668 on the submitted R1 (about 2%
below AAS). The text was rewritten more tightly: Sample and Methods merged, bullet lists folded into prose, the
gap-truncated-transit sweep sentence and the HD-number example dropped, table notes and acknowledgments shortened. A
number-by-number diff against tag `aas82207R1-submitted` shows only intentional removals (SNR ~20 population
remark, the Sector 91 range bound, the truncated-TIC example) and no new or changed numbers. New count: abstract 147,
AAS-counted text 1268 (about 1295 by AAS's count).

**Resubmitted 2026-10-13 after the length bounce:** text = commit 9bfccae (tag `aas82207R1-length-fix`). The AAS
merged PDF (3 pages) was checked: about 1316 words after the abstract excluding references, scaled from the
pdftotext count with the ratio calibrated on the bounce (AAS 1706 vs pdftotext 1753).

## Second revision request, 2026-10-15 (AAS82207R1 -> R2)

Editor: the start of Section 2 ("The event test") reads as something new. Fix: Section 1 now introduces the check by
name ("I first check that each catalogued event is present: ..."), and Section 2 opens "The check for the catalogued
event (Section 1) used pipeline light curves ...". No numbers changed. Count (`aas_wordcount.py`): abstract 147,
text 1275 (about 1325 by AAS's count). Upload file: `private/McBee_AAS82207_R2.tex`; response:
`private/aas82207_R2_response.txt`.

**Resubmitted 2026-10-15 as AAS82207R2** (received 13:25): text = commit 7b8647e (tag `aas82207R2-submitted`). AAS merged
PDF checked: 3 pages, both clarifications present, about 1323 words by the AAS count. Abstract field: plain text from
`private/aas82207_R1_abstract_plaintext.txt`.
