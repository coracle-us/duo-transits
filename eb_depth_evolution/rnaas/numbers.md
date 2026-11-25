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

## Fact-check fixes (2026-10-12, `eb_depth_evolution/review_fable/REPORT.txt`) and new provenance

**Word count:** abstract 136, body 932 (`rnaas/wordcount.py`). The "137 / 734" line above is superseded.

| Fix | New text and source |
|---|---|
| R1 rate | "0.31–0.36 deg/yr (2021–2026)" = linear fits to the three TESS seasons: blind 0.31, ours 0.358 (`scripts/29_note_figure.py` stdout). The difference is mostly the 2026 eclipse. "Constant excluded": a constant 0.31 predicts a ~1% K2 primary (blind REPORT; review_fable check) |
| R2 inclinations | Base values (`data/inclination_trend.csv`, variant base, ph14 0.3322): 85.74, 87.26, 87.91, 88.77, so 85.7, 87.3, 87.9, 88.8. The abstract gives the 2026 span of both analyses, 88.5–88.8 |
| R3/R4/R5 timing | "Interval from primary to secondary fell by 26 ± 2.5 min within one orbit (2021); 3–4 min/orbit on average 2021–23" (`data/sec_minus_pri_interval.csv`). The 2023 "~50 min" is dropped: it relied on the excluded partial S72 secondary |
| R6 totality | "2027–2031": ours 2027.2–2027.6 (threshold 89.25, i 88.77, rate 0.36–0.40); blind 2027.5–2031 |
| R7 turnaround | Now specific: "a quadratic trend through the four epochs, turning around in 2010" (`data/inclination_trend.csv`, quad columns) |
| R8 figure | New `rnaas/eb_history.png` (`scripts/29_note_figure.py`): measured values from both analyses, TESS-era linear fits, thresholds, totality band; no rejected model curves |
| R9 K2 secondary | "< 0.1%" (consistent with AAS82225; blind < 0.08%) |
| S3 duration | "about 8 h": model 8.1–8.5 h, blind 7.6–8.4 h (review_fable check_numbers_out.txt) |
| S6 apsidal motion | "GR + tidal + rotational < 0.01 deg/yr": GR 0.0035, tidal < 2e-4, rotational < 0.007 deg/yr (review_fable check_numbers_out.txt) |
| S7/S8 Gaia | "tertiary below ~0.3 Msun or an outer period near 100 d" (review_fable: 0.5 Msun at 200 d → 0.18 mas) |
| S1/S2 refs | Added Zasche+2023, Davenport+2021, Borkovits+2016, Luger+2016 (EVEREST), Vanderburg & Johnson 2014 (K2SFF), Butters+2010 (WASP). All verified in ADS (`scripts/27_verify_refs.py` → `refs_check.txt`); 20 journal refs plus 2 ASCL |
| S12 overlap | The intro says the companion note (J. McBee 2026, submitted) reported the deepening, and lists what's new here. §2 says "The depths repeat those of the companion note" |

**Deferred:**
- the exact DASCH NSF-grant acknowledgment sentence (check the DASCH site before submission);
- the NASA Exoplanet Archive acknowledgment wording.

**STATUS (superseded 2026-11-15, see below): content-final draft, HELD until the AAS82225 decision.**

## Standalone version, 2026-11-15 (hold lifted)

On 2026-11-15 the editor declined AAS82225 (the K2 note) and invited a revision with a literature search and a
blind transit search. The hold existed only to avoid overlap with AAS82225, so this note now goes first on its
own. A K2 revision should cite this note for EPIC 202126877. Changes:

- The "companion note (J. McBee 2026, submitted)" sentence and "The depths repeat those of the companion note" are
  removed. The intro now says the note shows the deepening itself. No numbers changed.
- **Predictions:** 2026 Dec 8 11:46, Dec 30 06:21, 2027 Jan 21 00:58 and Feb 11 19:35 UTC, from
  `scripts/28_predictions.py` with the start date set to 2026-11-20 (`predictions_2026_27.csv`, 8 epochs to 2027 May 9).
  Same ephemeris as before.
- **DASCH:** the acknowledgment sentence and the principal citation, Grindlay et al. 2012 (IAU Symp. 285, 29;
  2012IAUS..285...29G; doi:10.1017/S1743921312000166), are those asked for on https://dasch.cfa.harvard.edu/citing/
  (read 2026-11-15). This replaces Grindlay et al. 2009. The four authors come from arXiv:1211.1051. The bibcode
  comes from the DASCH page; this Mac has no ADS token, so `27_verify_refs.py` was not rerun.
- **NASA Exoplanet Archive:** the standard Archive acknowledgment, covering the KELT-North light curves.
- **AI statement:** now first person, matching the wording accepted in AAS82207.
- **No new TESS data:** TESScut and tess-point list only S43-45, S71-72 and S1751 (checked 2026-11-15), so the
  depths and the inclination history are current.
- **Novelty recheck:** an arXiv API search (2026-11-15) for 202126877, 86491361 and "TYC 1883-1721-1" found nothing.
  The ADS full-text search dates from 2026-10-12.
- **Word count** (`singles/rnaas/aas_wordcount.py`): abstract 142, AAS-counted text 1217 (limit 1350).
  `rnaas/wordcount.py`: abstract 136, body 907.
- **Build:** tectonic 0.17.0 with aastex701.cls 7.0.1 gives 4 pages with the figure and no undefined references.

**Submission form fields:**
- Title: as in `rnaas.tex`.
- Abstract: `rnaas/abstract_plaintext.txt`.
- UAT: Eclipsing binary stars 444, Trinary stars 1714, Orbital evolution 1178.
- Corridor: Stellar Physics and Stellar Types (pick the closest stellar corridor in the dropdown).
- Funding: none.
- MAST DOIs:
  - 10.17909/t9ws3r (K2);
  - 10.17909/t9-wpz1-8s54 (TESS-SPOC);
  - 10.17909/t9-r086-e880 (QLP);
  - 10.17909/0cp4-2j79 (FFIs).
- IPAC: yes. The note uses NASA Exoplanet Archive data (KELT).
- Upload: `rnaas.tex`, `eb_history.png` and aastex701.cls.

**STATUS: submitted 2026-11-24 as AAS84041 (see the end of this file).**


## Companion-note citation restored (2026-11-24)

The K2 revision went in this morning as AAS84033 and reports the same EPIC 202126877 depths (6% in 2014 to 36% in
2026; secondary 4% to 12%) without citing this note. To avoid two overlapping notes that don't mention each other,
Jon chose to restore the October wording before submitting: the intro again says the companion note on K2 single
events (J. McBee 2026, submitted) reported the deepening from the same TESS photometry and lists what is new here;
Section 2 again says "The depths repeat those of the companion note" and "The ground data are new". No numbers
changed. Count (`singles/rnaas/aas_wordcount.py`): abstract 142, AAS-counted text 1243 (was 1217; limit 1350).
Build: tectonic with aastex701.cls 7.0.1, 4 pages, no undefined references.

## Submitted 2026-11-24 as AAS84041

- Submitted through the AAS Peer Review System from Jon's signed-in Chrome, driven by Claude on Jon's instruction
  ("submit the EB note the same way" as AAS84040); the system answered "Manuscript Approved (AAS84041)".
- **Submitted text:** `rnaas.tex` at commit `0d0714c` (tag `aas84041-submitted`), i.e. with the companion-note
  citation restored that morning. Files: `rnaas.tex` (Manuscript file), `eb_history.png` (Figure file),
  `aastex701.cls` v7.0.1 (Other LaTeX file); all byte-identical to main at that commit.
- **Merged PDF:** built by AAS, 4 pages; checked page by page against our build (abstract, the companion-note
  sentence, Sections 1-4, Figure 1, AI statement, references). Approved.
- **Form fields:**
  - corridor: Stars and Stellar Physics;
  - author: Mr. Jon McBee, Independent researcher (from Jon's AAS profile);
  - title as in `rnaas.tex`; abstract: `abstract_plaintext.txt`;
  - ethics: "Yes, I agree";
  - AI usage: Yes, "Anthropic's Claude (Claude Code, with Claude Opus 5.5 and Claude Fable 5.2 models) wrote and ran
    the analysis code under my supervision, including a blind re-analysis and adversarial reviews by separate AI
    instances, and drafted and edited the note. I reviewed the text, checked every number against the saved outputs
    and every reference against ADS, and take full responsibility for the content. The figure was made by the
    analysis code from the data; no AI-generated figures.";
  - MAST DOIs: 10.17909/t9ws3r, 10.17909/t9-wpz1-8s54, 10.17909/0cp4-2j79, 10.17909/t9-r086-e880;
  - IPAC: Yes (KELT-North from the NASA Exoplanet Archive; the form asks for no IPAC DOI at submission, but its
    instructions say to cite an IPAC dataset DOI in a revised manuscript: worth adding if a revision is requested);
  - digital assets: none; funding: No;
  - UAT: Eclipsing binary stars (444), Trinary stars (1714), Orbital evolution (1178).

## AAS84041 bounced on length at the initial quality check (2026-11-24); shortened

AAS: "Your abstract is 143 words. Your main text is over the 1350 word limit by 300 words", counting text,
acknowledgments, figure/table captions **and references**. Our old counter left the references out (1243); with them
the submitted text is 1662 (`singles/rnaas/aas_wordcount.py`, which now prints both). Jon asked for 1,300 or fewer,
keeping every result, number and the AI disclosure, with no change to the analysis.

- **New count:** abstract 142 (unchanged); AAS-counted text with references 1298 (was 1662); without references 941.
  pdftotext of the new PDF gives 1181 words after the abstract including references (1667 for the submitted PDF).
- **What was cut:** wording only. Every number in the old text after the abstract is still there (checked by
  script); the abstract, title, figure and AI disclosure are unchanged.
  - Headings: "Introduction" dropped; "Data and depths" -> "Data"; Predictions merged into "Timing and a third body".
  - Intro: dropped "a few dozen Kepler and TESS systems" and "constrain the architecture of hierarchical triples";
    V907 Sco now cites Zasche et al. 2023 only (Lacy et al. 1999 removed from the references).
  - Data list as one sentence instead of bullets; fit description as prose instead of bullets.
  - Dropped: "The data cannot tell whether the rate is still increasing" (the predictions paragraph still says the
    coming timings will test it); "within reach of small telescopes"; "set by the timing variations" after
    +/-30 min.
  - Caption, from 92 to about 60 words; acknowledgments, from 85 to 39 words (same four MAST DOIs, same NSF
    grants for DASCH; standard long-form archive sentences shortened).
  - References: ASCL titles removed for Astrocut and Lightkurve; editors and publisher removed for Grindlay 2012.
- **Build:** tectonic, aastex701.cls 7.0.1, 4 pages, no undefined references; `rnaas.pdf` committed.
- **Note:** the AAS84041 form's IPAC answer is Yes. The form's instructions ask for an IPAC dataset DOI in a revised
  manuscript; none is in the text yet (adding one costs words).
- **Gate** (`scripts/rnaas_gate.py`, required from 2026-11-24): rnaas_gate PASS: eb_depth_evolution/rnaas/rnaas.tex: abstract 142/150; main text with references 1298/1300 (AAS limit 1350)

**Back to Author Approval, files uploaded (2026-11-24, evening).** A note went to AAS staff through the AAS84041
"Send Manuscript Correspondence" form (Recipient: Staff; no attachments; Jon solved the CAPTCHA; the system said
"1 Email sent."). Subject "AAS84041: please return to Author Approval for the shortened version"; it gave the new
count (1,298 with references, abstract 142 and unchanged) and asked staff to return the manuscript. It was then in
Awaiting Author Approval. Uploaded: `McBee_AAS84041_R1.tex` (= 4438ecf `rnaas.tex`) replacing the old manuscript
file; figure and class file unchanged (byte-identical). The portal did not rebuild its merged PDF, so the merged PDF
is our tectonic build of the same TeX (4 pages, 236,450 bytes), uploaded through the portal's "replace" on the
merged-PDF row. Re-reading the TeX reset two form fields, which were put back: the abstract (restored to the
143-word plain-text version already on file, identical to `abstract_plaintext.txt`) and the author (Find Person, the
account with "Independent researcher"). Stopped before "Approve Manuscript" for Jon. The stray new submission started
earlier (no manuscript number) was deleted with Jon's OK.

**Page-1 fix prepared (2026-11-25), not uploaded.** The 2026-11-25 audit (`rnaas/audit-2026-11-25.md`, branch
`rnaas-audit-2026-11-25`) found that the PDF uploaded from 4438ecf typesets the opening paragraph above the title:
the length cut removed `\section{Introduction}`, and AASTeX 7 sets the title block at the first `\section`. Fix on
branch `eb-aas84041-fix`: the heading is back, and "As in the companion note," is dropped (the K2 note no longer
lists these depths). No number changed. Gate (d385d67): abstract 142, main text 1294, AAS range 1285-1322; page 1
checked by eye. `rnaas.pdf` is the tectonic build of this TeX (4 pages). To reach AAS it needs the manuscript back
in Author Approval.

**Page-1 fix uploaded (2026-11-25), awaiting Jon's approval.** AAS staff returned AAS84041 to Author Approval after
our correspondence. Uploaded `McBee_AAS84041_R1.tex` = c5d040a `rnaas.tex`, replacing the manuscript file; figure
and class file unchanged. Gate (d385d67) on the upload set: abstract 142/150, main text with references 1294/1300,
AAS range 1285-1322, page 1 checked by eye. The merged PDF is our tectonic build of that TeX (4 pages, 236,791 bytes,
SHA-256 666b4bbd…69566), uploaded with "replace" on the merged-PDF row; the portal's stored copy has the same hash.
Abstract (143-word plain text) and author (Find Person, "Independent researcher") restored after the re-upload;
dataset DOIs match the four in the TeX. Stopped before "Submit Manuscript" for Jon. Open point: the Introduction still
cites "A companion note on K2 single events (J. McBee 2026, submitted)", and that note (AAS84033) was withdrawn on
2026-11-25.
