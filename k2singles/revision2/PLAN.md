# AAS82225 revision: plan and pre-registered rules (written 2026-11-15, before any search result)

## Why

On 2026-11-15 Chris Lintott declined AAS82225 and invited a revision that does two things:

1. a thorough search of the existing literature for every candidate;
2. a blind search for transits, rather than one for a particular set of properties.

The submitted note slid a box with the catalogued depth and duration along the data. It needed a depth of
0.5–2 times the catalogued value and searched only the 144 events with an expected TESS SNR ≥ 7.
"Blind" here means three things. No sample cut uses the catalogued depth, duration or expected SNR. No
detection uses the catalogued depth, duration or epoch. Detection thresholds come from an inverted-flux
control, not from the catalogued events.

The catalogued values are used only after detection, to ask whether the blind search recovered the
catalogued event and whether a blind period links to it.

## Sample

All 164 events on 163 stars in LaCourse & Jacobs (2018), `k2singles/data/k2_singles.csv`. There is no cut.
Stars without data outside the event campaign are reported as such.

## Data (`fetch_direct.py`)

- **TESS:** every SPOC 2-min, TESS-SPOC and QLP light curve. Our own TESScut photometry for sectors ≥ 80 that
  have no pipeline product, as in the submitted note. Sectors above 200 are skipped, as in the note's
  pipeline. The note took the 2026 January special pointing ("S1751") from TESScut only for EPIC 202126877.
- **K2:** the mission PDCSAP, EVEREST v2 and K2SFF light curves for every campaign that observed the star,
  not only the event campaign.

## Search 1: blind single-event search (`search_blind.py`)

- **Code and rules:** identical to the unconstrained dip search of the AAS82207 revision
  (`singles/revision_search/search_any.py`, imported, not copied):
  - local quadratic baseline;
  - red-noise-inflated SNR ≥ 8;
  - coverage ≥ 0.8 and n_in ≥ 3;
  - baseline on both sides;
  - ≥ 0.25 d from segment edges;
  - the v2 momentum-dump rule;
  - median depth ≥ 0.5 × mean depth;
  - depth ≥ 3 × the 1-h scatter;
  - overlapping detections merged.
- **Durations:** 1, 2, 4, 8, 16, 24 and 48 h. The catalogued duration is NOT added. Events longer than about
  48 h can be found only in K2, because a TESS orbit is too short for their baselines.
- **TESS:** the best product per sector is searched, and every detection is re-measured in each other product.
  This is unchanged from search_any.
- **K2:** each of EVEREST, K2SFF and PDCSAP is searched separately, because each one removes some real
  eclipses (PDC blanks eclipse cores, EVEREST can flag deep events as outliers, K2SFF has gaps).
  Detections are merged across products, and each is re-measured in the other two.
- **Control:** the identical search on 2 − flux for every star and product. False-alarm rates are quoted
  from it.
- **Common mode:** an event within 0.5 d of a detection on ≥ 2 other stars, in the same K2 campaign or
  TESS sector, is flagged as common-mode. It is reported but not used to link periods.

## Search 2: blind periodic search (`search_blind.py`, second stage)

**(a) Linking events.**
- **Families:** for each star, events with fitted depths within a factor 2 and durations within 1/2–2×
  are grouped into a family.
- **Candidate periods:** for every pair in a family, P = ΔT/n for P ≥ 1 d.
- **Support:** every family event must fall on the ephemeris within max(0.1 d, D/4).
- **Exclusion:** a candidate is excluded if any predicted epoch with coverage ≥ 0.8 is flat, meaning a
  measured depth below 0.5 × family depth − 3σ. This is the rule tested on 83 confirmed planets in
  `salinas-second-transits/results/round2/EPOCH_RULE_TEST_2026-11-01.md`.
- **Reported:** periods supported by ≥ 3 events (a periodicity found blindly) and, separately, two-event
  alias lists.

**(b) BLS** (astropy BoxLeastSquares).
- **Input:** the combined K2 + TESS light curve of each star, each segment divided by a running median of
  window 3 d.
- **Grid:** periods 1–100 d; durations 1, 2, 4, 8 and 16 h.
- **Detection:** SDE ≥ the 99th percentile of the maximum SDE from the same BLS on the inverted light curves
  of all 163 stars. The threshold is fixed from the control before the real peaks are looked at.
- **Iterations:** up to three signals per star, masking each in turn.
- **Changed 2026-11-15, before any real BLS result:** BLS runs per data block (data separated by gaps over
  60 d), with periods from 1 d to min(100 d, block span / 2). A single BLS over the combined 12-year baseline
  would need about 4e8 trial periods. Linking across blocks is left to Search 2a. The threshold is unchanged:
  the 99th percentile of the inverted blocks' maximum SDE.

## Comparison (`compare.py`)

- Was the catalogued event recovered blindly, and in which products?
- For each of the 14 systems in the submitted Table 1, does the blind search give the same period or alias
  list?
- Which periodic signals or linked periods are new, meaning absent from the submitted note's targeted
  search? Each new one goes through the note's vetting: a period fit, the products, TESS difference images,
  Gaia and an independent re-derivation before any claim.

## Literature (`literature.py`, needs an ADS API token in ADS_DEV_KEY)

- **Queries:** every star under every identifier (EPIC, TIC, Gaia DR3, 2MASS, HD/TYC/BD where present),
  with ADS full-text and object queries. SIMBAD bibliography counts and VizieR catalogue hits as before.
- **Classification:** each hit is marked as giving a period for this object, mentioning the object, or
  irrelevant.
- **Output:** a per-star table that the revised note summarises: how many of our periods are already
  published, and where.

## Things decided after first results

Any change goes below with its date and reason.

**2026-11-15, after a one-star test (EPIC 248811085) and before the full run: significance from the control.**

- **What the test showed:** in K2 the search found 12 events at SNR 8–37 besides the catalogued eclipse (SNR 529).
  The inverted control found 17 events at SNR 8–70. Most were 200–900 ppm and 15–30 h long, or sharp PDCSAP
  features. A fixed SNR ≥ 8 is therefore not a meaningful detection threshold in K2. The TESS events (SNR 222
  and 399) had no control counterpart.
- **Eligible events:** events and control events are eligible if they are confirmed in a second product of the
  same campaign or sector (`n_reductions ≥ 2`), or if no second product covers them (`n_alt_tested = 0`).
- **Threshold:** for each mission separately (K2, TESS), S* is the smallest SNR at which eligible control events
  number at most 0.05 per searched star, pooled over all stars.
- **Significant events:** an event is significant if it is eligible and has SNR ≥ S*. Only significant,
  non-common-mode events enter the period linking (Search 2a).
- **Reporting:** events with SNR ≥ 8 that fail this test are listed as low-significance and not used.
- **Unchanged:** SNR ≥ 8 and the other search_any rules stay as the first-stage cut; nothing about them changes.

**2026-11-15, AFTER the full pre-registered run: a post-hoc linked-period test (`link_calibrated.py`).**

- **What the full run showed:**
  - S* came out at SNR 53.3 for K2 and 60.0 for TESS.
  - 162 of 163 stars had at least one catalogued event recovered blindly; the miss, EPIC 211939692, has
    two catalogued events and one was found. Of these, 79 had it at a significant level.
  - 22 stars had a linked family and 35 a BLS detection. 11 of the 15 Table 1 stars got their table period
    among the blind allowed periods.
- **What set the high threshold:** the control events above SNR 45 are large brightenings, such as bad
  cadences, glints and scattered light, up to 150 times the stellar flux. Inverted, they become "dips". They
  have no dimming counterpart in real data, so the inverted control overstates single false dips.
- **What the blind first stage still found:** all four TOI-5560.01 transits (SNR 9–32), all three events of
  EPIC 211503363 and of Sanders 617, and the primaries of EPIC 211489484. Only the single-event threshold
  dropped them. For EPIC 211489484, the greedy family rule also merged primaries and secondaries.
- **Post-hoc test:** link all eligible, non-common-mode first-stage events pairwise. Calibrate the number of
  events on one ephemeris, `max_support`, against the same linking of the control events. N* is the smallest
  support reached by at most 0.05 control stars per star. This N* rule was fixed before the test ran.
- **How it is reported:** as post-hoc and secondary. The pre-registered result stays the primary one, and
  both go in the revision.

**2026-11-15, after the full run: the common-mode rule is replaced by a coincidence test.**

- **The problem:** the rule as written flags an event if ≥ 2 other stars have any first-stage event within
  0.5 d in the same campaign or sector. It came from the AAS82207 search, where events were sparse. Here a K2
  star has about 15 first-stage events per campaign, most of them shallow systematics, so chance coincidences
  are common. The rule flagged 39% of K2 events, 18% of TESS events, and 46 of the 163 recovered catalogued
  events. Those are real astrophysical events, including the C5 eclipse of Sanders 617 at SNR 277.
- **New rule:** count k, the number of other stars with an event within 0.5 d. The expected count is
  λ = Σ over other stars in the campaign or sector of n_events × 1 d / span, where span is that campaign's or
  sector's event time range. The event is common-mode if k ≥ 2 and the Poisson probability P(K ≥ k | λ) is
  below 0.001.
- **What it applies to:** both the pre-registered stage (significance) and the post-hoc linking.
- **Test before use:** the share of catalogued events flagged must be small. If more than 5% are flagged, the
  rule is reported as failed and not used.
- **Test result (2026-11-15):** the rule flags 1 of 163 recovered catalogued events, 0.2% of K2 events, 0% of TESS
  events and 36 control events. It passes.

## Run record, 2026-11-22

- **Saved stage-2 outputs are current.** `check_outputs_current.py` recomputes thresholds, common-mode flags and
  significant events from `data/per_star/`. All 163 stars match `data/periodic/` under the Poisson rule (3 common-mode
  events). The old rule would have flagged 522. Under it, 79 catalogued events were significant, the figure in the
  2026-11-15 entry above. Under the current rule it is 112, the figure in `compare_run.log`. `compare.py` reran with
  byte-identical output.
- **Post-hoc linked-period test, completed.** The 2026-11-15 run hung at 56 of 163 stars: the pool never replaced
  workers retired by `max_tasks_per_child`. That setting is removed (cc98328) and the test reran in 46 s.
  - N* = 4. 7 control stars reach it, against 8 allowed.
  - 11 stars are detected (`results/link_calibrated.csv`).
- **Table 1 (15 stars), primary and post-hoc together:**
  - The pre-registered search recovers 11.
  - The post-hoc test adds 211489484 (82.258 d, support 4, control 0) and 211924561 (TOI-5560.01, 140.49/70.25 d,
    support 4). For 211924561 the star's own control reaches 6, so treat it as weak.
  - Neither recovers 211503363 (support 3) or 246194159 (support 2).
- **Not yet done:** classifying the literature hits (`results/literature/hits.csv`; its `period_words` column is False
  for all 919 hits, which needs checking); vetting the new leads (`results/new_leads.csv`, 41 stars; post-hoc
  detections not in Table 1: 202085278, 204776782, 247795097, 215067200, 205966706, 247756662, 211598816, 212805678).
  Several have control support close to their real support (204776782 12 vs 11; 205966706 6 vs 10).

## Literature classification: rules (written 2026-11-22, before reading the rerun's snippets or the VizieR census)

**Fixes before the rerun:**
- `literature.py` looked up ADS highlights by bibcode, but ADS keys them by document id, so every snippet was empty (798e8b5).
- `vizier_all.py` had not been run. It now keeps Gaia DR3 variability (I/358), which the `I/` skip had dropped.
- The 2026-11-15 outputs are kept as `data/literature/*_2026-11-15_nosnippets*`.

**Recall check.** Before any classification, count how many of these known (star, paper) pairs each source finds. All of
them come from the submitted note:

| Pair | Paper |
|---|---|
| all 163 stars | LaCourse & Jacobs 2018 |
| 211411112 | Sandquist+2021 |
| 202126877, 202072917 | IJspeert+2024 |
| 202072917, 211489484 | Valizadegan+2025 |
| 202071902 | Kovalev+2024, Guo+2025 |
| 211503363 | Kovács 2020 |
| 211892898 | Libralato+2016 |

Report recall per source (ADS full text, SIMBAD, VizieR). The revised note states the search's known gaps; nothing is
tuned to raise the recall.

**Classes, one per (star, paper or catalogue):**

| Class | Meaning |
|---|---|
| PERIOD | Gives an orbital period, ephemeris or alias set for this object |
| ESTIMATE | A period estimate or limit without repeats, e.g. from transit shape and stellar density |
| ANALYSIS | Studies this object (spectra, RVs, characterisation) without a period |
| LISTED | Only a row in a survey or catalogue, with no period |
| NOMATCH | The identifier matched a different object |

**How a class is decided:**
- Title, abstract and ADS snippets first.
- Read the paper, or its table, if any of these hold:
  - the abstract or snippets mention a period, ephemeris, orbit, eclipse or transit for this object;
  - the title is about binaries, variables, planet candidates or periods;
  - a VizieR row has a period column.

  Then record the period value and its source.
- Large stellar-parameter and astrometric surveys go to LISTED without being read: APOGEE, GALAH, LAMOST, SDSS data
  releases, 2MASS, EPIC (Huber+2016), JSDC, HD, Tycho, UCAC and the Gaia main-catalogue papers. Gaia DR3 non-single-star
  and variability papers are NOT on this list: they can carry periods, so they are read.

**Per star:**

| Verdict | Meaning |
|---|---|
| published, agrees | A PERIOD source agrees with ours within 3σ, or our period is one of its aliases |
| published, differs | A PERIOD source gives a different period |
| estimate only | Only ESTIMATE sources |
| no prior period | No PERIOD or ESTIMATE source |

In the revised note a period is called new only if its star is "no prior period" or "estimate only".

**Independent check (Table 1 stars).** A second agent classifies the same hits from the same inputs. It does not see
these rules' outcomes or my classes. Disagreements are settled by reading the paper, and both classes are recorded.

## Literature search: run record, 2026-11-22

**ADS + SIMBAD rerun** (`literature.py`, snippets fixed): 922 (star, paper) rows on 97 of 163 stars (635 SIMBAD only,
219 ADS only, 68 both); 181 rows carry ADS snippets, 32 with period words. 66 stars have no paper at all.

**Recall check** (`recall_check.py` -> `data/literature/recall.csv`), ADS and SIMBAD:
- LaCourse & Jacobs 2018: 0 of 163 by either source. Its star list is not in ADS full text and SIMBAD does not link it.
- Other known pairs: 3 of 9 (Sandquist+2021, Kovács 2020, Libralato+2016, each found by both sources). Missed: IJspeert+2024
  (2), Valizadegan+2025 (2), Kovalev+2024, Guo+2025. All six name the star only in a machine-readable table.
- VizieR: pending (below).

**Fix: Gaia DR3 orbits were not covered.** `vizier_all.py` skips the `I/` section and keeps only I/358, so it never saw
I/357 (Gaia DR3 non-single stars), which the rules above say must be read. `gaia_nss.py` queries the Gaia archive directly
for all 163 stars (source_id from SIMBAD for 96, else the brightest Gaia DR3 source within 2" for 67). Orbits found:
211892898 SB1 72.534(29) d; 247756662 SB1 610.7(5.6) d; 229021605 SB1 114.838(66) d (plus a weaker 115.35 d validated
targeted-search orbit); 212026444 an acceleration solution (no period). No star has a `vari_eclipsing_binary` entry.

**Fix: the VizieR census cached failures as empty results.** The first run (4 threads, CfA mirror) failed outright for
42 of 112 stars (HTTP 504 at the mirror's 60 s gateway limit), and 66 of the 70 it cached were empty although each star is
in Gaia DR3 and 2MASS there (a targeted query finds them). That run is kept as `*_2026-11-22_emptyfail`. A response
without I/355 now counts as a failure; rerunning with 2 threads. The CDS host and the China, IUCAA, NAOJ and INASAN
mirrors did not answer; IDIA answered but timed out on the same query.

**Table 1 classification, two independent passes** (A = me; B = an agent that saw only `table1_hits.csv` and the class
definitions). 97 rows on the 10 Table 1 stars with any hit. `compare_classes.py`: 89 agree, 7 differ only LISTED vs
ANALYSIS (no period either way), 1 differs in class: Wang+2021 (RAA 21, 292) for 202071902. B was right: the online
Table 3 lists it as EA with P = 8.6973 d from AAVSO VSX; A had searched the PDF, which prints only a sample of the table
(`settle_classes.py`). Final: 84 LISTED, 5 ANALYSIS, 5 PERIOD, 2 ESTIMATE, 1 NOMATCH.

Also checked by reading: the Villanova K2 EB catalogue gives 211924561 and 211489484 "period 100, error -1", its
placeholder for single eclipses; Zhang+2019 copies the 100.0 d. Not a measurement (both passes agree).

**Per-star verdicts, Table 1, from what the search found** (VizieR still pending):

| EPIC | Verdict | Sources |
|---|---|---|
| 211411112 | published, agrees | Sandquist+2021 62.59484(3) (ours 62.594821(27), 0.5σ); Geller+2021 SB2 62.593(3) (0.6σ) |
| 211892898 | published, agrees | Gaia DR3 SB1 72.534(29) (ours 72.5089(14), 0.9σ); Libralato+2016 P > 50 d |
| 211503363 | estimate only | Kovács 2020: single transit, P ≳ 58 d (ours 122.905 or 61.453 d) |
| 202071902 | published, differs | Wang+2021 / VSX 8.6973 d, which the note already rules out. Guo+2025 and Kovalev+2024 (agree) were missed by the search |
| 202126877, 202072917 | no prior period found | IJspeert+2024 and Valizadegan+2025 (agree) were missed by the search |
| 211489484 | no prior period found | Kruse+2019 lists two eclipses, no period. The SPOC TCE (Valizadegan+2025, 106.94 d, which the note excludes) was missed |
| 211924561 | no prior period found | TOI-5560; the ExoFOP limit ≤ 702 d is not from the search |
| 246194159, 211390677 | no prior period | only survey rows (one NOMATCH) |
| 248811085, 246089124, 246696483, 247967714, 212012030 | no paper found | |

So the search alone recovers the prior period of 2 of the 5 Table 1 stars known to have one (211411112, 211892898), and
misses the 3 C0 binaries' catalogue periods. VizieR is the source that should find those.
