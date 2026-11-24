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

**2026-11-22, after the VizieR census failed: a fixed catalogue list (Jon chose "both").**

- **The problem:** the pre-registered all-catalogue census (`vizier_all.py`) cannot run. The CfA mirror's catalogue
  database (`metaviz@CFA`) was down most of 2026-11-22/23. It answers HTTP 200 with "The database is not currently
  reachable", which astroquery returns as an empty result. When the database is up, one star's all-catalogue reply
  takes several minutes. The CDS host and the other mirrors did not answer.
- **Decision:** `vizier_targeted.py` queries a fixed list of catalogues, chosen before any was queried:
  - the recall papers' catalogues, plus Wang+2021;
  - variable-star and period catalogues: VSX, Gaia DR3 variability and NSS, ASAS-SN, ATLAS, ZTF, WISE, TESS EBs;
  - 11 K2 variability, EB, rotation and planet-candidate catalogues;
  - for Table 1 stars only, the VizieR catalogue of every paper the ADS/SIMBAD search found.

  Same 3" radius. `vizier_all.py` keeps running serially; if it finishes, it is the primary result.
- **Consequence for recall:** the list holds the recall papers' catalogues by construction, so finding them measures
  only whether the query works, not how a blind search would do. The revised note must say so.

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

**Per-star verdicts, Table 1, from what the search found** (before VizieR; superseded below):

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

## Literature search: VizieR, 2026-11-23

**All-catalogue census** (`vizier_all.py`):
- Second run: 160 of 163 failed the I/355 check (database down).
- Third run, serial, 8 tries per star: 12 of 12 failed in about 12 h. Still running.
- Usable results: 3 stars. VizieR recall from the census is not measured.

**Fixed-list query** (`vizier_targeted.py`, post-hoc, see "Things decided after first results"):
- Catalogue IDs: 93 checked. 78 exist on the mirror, 13 are not found, and 2 list no tables in a cone search (no
  positions).
- Not found, among those that matter:
  - recall papers: Sandquist+2021 (J/AJ/161/59) and Kovács 2020 (J/A+A/643/A169);
  - general list: Zink+2021 (J/AJ/161/259), Pope+2016 (J/MNRAS/461/3399) and Armstrong+2016 (J/MNRAS/456/2260).
    They were not searched.
- Result: all 163 stars queried; 296 table matches on 95 stars (`data/literature/vizier_targeted.csv`).

**Fixes during the run:**
- A Postgres error in a "Table or Catalog not found" reply had been read as an outage.
- `-out.all=` (empty) returns default columns only; it is now `-out.all=1`.
- J/A+AS/100/243 answers with no tables; that now has its own status.
- EPIC 211498244 moves 346 mas/yr, and its Gaia DR3 position is 5.2" from J2000. So a reply now passes on Gaia DR3 *or*
  2MASS being present, in both scripts.

**Known gaps:**
- **Dropped tables.** The mirror sometimes drops single catalogues from an otherwise complete reply. One query for
  202071902 returned nothing; on the Table 1 refetch, I/357/tbosb1, a VSX row and II/361 were missing until retried.
  The 163-star run checks only that the Gaia DR3 or 2MASS control row is present, so it may have missed matches.
  `table1_vizier_rows.py` retries until every table the main run found comes back, but it cannot recover a table
  that the main run itself lost.
- **High proper motion.** A 3" cone at J2000 can miss catalogues that give positions at Gaia's epoch for fast-moving
  stars. Not yet counted: the Gaia archive was in downtime on 2026-11-23.

**Fixed-list "recall"** (not recall; `recall.csv`, column `vizier_fixed_list`):
- Found: 6 of the 7 pairs whose catalogue is on VizieR (IJspeert ×2, Valizadegan ×2, Kovalev, Guo).
- Missed: Libralato+2016 for 211892898. Its table3 lists the star at 4.35" (blank 2MASS magnitudes, an instrumental
  Kp; apparently their position for a saturated star), outside the pre-registered 3".

**Table 1 classification, two independent passes on the VizieR rows:**
- Input: `table1_vizier_rows.py` -> `table1_vizier_rows.txt` (130 rows, 106 star–table pairs, 11 stars).
- B saw only that file and the class definitions.
- `compare_classes.py --vizier`: 105 agree, including every PERIOD pair and its values. 1 differs LISTED vs ANALYSIS
  (Douglas+2017 table3 for 211892898; no period either way).
- Final (A): 83 LISTED, 15 PERIOD, 8 ANALYSIS (`table1_vizier_classes_final.csv`).
- Both passes class SPOC TCE periods (Valizadegan+2025) as PERIOD, kind TCE.

**Per-star verdicts, Table 1, search plus fixed-list VizieR:**

| EPIC | Ours (d) | Verdict | Sources |
|---|---|---|---|
| 211411112 | 62.594821(27) | published, agrees | Sandquist+2021 (0.5σ); Geller+2021 SB2 62.593(3) (0.6σ) |
| 211892898 | 72.5089(14) | published, agrees | Gaia DR3 SB1 72.534(29) (0.9σ). Also: Rizzuto+2017 Prot 8.51 (rotation; cluster parallax assumed, the star is a background giant) |
| 202071902 | 22.904128(3) | published, agrees | Guo+2025 22.89(40) (0.04σ); Kovalev+2024 Keplerian 22.859(15) (3.0σ). Differs: Wang+2021/VSX 8.6973 (EA), which the note rules out |
| 202072917 | 34.745586(7) | published, agrees | Valizadegan+2025 TCE 34.74574 (S1-46); IJspeert+2024 QLP 34.74438. No σ published; differences 1.5e-4 and 1.2e-3 d. Differs: VSX 3.842, IJspeert SPOC 23.255, other TCEs (34.742, 1007.6, 15.17, 14.20, 5.579) |
| 202126877 | 21.77465(22) | published, agrees | IJspeert+2024 SPOC 21.7808, QLP 21.7755 (no σ published). Differs: VSX 3.877 |
| 211503363 | 122.905 or 61.453 | estimate only | Kovács 2020 P ≳ 58 d. Also: Zhang+2020 Prot 14.05 and 12.47 (rotation) |
| 211489484 | 82.258438(52) | published, differs | Valizadegan+2025 TCEs 106.936 (the note excludes it) and 0.2405, both ExoMiner FP |
| 211924561 | 140.49392(14) or 70.24696(7) | published, differs | Valizadegan+2025 SPOC S42-46 TCEs 58.290 (MES 11.5) and 31.269 (MES 7.4), both ExoMiner FP. See below |
| 212012030 | 220.32 | no prior period | ATLAS lists a "dubious" variable, periodogram peaks 0.721 / 1.258 / 6.794 d |
| 246194159, 211390677 | | no prior period | survey rows only |
| 248811085, 246089124, 246696483, 247967714 | | no paper or catalogue found | |

**211924561 (TOI-5560) TCEs.** Checked 2026-11-23 against the SPOC DV reports (XML, from MAST with Jon's OK:
`tess2021233042500-s0042-s0046-…-00550_dvr.xml`, `tess2023315161428-s0072-s0072-…-00827_dvr.xml`,
`tess2018235142541-s0002-s0072-…-00829_dvr.xml`; `dv_tce_check.py`, output `results/literature/dv_tce_check.txt`).
SPOC has four TCEs on this star, not two. The S2–S72 and S72 runs were not in the literature search: they are MAST
data products, not papers or catalogues.

| SPOC run | P (d) | MES | Its transit times in the run's data |
|---|---|---|---|
| S42–46 | 58.283 | 11.5 | S44 2510.01 (rejected scattered-light step); S46 2568.30 (our transit) |
| S42–46 | 31.253 | 7.4 | S44 2509.80 (the same step); S45 2541.06 (vetted, < 0.5 ppt on a background ramp); S46 2572.31 (nothing in our search) |
| S2–S72 | 58.537 | 14.6 | S44 2509.77 (the step); our S46 and S72 transits. 58.54 ≈ 702.47/12, the S46–S72 gap |
| S72 | 10.285 | 7.8 | 3260.48 (0.3 d after the first S72 cadence; nothing in our search); our S72 transit; 3281.05 (not seen by SPOC either) |

- None of the four periods fits all four of our transits (C5, C16, S46, S72). The C5 transit falls 1.2–14.1 d from
  each TCE's ephemeris.
- Every TCE leans on the rejected S44 step or on a time where our search found nothing. Three of the four also use
  one or both of our real TESS transits.
- Correction to the 2026-11-23 note above: the 58.29 d TCE does not pair two rejected events. It pairs the rejected
  S44 step with our real S46 transit. The 31.25 d TCE uses no transit of ours.

**211489484 TCEs.** Checked 2026-11-23 against its SPOC DV report (XML, from MAST with Jon's OK:
`tess2018206190142-s0001-s0046-0000000006951741-00555_dvr.xml`; same script and output). One run (S1–S46), two TCEs:

| SPOC run | P (d) | MES | Its transit times |
|---|---|---|---|
| S1–S46 | 106.936 | 31.8 | our S7 and S46 primary eclipses, ten cycles apart (1069.36 d = 10 × 106.936 = 13 × 82.258) |
| S1–S46 | 0.2405 | 8.5 | 363 of 4520 expected events, 0.5 h pulse |

- 106.936 d is the 10-cycle alias of the S7–S46 pair. Our C5, C18, S34 and S71 primaries fall 8.2–32.9 d from its
  ephemeris, so it is excluded.
- 0.2405 d cannot be the eclipse period: our 4-h eclipses would then recur every 5.8 h, and our searches find
  isolated events only.
- **Cause of the 0.2405 d TCE** (`tce024_check.py`, output `results/literature/tce024_check.txt`; SPOC 2-min light
  curves of S7, S44, S45, S46 downloaded with Jon's OK; S34 was not in this SPOC run):
  - The target's light curve, eclipses masked, has a δ Scuti-type pulsation at 20.786 c/d (P = 0.048109 d,
    semi-amplitude 3.2–7.8 ppt) in every sector. 0.2405 d is five pulsation cycles (5 × f_TCE = 20.789 c/d): the TPS
    fold stacks every fifth minimum. Folded at 0.2405 d the light curve shows five ripples, not one dip.
  - The DV difference images put the source 24.5 ± 3.5" from the target (S44–46), at RA +2.6", Dec −24.4".
    Gaia DR3 651067327828913280 is at +2.4", −26.0" (26.3"), G = 15.81, BP−RP = 0.42. Gaia DR3 classifies it
    DSCT|GDOR|SXPHE (score 0.89) and VSX lists it the same. ATLAS gives 0.0506 d, not the TESS 0.0481 d (ATLAS
    sampling is sparse; not resolved).
  - So the 0.2405 d TCE is a subharmonic of a pulsating neighbour 26" away, unrelated to the eclipses.

**Decided 2026-11-23 (Jon): TCE periods do block "new".** The rules call a period new only if its star is "no prior
period" or "estimate only", and both classification passes class the TCEs as PERIOD. So the rule stays as
pre-registered, and the note does not call the periods of 211489484 and 211924561 new.
- The note presents each as a resolved period that conflicts with the SPOC TCEs, and says why each TCE is excluded.
- 211489484, 106.936 d: the 10-cycle alias of our S7 and S46 primaries, excluded by the other four (table above).
- 211489484, 0.2405 d: a sub-day signal, not a competing eclipse period. One line saying so.
- 211924561: all four SPOC TCEs (58.28, 31.25, 58.54, 10.28 d) are checked against the DV epochs (table above).
  Each combines our TESS transits with a rejected event or an empty time, and none fits the K2 transits.
