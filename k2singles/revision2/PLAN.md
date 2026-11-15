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

None yet. Any change goes below with its date and reason, and detection thresholds are not changed.
