# Hand vetting of three "repeats-inconsistent" stars with a hidden period (2026-10-12)

**Why these three.** `singles/results/period_recheck/` found 6 of the 23 repeats-inconsistent stars
have a period that explains all but one repeat. These three were chosen for vetting because they are
on target (or inconclusive), Gaia-clean, and have pipeline-product repeats. The question per star: is
the left-out repeat spurious, so that the rest share one ephemeris?

**Script.** `singles/scripts/vet_hidden_periods.py`. It wraps `scripts/vet_candidate.py` (joint fit,
per-product depths, centroids, dumps/edges, Gaia neighbours) and adds:
- per-event trapezoid fits in every product;
- our own TESScut photometry for every event, including a version that keeps stray-light-flagged
  cadences (`OWNSL`; its background estimate is independent of SPOC's);
- quality flags at each event;
- cross-star timing coincidences;
- predicted transits.

Own photometry is regenerated into `own/`, which is not committed: it is rebuilt from the cached
cutouts.

**Outputs.**
- Per-event fits: `events.csv`.
- Flags: `quality_flags.csv`.
- Joint fits: `summary.csv` and `<TIC>_<label>.{json,png,_epochs.csv}`.
- Wide-window plots: `<TIC>_events_wide.png`.
- Predictions: `140215502_P297_predictions.csv`.

No repeat of these stars coincides within 0.25 d with a repeat on any other singles star
(`coincidences.csv`).

## Verdicts

| TIC | Left-out repeat | Verdict | Period |
|---|---|---|---|
| 140215502 | S5 (TESS-SPOC) | **Tentative lead.** S32 + S87 + S107 lie on one ephemeris (O−C 10 min; chance alignment ≈ 0.3–0.6%), with consistent depth and duration. S5 is compromised and not on it | 297.32 / m, m = 1, 2, 3, 4, 6: **297.3, 148.7, 99.1, 74.3, 49.6 d** |
| 162031311 | S83 (QLP) | **Duotransit, S83 likely spurious.** S16 matches the event in depth in three reductions. S83 is 1.8× too deep, starts at a dump after a gap, and has a similar dip 25 h earlier | 1084.31 / n: **24 aliases, 24.1–1084.3 d** |
| 76140400 | — | **No credible signal.** The catalogued S26 event appears only in SPOC PDC. The S40 dips are stray-light-flagged, V-shaped, show centroid shifts, and their depth doubles in our photometry | none |

None of this changes a pipeline class. MANUAL notes recording these verdicts are in
`classify_final.py`; MANUAL notes never change a class.

---

## TIC 140215502 (CTOI 140215502.01, steuer 2021; Tmag 11.2, Teff 5528 K, R* 0.93 R☉)

**Events.** Trapezoid fits, `events.csv`:

| Event | Product | Depth (ppt) | T14 (h) | Flags in transit |
|---|---|---|---|---|
| S32 (CTOI event) | SPOC 2-min / TESS-SPOC / QLP / own | 4.45 / 4.35 / 4.11 / 4.31 (± 0.15–0.17) | 7.22–7.30 | none |
| S87 | QLP / own | 3.76 ± 0.18 / 4.04 ± 0.16 | 8.2 / 7.6 | none |
| S107 | own only | 4.99 ± 0.19 | 7.5 | none |
| **S5** | TESS-SPOC / SAP / own | **2.79 ± 0.19 / 2.99 / 3.25 ± 0.16** | **6.2–6.7** | **momentum dump at +1.1 h (mid-transit); QLP flags 9 of 13 in-transit cadences (scattered light, bit 4096)** |

- **S32, S87 and S107 share one ephemeris.** Times from S32 and S87 predict S107 to within
  10 min (P = 297.3249 d joint fit).
- **That alignment is unlikely by chance.** Without S107/S108 data, 22 of 297 S32+S87 aliases
  remain allowed. Only 0.3% of S107+S108 cadences lie within ±15 min of one of their epochs (0.6%
  within ±30 min).
- **Depths and durations agree** within ~1.2× and 0.5 h, comparing like products: QLP reads ~0.92×
  SPOC on the S32 event.
- **S5 is not on this ephemeris** for any surviving alias. It is 0.64–0.75× the S32 depth (≥ 4.8σ in
  own photometry), slightly shorter, and has a dump mid-transit. It is either a systematic or a
  second body; we cannot tell which.
- **Source and host are clean:**
  - difference images on target: S32 2.1″, S87 1.5″, S107 5.4″ ± 1.6, S5 3.2″;
  - no Gaia neighbour within 60″ could host the signal;
  - Gaia RUWE 1.10, not RV-variable, no NSS solution.
- **Rp ≈ 6.7 R⊕.** Central-transit durations for this star run from ≈ 6.3 h (49.6 d) to ≈ 11.4 h
  (297 d). The observed 7.3 h fits every alias for some impact parameter.

**Prior follow-up: TFOP has worked this star, on a different pairing.** ExoFOP shows TESS Transit
Finder targets (Ulmer-Moll) treating it as an S5 + S32 duotransit with a 146.78 d alias
(= 733.89 / 5). A TFOP SG1 campaign on 2024-12-17/18 (15 observations, BJD 2460661.6–2460662.8) gave
mutually inconsistent results:
- three weather-truncated CTIO runs report a "possible" ingress near 2460661.67;
- PEST reports a "likely" ingress 105 min late (Tc 2460662.28);
- one CTIO run reports a "possible" egress at 2460662.775;
- the rest are non-detections.

No period has been posted. Our S32+S87+S107 ephemeris predicts no transit in that window (nearest
epochs BTJD 3631.5 and 3681.06). TESS S87 started 1.4 d after it.

**Discriminating epoch.** On **2026-12-21 14:21 UTC (BTJD 4396.10)**, every S5 + S32 alias (733.9/n,
n = 1, 3, 5, 7) predicts a transit and no S32+S87+S107 alias does.

**Next S32+S87+S107 transits** (`140215502_P297_predictions.csv`). Formal errors are ±5 min; inflate
them about 3×.

| Date and time (UTC) | Aliases predicting a transit |
|---|---|
| 2026-11-05 13:02 | 74.3 d |
| 2026-11-30 07:41 | 99.1 and 49.6 d |
| 2027-01-18 20:59 | 148.7, 74.3 and 49.6 d |
| **2027-06-16 12:53** | **all five** |

**Status: tentative.** The third transit is only in our own S107 photometry. MAST has no S107
pipeline product yet; checked 2026-10-12: S87 and S98 QLP are the latest. Per the S106–108 lesson,
it stays tentative until one appears.

Before any outreach:
- an independent (blind) re-derivation (CLAUDE.md §5);
- a literature check.

The natural recipients are the TTF/TFOP coordinators, whose 146.78 d hypothesis this would supersede.

## TIC 162031311 (CTOI 162031311.01, UBC planet search 2025, "in prep"; Tmag 11.2, Teff 6160 K, R* 1.57 R☉)

| Event | Product | Depth (ppt) | T14 (h) | Flags in transit |
|---|---|---|---|---|
| S56 (CTOI event) | QLP / own | 2.68 ± 0.20 / 2.31 ± 0.15 | 6.75 | none |
| S16 | TESS-SPOC / SAP / own with stray light kept | 2.75 ± 0.24 / 2.50 / 2.76 ± 0.28 | 8.6–9.2 (± 0.9) | all 13 FFI cadences stray-light (2048); QLP and default own photometry drop them |
| **S83** | QLP / own | **4.69 ± 0.19 / 4.15 ± 0.13** | 7.4–8.0 | none, but data resume after a gap with momentum dumps at −3.6 h (ingress unobserved) |

**S16 is real and matches the event.**
- Joint-fit scales are 1.00 and 1.00 (χ² 0.0/1).
- The dip appears in three reductions. One of them, our own with an independent background, rules
  out a SPOC background-subtraction artefact.
- Centroid shifts are 0.03σ and 0.08σ.
- It is about 2σ longer than the event.

**S83 is probably not a transit.**
- It is 1.8× the event depth (χ² 84/1).
- It begins right after a gap at a dump.
- QLP and own photometry both show a similar unflagged dip about 25 h earlier, which points to
  systematics or variability in that stretch.
- Its difference image is on target (3.7″, SNR 10.8), so the flux deficit is at the star. That fits a
  systematic as well as a real dip.

**Result.** S16 + S56 give **P = 1084.31/n, 24 aliases from 24.1 to 1084.3 d** (n = 1–12, 14, 15,
17–21, 24, 27, 29, 30, 45).
- The 6.8 h duration on this low-density star (0.43 g cm⁻³) favours P ≳ 50 d for a circular orbit.
  This is soft.
- **Rp ≈ 9 R⊕ (0.8 R_J).**
- **Source and host:**
  - S56 difference image on target (1.9″ ± 2.1, SNR 6.6). It excludes the G = 13.7 star at 22.8″,
    which would need a 2.1% eclipse.
  - S16 has no usable image.
  - Gaia RUWE 1.20, not RV-variable.

**Caveats.**
- S16 sits in a stray-light interval.
- The event is FFI-only.
- The submitting group (UBC) says a paper is in preparation and may already have S16; check before
  any outreach.

## TIC 76140400 (CTOI 76140400.01, steuer 2020; Tmag 10.6, Teff 5059 K)

**The catalogued S26 event is not confirmed.**

| Reduction | Depth (ppt) |
|---|---|
| SPOC 2-min PDC | 1.13 ± 0.13 |
| SPOC SAP | 0.19 ± 0.13 |
| TESS-SPOC | 0.31 |
| Own, stray light kept | not detected |

All 28 FFI cadences in transit carry the stray-light flag, and QLP drops them all.

**The two S40 dips (BTJD 2398.18, 2412.73) are not a planet signal:**
- they are in stray-light-flagged FFI intervals;
- they are V-shaped (ingress fraction 0.43–0.50, T14 11–15 h);
- they show 7.2σ and 8.1σ shifts in SPOC 2-min `mom_centr2`;
- their depth is 2.7–3.0 ppt in SPOC/TESS-SPOC but 5.3–5.4 ppt in our own photometry, i.e. it depends
  on the background or aperture;
- they are not periodic at 14.74 d: the predicted S53/S54 epochs are flat (scale −0.05 ± 0.21 and
  0.15 ± 0.20, against 3.3–3.9 for the dips).

**Verdict: scattered light or a blend.** There is no credible on-target transit, so this star should
be treated as event-not-recovered.

## Effect on the RNAAS note (AAS82207)

**None as submitted.** The 23 "repeats that no single period explains" is a pipeline-rule count, and
it is still true. A revision, if the editor asks for one, could add: hand vetting finds that for 2 of
the 23 (TIC 140215502, tentative; TIC 162031311) one repeat is compromised and the rest fit one
ephemeris. That is a proposal only.
