# Candidate-period re-check with `robust_periods` (2026-10-12)

**Why.** The Fable 5.1 review of the K2 project (finding C2, `k2singles/README.md`) showed that
`pilot.periods_from_repeats` proposes trial periods only from the repeat nearest the event and needs
every repeat on one ephemeris, so one spurious repeat can delete the true period. The singles
project used it for `candidate_periods` in every result JSON, and `classify_final.py` classes
depend on the count. This re-check replaces it with `k2lib.robust_periods` (every repeat anchors
trial periods; each period reports how many repeats it explains).

**Script.** `singles/scripts/recheck_periods_robust.py`. Results JSONs, scripts and the submitted
note are untouched. For each of the 102 JSONs with ≥ 1 repeat (`pilot`, `batch2`, `batch3`,
`full`, `rerun_2026-10-11`), it rebuilds the original grid (`pilot.grid_windows`, OWN products
included, the JSON's `event_T`, catalogue D and depth). It uses only light-curve files that already
existed when the JSON was first committed, so later downloads and later own photometry (e.g.
TIC 147035057 S108) do not leak in.

**Outputs.** `period_recheck.csv` has one row per JSON; `in_frozen_note` marks the 97 JSONs behind
`final_classification_2026-10-11_rnaas.csv`. `stdout.txt` is the run log.

## Reproduction checks (all pass)

| Check | Result |
|---|---|
| Old `periods_from_repeats` on the rebuilt grid = stored `candidate_periods` | 102/102 |
| Repeat search on the rebuilt grid = stored `repeats` | 102/102 |
| `classify_final.py` ladder with the stored count = frozen class | 97/97 |

## Results

Rules fixed before the run (script docstring). The like-for-like replacement for the old list is the
set of robust periods that explain **every** repeat.

**Classes: no change.** 0 of 102 JSONs change class, so the frozen counts stand:
- repeats-inconsistent 23;
- EB on target (Gaia) 21;
- off-target EB 18;
- hand-vetted 17 (14 duotransit/narrowed, 2 periodic on target, 1 periodic possibly off-target).

**No period was gained anywhere.** The failure mode where the nearest repeat is spurious and the true
period is never proposed did not occur in this sample.

**Period sets: 3 change, all losses, none affecting a class or a claim.**

| TIC | Class (frozen) | Stored → robust | Mechanism |
|---|---|---|---|
| 51079186 | duotransit/narrowed; rejected by hand (depths inconsistent) | 4 → 3; loses 5.3066 d | Old joint fit split the residuals (±37 min < D/3 = 56 min). Robust admits a repeat only within D/3 of the event + anchor ephemeris |
| 22642772 | off-target EB (firm) | 1 → 0; loses 22.2104 d | Robust checks flat windows at the event epoch: one predicted epoch lands on flat data. The old fit shifted T0 by 22 min and missed it |
| 48882295 | off-target EB (firm) | 2 → 0; loses 63.31 and 31.65 d | Same as 51079186: joint-fit residuals up to 108 min (D/3 = 135 min) |

These are differences in tolerance handling, not the C2 failure. Stars with one repeat are
unaffected by construction, since there is only one anchor. That covers 27 of the 50 JSONs that have stored periods.

**Named checks (published values come from separate hand analyses):**

| TIC | Stored | Robust | Status |
|---|---|---|---|
| 159614322 | 5 aliases, 719.657/n, n = 1, 2, 4, 5, 8 | identical | Unaffected. The hand values 89.9651 / 179.9302 d (n = 8, 4) are in the set |
| 147035057 | 15 aliases, 708.976/n, incl. n = 19 (37.3145 d) | identical | Unaffected. 37.31497 d comes from the hand analysis with own S108, which is not in the pilot grid |
| 229055790 | 38 aliases, 1838.945/n | identical | Unaffected |

## The 23 repeats-inconsistent stars: periods hidden by one repeat?

For ≥ 2 repeats, the CSV lists periods that explain all but one repeat (`n_robust_allbut1`,
`allbut1_left_out_T`) and the best support (`max_support`).

**6 of the 23 have such periods.** None changes class under the current rule. Each would move to
"duotransit / narrowed" only if hand vetting showed the left-out repeat to be spurious. That vetting
has **not** been done.

| TIC | Source, Gaia | Repeats (T, sector, product, depth ppm, SNR) | Survives if one repeat is dropped |
|---|---|---|---|
| 140215502 | on-target, clean | 1460.54 S5 TESS-SPOC 2786 16.0; 3681.07 S87 QLP 3845 14.7; 4275.76 S107 **OWN** 4671 11.4 | Drop S5: 297.33/n, n = 1–4, 6 (297.3, 148.7, 99.1, 74.3, 49.6 d). Relies on an OWN S107 repeat (tentative per the S106–108 lesson) |
| 162031311 | on-target, clean | 1745.69 S16 TESS-SPOC 2588 12.1; 3578.21 S83 QLP 4680 21.8 | Either one: 40 aliases. The two depths differ 1.8× |
| 76140400 | inconclusive, clean | 2398.10 S40 SPOC 1705 18.8; 2412.84 S40 SPOC 1515 17.3 | Either one: 8 aliases (30.8–400 d). Two strong dips 14.74 d apart in one sector that do not phase with the event: worth a look |
| 407663579 | mixed, clean | 3914.70 S96 QLP 3182 8.4; 3954.63 S97 QLP 4759 8.9 | Either one: 13 aliases. Both near threshold |
| 441056242 | off-target (possible), clean | 1372.62 S2 SPOC 588 8.1; 3914.54 S96 SPOC 647 8.5 | Drop S2 (3.46 d before the event): 64 aliases 2538.46/n |
| 248960573 | on-target, RUWE 1.77 | 4304.93 S108 **OWN** 3523 15.3; 4311.48 S108 **OWN** 1275 11.6 | Either one: 109 aliases. Both OWN, in a scattered-light sector; most likely both systematics |

**The other 17 have no all-but-one period:**
- **3 have a single repeat whose aliases all fall on flat data** (239190192, 283222856, 293845053).
  Under the old function, these single-repeat stars also returned nothing.
- **14 have 3–39 repeats with support ≤ 10.** These are dense dip trains: 53847727 is the known
  3.01-d EB, and 83384234 and 76804724 have spacings of about 4.3 and 3.7 d. A P ≥ 5 d search cannot
  return them by design.

Outside the repeats-inconsistent class, 11 rows have all-but-one periods but no full-support period (6
event-not-recovered, 3 Gaia EBs including one re-run JSON, 2 off-target EBs). They cannot change class, because the period count
is not used for them.

## Effect on the RNAAS note (AAS82207, under review)

**No number in `singles/rnaas/rnaas.tex` changes:**
- per-class counts: 306 / 34 / 193 / 79 / 21 / 18 / 23 / 17;
- the 17 hand-vetted stars, as a set;
- the 2 / 3 / 5 / 7 split;
- the two leads' periods and alias counts;
- the Table 1 periods.

The Methods sentence "Candidate periods come from linear-ephemeris fits to all events" still
describes what was done.

**Optional revision text,** for use only if the editor asks for changes. It is not applied:

> A re-check in which every repeat may anchor trial periods changes no classification; for 6 of the
> 23 stars without a single period, a period explains all but one repeat.

**Open follow-up.** Hand-vet the left-out repeats of 140215502, 162031311 and 76140400 (on-target or
inconclusive, Gaia-clean, pipeline-product repeats). If one proves spurious, the 23 / 17 split
would shift in a future version, not in AAS82207 as submitted.
