# Repeat recovery: where it is lost

**Analysis made after the outcome was known. Not part of the pre-registered test.**
From the independent rows only (`ind_a_rows.csv`, `ind_b_rows.csv`) and the defs.
Denominator: trials with at least one true repeat in clean data.

## A2 one-TOI

Round 3 (R4 and `sec_box_snr >= 9`): 49 of 90 (54.4%). Without the round-3 veto (R4 only): 53 of 90 (58.9%).

| fate of the true repeat | trials | share |
|---|---|---|
| detected | 49 | 54.4 |
| fails R4 | 35 | 38.9 |
| lost to T5 only (passes R4, sec_box_snr < 9 or undefined) | 4 | 4.4 |
| no event at the truth among the top 3 of its sector | 2 | 2.2 |

R4 conditions failed by the strongest event at the truth (a trial can fail several):

| failed | trials |
|---|---|
| snr | 14 |
| snr+anti | 10 |
| depth+clean | 2 |
| snr+anti+clean | 2 |
| snr+depth+anti | 2 |
| anti | 1 |
| snr+depth | 1 |
| anti+clean | 1 |
| snr+clean | 1 |
| snr+depth+anti+clean | 1 |

Lost to T5 only: best `sec_box_snr` of an R4-passing event at the truth - undefined in 0, otherwise min 6.4, median 8.6, max 8.9; event SNR median 9.0 (min 8.6, max 10.6).

By reference SNR:

| snr_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| unmeasured | 26 | 0 | 0 | 0.0 | 0.0 | 0 |
| 5-8 | 20 | 6 | 9 | 30.0 | 45.0 | 3 |
| 8-12 | 10 | 9 | 10 | 90.0 | 100.0 | 1 |
| 12-20 | 8 | 8 | 8 | 100.0 | 100.0 | 0 |
| over 20 | 26 | 26 | 26 | 100.0 | 100.0 | 0 |

By true depth (catalogued for arm A, injected for arm B):

| depth_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 1000 ppm | 26 | 8 | 9 | 30.8 | 34.6 | 1 |
| 1000-2000 | 27 | 10 | 12 | 37.0 | 44.4 | 2 |
| 2000-5000 | 1 | 1 | 1 | 100.0 | 100.0 | 0 |
| over 5000 | 36 | 30 | 31 | 83.3 | 86.1 | 1 |

By true duration:

| dur_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 3 h | 36 | 14 | 15 | 38.9 | 41.7 | 1 |
| 3-6 h | 37 | 27 | 27 | 73.0 | 73.0 | 0 |
| 6-12 h | 13 | 6 | 9 | 46.2 | 69.2 | 3 |
| over 12 h | 4 | 2 | 2 | 50.0 | 50.0 | 0 |

By expected SNR of one repeat (true depth / duration-scale noise of the merged data):

| exp_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 5 | 26 | 0 | 0 | 0.0 | 0.0 | 0 |
| 5-8 | 16 | 3 | 5 | 18.8 | 31.2 | 2 |
| 8-12 | 10 | 8 | 10 | 80.0 | 100.0 | 2 |
| 12-20 | 12 | 12 | 12 | 100.0 | 100.0 | 0 |
| over 20 | 26 | 26 | 26 | 100.0 | 100.0 | 0 |

By number of clean repeats:

| rep_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| 2 | 4 | 2 | 3 | 50.0 | 75.0 | 1 |
| 3 or more | 86 | 47 | 50 | 54.7 | 58.1 | 3 |

By product and cadence (min) of the sectors holding the clean repeats:

| repeat_products | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| SPOC:10 | 61 | 30 | 33 | 49.2 | 54.1 | 3 |
| SPOC:10 TESS-SPOC:30 | 13 | 9 | 9 | 69.2 | 69.2 | 0 |
| QLP:30 SPOC:10 | 5 | 2 | 3 | 40.0 | 60.0 | 1 |
| SPOC:10 TESS-SPOC:10 TESS-SPOC:30 | 4 | 3 | 3 | 75.0 | 75.0 | 0 |
| QLP:10 SPOC:10 | 3 | 2 | 2 | 66.7 | 66.7 | 0 |
| QLP:10 QLP:30 SPOC:10 | 1 | 0 | 0 | 0.0 | 0.0 | 0 |
| QLP:30 SPOC:10 TESS-SPOC:10 | 1 | 1 | 1 | 100.0 | 100.0 | 0 |
| QLP:30 SPOC:10 TESS-SPOC:30 | 1 | 1 | 1 | 100.0 | 100.0 | 0 |
| SPOC:10 TESS-SPOC:10 | 1 | 1 | 1 | 100.0 | 100.0 | 0 |

By product and cadence of the reference sector:

| ref_prod | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| SPOC:10 | 62 | 31 | 34 | 50.0 | 54.8 | 3 |
| TESS-SPOC:30 | 18 | 13 | 13 | 72.2 | 72.2 | 0 |
| QLP:30 | 8 | 4 | 5 | 50.0 | 62.5 | 1 |
| QLP:10 | 1 | 0 | 0 | 0.0 | 0.0 | 0 |
| TESS-SPOC:10 | 1 | 1 | 1 | 100.0 | 100.0 | 0 |

## B

Round 3 (R4 and `sec_box_snr >= 9`): 203 of 365 (55.6%). Without the round-3 veto (R4 only): 214 of 365 (58.6%).

| fate of the true repeat | trials | share |
|---|---|---|
| detected | 203 | 55.6 |
| fails R4 | 95 | 26.0 |
| no event at the truth among the top 3 of its sector | 56 | 15.3 |
| lost to T5 only (passes R4, sec_box_snr < 9 or undefined) | 11 | 3.0 |

R4 conditions failed by the strongest event at the truth (a trial can fail several):

| failed | trials |
|---|---|
| snr+anti | 31 |
| snr | 27 |
| depth | 12 |
| snr+depth+anti | 8 |
| anti | 5 |
| snr+anti+clean | 4 |
| clean | 2 |
| snr+depth | 2 |
| depth+anti | 1 |
| snr+clean | 1 |
| depth+clean | 1 |
| snr+depth+anti+clean | 1 |

Lost to T5 only: best `sec_box_snr` of an R4-passing event at the truth - undefined in 0, otherwise min 5.5, median 7.5, max 8.2; event SNR median 9.2 (min 8.2, max 15.4).

By reference SNR:

| snr_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| unmeasured | 103 | 3 | 3 | 2.9 | 2.9 | 0 |
| 5-8 | 40 | 7 | 10 | 17.5 | 25.0 | 3 |
| 8-12 | 32 | 16 | 23 | 50.0 | 71.9 | 7 |
| 12-20 | 47 | 41 | 42 | 87.2 | 89.4 | 1 |
| over 20 | 143 | 136 | 136 | 95.1 | 95.1 | 0 |

By true depth (catalogued for arm A, injected for arm B):

| depth_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 1000 ppm | 75 | 16 | 18 | 21.3 | 24.0 | 2 |
| 1000-2000 | 73 | 32 | 36 | 43.8 | 49.3 | 4 |
| 2000-5000 | 82 | 46 | 49 | 56.1 | 59.8 | 3 |
| over 5000 | 135 | 109 | 111 | 80.7 | 82.2 | 2 |

By true duration:

| dur_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 3 h | 11 | 2 | 2 | 18.2 | 18.2 | 0 |
| 3-6 h | 117 | 60 | 62 | 51.3 | 53.0 | 2 |
| 6-12 h | 167 | 93 | 100 | 55.7 | 59.9 | 7 |
| over 12 h | 70 | 48 | 50 | 68.6 | 71.4 | 2 |

By expected SNR of one repeat (true depth / duration-scale noise of the merged data):

| exp_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 5 | 88 | 0 | 0 | 0.0 | 0.0 | 0 |
| 5-8 | 45 | 2 | 6 | 4.4 | 13.3 | 4 |
| 8-12 | 33 | 14 | 20 | 42.4 | 60.6 | 6 |
| 12-20 | 47 | 42 | 43 | 89.4 | 91.5 | 1 |
| over 20 | 152 | 145 | 145 | 95.4 | 95.4 | 0 |

By number of clean repeats:

| rep_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| 1 | 168 | 85 | 91 | 50.6 | 54.2 | 6 |
| 2 | 72 | 37 | 38 | 51.4 | 52.8 | 1 |
| 3 or more | 125 | 81 | 85 | 64.8 | 68.0 | 4 |

By product and cadence (min) of the sectors holding the clean repeats:

| repeat_products | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| SPOC:10 | 134 | 95 | 100 | 70.9 | 74.6 | 5 |
| QLP:10 | 112 | 35 | 37 | 31.2 | 33.0 | 2 |
| TESS-SPOC:10 | 29 | 18 | 19 | 62.1 | 65.5 | 1 |
| QLP:10 SPOC:10 | 26 | 17 | 17 | 65.4 | 65.4 | 0 |
| QLP:10 TESS-SPOC:10 | 17 | 10 | 10 | 58.8 | 58.8 | 0 |
| SPOC:10 TESS-SPOC:10 | 11 | 6 | 7 | 54.5 | 63.6 | 1 |
| QLP:10 SPOC:10 TESS-SPOC:10 | 8 | 6 | 6 | 75.0 | 75.0 | 0 |
| QLP:10 QLP:30 | 6 | 2 | 3 | 33.3 | 50.0 | 1 |
| QLP:10 TESS-SPOC:10 TESS-SPOC:30 | 6 | 4 | 4 | 66.7 | 66.7 | 0 |
| SPOC:10 TESS-SPOC:30 | 3 | 2 | 2 | 66.7 | 66.7 | 0 |
| TESS-SPOC:10 TESS-SPOC:30 | 3 | 2 | 2 | 66.7 | 66.7 | 0 |
| QLP:10 QLP:30 TESS-SPOC:10 TESS-SPOC:30 | 2 | 1 | 2 | 50.0 | 100.0 | 1 |
| QLP:10 SPOC:10 TESS-SPOC:10 TESS-SPOC:30 | 2 | 2 | 2 | 100.0 | 100.0 | 0 |
| QLP:10 SPOC:10 TESS-SPOC:30 | 2 | 2 | 2 | 100.0 | 100.0 | 0 |
| QLP:30 | 2 | 0 | 0 | 0.0 | 0.0 | 0 |
| SPOC:10 TESS-SPOC:10 TESS-SPOC:30 | 2 | 1 | 1 | 50.0 | 50.0 | 0 |

By product and cadence of the reference sector:

| ref_prod | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| SPOC:10 | 166 | 108 | 114 | 65.1 | 68.7 | 6 |
| TESS-SPOC:30 | 81 | 58 | 59 | 71.6 | 72.8 | 1 |
| QLP:30 | 47 | 14 | 15 | 29.8 | 31.9 | 1 |
| TESS-SPOC:10 | 36 | 16 | 17 | 44.4 | 47.2 | 1 |
| QLP:10 | 35 | 7 | 9 | 20.0 | 25.7 | 2 |

By impact parameter:

| b_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| 0-0.3 | 118 | 62 | 70 | 52.5 | 59.3 | 8 |
| 0.3-0.6 | 117 | 69 | 70 | 59.0 | 59.8 | 1 |
| 0.6-0.9 | 130 | 72 | 74 | 55.4 | 56.9 | 2 |

By catalogued / injected depth (the 10% catalogue error):

| cd_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 0.9 | 60 | 34 | 37 | 56.7 | 61.7 | 3 |
| 0.9-1.1 | 249 | 139 | 147 | 55.8 | 59.0 | 8 |
| over 1.1 | 56 | 30 | 30 | 53.6 | 53.6 | 0 |

By catalogued / injected duration:

| cu_bin | n | detected | r4_only | recovery | without_T5 | lost_to_T5 |
|---|---|---|---|---|---|---|
| under 0.9 | 49 | 22 | 23 | 44.9 | 46.9 | 1 |
| 0.9-1.1 | 255 | 144 | 153 | 56.5 | 60.0 | 9 |
| over 1.1 | 61 | 37 | 38 | 60.7 | 62.3 | 1 |

## What the round-3 veto removed on the false-event side (same rows)

- A1 one-TOI stars with a passing event: R4 only 4 of 90, round 3 3 of 90.
- B trials with a false event, own events removed: R4 only 10 of 1072, round 3 10 of 1072.
- B trials with a false event, all events: R4 only 114 of 1072, round 3 93 of 1072.

