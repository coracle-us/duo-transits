# Pilot: wide single-transit search, Sector 104 (QLP)

Run under `PLAN.md` v1.0 (committed before any download). Counts only: files naming candidate stars stay in
`tess_wide/data/` (outside git) and are not public. Scripts: `s1_starcut.py` … `s7_report.py`.

## Stars

| | Stars |
|---|---|
| QLP light curves in S104 | 1,152,719 |
| fail T ≤ 12 | 839,464 |
| fail R★ ≤ 2 R☉ (or no radius) | 220,163 |
| fail TIC disposition (split/artefact/duplicate) | 17,836 |
| fail contratio ≤ 1 | 6,475 |
| fail 3000–7500 K | 2,900 |
| fail log g ≥ 3.5 | 0 |
| **pass the cut** | **65,881** |
| downloaded (failed) | 65,881 (0) |
| searched | 65,881 |

## Events at each stage (real search, and the identical search on 2 − flux)

| Stage | Real | Control |
|---|---|---|
| detections (SNR ≥ 8 and the other `search_any` rules) | 11,305 on 5,060 stars | 7,906 on 4,512 stars |
| common-mode bins flagged | 29 | 35 |
| after common-mode veto | 8,618 | 6,077 |
| in stars with exactly one event (single) | 1,532 | 1,898 |
| (events in "repeating" / "multi-event" stars) | 830 / 6,256 | 806 / 3,373 |
| after size test (Rp ≤ 2 R_J) | 1,412 | 1,739 |
| CNN pass | 225 | 300 |
| CNN fail / out of range / error | 1148 / 39 / 0 | 1391 / 48 / 0 |
| CNN pass, already known (cross-match) | 22 | – |
| **CNN pass, not known: new single-event candidates** | **203** | – |
| out of CNN range, not known (counted separately) | 32 | – |

Cross-match hits among CNN passes, by catalogue: {"toi": 8, "ctoi": 7, "confirmed": 1, "hawthorn": 1, "salinas": 0, "prsa2022": 0, "gaia_veb": 9, "vsx": 10}. VizieR queries that failed:
0. Catalogue snapshots: 2026-11-25.

**Control rule (PLAN 4):** final-stage control count 300 against real 225:
**control-level** — by the rule fixed in PLAN.md, the search is not yet separating signal from systematics in this sector.

**Flags among CNN passes (PLAN 5.4; recorded, not used to cut):**

| | Real | Control |
|---|---|---|
| CNN passes | 225 | 300 |
| in-event SAP_BKG change > 3σ | 119 | 104 |
| SAP_X or SAP_Y shift > 3σ | 162 | 196 |
| either flag | 181 | 212 |
| by camera (1/2/3/4) | 94/55/42/34 | 137/74/56/33 |

## Completeness (injection–recovery, 3,000 injections)

Overall: detection 38 ± 1 %; end to end (through the CNN pass,
true Rp ≤ 2 R_J, 2,956 injections) 23 ± 1 %
(24 % counting CNN out-of-range as kept). Injections falling in data gaps count as missed.

**At detection** (bins of true depth and injected duration; ± is the 1σ Wilson interval):

| depth (ppm) \ duration (h) | 2.0–3.2 | 3.2–5.1 | 5.1–8.1 | 8.1–13.0 |
|---|---|---|---|---|
| 300–656 | 1 ± 1 % (1/108) | 2 ± 1 % (2/119) | 3 ± 2 % (4/117) | 6 ± 2 % (7/126) |
| 656–1,435 | 4 ± 2 % (5/132) | 7 ± 2 % (12/165) | 14 ± 3 % (23/160) | 23 ± 3 % (40/176) |
| 1,435–3,137 | 30 ± 4 % (48/158) | 36 ± 4 % (55/154) | 57 ± 4 % (84/147) | 48 ± 4 % (79/166) |
| 3,137–6,860 | 56 ± 4 % (95/169) | 68 ± 4 % (93/136) | 59 ± 4 % (97/165) | 56 ± 4 % (84/150) |
| 6,860–15,000 | 64 ± 4 % (89/139) | 69 ± 4 % (107/154) | 67 ± 5 % (69/103) | 67 ± 4 % (99/147) |

**End to end:**

| depth (ppm) \ duration (h) | 2.0–3.2 | 3.2–5.1 | 5.1–8.1 | 8.1–13.0 |
|---|---|---|---|---|
| 300–656 | 0 ± 0 % (0/108) | 2 ± 1 % (2/119) | 2 ± 1 % (2/117) | 2 ± 1 % (2/126) |
| 656–1,435 | 2 ± 1 % (3/132) | 4 ± 1 % (6/165) | 11 ± 3 % (18/160) | 14 ± 3 % (25/176) |
| 1,435–3,137 | 18 ± 3 % (28/158) | 21 ± 3 % (32/154) | 35 ± 4 % (52/147) | 34 ± 4 % (56/166) |
| 3,137–6,860 | 28 ± 3 % (47/169) | 49 ± 4 % (66/136) | 42 ± 4 % (70/165) | 40 ± 4 % (60/150) |
| 6,860–15,000 | 35 ± 4 % (47/136) | 38 ± 4 % (56/146) | 44 ± 5 % (44/101) | 48 ± 4 % (66/137) |

## Runtime and disk (this MacBook: 18 cores, 128 GB)

| Step | Time | Notes |
|---|---|---|
| TIC rows for 1,152,719 stars (MAST TAP, 300 ids per query) | 20 min | the TAP server slowed for ~10 min mid-run |
| download + extract 65,681 stars (two runs, started on a partial star list; 200 more in a timing test) | 62 min | 44.8 GB in; FITS deleted after extraction; 16 parallel requests |
| search (real + control) | 1.5 min on 16 cores | 22 ms CPU per star |
| injections (3,000) | 0.1 min | |
| CNN (real / control / injections) | 51 / 62 / 36 s | CPU, weights loaded once |
| disk kept (compressed arrays) | 10.8 GB | 164 kB per star |

## Projection for all QLP sectors (S1–S104) on this MacBook

Assumes every sector passes the same fraction of its QLP stars as S104 (5.7 %), file sizes by cadence era
as measured on MAST, and the pilot's download rate (17.7 stars/s) and per-star search cost.
- stars (star-sectors): 5.85 million over 98 sectors (MAST has no QLP bulk
  script for S74, S75, S76, S77, S78, S79; they are left out);
- download: 2.49 TB in, 599 GB kept, about 3.8 days at the pilot's rate;
- search: 36 CPU-hours (2.2 h on 16 cores), plus injections;
- the CNN cost is per surviving event and negligible next to downloads.

The download dominates. The pass fraction is the weakest assumption: it depends on galactic latitude and on how
deep QLP went in each sector.
