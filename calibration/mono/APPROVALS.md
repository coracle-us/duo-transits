# Pre-registration approvals (single-transit calibration)

| Version | Approved by | Date | SHA-256 of `prereg_mono_calibration.yaml` | Note |
|---|---|---|---|---|
| 1.0 | Jon McBee | 2026-11-22 | `9a5676efa7a117ba6eca68dac578ce93bc90b85137c4a5cddf6a1bfd2669712c` | Approved as drafted ("Yes"), in the Claude Code session, before any light curve was fetched or any trial run. Three open items listed in the file were open at approval. |
| 1.1 | Jon McBee | 2026-11-22 | `627ca32bb66a5f689ff3dfb3c53af0820a51246c48827e8da1599af930fe502e` | "Approved", in the Claude Code session, after the six gaps in 1.0 and the recommendation to use the Claret table were put to him. No trial had been scored. The `implementation_notes` in the file were not individually approved. |

## Round 2 (`prereg_round2.yaml`)

| Version | Approved by | Date | SHA-256 of `prereg_round2.yaml` | Note |
|---|---|---|---|---|
| 1.0 | Jon McBee | 2026-11-22 | `96ff8d7c0031abbefa8821d651114955e32ceec4c7aee1f6711d9f73ca91845e` | Step 1 of 2. "Approved", in the Claude Code session, after the draft and six changes from round 1 were put to him. Covers the tuning plan, the test design and the pass criteria. No tuning had been run and no test-set star loaded. Step 2 (the frozen rule) is a separate approval. |
| 1.1 | Jon McBee | 2026-11-22 | `376532f74937f5877847c1a0be4b507c2e875324e6513c12bd5b785d9a42f004` | Still step 1 of 2. "Approved", in the Claude Code session, after `round2/TUNING_LOG.md` showed the event criteria could not be met as written. Three changes, listed in the file under `changes_from_1_0`. Tuning-set numbers had been seen; no rule was frozen and no test-set star loaded. Version 1.0 is commit 0be82f0. |
| 1.1, step 2 | Jon McBee | 2026-11-23 | `376532f74937f5877847c1a0be4b507c2e875324e6513c12bd5b785d9a42f004` | Step 2 of 2 (the freeze). "Approved", in the Claude Code session, after the proposed rule, its tuning-set numbers, the 60-of-60 blind cross-check and its limits, and the known weak points were put to him. Rule: `Rule(floor_sig=3.0, min_thr_sig=-inf, weak_ref_excludes=False, w_dur=0.5, w_dur_multi=1.0, snr_min=8.0, depth_ratio_max=1.7, dur_ratio_max=2.0, anti_ratio=1.5, local_anti_max=inf, require_clean=True)`, where `min_thr_sig=-inf` and `local_anti_max=inf` mean switched off, with the R5 table of `round2/FROZEN_RULE.md` (SHA-256 `abdd204651d95378bbe599eef4c30697c755cab0105743f7c90b139ca51defb3`). Code: commit 0deedb1; file hashes in `round2/FREEZE.json`. No test-set light curve had been loaded; the samples were drawn from catalogues only. |

## Round 3 (`prereg_round3.yaml`)

| Version | Approved by | Date | SHA-256 of `prereg_round3.yaml` | Note |
|---|---|---|---|---|
| 1.0 | Jon McBee | 2026-11-23 | `9d17bf0bf972557675655a08ca908be23291d4dfce748ae54ce9ef3e64409e4c` | Step 1 of 2. "Approved", in the Claude Code session, after the draft (0.1, commit 9206386) was put to him with the choice of adding a planet-candidate arm; taken as approval as drafted, without that arm. Covers the tuning plan, the test design and the pass criteria. Round 2's results were known; no round-3 tuning had been run and no round-3 test star loaded. The freeze is a separate approval. |
