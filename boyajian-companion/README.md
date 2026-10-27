# Boyajian's star companion: RNAAS note workspace (Oct 27, 2026)

Draft note: Claude Doc "RNAAS draft: Boyajian's star companion period".

## Results so far (cloud run, Kepler + TESS from MAST's public AWS bucket)
- TESS S15 trapezoid fit: tc BTJD 1730.154 ± 0.003, depth 1.16%, T14 20.3 h, T23 9.5 h
  (MF26 T0 = BJD 2458730.164 ± 0.004 adopted).
- Additive-deficit test, 770-845 d: survivors 776.17-776.60, 790.50-790.66, 792.08-792.92,
  793.25-793.40 d (kepler_alias_survivors.csv). 781.6 and 800.9 d excluded.
- Same test recovers MF26 windows 1, 2, 3, 5; window 4 (1345.648 d) is excluded: its n=-2
  epoch (BKJD 1205.9) sits on a 0.45%-deep dip.
- Window 3 (2026 Oct 18-24) not yet tested: AAVSO, ASAS-SN, IRSA unreachable from cloud.

## Files
- kepler_alias_scan.py: standalone, reproduces the survivor table.
- ground_check.py + CLAUDE_CODE_BRIEF.md: run locally for window 3 and the 2021/2023 epochs.
- boyajian_alias_figure.png: Figure 1. make_figure_cloud.py made it (needs a quiet.npy cache
  from the scan, cloud-only helper; kept for reference).

## Local run (2026-10-27)
See RESULTS_LOCAL.md: Kepler survivors reproduced; no TESS sector covers window 3; ASAS-SN
(maintenance) and IRSA ZTF (outage) were unreachable; AAVSO download still to do.
