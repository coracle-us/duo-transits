# Blind period re-derivation: four Salinas et al. single-transit stars

Independent analysis done 2026-11-22. I did not look at the first analysis, any `coracle/` folder other
than this one, Drive, ExoFOP, TOI/CTOI tables, or any literature period for these TICs. The shared notes
folder was empty. I did not come across any period from another source.

All times are BTJD (BJD_TDB − 2457000). "E" is the epoch number, with E = 0 at the paper's event.

## Summary

| TIC | Repeats? | N transits in fit | P (d) | T0 at E=0 (BTJD) | χ²ᵣ (dof) | >3σ outliers | Aliases ruled out |
|---|---|---:|---|---|---|---|---|
| 275647095 | **yes** | 10 | 17.064042 ± 0.000028 ‡ | 2315.51254 ± 0.00200 ‡ | 1.47 (8) | none | P/2 … P/149 (all P ≥ 1 d), ≥ 24.7σ each |
| 241550976 | **yes** | 7 | 21.196873 ± 0.000018 | 2339.44506 ± 0.00105 | 1.02 (5) | none | P/2 … P/85, ≥ 29.9σ each |
| 107113048 | **yes** | 4 | 49.70148 ± 0.00013 | 2312.94325 ± 0.00302 | 0.03 (2) | none | P/2 … P/49, ≥ 15σ each (P/2 by 2 epochs at 17σ) |
| 230133616 | **yes, 37 transits** | 37 | 23.752987 ± 0.000023 | 2396.12748 ± 0.00070 | 0.71 (35) | none | P/2 … P/80, ≥ 20.9σ each |

‡ Errors scaled by √χ²ᵣ = 1.21. Unscaled values are given below. For the other stars χ²ᵣ ≤ 1.02, so the
errors are left as fitted, not shrunk.

The longer aliases (2P, 3P, …) are all excluded by the data. Every star has transits at consecutive
epochs, or epochs whose differences have no common factor (TIC 107113048: E = 0, 15, 29, 36). So no
period longer than the adopted one fits every detected event.

## Data

- **Light curves:** all TESS-SPOC and QLP HLSP light curves at MAST for each TIC, from `lightkurve`,
  stored in `data/lk/`. There are no 2-min SPOC light curves for these stars.
- **TESScut:** 15×15 cutouts for every sector with no light curve, and for every sector from S94 on.
  These are in `data/tesscut/`. Photometry is in `data/tesscut_phot/`, made by `s03`.
  - The cutouts for TIC 241550976 S11, TIC 107113048 S63 and TIC 230133616 S14, S15, S55 and S85 have
    the target off the science pixels (target pixels ≈ 0 e/s). So **those sectors have no usable data**.
  - TESScut times were re-barycentred to the target direction: t_new = t + (bc_target − TIMECORR).
    The geocentric correction comes from astropy. The change was −0.86 to +0.60 min, depending on
    sector. I checked the method against TESS-SPOC's own TIMECORR for TIC 275647095 S37: they agree to
    0.03 s on average, 0.4 s at most.
- **QLP from S94 on is never used for timing.** In those sectors, TESScut supplies the mid-times. Before
  S94, the timing source for each transit is TESS-SPOC if it covers the transit, otherwise QLP. Other
  sources were fitted for comparison only.

## Method (scripts in `scripts/`, run in order)

1. `s01_search.py`, `s02_download.py`: query MAST and download the data.
2. `s03_tesscut_phot.py`: TESScut aperture photometry. Circular aperture of radius 1–3 px, chosen for
   minimum 1-h scatter. Background from the faintest pixels. Barycentric re-correction as above.
3. `s04_search.py`: blind single-event search over every sector. Detrending with a wotan biweight
   (window 3× the paper's duration), then a sliding box of the paper's duration. Events with SNR > 5 are
   listed with depth relative to the paper (`out/search_<tic>.txt` and `.png`). The repeats were found
   here without assuming any period.
4. `s05_fit_times.py`: per-transit mid-times.
   - Each sector is detrended with the predicted transits masked.
   - A joint batman fit sets one transit shape per star (Rp/R*, a/R*, b; quadratic LD 0.40/0.25).
   - Each transit is then fitted for Tc, depth scale and a linear baseline over ±1.5 durations.
   - Tc error = max(covariance error, Δχ² = 1 profile half-width) × β (time-correlated-noise factor),
     with flux errors rescaled to the local out-of-transit scatter.
   - The fit was iterated once with the refined ephemeris.
   - Output: `out/times_<tic>.csv`, `out/transits_<tic>.png`.
5. `s06_ephemeris.py`: weighted linear ephemeris. It uses timing-source transits with in-transit
   coverage ≥ 0.4 and ingress or egress covered (`out/ephem_<tic>.txt/.json`). Fitting only fully covered
   transits (`--full-only`, `out/ephem_<tic>_full.txt`) gives the same P within 0.3σ.
6. `s07_aliases.py`: alias scan, independent of detrending.
   - Every period P_N = (T_last − T_first)/N ≥ 1 d that puts all detected transits within 45 min of
     integer epochs is a candidate.
   - For each candidate, every other predicted epoch in data with ≥ 50% coverage is tested. The fixed
     transit shape is fitted to the *raw* normalised flux with a quadratic baseline; the depth scale s
     should be 1 if the period is real.
   - A candidate is ruled out if any such epoch has s below 1 by more than 5σ.
   - Output: `out/aliases_<tic>.txt`, with visual check in `out/aliases_<tic>.png` (P, P/2, P/3 epochs).
7. `s09_checks.py`: odd/even depths, depth consistency, eclipse at phase 0.5, stellar density, and
   QLP − TESScut timing comparison (`out/checks_<tic>.txt`).
8. `s10_tables.py`: generates the tables below from the JSON output.

Environment: `uv venv --python 3.12 .venv`, plus lightkurve 2.6.0, batman 2.5.1, wotan, scipy and
astroquery. `setuptools` is needed for batman on Python 3.12.

## QLP clock-offset check (S94+)

QLP − TESScut mid-time differences for transits seen in both:

| TIC | Sector/E | QLP − TESScut (min) |
|---|---|---|
| 275647095 | S100 E104 | −1.7 ± 6.2 |
| 275647095 | S101 E106 | +0.9 ± 5.6 |
| 241550976 | S100 E82 | +0.9 ± 3.3 |
| 241550976 | S102 E85 | **−6.4 ± 2.9** |
| 107113048 | S101 E36 | −1.8 ± 8.1 |

This is consistent with the stated caveat. S102 for TIC 241550976 shows a ≈ −6 min QLP offset at about
2σ. QLP's S102 time for that star would have been a −3.5σ outlier against the ephemeris. It was not used.

## Per-star results

### TIC 275647095: repeated transits, P = 17.0640 d

- **Repeats:** yes. In the blind search, events of the same depth and duration show up in S10, S37,
  S64 (2), S90 (2), S99, S100 and S101 (2). Consecutive transits 17.06 d apart appear in S64, S90 and
  S100–S101. The second S37 transit (E1, 2332.59) is half-covered at the end of the sector, which is
  probably why it was reported as a single.
- **Ephemeris (unscaled):** P = 17.064042 ± 0.000023 d; T0(E=0) = 2315.51254 ± 0.00165;
  T0(E=56) = 3271.09889 ± 0.00101. χ² = 11.78 / 8 dof (χ²ᵣ = 1.47), so the summary quotes errors × 1.21.
  Using only fully covered transits: P = 17.064035 ± 0.000027, χ²ᵣ = 0.80.
- **Aliases:**
  - 2P and longer are ruled out by the consecutive-epoch pairs.
  - P/2 = 8.532 d is ruled out by 11 flat predicted epochs, including S10 (1573.23, 1590.29), S37
    (2324.05), S64 (3057.80), S99 (4047.51 and 4064.58) and S100/S101. Each gives s ≈ 0 ± 0.05, and the
    weakest is 24.7σ below 1.
  - All P/n down to 1 d are excluded the same way.
- **Predicted epochs in good data:** all show the transit (`out/aliases_275647095.txt`, detail list).
- **Concerns:**
  - **Transit is long for this star.** A circular orbit with the fitted a/R* = 16.9 and b = 0 (at its
    bound) implies ρ* ≈ 0.22 ρ☉. The TIC gives 0.63 ρ☉, so the transit is about 1.4× longer than a
    circular orbit allows. Possible explanations: eccentricity (transit near apastron), a larger or
    evolved star than in the TIC, or dilution. TIC contamination ratio is 0.21, the highest of the four.
    This is worth a look before calling it a planet.
  - Odd/even depth: odd 1.036 ± 0.024 vs even 0.979 ± 0.017, a 1.9σ difference. This is not significant
    but worth noting. No eclipse at phase 0.5: −29 ± 153 ppm.
  - Depth scatter between transits is larger than the errors (QLP χ² 16.8/4). The TESScut depths run
    about 7% shallower than QLP and SPOC, which is an aperture/dilution difference. Partial transits
    (E1, E103) fit shallower.
  - E1 (S37, partial, egress missing) is +20 min, +2.2σ. E−43 (S10, 30-min cadence, egress in a gap)
    and E103 (S99, TESScut-only; QLP drops those cadences) are partials too. None is more than 3σ off.
  - The SPOC mid-time of the paper's own event (2315.5113 ± 0.0022) is 8 min earlier than the paper's
    2315.5171. QLP gives 2315.5149.

Ephemeris fit inputs (generated by `s10_tables.py`):

- P = 17.064042 ± 0.000023 d (× √χ²ᵣ: ± 0.000028)
- T0 (E=0, paper event) = 2315.51254 ± 0.00165 BTJD (scaled ± 0.00200)
- T0 (E=56, least-correlated) = 3271.09889 ± 0.00101 BTJD (scaled ± 0.00123)
- χ² = 11.78 for 8 dof, reduced χ² = 1.47; N = 10

| E | Sector | Source | Tc (BTJD) | ± (min) | O−C (min) | O−C/σ | coverage in/ingress/egress |
|---:|---:|---|---|---:|---:|---:|---|
| -43 | 10 | QLP | 1581.75418 | 8.4 | -6.6 | -0.78 | 0.41/0.98/0.00 |
| 0 | 37 | TESS-SPOC | 2315.51129 | 3.2 | -1.8 | -0.57 | 1.00/0.98/0.98 |
| 1 | 37 | TESS-SPOC | 2332.59046 | 9.0 | +20.0 | +2.23 | 0.51/1.05/0.00 |
| 43 | 64 | QLP | 3049.26613 | 3.5 | -0.3 | -0.09 | 1.00/1.00/1.00 |
| 44 | 64 | QLP | 3066.33381 | 5.5 | +4.9 | +0.89 | 1.00/0.98/0.98 |
| 84 | 90 | QLP | 3748.89090 | 4.3 | -1.7 | -0.39 | 1.00/1.00/1.00 |
| 85 | 90 | QLP | 3765.95580 | 4.5 | -0.5 | -0.10 | 1.00/0.98/0.98 |
| 103 | 99 | TESScut | 4073.11403 | 5.6 | +7.4 | +1.33 | 0.59/0.85/0.20 |
| 104 | 100 | TESScut | 4090.17493 | 4.7 | +2.9 | +0.62 | 1.00/1.00/1.00 |
| 106 | 101 | TESScut | 4124.29580 | 4.5 | -7.5 | -1.66 | 1.00/1.00/1.00 |

### TIC 241550976: repeated transits, P = 21.1969 d

- **Repeats:** yes. Transits in S38, S64, S65, S100 (2: E82 and E83, the latter only in TESScut), S101
  (E84, partial) and S102. S64/S65 (E34/E35) and S100–S102 (E82–E85) each show consecutive transits 21.19 d apart.
- **Ephemeris:** P = 21.196873 ± 0.000018 d; T0(E=0) = 2339.44506 ± 0.00105;
  T0(E=47) = 3335.69810 ± 0.00066. χ² = 5.08 / 5 dof (χ²ᵣ = 1.02). Using only full transits:
  P = 21.196872 ± 0.000018.
- **Aliases:**
  - 2P and longer are ruled out by consecutive epochs (E34/35, E82–85).
  - P/2 = 10.598 d is ruled out by 8 flat epochs: S38 2350.04, S64 3049.54, S65 3070.74 and 3091.93,
    and S99–S102. Each has s ≈ 0 ± 0.04, and the weakest is 29.9σ below 1.
  - All P/n down to 1 d are also excluded.
- **Predicted epochs in good data:** all show the transit.
- **Concerns:**
  - E84 (S101) sits at the end of a data segment, followed by a gap and a flux jump. It is a partial
    with depth scale 0.76. It is in the fit (O−C +0.5σ), and leaving it out changes nothing.
  - Depths differ by source: SPOC 1.07, QLP 0.88–1.06, TESScut 0.97–1.04 of the joint model. This
    looks like pipeline dilution rather than astrophysics.
  - Odd/even agree (0.55σ). No eclipse at phase 0.5: −105 ± 256 ppm. Fitted density 0.94 ρ☉ vs TIC
    1.28 ρ☉ for a circular orbit, which is fine given the b–ρ degeneracy.
  - QLP's S102 time is about 6 min early (see the QLP offset table above).

Ephemeris fit inputs (generated by `s10_tables.py`):

- P = 21.196873 ± 0.000018 d (× √χ²ᵣ: ± 0.000018)
- T0 (E=0, paper event) = 2339.44506 ± 0.00105 BTJD (scaled ± 0.00106)
- T0 (E=47, least-correlated) = 3335.69810 ± 0.00066 BTJD (scaled ± 0.00066)
- χ² = 5.08 for 5 dof, reduced χ² = 1.02; N = 7

| E | Sector | Source | Tc (BTJD) | ± (min) | O−C (min) | O−C/σ | coverage in/ingress/egress |
|---:|---:|---|---|---:|---:|---:|---|
| 0 | 38 | TESS-SPOC | 2339.44452 | 1.7 | -0.8 | -0.47 | 1.00/0.96/0.96 |
| 34 | 64 | QLP | 3060.14180 | 3.1 | +4.4 | +1.43 | 0.98/0.99/0.99 |
| 35 | 65 | QLP | 3081.33580 | 3.0 | +0.3 | +0.08 | 0.97/0.86/0.99 |
| 82 | 100 | TESScut | 4077.58629 | 2.5 | -3.4 | -1.35 | 0.99/1.02/0.99 |
| 83 | 100 | TESScut | 4098.78706 | 2.6 | +2.2 | +0.84 | 1.00/1.02/0.99 |
| 84 | 101 | TESScut | 4119.98460 | 6.4 | +3.1 | +0.50 | 0.76/0.92/0.00 |
| 85 | 102 | TESScut | 4141.17903 | 2.2 | -0.4 | -0.17 | 1.00/0.99/1.02 |

### TIC 107113048: repeated transits, P = 49.7015 d

- **Repeats:** yes, 4 transits: S37 (E0), S64 (E15), S90 (E29) and S101 (E36, at the very start of the
  sector). There is only one transit per sector, so the period rests on the alias analysis. The epoch
  differences 15, 14 and 7 have no common factor, so 49.70 d is the *longest* period consistent with all
  four.
- **Ephemeris:** P = 49.701483 ± 0.000129 d; T0(E=0) = 2312.94325 ± 0.00302;
  T0(E=18) = 3207.56995 ± 0.00187. χ² = 0.06 / 2 dof (χ²ᵣ = 0.03), so the residuals are much smaller
  than the errors. With 2 dof this happens by chance about 3% of the time. The errors are probably
  somewhat conservative. I did not shrink them.
- **Aliases:**
  - **P/2 = 24.851 d** is ruled out by two well-covered flat epochs. One is S10 1592.27 (TESS-SPOC
    30-min, coverage 0.93, s = 0.02 ± 0.06). The other is S100 4077.35 (TESScut, coverage 1.0,
    s = 0.00 ± 0.07). They are 17σ and 14σ below a transit. See `out/aliases_107113048.png`.
  - P/3 = 16.567 d is ruled out by 4 epochs (S37, S90, S100 and S101).
  - P/4 and all shorter P/n down to 1 d are ruled out (48 candidates; the weakest is 15σ).
  - These are all the periods consistent with the four events.
- **Predicted epochs in good data:** E0, 15, 29 and 36 all show the transit. No other P epoch falls in
  data. (S10 covers no P epoch, and S63 has no usable data.)
- **Concerns:**
  - Only 4 transits, and the P/2 rejection rests on 2 epochs. Both are clean, but this is the
    least-redundant of the four results.
  - E36 is at the start of S101. QLP and TESScut both cover it fully, and the depth is a bit shallow
    (TESScut 0.90, QLP 0.74–0.83).
  - Odd/even: 1.069 ± 0.044 vs 0.940 ± 0.037, a 2.2σ difference from 2 vs 2 transits that mixes sources
    (the odd epochs are both QLP). Not significant. No eclipse at phase 0.5: 142 ± 382 ppm, 2 epochs.
  - Density 0.70 ρ☉ vs TIC 0.91 ρ☉. Consistent.

Ephemeris fit inputs (generated by `s10_tables.py`):

- P = 49.701483 ± 0.000129 d
- T0 (E=0, paper event) = 2312.94325 ± 0.00302 BTJD
- T0 (E=18, least-correlated) = 3207.56995 ± 0.00187 BTJD
- χ² = 0.06 for 2 dof, reduced χ² = 0.03; N = 4

| E | Sector | Source | Tc (BTJD) | ± (min) | O−C (min) | O−C/σ | coverage in/ingress/egress |
|---:|---:|---|---|---:|---:|---:|---|
| 0 | 37 | TESS-SPOC | 2312.94309 | 4.6 | -0.2 | -0.05 | 0.91/0.96/0.96 |
| 15 | 64 | QLP | 3058.46641 | 6.9 | +1.3 | +0.19 | 1.01/1.01/1.01 |
| 29 | 90 | QLP | 3754.28588 | 4.7 | -0.6 | -0.12 | 1.00/0.99/0.99 |
| 36 | 101 | TESScut | 4102.19688 | 6.3 | +0.3 | +0.05 | 1.01/1.01/1.01 |

### TIC 230133616: many repeated transits, P = 23.7530 d

- **Repeats:** yes, and strongly so. This star is in the northern CVZ, with 40 usable sectors (S16–S86).
  The blind search finds **37 transits** of matching depth (≈ 6300 ppm) and duration (≈ 9 h), one every
  23.75 d. That includes the paper's S40 event, which is E0 here.
  - 33 transits are timed from TESS-SPOC and 4 from QLP: E9 (S48), E41 (S76) and E46 (S80), where
    SPOC's cadence mask drops the transit, and E21 (S58, which has no SPOC).
  - It is surprising that this was reported as a single transit. The match to the paper's event is
    exact (SPOC Tc 2396.1325 vs paper 2396.1307), so either the paper looked at limited data, or the
    "single transit" label refers to something else. **Worth checking with the paper and the other
    analysis.**
- **Ephemeris:** P = 23.752987 ± 0.000023 d; T0(E=0) = 2396.12748 ± 0.00070;
  T0(E=14) = 2728.66929 ± 0.00063. χ² = 24.81 / 35 dof (χ²ᵣ = 0.71). No transit is more than 2σ off.
- **Aliases:**
  - 2P is ruled out by many consecutive epochs.
  - P/2 = 11.877 d is ruled out by 37 flat half-phase epochs (weakest 20.9σ below a transit).
  - All P/n down to 1 d are also excluded.
- **Predicted epochs in good data:** every one with ≥ 50% coverage is a detected transit (the alias scan
  finds 0 "extra" epochs for P). The partials are E38 (S73, egress missing, s = 0.75) and E40 (S75, QLP
  only, 42% coverage, unusable fit, not in the ephemeris).
- **Concerns:**
  - Odd/even agree: 1.021 ± 0.011 vs 1.018 ± 0.012 (0.2σ). Phase-0.5 eclipse: −113 ± 74 ppm over 37
    epochs, so no secondary eclipse.
  - The fitted shape (b → 0, a/R* = 21) gives ρ* = 0.22 ρ☉ for a circular orbit, while the TIC says
    0.07 ρ☉ (R* = 2.5 R☉). The transit is *shorter* than central for this star, which a non-zero b would
    explain. The joint fit's b is poorly constrained at 30-min and 10-min cadence, so this is not alarming.
  - SPOC depth χ² = 91/32. Some extra scatter in transit depth (0.72–1.15) is likely from sector-to-sector
    dilution/background differences in the FFI light curves.
  - E9 (S48, QLP) has a baseline step, and its error (12.4 min) reflects that.

Ephemeris fit inputs (generated by `s10_tables.py`):

- P = 23.752987 ± 0.000023 d
- T0 (E=0, paper event) = 2396.12748 ± 0.00070 BTJD
- T0 (E=14, least-correlated) = 2728.66929 ± 0.00063 BTJD
- χ² = 24.81 for 35 dof, reduced χ² = 0.71; N = 37

| E | Sector | Source | Tc (BTJD) | ± (min) | O−C (min) | O−C/σ | coverage in/ingress/egress |
|---:|---:|---|---|---:|---:|---:|---|
| -27 | 16 | TESS-SPOC | 1754.79661 | 4.6 | -0.3 | -0.07 | 0.98/1.09/0.91 |
| -26 | 17 | TESS-SPOC | 1778.54848 | 5.3 | -1.9 | -0.37 | 1.03/0.91/0.91 |
| -24 | 19 | TESS-SPOC | 1826.06146 | 6.1 | +8.1 | +1.34 | 0.98/0.91/1.09 |
| -23 | 20 | TESS-SPOC | 1849.81034 | 3.5 | +2.2 | +0.64 | 0.98/0.91/1.09 |
| -22 | 21 | TESS-SPOC | 1873.56797 | 5.4 | +8.9 | +1.65 | 0.98/0.91/1.09 |
| -21 | 21 | TESS-SPOC | 1897.31711 | 5.1 | +3.4 | +0.66 | 0.98/1.09/1.09 |
| -20 | 22 | TESS-SPOC | 1921.06637 | 5.6 | -2.0 | -0.36 | 0.87/1.09/0.91 |
| -19 | 23 | TESS-SPOC | 1944.81564 | 5.4 | -7.3 | -1.35 | 1.03/1.09/0.91 |
| -17 | 25 | TESS-SPOC | 1992.31631 | 10.2 | -15.0 | -1.47 | 1.03/0.91/0.72 |
| -16 | 26 | TESS-SPOC | 2016.07830 | 5.8 | -2.0 | -0.34 | 0.92/0.72/0.91 |
| 0 | 40 | TESS-SPOC | 2396.13253 | 4.9 | +7.3 | +1.47 | 1.00/0.97/0.97 |
| 8 | 47 | TESS-SPOC | 2586.15586 | 8.9 | +6.5 | +0.73 | 0.94/0.97/0.79 |
| 9 | 48 | QLP | 2609.90097 | 12.4 | -4.9 | -0.39 | 1.00/1.03/1.03 |
| 10 | 48 | TESS-SPOC | 2633.65093 | 4.8 | -9.2 | -1.91 | 1.00/0.97/0.97 |
| 11 | 49 | TESS-SPOC | 2657.40602 | 5.9 | -6.2 | -1.06 | 1.01/1.03/1.03 |
| 14 | 52 | TESS-SPOC | 2728.66758 | 5.2 | -2.5 | -0.47 | 1.00/0.97/0.97 |
| 15 | 53 | TESS-SPOC | 2752.41861 | 4.5 | -5.3 | -1.17 | 1.00/1.03/1.03 |
| 16 | 54 | TESS-SPOC | 2776.17803 | 10.6 | +4.0 | +0.38 | 0.92/1.03/0.97 |
| 19 | 56 | TESS-SPOC | 2847.43515 | 4.5 | +1.3 | +0.30 | 1.00/1.01/0.99 |
| 20 | 57 | TESS-SPOC | 2871.18411 | 4.4 | -4.5 | -1.01 | 1.00/1.01/1.01 |
| 21 | 58 | QLP | 2894.93958 | 5.2 | -0.9 | -0.17 | 1.00/0.99/1.01 |
| 22 | 59 | TESS-SPOC | 2918.69436 | 4.2 | +1.7 | +0.40 | 0.99/0.99/1.01 |
| 23 | 60 | TESS-SPOC | 2942.44666 | 4.9 | +0.7 | +0.15 | 1.00/1.01/0.99 |
| 38 | 73 | TESS-SPOC | 3298.74399 | 12.1 | +4.3 | +0.36 | 1.00/1.01/0.56 |
| 39 | 74 | TESS-SPOC | 3322.49554 | 4.6 | +2.3 | +0.49 | 1.00/0.99/1.01 |
| 41 | 76 | QLP | 3369.99554 | 8.8 | -6.3 | -0.72 | 1.00/1.01/0.99 |
| 42 | 76 | TESS-SPOC | 3393.75294 | 4.7 | +0.0 | +0.01 | 1.00/0.99/1.01 |
| 43 | 77 | TESS-SPOC | 3417.50828 | 6.9 | +3.4 | +0.50 | 1.00/0.95/1.01 |
| 44 | 78 | TESS-SPOC | 3441.25701 | 5.3 | -2.7 | -0.51 | 1.00/1.01/0.99 |
| 45 | 79 | TESS-SPOC | 3465.01124 | 5.2 | -0.9 | -0.17 | 1.00/1.01/1.01 |
| 46 | 80 | QLP | 3488.76244 | 5.8 | -3.5 | -0.60 | 1.00/0.99/1.01 |
| 47 | 81 | TESS-SPOC | 3512.51879 | 4.7 | +1.4 | +0.29 | 1.00/1.01/0.99 |
| 48 | 82 | TESS-SPOC | 3536.27087 | 6.8 | +0.1 | +0.01 | 1.00/1.01/1.01 |
| 49 | 83 | TESS-SPOC | 3560.02936 | 6.6 | +8.0 | +1.21 | 1.00/0.99/1.01 |
| 50 | 83 | TESS-SPOC | 3583.77628 | 11.9 | -0.8 | -0.06 | 1.00/1.01/0.99 |
| 51 | 84 | TESS-SPOC | 3607.53143 | 5.2 | +2.4 | +0.46 | 1.00/1.01/0.99 |
| 53 | 86 | TESS-SPOC | 3655.03962 | 6.1 | +5.5 | +0.91 | 1.00/0.99/1.01 |

## Things that would change the answer

- For TIC 107113048, an unseen transit at a P/2 epoch would need both S10 1592.27 and S100 4077.35 to be
  bad data. Both look clean in `out/aliases_107113048.png`.
- None of the four periods depends on QLP S94+ times. Dropping every QLP mid-time would leave
  TIC 275647095 with fewer epochs (E−43, 43, 44, 84 and 85 are QLP-only), but the period is still fixed
  by SPOC E0/E1 and TESScut E103–E106, and the alias scan (which uses all sources for depth) would not
  change.
