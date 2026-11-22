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
9. Follow-up checks (stellar size and source of the dip), each taking the TIC as its first argument;
   run for TIC 275647095, 241550976, 107113048 and 230133616:
   - `s11_dl_cutouts.py TIC SECTOR…`: 21×21 TESScut cutouts for difference imaging (`data/tesscut_diff/`).
   - `s12_gaia.py`: Gaia DR3 parameters, and neighbours within 200″ with the eclipse depth each would
     need to mimic the dip (`out/gaia_<tic>.txt/.json`).
   - `s14_radius.py`: radius and density from parallax, G and the TIC Teff/E(B−V) (`out/radius_<tic>.txt`).
   - `s13_diffimg.py`: difference images, scene fit, source position and own-pixel test
     (`out/diffimg_<tic>.txt/.png`). Off-silicon cutouts are skipped with the same rule as `s03`.
     `--exclude E,…` reports an epoch but leaves it out of the combined results.
   - `s15_rho_fixed.py`: joint shape fit with a/R* free vs fixed to the Gaia density for a circular orbit
     (`out/rhofixed_<tic>.txt`).
   - `s16_eb_neighbour_107113048.py` (no argument): the eclipsing-binary neighbour of TIC 107113048
     (`out/eb_neighbour_107113048.txt/.png`).

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
    circular orbit allows. TIC contamination ratio is 0.21, the highest of the four.
    **Follow-up (below): the star is not larger than the TIC says, and the dip comes from the target's
    position, so eccentricity (e ≳ 0.34) is now the leading explanation.**
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
    **Follow-up (below): confirmed. Gaia gives 0.96 R☉, the dip is on the target, and a circular orbit
    at the Gaia density fits.**
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
  - **Follow-up (below):** the dip is on the target and a circular orbit at the Gaia density fits. A
    G = 12.74 eclipsing binary 97″ away has an eclipse that overlaps E29, and E29 is the only transit
    that rules out 74.55 d and 149.10 d. The timing and leakage checks show E29 is a real transit
    of the target.

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
  - **Follow-up (below):** the dip is on the target, but Gaia gives **RUWE = 3.8**, so there is
    probably an unresolved companion. The transit shape rejects a circular orbit around the 2.45 R☉
    star Gaia implies (Δχ² ≈ 49), and the companion is very large (≈ 19 R⊕ ≈ 1.7 R_J). The first pass
    called the density mismatch "not alarming"; it is now this star's main open question.

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

## Follow-up: TIC 275647095 stellar size and source of the dip

These checks were added after the period work. They use only Gaia DR3 and my own TESScut cutouts,
and do not change any period.

Scripts: `s11_dl_cutouts.py 275647095 37 64 90`, then `s12_gaia.py`, `s14_radius.py`, `s13_diffimg.py`
and `s15_rho_fixed.py` with argument 275647095 (method step 9).

Outputs: `out/gaia_275647095.txt/.json`, `out/radius_275647095.txt`, `out/diffimg_275647095.txt/.png`,
`out/rhofixed_275647095.txt`.

**Gaia DR3 for the target.**
- G = 12.74, BP−RP = 0.98, parallax 2.554 ± 0.015 mas (≈ 390 pc).
- RUWE = 1.03, `ipd_frac_multi_peak` = 0, `non_single_star` = 0. There is no astrometric or imaging hint
  of a close companion.
- DR3 gives no FLAME or GSP-Phot parameters for this star.

**Radius.** I derived the radius from the parallax, G, the DR3 bolometric correction, the TIC Teff
(5489 ± 131 K) and A_G = 2.74 × E(B−V) (TIC E(B−V) = 0.058).
- **R★ = 1.14 ± 0.06 R☉** and L = 1.07 L☉. With the TIC mass (0.96 M☉ ± 10%), **ρ★ = 0.64 ± 0.13 ρ☉**,
  which agrees with the TIC.
- The circular-orbit transit fit needs ρ★ = 0.22 ρ☉, which would mean R★ ≈ 1.63 R☉. None of the
  20,000 Monte Carlo draws reaches it. **So the star's size does not explain the long transit.**
- The duration is 1.42× the circular expectation. That needs **e ≥ 0.34**, the minimum for a transit
  near apoapsis. b = 0 was preferred by the fit, and a larger b would only raise that minimum.
- Direct test (`s15_rho_fixed.py`): holding a/R* at the circular value for the Gaia density (24.0
  instead of the fitted 16.9) makes the joint fit worse by Δχ² = 1806, or 970 after scaling by
  χ²/dof = 1.86. A circular orbit around this star does not fit the transit shape.

**Neighbours.** Gaia has no other source within 10.3″. The 10.3″ star (G = 17.9) is too faint to make
the dip even if fully eclipsed: it would need a 114% eclipse. Within 200″, 41 stars are bright enough
that an eclipse of < 50% could mimic the dip. The ones that matter most:
- **G = 14.31 at 19.6″** (under one TESS pixel) needs a 5% eclipse.
- **G = 11.37 at 132″** needs a 1.2% eclipse.
- **G = 12.08 at 190″** needs a 1.5% eclipse.

The G = 14.31 star also accounts for the TIC's 0.21 contamination ratio (flux ratio 0.24).

**Difference images.** I used 7 full transits: S37 E0, S64 E43/E44, S90 E84/E85, S100 E104 and S101
E106. The two partial transits (E103, E105) were excluded. Each image is the out-of-transit flux minus
the in-transit flux, after fitting a linear trend in time to each pixel.
- **Location of the dip.** A Gaussian-PSF scene model with Gaia positions, refined by a global shift
  fitted to the out-of-transit image, puts the source **0.1–1.6″ from the target in 5 of the 7
  transits**. S64 gives 4–5″, with worse fits (χ² 110 and 143 for 61 dof) and a visible jitter residual
  on the bright G = 11.4 star in the E44 image.
- **Combined offset:** dRA = 0.0 ± 0.6″, dDec = −1.5 ± 0.6″ (error from the transit-to-transit
  scatter). This is within the systematic error of a Gaussian PSF on TESS data.
- **19.6″ neighbour:** the "neighbour is the source" model is worse than "target is the source" in every
  transit. The total Δχ² is 2989, and the neighbour's position is 42σ from the measured source position.
- **Resolved neighbours, own-pixel test.** If a star were the source, its own 3×3 pixels would lose
  (total lost flux) × (its PSF fraction there). The measured ratio is ≈ 0 for all 41 candidates, against
  1 for the true source; each is excluded at ≥ 42σ.
  - G = 11.37 at 132″: −0.01 ± 0.02.
  - G = 12.08 at 190″: −0.04 ± 0.03.
  - The largest ratio is 0.14, for a star 52″ (2.5 px) away, which picks up the target's PSF wing.

**Conclusion for TIC 275647095.** The dip comes from within about 2″ of the target, and no Gaia star
there could make it. The star's size agrees with the TIC. What remains:
- an **eccentric (e ≳ 0.34) giant planet**, with Rp/R★ = 0.089, so Rp ≈ 11 R⊕ for R★ = 1.14 R☉; or
- an eclipsing system blended within ~2″ and below Gaia's detection or resolution. RUWE 1.03 and the
  flat-bottomed transit shape argue against this.

Radial velocities, or a ground-based transit with high-resolution imaging, would separate the two.
The period result is unchanged.

**How the follow-up changes the first-pass concerns.**

| First pass | Follow-up finding | Status |
|---|---|---|
| **Transit too long:** circular fit needs ρ★ = 0.22 ρ☉; TIC gives 0.63 ρ☉. Possible causes: bigger star, blend or eccentric orbit. | Gaia radius is 1.14 ± 0.06 R☉ (ρ★ = 0.64 ρ☉). A circular orbit would need 1.63 R☉, reached by 0 of 20,000 Monte Carlo draws. | **Bigger star ruled out.** The duration is 1.42× circular, so e ≥ 0.34. b = 0 is the most favourable case; a larger b raises the minimum. |
| **Contamination 0.21**, highest of the four; source unknown. | It comes from the G = 14.3 star at 19.6″ (flux ratio 0.24). That star is disfavoured as the source in all 7 transits (Δχ² = 2989) and lies 42σ from the measured dip position. | **Explained; not the source.** |
| **Blend** (nearby eclipsing pair): not tested. | The dip is within about 2″ of the target. All 41 Gaia stars within 200″ bright enough to mimic it fail the own-pixel test, including G = 11.4 and G = 12.1. Gaia has no other source within 10″, and RUWE is 1.03. | **Narrowed** to an eclipsing pair within about 2″ that Gaia does not detect. |
| **TESScut depths about 7% shallower** than QLP and SPOC, attributed to aperture dilution. | The 19.6″ star sits in the TESScut aperture and the TESScut photometry has no crowding correction (SPOC has one). Fully in the aperture it would make transits about 20% shallower; partly in, about 7%. | **Likely cause found.** Direction and rough size agree; not computed exactly. |
| **Odd/even depths differ by 1.9σ;** no phase-0.5 eclipse (−29 ± 153 ppm). | Not tested by the follow-up. | **Open but weak.** At 1.9σ this is not evidence of an eclipsing binary. |
| **Depth scatter exceeds the errors** (QLP χ² 16.8/4). | Not tested by the follow-up. | **Open.** Most likely detrending or partial-coverage systematics. |
| **Period:** P = 17.064042 ± 0.000028 d, aliases ruled out at ≥ 24.7σ. | Unchanged. | — |

## Follow-up: TIC 241550976 stellar size and source of the dip

The same checks as for TIC 275647095 (method step 9). They use only Gaia DR3 and my own TESScut
cutouts, and do not change the period.

Scripts: `s11_dl_cutouts.py 241550976 38 64 65`, then `s12_gaia.py`, `s14_radius.py`, `s13_diffimg.py`
and `s15_rho_fixed.py` with argument 241550976. Outputs: `out/gaia_241550976.txt/.json`,
`out/radius_241550976.txt`, `out/diffimg_241550976.txt/.png`, `out/rhofixed_241550976.txt`.

**Gaia DR3 for the target.**
- G = 12.63, BP−RP = 0.90, parallax 2.853 ± 0.014 mas (≈ 350 pc).
- RUWE = 1.02, `ipd_frac_multi_peak` = 0, `non_single_star` = 0. There is no astrometric or imaging hint
  of a close companion.
- Unlike TIC 275647095, DR3 has FLAME and GSP-Phot parameters: R = 0.96 (0.94–0.98) R☉,
  M = 1.06 (1.02–1.10) M☉, L = 0.98 L☉, evolutionary stage 115 (main sequence), GSP-Phot
  Teff = 5868 K (TIC: 6034 ± 135 K).

**Radius.** From the parallax, G and the TIC Teff and E(B−V) = 0.128: **R★ = 0.96 ± 0.05 R☉**,
L = 1.09 L☉. That matches FLAME and the TIC (0.956 R☉). With the TIC mass (1.12 M☉ ± 10%),
**ρ★ = 1.27 +0.24/−0.21 ρ☉**.
- The free shape fit gives ρ★ = 0.94 ρ☉ for a circular orbit (b = 0.45, a/R* = 31.6). 5% of the Gaia
  posterior lies at or below that, so this is about a 1.6σ difference.
- Direct test (`s15_rho_fixed.py`): with a/R* fixed at 34.9 (the circular value for the Gaia density),
  the fit moves to b = 0.16 and Rp/R* = 0.122. It is worse by only Δχ² = 3.2, or 1.9 after scaling by
  χ²/dof = 1.65. **A circular orbit fits; there is no sign of eccentricity**, and the first pass's
  "b–ρ degeneracy" reading holds.
- Rp/R* = 0.122–0.125 gives **Rp ≈ 12.7–13.0 R⊕** (about 1.15 R_J) for R★ = 0.96 R☉.

**Neighbours.** The closest Gaia sources (4.5″, 11.0″, 11.8″, …) are G = 18.8–20.8. The brightest star
within 30″ is G = 17.2 at 17.7″, which would need a 122% eclipse. None can make the dip. Within 200″
this is a crowded field (946 Gaia sources), and 46 are bright enough that an eclipse of < 50% could
mimic the dip. The nearest is **G = 15.13 at 36.4″** (1.7 px; needs 19.6%), and the brightest are
**G = 10.89 at 113″** (2.1%) and **G = 11.75 at 109″** (2.6%). Together, the G = 15.1 star at 36″ and the
G = 15.3 star at 43″ account for most of the TIC's 0.145 contamination ratio.

**Difference images.** I used 7 full transits: S38 E0, S64 E34, S65 E35, S100 E82/E83, S101 E84 and S102
E85. S11 was skipped because the target is off-silicon in that cutout (no flux). E84 passes the
full-transit rule but sits before a gap and flux jump, and fits poorly (χ² 302 for 63 dof).
- **Location of the dip.** The source is 0.3–2.3″ from the target in 4 transits, and 3.4–5.2″ in S102,
  S64 (χ² 283) and S101 (χ² 302).
- **Combined offset:** weighted dRA = −1.9 ± 0.3″, dDec = 0.0 ± 0.3″. From the transit-to-transit
  scatter: dRA = −2.2 ± 0.5″, dDec = −0.1 ± 0.8″. The westward shift is consistent between transits,
  but it is within the few-arcsec systematic of a Gaussian PSF on TESS pixels in a crowded field. Gaia
  has no star within 4.5″ bright enough to be the source.
- **G = 15.13 neighbour at 36.4″:** "target is the source" beats "neighbour is the source" in every
  transit, with total Δχ² = 7598. The neighbour is 95σ from the measured source position.
- **Own-pixel test:** every one of the 46 candidates is excluded (≥ 52σ). Ratios are 0.43 and 0.28 for
  the stars at 36″ and 43″, which pick up the target's PSF wing (1.7–2 px away); beyond 55″ they are
  all < 0.09. The two brightest give −0.02 ± 0.02 (G = 10.89) and 0.04 ± 0.02 (G = 11.75).
- Unrelated: S64 shows a strong variable signal at the bottom edge of the cutout, about 10 px
  (210″) south of the target, beyond the 200″ search radius. It is not seen in S38 or S65 (the S100–S102
  cutouts are 15×15 and do not reach it) and does not affect the target's pixels.

**Conclusion for TIC 241550976.** The dip comes from within about 2″ of the target, no resolved Gaia star
can make it, and the star is a main-sequence G dwarf of 0.96 R☉ whose density fits a circular orbit.
The leading explanation is a **giant planet of about 13 R⊕ on a 21.2-day circular orbit**. The fallback,
as for TIC 275647095, is an eclipsing pair within about 2″ that Gaia does not detect; RUWE 1.02 and
the clean odd/even and phase-0.5 checks argue against it. The period result is unchanged.

**How the follow-up changes the first-pass concerns.**

| First pass | Follow-up finding | Status |
|---|---|---|
| **Density:** fit gives 0.94 ρ☉ vs TIC 1.28 ρ☉; judged fine given the b–ρ degeneracy. | Gaia gives R★ = 0.96 R☉ and ρ★ = 1.27 ρ☉. With a/R* fixed to that density, the fit is only Δχ² = 1.9 (scaled) worse, at b = 0.16. | **Confirmed fine.** A circular orbit fits; no eccentricity needed. |
| **Depths differ by source** (SPOC 1.07, QLP 0.88–1.06, TESScut 0.97–1.04), attributed to pipeline dilution. | The TIC contamination is 0.145, mostly from the G = 15.1 and 15.3 stars at 36–43″. TESScut has no crowding correction; SPOC does. A contamination of 0.145 would make undiluted-aperture depths up to about 13% shallower; TESScut is about 9% below SPOC. | **Consistent with dilution.** Not computed exactly. |
| **Blend:** not tested. | The dip is within about 2″ of the target (a consistent −2″ RA shift, within PSF systematics). All 46 resolved stars that could mimic it fail the own-pixel test. Nothing in Gaia within 4.5″ could make it. | **Narrowed** to an eclipsing pair within about 2″ that Gaia does not detect. |
| **E84 partial** before a gap and flux jump (depth 0.76). | It also gives the worst difference-image fit (χ² 302/63), but still points at the target. | **Unchanged:** the data there are poor; it does not move the period. |
| **Odd/even agree** (0.55σ), no phase-0.5 eclipse. | Not tested by the follow-up. | **No concern.** |
| **QLP S102 about 6 min early.** | Not tested by the follow-up. | **Unchanged** (pipeline timing, see the QLP offset table). |
| **Period:** P = 21.196873 ± 0.000018 d, aliases ruled out at ≥ 29.9σ. | Unchanged. | — |

## Follow-up: TIC 107113048 stellar size, source of the dip, and an eclipsing neighbour

The same checks (method step 9), plus `s16_eb_neighbour_107113048.py`. Only Gaia DR3 and my own TESScut
cutouts are used.

Scripts: `s11_dl_cutouts.py 107113048 37 64 90`, then `s12_gaia.py`, `s14_radius.py`, `s15_rho_fixed.py`
with argument 107113048, `s13_diffimg.py 107113048 --exclude 29`, and `s16_eb_neighbour_107113048.py`.
Outputs: `out/gaia_107113048.txt/.json`, `out/radius_107113048.txt`, `out/rhofixed_107113048.txt`,
`out/diffimg_107113048.txt/.png`, `out/eb_neighbour_107113048.txt/.png`.

**Gaia DR3 for the target.**
- G = 12.89, BP−RP = 0.87, parallax 2.450 ± 0.020 mas (≈ 410 pc).
- RUWE = 1.09, `ipd_frac_multi_peak` = 0, `non_single_star` = 0. There is no hint of a close companion.
- FLAME: R = 1.00 (0.99–1.01) R☉, M = 0.97 (0.93–1.01) M☉, L = 0.95 L☉, age 5.9 Gyr, evolutionary
  stage 286 (main sequence). GSP-Phot Teff = 5541 K (TIC: 5694 ± 119 K).

**Radius.** From the parallax: **R★ = 1.01 ± 0.05 R☉**, L = 0.98 L☉, and with the TIC mass (1.01 M☉)
**ρ★ = 0.96 +0.18/−0.16 ρ☉**. All of this matches FLAME and the TIC.
- The free shape fit gives 0.70 ρ☉ for a circular orbit (b = 0.48). 4% of the Gaia posterior is at or
  below that.
- With a/R* fixed at 56.1 (circular orbit, Gaia density), the fit moves to b = 0.24 and Rp/R* = 0.085.
  It is worse by only Δχ² = 0.5, or 0.3 after scaling. **A circular orbit fits.**
- Rp/R* = 0.085–0.087 gives **Rp ≈ 9.4–9.6 R⊕** (about 0.85 R_J).

**Neighbours.** This is a sparse field (123 Gaia sources within 200″), and Gaia has **nothing within
28″**. The closest sources are G = 20.2 at 28.2″ and G = 19.7 at 39.5″, which would need eclipses of
739% and 463%. 21 stars within 200″ could mimic the dip with an eclipse of < 50%. The nearest is
G = 17.03 at 47.2″ (needs 40%); the brightest are **G = 12.74 at 97″** (1.6%) and **G = 13.15 at 62″**
(2.0%). The TIC contamination ratio is only 0.017.

**Difference images.** There are 4 full transits: S37 E0, S64 E15, S90 E29 and S101 E36. S63 was skipped
(target off-silicon).
- **E0, E15 and E36:** the source is 1.0″, 1.0″ and 2.1″ from the target. Combined: dRA = −0.2 ± 0.5″,
  dDec = +0.7 ± 0.5″ (scatter-based: 0.0 ± 0.5″, +0.4 ± 0.8″). The G = 17.0 star at 47″ is 63σ away,
  and "target is the source" wins by Δχ² = 1293 in total.
- **Own-pixel test (3 transits):** all 21 candidates are excluded at ≥ 27σ. G = 12.74 gives
  0.02 ± 0.03 and G = 13.15 gives −0.01 ± 0.03.
- **E29 (S90) is excluded from the combined results.** Its difference image is dominated by the G = 12.74
  star 97″ (4.6 px) away. A Gaussian PSF cannot model that residual, so the fitted source position is
  meaningless there (χ² 13162 for 61 dof). See the next paragraph for why.

**The G = 12.74 star at 97″ is an eclipsing binary, and one of its eclipses overlaps E29.** Its own 3×3
pixels show eclipses at 2332.47 (S37, 15.7% deep, at the end of the sector), 3043.40 (S64, 17.5%,
FWHM 7.7 h), 3754.26 (S90, 18.0%, FWHM 7.8 h) and 4118.83 (S101, 7.3%). Two broad 25-hour dips in S101
appear in both stars' pixels, so they are scattered light, not eclipses. This matters because, as the
alias scan above shows, E29 is what rules out 149.10, 74.55 and 29.82 d. Without it those periods fit
E0/E15/E36 with no other covered epoch to test them. Three tests show **E29 is a real transit of the
target**:
- **Timing.** E0, E15 and E36 alone predict E29 at 3754.28665 ± 4.8 min. E29 is observed at
  3754.28588 ± 4.7 min (−1.1 min, −0.2σ). The binary's eclipse centre is 41 min earlier, at 3754.2581.
  Light leaking from the eclipse would be centred on the eclipse, not on the target's ephemeris.
- **No leakage elsewhere.** During the binary's other eclipses the target's 3×3 pixels do not dim. In
  S37 (15.7% eclipse) they change by −25 ppm against a 1-h scatter of 219 ppm. In S64 (17.5% eclipse)
  they change by +2576 ppm, i.e. brighter, not fainter. In S101 (7.3%) the change is +1165 ± 1468 ppm.
  A transit-like 0.86% dip from leakage would need far more light to spill across 4.6 px than any
  of these show.
- **Depth.** E29's depth scale is 1.11 ± 0.06, consistent with the other transits (0.90–1.00). A small
  leakage contribution can't be excluded, but it is not what makes the dip.

**Conclusion for TIC 107113048.** The dip comes from within about 2″ of the target, there is no Gaia
source within 28″, and the star is a 1.0 R☉ main-sequence G dwarf whose density fits a circular orbit.
The leading explanation is a **giant planet of about 9.5 R⊕ on a 49.70-day circular orbit**. The
fallback is an eclipsing pair within about 2″ that Gaia misses; RUWE 1.09 argues against it. The
period stands, but **E29 carries more weight than the first pass said**: it is the only transit that
rules out 74.55 and 149.10 d. A transit at the next E42 epoch (2026-12-25 21:43 UTC), or a transit
recovered at any other epoch, would settle this independently of S90.

**How the follow-up changes the first-pass concerns.**

| First pass | Follow-up finding | Status |
|---|---|---|
| **Only 4 transits;** P/2 rejected by 2 flat epochs. "Least-redundant of the four." | Still true. In addition, **E29 alone** rules out 3P = 149.10, 3P/2 = 74.55 and 3P/5 = 29.82 d, and E29 overlaps an eclipse of a G = 12.74 binary 97″ away. | **Sharper concern, but E29 holds up.** Its timing matches the other three to 1.1 min (the binary's eclipse is 41 min off), and the binary's other eclipses do not dim the target's pixels. |
| **E36 depth a bit shallow** (TESScut 0.90, QLP 0.74–0.83). | The difference image for E36 puts the source 2.1″ from the target with a good fit (χ² 61/63). The two broad scattered-light dips in S101 show this sector's photometry is noisy. | **No concern for the source.** The shallow depth is most likely S101 systematics. |
| **Odd/even 2.2σ,** from 2 vs 2 transits; both odd epochs are QLP (E15, E29). | E29's depth (1.11) may include a little leakage from the binary's eclipse, which would push the odd average up. | **Weaker.** A small leak at E29 is a plausible cause. |
| **Density 0.70 vs 0.91 ρ☉.** | Gaia gives 0.96 ρ☉; with a/R* fixed to it the fit is worse by Δχ² = 0.3 (scaled). | **Confirmed fine.** A circular orbit fits. |
| **Blend:** not tested. | The dip is within about 2″ of the target, no Gaia source within 28″, and all 21 candidates fail the own-pixel test. | **Narrowed** to an undetected pair within about 2″. |
| **Period:** P = 49.701483 ± 0.000129 d. | Unchanged. E0/E15/E36 alone give 49.701496 ± 0.000151 d. | — |

## Follow-up: TIC 230133616 stellar size, source of the dip, and a likely companion

The same checks (method step 9). Only Gaia DR3 and my own TESScut cutouts are used.

Scripts: `s11_dl_cutouts.py 230133616 40 49 57 74 79 84`, then `s12_gaia.py`, `s14_radius.py`,
`s15_rho_fixed.py` and `s13_diffimg.py` with argument 230133616. Outputs: `out/gaia_230133616.txt/.json`,
`out/radius_230133616.txt`, `out/rhofixed_230133616.txt`, `out/diffimg_230133616.txt/.png`.

**Gaia DR3 for the target.**
- G = 12.59, BP−RP = 0.74, parallax 0.953 ± 0.040 mas (≈ 1.05 kpc).
- **RUWE = 3.79.** This is well above the ≈ 1.4 usually taken as the sign of a poor single-star
  astrometric fit, and the commonest cause is an unresolved companion. `ipd_frac_multi_peak` = 0 and
  `ipd_gof_harmonic_amplitude` = 0.02, so the scan images do not show a resolved double; the companion
  would be close (≲ 0.5″) or the RUWE comes from orbital motion. `non_single_star` = 0 (no orbital
  solution published).
- DR3 gives no FLAME or GSP-Phot parameters for this star.

**Radius.** From the parallax, G and the TIC Teff (6162 ± 134 K) and E(B−V) = 0.034: **R★ = 2.45 ± 0.12
R☉**, L = 7.8 L☉ (a subgiant; TIC 2.53 R☉, log g 3.70). With the TIC mass (1.17 M☉), **ρ★ = 0.08 ± 0.01
ρ☉**. With RUWE 3.8 the quoted parallax error is probably too small, so the real uncertainty is larger.

**Transit shape vs this star.** The free fit gives a/R* = 21.0 and b = 0, i.e. 0.22 ρ☉ for a circular
orbit. With a/R* fixed at 15.0 (circular orbit, Gaia density), the fit needs b = 0.70 to shorten the
transit, and it is worse by **Δχ² = 47 (49 after scaling)**, across 34 well-sampled transits. The
first pass guessed that a higher b would reconcile the two, but the transit's ingress and egress are too
sharp for that. The transit is about 0.71× the duration a central transit of this star would last. Three
ways to explain it:
- **The star is smaller than Gaia's parallax implies.** A circular orbit needs about 1.74 R☉. That
  happens if the target supplies only about half the light, e.g. a near-equal companion, which fits the
  high RUWE. It would also happen if the parallax were about 1.34 mas instead of 0.95, but that needs a
  0.39 mas error, about 10× the quoted one.
- **An eccentric orbit with the transit near periastron**, e ≥ 0.33.
- **An eclipsing pair blended with the target** (see below).

**Size of the transiting object.** Rp/R* = 0.072–0.078 with R★ = 2.45 R☉ gives **Rp ≈ 19–21 R⊕
(1.7–1.9 R_J)**. In the companion case (target gives half the light, R★ ≈ 1.73 R☉), the undiluted
Rp/R* is about 0.10, and Rp comes out the same, about 19 R⊕. That is large for a planet on a 23.75 d
orbit (a ≈ 0.17 AU, Teq ≈ 1100 K, so only mild inflation is expected). Brown dwarfs are about 1 R_J, so
an object this size is more likely an inflated giant planet at the top of the known range, or a
low-mass star of about 0.17 R☉. The flat phase-0.5 light curve (−113 ± 74 ppm) rules out a bright
secondary eclipse but not an M dwarf. Only radial velocities can tell these apart.

**Neighbours.** Gaia has **nothing within 53″**. The closest are G = 18.1 at 53.5″, which would need a
102% eclipse, and G = 14.99 at 53.6″ (6.4%). 18 stars within 200″ could mimic the dip with an eclipse
of < 50%; the brightest are G = 12.52 at 90″ (1.2%) and G = 12.78 at 195″ (1.4%). TIC contamination
is 0.0075.

**Difference images.** I used 6 full transits spread over four years: S40 E0, S49 E11, S57 E20, S74 E39,
S79 E45 and S84 E51. The older 15×15 cutouts (S14, S55, S85) have the target off-silicon, and S15 has
no good cadences; all four were skipped.
- The source is 0.7–3.7″ from the target in each transit. Combined: dRA = −0.3 ± 0.4″,
  dDec = −0.4 ± 0.4″ (scatter-based: −0.2 ± 0.3″, −0.5 ± 0.8″). This is the tightest of the four stars.
- The G = 15.0 star at 53.6″ is 86σ away, and "target is the source" wins by Δχ² = 2115 in total.
- All 18 own-pixel candidates are excluded at ≥ 50σ. G = 12.52 at 90″ gives 0.01 ± 0.02.

**Conclusion for TIC 230133616.** The dip comes from the target to within about 1″, and no resolved
star can make it. But this is the least clear-cut of the four:
- Gaia's RUWE of 3.8 says the target probably has a close companion.
- The transit does not fit a circular orbit around a star as large as Gaia's parallax implies.
- The transiting object is 1.7–1.9 R_J, large for a planet at 23.75 d.

The possibilities are a giant planet or a low-mass star transiting a star that is either smaller and
diluted by a companion, or orbiting eccentrically. An eclipsing pair hidden in the unresolved companion
also fits everything here, and RUWE makes that more plausible for this star than for the other three.
High-resolution imaging (to find the companion) and radial velocities are the next steps. The period
(37 transits) is unaffected.

**How the follow-up changes the first-pass concerns.**

| First pass | Follow-up finding | Status |
|---|---|---|
| **Density:** fit gives 0.22 ρ☉ vs TIC 0.07 ρ☉; "the transit is shorter than central, which a non-zero b would explain … not alarming." | Gaia gives R★ = 2.45 R☉, ρ★ = 0.08 ρ☉. With a/R* fixed to that density the fit needs b = 0.70 and is worse by Δχ² ≈ 49 over 34 transits. | **Now the main concern.** b alone does not explain it. The options are a diluting companion, e ≥ 0.33, or a blend. |
| (Not raised) | **RUWE = 3.79:** a probable unresolved companion. | **New.** It favours the diluting-companion and hidden-blend options. |
| (Not raised) | Rp ≈ 19–21 R⊕ (1.7–1.9 R_J), under either radius scenario. | **New.** Large for a planet; a low-mass star is possible. |
| **Odd/even agree** (0.2σ), no secondary (−113 ± 74 ppm). | Not tested by the follow-up. | **Unchanged.** Rules out a near-equal eclipsing pair on the target, not an M dwarf or a faint blend. |
| **SPOC depth scatter** (χ² 91/32). | Not tested by the follow-up. | **Unchanged.** |
| **Blend:** not tested. | The dip is within about 1″ of the target, with nothing in Gaia within 53″ and all 18 candidates excluded. | **Resolved stars ruled out.** The unresolved companion is untested. |
| **Single-transit label in the paper.** | Not affected. | **Unchanged:** 37 transits. |
| **Period:** P = 23.752987 ± 0.000023 d. | Unchanged. | — |

## Things that would change the answer

- For TIC 107113048, an unseen transit at a P/2 epoch would need both S10 1592.27 and S100 4077.35 to be
  bad data. Both look clean in `out/aliases_107113048.png`.
- For TIC 107113048, if E29 (S90) were not a transit of the target, 149.10, 74.55 and 29.82 d would fit
  E0, E15 and E36 with no other covered epoch to test them. E29 overlaps an eclipse of a neighbouring
  eclipsing binary, but it falls 1.1 min from the time E0/E15/E36 predict and 41 min after the binary's
  eclipse centre, and the binary's other eclipses do not dim the target's pixels (see the follow-up).
- None of the four periods depends on QLP S94+ times. Dropping every QLP mid-time would leave
  TIC 275647095 with fewer epochs (E−43, 43, 44, 84 and 85 are QLP-only), but the period is still fixed
  by SPOC E0/E1 and TESScut E103–E106, and the alias scan (which uses all sources for depth) would not
  change.
