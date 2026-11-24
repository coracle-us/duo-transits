# Step 0: what Lee et al. (arXiv:2611.17590) already did

Paper: "The Eclipsing δ Scuti Binary HIP 7666: Absolute Dimensions and Pulsation Properties",
J. W. Lee et al., arXiv:2611.17590. Full text read on 2026-11-24 from the arXiv PDF: body,
figure captions and Tables 1–5. The paper has no appendix.

## Verdict

**They did no third-body, light-travel-time, quadratic-ephemeris or pulsation-timing analysis.
The timing they did use is one TESS sector (S58) plus historical minima, fit with a linear
ephemeris.** Step 0 therefore does not stop the project. The gap the note fills is all the other
TESS sectors, timed by both eclipses and pulsation phases.

## Data they used

- **TESS:** Sector 58 only, SPOC 120-s, SAP flux detrended and converted to magnitudes. It runs
  2022 Oct 29 – Nov 26 (BJD 2459882.33–2459910.06), 18,631 points, CROWDSAP 0.961. The light-curve
  DOI is 10.17909/d6tf-9118.
- **Spectra:** 63 echelle spectra, 60 from BOAO/BOES and 3 from TNO/MRES.
- **Historical minima:** Escolà-Sirisi et al. 2005 and Feng et al. 2021, converted from HJD to
  BJD_TDB.

## Their timing work (§2, Table 1, Figure 2)

- **S58 eclipse times:** 23 mid-eclipse times (11 primary, 12 secondary) by Kwee & van Woerden
  (1956), with errors of 0.00004–0.00015 d.
- **Linear ephemeris** from an OLS fit to the TESS and historical minima, their Eq. 1:
  `C = BJD 2451769.56299(50) + 2.37222188(18) × E`.
- **O−C:**
  - The plot (Fig. 2) spans epochs 0 to about 4,600 (2000–2022), with scatter of a few × 10⁻³ d
    in the historical points.
  - In S58, the secondaries sit about 0.00019 d late relative to the primaries. They attribute
    this to the primary's pulsations, not apsidal motion.
  - They did not fit a quadratic or LTTE model and say nothing about a third body.
- **W-D ephemeris (Table 3):** T0 = BJD 2459892.050858(82), P = 2.37222122(33) d.
- **Third light:** l3 = 0.0375 ± 0.0024, matching the TESS contamination of 3.9%. No third body
  is discussed.

## Binary (Tables 3–4)

| | value |
|---|---|
| e | 0 (RV Keplerian fit 0.003 ± 0.002) |
| i | 81.366 ± 0.053° |
| a | 10.872 ± 0.054 R☉ |
| q | 0.6842 ± 0.0068 |
| M_A, M_B | 1.822 ± 0.029, 1.246 ± 0.020 M☉ (total 3.068 M☉) |
| R_A, R_B | 2.154 ± 0.017, 1.170 ± 0.015 R☉ |
| distance | 316 ± 8 pc (Gaia DR3: 313 ± 1 pc) |

## Their 15 frequencies (Table 5; d⁻¹, amplitude in mmag)

The frequencies come from S58 out-of-eclipse residuals (phases 0.044–0.456 and 0.544–0.956)
after cycle-by-cycle W-D subtraction, found with PERIOD04 and kept at S/N > 5.

| | f | A | S/N | remark |
|---|---|---|---|---|
| f1 | 25.4656(5) | 0.865(66) | 22.43 | |
| f2 | 0.4203(6) | 0.741(63) | 20.17 | f_orb |
| f3 | 21.1939(7) | 0.536(51) | 17.84 | f4 − 3f_orb (?) |
| f4 | 22.4278(7) | 0.601(58) | 17.78 | f3 + 3f_orb (?) |
| f5 | 21.5836(7) | 0.519(53) | 16.64 | f3 + f_orb |
| f6 | 20.3533(9) | 0.363(47) | 13.28 | f3 − 2f_orb |
| f7 | 24.6340(16) | 0.297(65) | 7.83 | f1 − 2f_orb |
| f8 | 28.0866(18) | 0.270(69) | 6.74 | 2f7 − f3 |
| f9 | 21.8091(12) | 0.328(54) | 10.38 | |
| f10 | 27.1883(15) | 0.334(68) | 8.34 | f1 + 4f_orb |
| f11 | 22.6767(14) | 0.300(59) | 8.70 | f9 + 2f_orb |
| f12 | 28.0289(23) | 0.218(68) | 5.46 | f10 + 2f_orb |
| f13 | 21.0911(19) | 0.194(50) | 6.64 | |
| f14 | 20.9630(19) | 0.188(50) | 6.47 | f9 − 2f_orb |
| f15 | 37.8159(17) | 0.123(28) | 7.41 | |

The modes are weak, below 0.9 mmag. Pulsation timing will therefore be much less precise than
eclipse timing. That is expected and will be stated.

## Other prior work they cite

- **Pawar et al. 2024:** the same S58 data plus the Feng RVs. They found 7 frequencies and
  e = 0.027. Lee et al. do not credit them with any timing analysis; I still need to check
  that paper directly in Step 4.
- **Feng et al. 2021:** ground-based BVR photometry, with minima used in the O−C. Whether they fit
  more than a linear ephemeris is also to be checked in Step 4.
