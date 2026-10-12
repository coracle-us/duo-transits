# Hand vetting of the full-run leads (2026-10-11)

Leads are the Gaia-clean, not-off-target stars with repeat events, as printed by
`singles/scripts/classify_final.py`. Each was run through `scripts/vet_candidate.py`, which does:
- a batman joint fit with a per-epoch depth scale;
- a per-product depth check;
- centroids, dumps/edges and Gaia neighbours.

Plots and numbers are in `singles/results/vetting/<TIC>_n<n>.{png,json}`.

**Independent re-check.** Each event was re-measured with the loader from `scripts/independent_bls.py` (no shared code):
- 1-day rolling-median detrend, which suppresses long transits, so absolute depths read low;
- local quadratic baseline;
- in-transit mean.

That loader does not read our own TESScut photometry. So repeats seen only in S≥100 own photometry cannot be re-checked independently yet.

## Leads (on target, consistent, independently re-detected)

| TIC | Events | Depth scale per epoch | Independent depths (ppm, SNR) | Source | Rp | Surviving periods |
|---|---|---|---|---|---|---|
| 159614322 | S14 (TESS-SPOC) + S41 (SPOC 2-min; the CTOI event) | 0.99 ± 0.05, 1.00 ± 0.03 | 2863 (19), 2640 (19) | S41 diff image on target (1.2″, SNR 15); S14 image too weak | ≈ 0.96 R_J (R* = 1.75) | 719.657 / n for n = 1, 2, 4, 5, 8: 719.66, 359.83, 179.91, 143.93, 89.96 d |
| 229055790 | S1 + S68 (SPOC 2-min) | 1.04 ± 0.05, 0.97 ± 0.05 | 819 (12.6), 774 (11.3) | S68 diff image on target (3.6″, 1.1σ, SNR 5.5) | ≈ 3.8 R⊕ | 38 aliases 1838.945 / n, from 26.7 d to 1838.9 d |

**TIC 159614322** (Tmag 11.2, Teff 5650 K, R* 1.75 R☉):
- Every product agrees, scale 0.99–1.07 across SPOC PDC/SAP, TESS-SPOC PDC/SAP and QLP.
- Centroid shifts 0.07σ and 0.2σ.
- No dumps or edges.
- No Gaia neighbour within 39″ could host the signal.
- Gaia RUWE 1.51, mildly elevated; not RV-variable (22 transits). Worth a note, not a rejection.
- ExoFOP: CTOI by vogl (2021) with no period; no TOI, imaging or spectra.
- No literature hit in a web search.

**TIC 229055790** (Tmag 9.1, Teff 5865 K, R* 1.06 R☉):
- A Planet Hunters TESS candidate (listed on the PHT results page; CTOI by eisner, 2020).
- The S1 epoch lies near a momentum dump, but its depth matches S68 in both our fit and the independent re-check.
- ExoFOP lists SMARTS/CHIRON spectra (PI Johnson, Dec 2022), so someone is already following it up. Our contribution is the second transit and the alias list; cite the PHT work and do not claim discovery.

## Tentative (repeat seen only in our own S≥100 photometry) → recheck list

- **139726266** (K dwarf, R* 0.75, so Rp ≈ 2.8 R⊕):
  - S32 SPOC event plus an own-photometry S107 event at the same depth (scale 0.95 vs 1.07).
  - Diff images inconclusive (SNR < 5).
  - A G = 16.1 star at 6.7″ would need a 23% eclipse.
  - 27 aliases of 2109.0 / n.
- **168702709:**
  - S31 SPOC event is clean;
  - the S107 own-photometry event is partial: it has a gap through mid-transit, then the egress;
  - 2-min centroid shift of 3.9σ.
- **268534931:** already on the weekly recheck.

## Not confirmed / rejected

- **69745398:** the S96 repeat is absent in the independent re-check (0.04×, 0.7σ); the S28 dip is shorter than the catalogued 17 h.
- **337665591:**
  - depths differ (χ² 12.5/1);
  - the field is crowded (submitter's note);
  - the S66 diff image is 8.4″ off at SNR 11.
- **365250897:** depths differ (χ² 40/1); QLP shows no dip.
- **51079186:** depths inconsistent across epochs and products (χ² 22/4).
- **58946313, 64549889, 67884964:** see the MANUAL notes in `classify_final.py` (earlier vetting).
