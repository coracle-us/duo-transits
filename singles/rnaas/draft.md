# RNAAS draft (working). Numbers filled from singles/results/final_classification.csv (2026-10-11).

**Target:** Research Notes of the AAS, at most ~1,000–1,500 words, one table or one figure.
AASTeX version to follow once numbers are final.

**AAS policy (checked 2026-10-11, https://journals.aas.org/author-llm-guidelines):**
- AI cannot be an author.
- Disclose the tool and its purpose in a statement at submission.
- Supervised drafting is allowed; authors must check all text.
- Figure captions name the tool if one generated the figure.
- Every reference must be verified to exist (via ADS).

---

## Title

Revisiting stale TESS single-transit community candidates: mostly eclipsing binaries, and two new duotransits

## Author

Jon McBee (Independent researcher)

## Abstract

Community TESS Objects of Interest (CTOIs) flagged as single transits are rarely revisited
once uploaded. We re-examined 309 such candidates on 306 stars, most uploaded in 2020–2022,
against all TESS data through Sector 108, using a validated repeat-event search, our own
photometry for sectors not yet processed by the TESS pipelines, difference-image source
localisation and Gaia DR3 binarity diagnostics. The catalogued event is recovered for 272
stars, and 79 of these show repeat events. Of the 79, 21 are on stars that Gaia DR3 flags as
binaries (three Gaia spectroscopic orbits independently reproduce our periods), 18 are blends
with nearby stars, and 23 have repeats no single period explains. For CTOI 147035057.01 a third
eclipse fixes the period at 37.31497 ± 0.00003 d, but the eclipse is grazing and the host is
RV-variable, so the companion is most likely a star. Two candidates survive vetting as
on-target duotransits on Gaia RV-stable stars. For TIC 159614322 (≈0.9 R_J), a partial third
transit in Sector 21 narrows the period to 89.965 or 179.930 d. TIC 229055790 (≈3.8 R⊕) has 38
allowed periods ≥ 27 d. 193
candidates remain single, with a median 52% of log-period space (5–1500 d) now excluded.

## 1. Sample

- **Starting set:** 450 open CTOIs flagged as single or duo transits, from groups not currently
  working their candidates. Source: ExoFOP, downloaded 2026-10-11.
- **Removed:** 46 with missing parameters, 8 on stars that host confirmed planets, 4 on stars
  with a TOI, and 24 on catalogued TESS eclipsing binaries (Prša et al. 2022). This leaves
  368 CTOIs on 359 stars.
- **Not re-analysed:** 32 stars whose ExoFOP entries already carry a community period.
- **Stale entries:** for example, CTOI 1042868.01 is the confirmed planet HD 88986 b
  (Heidari et al. 2024).

## 2. Methods

- **Light curves:**
  - SPOC 2-min, TESS-SPOC and QLP products from MAST.
  - For sectors ≥ 100 without pipeline products, our own aperture photometry from TESScut
    cutouts. Validated by recovering a known eclipse to 2%.
- **Repeat search:** a box of the catalogued depth and duration on a D/4 grid, with
  segment-aware local baselines. Detections need SNR ≥ 8, depth 0.5–2× the catalogued value,
  ≥ 80% coverage, and baseline on both sides. Candidate periods come from linear-ephemeris fits
  to all events, excluding any period whose predicted epochs fall in flat windows.
- **Validation:** applied blind from the first transit only to 26 TESS duotransit hosts with
  known periods. The true period was kept in 26/26 and recovered uniquely in 19/26, with no
  wrong unique answers.
- **Source location:** difference images (out-of-transit minus in-transit) from TESScut.
  Source offsets are measured from the target's Gaia position. On confirmed planets the
  controls land within 1–3″. Off-target means ≥ 6″ and ≥ 3σ.
- **Binarity:** Gaia DR3 RUWE, non_single_star, NSS orbits, and RV variability
  (rv_chisq_pvalue ≤ 0.01, rv_renormalised_gof > 4; Katz et al. 2023).

## 3. Results

**Counts (306 stars).**

| Class | Stars |
|---|---|
| Event not recovered | 34 |
| Single (no repeat) | 193 |
| Repeats, Gaia binarity flag (EB on target) | 21 |
| Repeats, firmly off-target (≥ 2 agreeing difference images) | 18 |
| Repeats, no single period explains them | 23 |
| Periodic or duotransit, hand-vetted | 17 |

**The 17 hand-vetted stars:**
- 2 leads (below);
- 3 tentative, where the repeat is so far seen only in our own S≥100 photometry;
- 5 possibly off-target (one difference image);
- 7 rejected (depth inconsistency, absent in independent reduction, or host too large).

**New duotransits.** Both have consistent depths across epochs and data products, are on
target in difference images, and are re-detected by an independent reduction.
- **TIC 159614322:**
  - S14 and S41 are full transits; flat TESS data leave P = 719.72/n with n ∈ {1, 2, 4, 5, 8}.
  - A partial transit in S21 (ingress plus ~9 h; egress lost in a downlink gap) is seen in
    TESS-SPOC, QLP and our own photometry, and was found again in a blind re-analysis. It
    selects n = 8 or n = 4.
  - **P = 89.9651 ± 0.0006 d or 179.9302 ± 0.0012 d**, from the two full transits; S21 O−C is
    +13 to +26 min.
  - Rp ≈ 0.9 R_J; Gaia RUWE 1.51.
  - Next transit, predicted by both periods: 2026-12-28 06:42 UTC ± 45 min.
  - Source: `singles/writeups/tic159614322/`.
- **TIC 229055790:** S1 + S68; a Planet Hunters TESS candidate with CHIRON spectra on ExoFOP;
  38 aliases, P = 1838.945/n ≥ 26.7 d; Rp ≈ 3.8 R⊕.

**Independent check from Gaia DR3 SB1 orbits.** Our periods match Gaia's for three targets:

| TIC | TESS period | Gaia SB1 period |
|---|---|---|
| 231969683 | 21.978 d | 21.969 ± 0.005 d |
| 230130556 | 22.259 d | 22.245 ± 0.024 d |
| 458419988 | an alias of ΔT = 734.6 d at n = 40 | 18.355 d |

These systems are eccentric eclipsing binaries.

**CTOI 147035057.01.**
- **Eclipses:** S7, S33 and S108. The S108 eclipse is from our own photometry, at the time
  predicted by the n = 19 alias.
- **Ephemeris:** P = 37.31497 ± 0.00003 d, T0 = BJD_TDB 2458515.0661 ± 0.0011.
- **Checks:** on target, depth consistent, no secondary eclipse at phase 0.5.
- **Shape:** grazing. With b and k free, ΔBIC ≈ 10 against the best full-transit model, and
  R₂/R₁ = 0.34 (95%: 0.09–0.75), so R₂ > 1.7 R_J.
  Source: `singles/writeups/tic147035057/shape_fit.json`.
- **Host:** RV-variable in Gaia DR3 (amplitude ≈ 50 km/s). That implies M₂ ≳ 0.6 M☉, so this
  is likely an eccentric EB. RUWE is a low 1.12, which leaves a near-equal SB2 pair possible.
  Gaia DR4 epoch RVs will decide.
- **Next eclipse:** 2027-01-20 03:11 UTC.

**Off-target cases.** 18 firm and 36 possible off-target sources across the sample; for each, we give the eclipse depth the neighbour would need.

**Table 1.** Per-CTOI classification: EB on target / blend / period / still single with
P_min allowed.

## 4. Conclusion

Stale single-transit CTOIs with repeat signals are mostly stellar: cheap Gaia binarity checks
and difference imaging dispose of most of them. A small residue survives: two new on-target
duotransits whose alias lists are now short enough to target, and three tentative repeats
that pipeline products for Sectors ≥ 100 will confirm or reject.

## Acknowledgments and AI disclosure

This work uses data from TESS (MAST; SPOC, TESS-SPOC, QLP, TESScut), ESA Gaia DR3, ExoFOP and
the NASA Exoplanet Archive. The author used an AI assistant (Anthropic's Claude) to write the
analysis code and to draft and edit this text. The author has checked the text, figures and
references.

## References (to verify on ADS before submission)

- Prša et al. 2022, ApJS 258, 16
- Heidari et al. 2024, A&A 681, A55
- Katz et al. 2023, A&A 674, A5 (Gaia DR3 RVs)
- Gaia Collaboration 2023 (DR3)
- Hawthorn et al. 2024, MNRAS 528, 1841
- Eisner et al. (Planet Hunters TESS catalogue paper; identify exact reference for TIC 229055790)
- TESScut: Brasseur et al. 2019
- lightkurve, astroquery, batman (Kreidberg 2015)
