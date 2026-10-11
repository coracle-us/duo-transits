# RNAAS draft (working). Numbers marked {{...}} are filled in from the full run.

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

Most repeating signals among stale TESS single-transit community candidates are eclipsing binaries

## Author

Jon McBee (Independent researcher)

## Abstract

Community TESS Objects of Interest (CTOIs) flagged as single transits are rarely revisited
once uploaded. We re-examined {{N_ANALYSED}} such candidates on {{N_STARS}} stars, most uploaded
in 2020–2022, against all TESS data through Sector 108. Our methods were:
- a validated search for repeat events;
- our own photometry for sectors not yet processed by the TESS pipelines;
- difference-image source localisation;
- Gaia DR3 binarity diagnostics.

{{N_REPEAT}} candidates show repeat events. Every repeating signal we could classify is either an
eclipsing binary on the target (Gaia DR3 spectroscopic orbits independently reproduce three of
our periods) or a blend with a nearby star. For CTOI 147035057.01 (TYC 5963-01913-1), a third
eclipse in Sector 108 fixes the period at 37.31497 ± 0.00003 d, and Gaia flags the host as
RV-variable. {{N_SINGLE}} candidates remain single, with {{FRAC}} of log-period space (5–1500 d)
now excluded. We list updated classifications to help retire stale entries.

## 1. Sample

- **Starting set:** 450 open CTOIs flagged as single or duo transits, from groups not currently
  working their candidates. Source: ExoFOP, downloaded 2026-10-11.
- **Removed:** 46 with missing parameters, 8 on stars that host confirmed planets, 4 on stars
  with a TOI, and 24 on catalogued TESS eclipsing binaries (Prša et al. 2022). This leaves
  368 CTOIs on 359 stars.
- **Not re-analysed:** {{N_SKIP}} that already carry a period on ExoFOP.
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

{{RESULTS PARAGRAPH: counts by class from the full run}}

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
- **Host:** RV-variable in Gaia DR3 (amplitude ≈ 50 km/s), so likely an eccentric EB.
- **Next eclipse:** 2027-01-20 03:11 UTC.

**Off-target cases.** Blends with nearby EBs dominate the remaining repeating signals.

**Table 1.** Per-CTOI classification: EB on target / blend / period / still single with
P_min allowed.

## 4. Conclusion

Stale single-transit CTOIs with repeat signals are overwhelmingly stellar. Cheap Gaia
binarity checks and difference imaging dispose of most of them. The surviving singles, on
Gaia-clean stars with on-target events, are the residual long-period planet candidates;
{{N_RESIDUAL}} remain.

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
- TESScut: Brasseur et al. 2019
- lightkurve, astroquery, batman (Kreidberg 2015)
