# Be star pulsations vs WISE disk phase: literature check (2026-11-22, extended 2026-11-23)

**Question.** Has anyone already shown that, within a Be star, TESS variability changes with
the disk state? In particular, has anyone shown the follow-up result
(`RESULTS_FOLLOWUP.md` §1)?
- Low-frequency stochastic power below the star's dominant frequency grows with W1 excess,
  at about +8% per 0.1 mag.
- All frequencies grow at about +3% per 0.1 mag.
- The dominant pulsation group sits at that broadband level.

**Answer: no prior study found.** The test looks new. It is consistent with the day-scale
outburst results of Labadie-Bartz et al. The disagreement with Niwano et al. 2024 mostly
goes away once single-star scatter is taken into account (§4).

## 1. How the check was done

**First pass (2026-11-22)**
- 8 web searches.
- Papers read:
  - Jia et al. 2026: HTML.
  - Niwano et al. 2024: HTML, partly as a summary.
  - Labadie-Bartz et al. 2022 and Balona & Ozuyar 2020: PDF, keyword passages.
  - Labadie-Bartz et al. 2025: HTML, partly as a summary.
  - Balona & Ozuyar 2019: abstract.

**Second pass (2026-11-23)**
- **Forward citations**, all read-only and anonymous:
  - OpenAlex: works citing Labadie-Bartz et al. 2022 (80) and 2025 (19).
  - OpenAlex and Semantic Scholar: works citing Niwano et al. 2024. OpenAlex counts 1 citing
    work but returned none; Semantic Scholar returned 0.
  - This gives 87 unique citing works. I screened titles and abstracts by keyword and read
    the abstracts of the 10 closest.
  - Jia et al. 2026 is three weeks old, with no citations yet.
- **arXiv abstract search**, 5 queries combining "Be star(s)" with: TESS + WISE;
  TESS + infrared (excess); pulsation + disk + infrared; low-frequency + disk + TESS.
- **Full reads:**
  - Labadie-Bartz et al. 2022 §5.2.4 (flickers and frequency groups), plus the §6.1 and §7
    passages.
  - Labadie-Bartz et al. 2025 §5.4 (stellar vs circumstellar photometric signals).
- **Not done:**
  - an ADS search; ADS needs an account token;
  - a SIMBAD/VizieR query: CDS returned "page unavailable" all day on 2026-11-23.

**What the searches found.** No citing paper and no arXiv hit tests TESS variability
against a multi-year infrared disk record within stars across a sample. The closest are:
- NGC 7419 (MNRAS 2025): ZTF + NEOWISE variability classes for about 49 Be stars, with no
  within-star coupling test.
- PY Gem / HD 253659 (MNRAS 2026): two stars, TESS flickers plus colour loops.
- "Exploring the variability of three Be stars with TESS" (Universe 2025): Hα EW
  anticorrelates with brightness in one star.
- OGLE Be demographics (ApJ 2025): disk duty cycles from optical light curves only.

## 2. Closest prior work

| Paper | Sample, data | Timescale | Finding relevant to ours |
|---|---|---|---|
| **Niwano et al. 2024** (MNRAS, doi:10.1093/mnras/stae2160; arXiv:2409.09581) | 5 Be X-ray binaries. TESS FFI (difference imaging), 2–3 sectors each. Gattini-IR J band as the disk proxy. | Sector to sector | Main coherent peaks (2–4 c/d) **anticorrelate** with J excess in 4 of 5. Red noise is removed, not tested. No significance: "too small a sample size". |
| **Labadie-Bartz et al. 2022** (AJ 163, 226) | 432 Be stars, TESS southern ecliptic | Days | 83% (62/75) of flickering stars show groups emerging or strengthening at the flicker. Frequencies about 10% below the dominant one (Štefl frequencies) are tied to recently ejected circumstellar material. Circumstellar signals are expected to be stronger at high inclination (Baade et al. 2016). |
| **Labadie-Bartz et al. 2025** (A&A 699, A82; "Birth of Be star disks I") | About 30 outbursts in 13 stars, TESS plus spectroscopy | Days | Extra power mostly on the low-frequency side of groups during outbursts. Circumstellar material may contribute but "cannot be the sole explanation". ω CMa has a circumstellar 0.67 c/d signal that stays coherent during high states, when the inner disk is densest. |
| **Huat et al. 2009** (HD 49330, CoRoT; via reviews) | 1 star | Days | p-mode amplitudes drop and g-mode groups rise through an outburst. |
| **Balona & Ozuyar 2019, 2020** | 57 and 441 Be stars, TESS | Per sector | The drop in β Cep amplitude in outbursts is read as obscuration. No disk-strength test. |
| **Jia et al. 2026** (arXiv:2611.02450; the source catalogue) | 3,300 WISE Be transition stars | — | TESS is used only for eclipsing binaries. |

## 3. How our result sits against it

1. **The test is probably new.** It is many stars (1,084), within-star, over years of
   sectors, against the multi-year WISE record, with cross-star nulls. Nothing in the
   citation trees or arXiv searches does this.
2. **It is consistent with the day-scale outburst picture, and adds a frequency detail.**
   - The low-frequency excess below f1 matches the "low-frequency side" enhancement during
     outbursts (Labadie-Bartz et al. 2025).
   - In our relative-frequency curve, β is 0.53 at f/f1 ≈ 0.9 against 0.35 at f/f1 ≈ 1.1. That
     points the same way as Štefl frequencies, which lie about 10% below the dominant
     frequency and are tied to fresh ejecta. It is suggestive only: the bins are 0.1 dex
     wide, and the analysis is post hoc.
   - Simplest reading: a star with a growing or dense disk has more ejection episodes per
     sector, and our sector-averaged spectra pick these up.
   - Framing for a note: an extension of the event-scale results to the disk-state, multi-year
     timescale across a large sample. Not a new mechanism.
3. **Inclination.** Circumstellar signals should be stronger at high inclination (Baade et
   al. 2016). Labadie-Bartz et al. 2025 also see the same behaviour in pole-on stars. Our v
   sin i split finds no difference, but it is weak (±0.3 per mag). So it neither supports
   nor rules out a circumstellar share.

## 4. Be X-ray binaries in our sample (`src/s7_bexrb_xmatch.py`, `src/s7b_bexrb_slopes.py`)

**Cross-match**
- Inputs: Niwano et al.'s 5 sources plus an ad hoc list of 33 well-known Galactic Be X-ray
  binaries.
  - The list is from memory and not complete. With SIMBAD down there was no full list.
  - Each name was resolved with the MAST name lookup (36 of 38 resolved) and cone-matched
    within 5″ to the Jia et al. catalogue.
  - Some resolved positions are X-ray positions (arcmin errors), so a few true members could
    be missed.
- **4 matches**, all passing the sample cuts:

| Source | TIC | Main sample | Sectors | W1-excess range (mag) | A1 per-star slope (per mag) |
|---|---|---|---|---|---|
| **RX J0440.9+4431** (one of Niwano's 5) | 410336237 | yes | 3 (S19, S59, S73; TESS-SPOC) | 0.32 | **−0.86** |
| SAX J2239.3+6116 | 455949380 | no (crowding) | 6, QLP | 0.17 | −7.1 |
| GRO J1008−57 | 462886452 | no (crowding) | 6, QLP | 0.04 (too small to mean anything) | −18.8 |
| Swift J0243.6+6124 | 50518138 | no (crowding) | 2, QLP | 0.17 | +10.5 |

**RX J0440.9+4431**
- A1 falls from 8.8 to 6.3–7.0 ppt as W1 excess rises from 0.50 to 0.82 mag. That is the
  same sign as Niwano et al., and likely partly the same TESS sectors.
- The stage-4 dilution correction adds 0.921·k per mag to the slope. On the stage-4 grid
  (k ≤ 0.5) the slope stays negative (−0.40 at k = 0.5). It would only flip at k ≥ 0.93.

**Context from single-star scatter**
- In the main sample, **48% of stars have a negative per-star A1 slope**.
  - The median per-star slope is 0.17. The pooled β is 0.40.
  - The pooled estimate weights stars by their disk swing, so stars with large swings carry
    it.
- Against the main-sample distribution:
  - RX J0440 sits at the 38th percentile;
  - SAX J2239 at the 15th;
  - Swift J0243 at the 90th;
  - GRO J1008, at the 7th, has too little disk change to count.
- For individual stars with 2–3 sectors, a negative slope is unremarkable.
  - If each star is negative with probability 0.48, 4 or more of 5 negative happens about 16%
    of the time.
  - So Niwano et al.'s 4-of-5 is consistent with our positive population mean.
- **The disagreement mostly goes away.** Their result and ours can both hold.
  - A within-star population signal does not imply every star follows it.
  - Niwano et al.'s small sample could not test the population mean.

**Robustness.** Removing the matched star from the main sample leaves A1 unchanged:
β = 0.41 [0.31, 0.50], p < 10⁻⁴ under both nulls, 1,083 stars.

## 5. Verdict

- **Novelty holds** after the citation pass. A short RNAAS note is defensible, framed as:
  - the first population-scale, within-star link between TESS variability and the WISE disk
    record;
  - concentrated at low frequency, below the dominant pulsation;
  - consistent with the day-scale outburst results of Labadie-Bartz et al.;
  - compatible with the Be X-ray binary anticorrelations of Niwano et al. once single-star
    scatter is allowed for.
- **The note should say:**
  - the per-frequency analysis is post hoc;
  - QLP dominates, and the SPOC-only A1 is not significant (p about 0.1);
  - the broadband rise is not separated from optical normalisation;
  - about half of individual stars go the other way.
- **Still open:**
  - an ADS search;
  - a complete X-ray-binary cross-match once SIMBAD is back.
  
  Neither is likely to change the verdict.

## Sources

- Jia et al. 2026: https://arxiv.org/abs/2611.02450
- Niwano et al. 2024: https://arxiv.org/abs/2409.09581 (MNRAS, doi:10.1093/mnras/stae2160)
- Labadie-Bartz et al. 2022: https://arxiv.org/abs/2010.13905 (AJ 163, 226)
- Labadie-Bartz et al. 2025: https://arxiv.org/abs/2504.07571 (A&A 699, A82)
- Balona & Ozuyar 2019: https://arxiv.org/abs/1911.03068
- Balona & Ozuyar 2020: https://arxiv.org/abs/2008.06288
- Huat et al. 2009: https://www.aanda.org/articles/aa/abs/2009/40/aa11928-09/aa11928-09.html
- NGC 7419 (2025): https://doi.org/10.1093/mnras/staf1104
- PY Gem / HD 253659 (2026): https://arxiv.org/abs/2609.28889
- Three Be stars with TESS (2025): https://doi.org/10.3390/universe11020071
- OGLE Be demographics (2025): https://doi.org/10.3847/1538-4357/ae0612
- Citation lists: OpenAlex API (works citing W3096488769, W4409026914, W4402679443).
