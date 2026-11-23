# Be star pulsations vs WISE disk phase: literature check (2026-11-22)

**Question.** Has anyone already shown that, within a Be star, TESS variability changes with
the disk state? In particular, has anyone shown the follow-up result
(`RESULTS_FOLLOWUP.md` §1)?
- Low-frequency stochastic power below the star's dominant frequency grows with W1 excess,
  at about +8% per 0.1 mag.
- All frequencies grow at about +3% per 0.1 mag.
- The dominant pulsation group sits at that broadband level.

**How the check was done**
- Web searches (8).
- Papers read:
  - Jia et al. 2026: HTML.
  - Niwano et al. 2024: HTML, partly as a summary.
  - Labadie-Bartz et al. 2022: PDF, keyword passages.
  - Balona & Ozuyar 2020: PDF, keyword passages.
  - Labadie-Bartz et al. 2025: HTML, first ~40k characters, the rest as a summary.
  - Balona & Ozuyar 2019: abstract.
  - Huat et al. 2009: via reviews, not read directly.
- There was no ADS citation search. This is a first pass, not an exhaustive review.

## Closest prior work

| Paper | Sample, data | Timescale | Finding relevant to ours |
|---|---|---|---|
| **Niwano et al. 2024** (arXiv:2409.09581, MNRAS) | 5 Be X-ray binaries. TESS FFI (difference imaging), 2–3 sectors each. Palomar Gattini-IR J band as the disk proxy. | Sector to sector | Amplitudes of the **main coherent peaks** (2–4 c/d) **anticorrelate** with J excess in 4 of 5. Max/min amplitude ratio is typically 5–6. Red noise is modelled and removed; it is not tested against the disk. No significance: "too small a sample size". |
| **Labadie-Bartz et al. 2022** (AJ 163, 226; arXiv:2010.13905) | 432 Be stars, TESS southern ecliptic. | Days, within a sector | Flickers (short mass-ejection events) occur in 30% of early types. In 82% of flickering stars, the frequency-group amplitude rises during the flicker. Future work they list includes dependence on inclination and on the strength [of the disk]. No WISE. |
| **Labadie-Bartz et al. 2025** (arXiv:2504.07571, "Birth of Be star disks I") | About 30 outbursts in 13 Be stars, TESS plus spectroscopy. | Days | Changes to frequency groups happen during or near outbursts, with extra power **mostly on the low-frequency side** of a group, lasting days. They say TESS is "poorly suited" to events of a month or longer. |
| **Huat et al. 2009** (HD 49330, CoRoT; via reviews) | 1 star | Days | During an outburst, the β Cep p-mode amplitudes drop, and g-mode groups appear and peak. |
| **Balona & Ozuyar 2019, 2020** (arXiv:1911.03068, 2008.06288) | 57 and 441 Be stars, TESS | Per sector | Groups are read as rotational or non-coherent. The drop in β Cep amplitude during outbursts is read as obscuration. No disk-strength test. |
| **Jia et al. 2026** (arXiv:2611.02450; the source catalogue) | 3,300 WISE Be transition stars | — | TESS is used **only to find eclipsing binaries** (63). There is no pulsation or variability analysis and no TESS-versus-disk-phase comparison. |

Searches for a population, within-star study of TESS variability against a multi-year IR
disk measure (WISE/NEOWISE or other) found none.

## How our result sits against it

1. **Novelty of the test: probably new.** No study found does what this one does:
   - many stars (1,084), within-star;
   - years of sectors against a multi-year infrared disk record;
   - with cross-star nulls.
   
   The source catalogue does not attempt it. Labadie-Bartz et al. 2022 lists disk-strength
   dependence as future work.
2. **Consistent with the event-scale picture.**
   - Labadie-Bartz et al. 2022 and 2025 find low-frequency power and group power rising
     during mass-ejection events lasting days.
   - Our low-frequency excess below f1 could be the long-term average of those events.
     Stars with growing disks should have more flickers per sector.
   - That would make our result an extension to disk-state timescales, not new physics.
   - A note should say this, not claim discovery of the mechanism.
3. **In tension with Niwano et al. 2024.**
   - They find coherent mode amplitudes falling as the disk grows.
   - We find A1 rising with W1 excess (β = 0.40 per mag, main sample).
   - Possible reasons:
     - The populations differ: their stars are X-ray binaries in giant outburst.
     - The normalisation differs: they use a fixed flux from Tmag; we use each sector's median.
     - The IR band differs: J against W1.
     - Their sample is 5 stars without significance tests.
     - Our cleaner SPOC-only A1 is itself not significant (p about 0.1).
   - The tension is worth one sentence, and possibly a check: whether any of their 5 stars,
     or other Be X-ray binaries, are in our sample.
4. **"Coherent pulsation not singled out" fits Balona & Ozuyar and Huat et al.**
   - Balona & Ozuyar read β Cep amplitude drops as obscuration.
   - Huat et al. saw p modes drop while the g groups rose.
   - Neither contradicts ours.

## Verdict

- **A short note is defensible on novelty**, if it is framed this way:
  - the first population-scale, within-star link between TESS variability and the WISE disk
    record;
  - consistent with the event-scale results of Labadie-Bartz et al.;
  - in tension with Niwano et al.'s small Be X-ray binary sample.
- **Before drafting:**
  1. An ADS forward-citation pass on Labadie-Bartz et al. 2022 and Niwano et al. 2024 for
     anything from 2025–26 that this first pass missed.
  2. Read Labadie-Bartz et al. 2022 §6 and Labadie-Bartz et al. 2025 §5 in full; only
     keyword passages and summaries were read here.
  3. Optional: cross-match our sample with known Be X-ray binaries and Niwano's 5 stars.
- **Unchanged from the follow-up write-up:**
  - The per-frequency analysis is post hoc.
  - QLP dominates the sample.
  - The rise at all frequencies is not separated from optical normalisation.

## Sources

- Jia et al. 2026: https://arxiv.org/abs/2611.02450
- Niwano et al. 2024: https://arxiv.org/abs/2409.09581
- Labadie-Bartz et al. 2022: https://arxiv.org/abs/2010.13905 (AJ 163, 226)
- Labadie-Bartz et al. 2025: https://arxiv.org/abs/2504.07571
- Balona & Ozuyar 2019: https://arxiv.org/abs/1911.03068
- Balona & Ozuyar 2020: https://arxiv.org/abs/2008.06288
- Huat et al. 2009: https://www.aanda.org/articles/aa/abs/2009/40/aa11928-09/aa11928-09.html
