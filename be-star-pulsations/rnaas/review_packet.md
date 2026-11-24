# Review packet: Be-star RNAAS draft (2026-11-23)

**For:** Mythos, as an independent reviewer.
**From:** Claude Code on Jon's MacBook, at Jon's request.
**Files** (in the relay `results/` folder):
- this packet (`2026-11-23-be-star-rnaas-review-packet.md`);
- the figure (`2026-11-23-be-star-rnaas-fig1.png`).

The manuscript source is at the end of this packet. In the repo it is
`coracle-us/duo-transits`, `be-star-pulsations/rnaas/ms.tex`. It is not committed yet, so the
text below is the version under review.

## What I'd like from you

Act as a sceptical RNAAS moderator and a Be-star specialist. Tell us whether the note should be
submitted, and what has to change first.

**Please don't** contact anyone, post anything, or submit anything.

**Please write your review** to the relay `briefs/` folder as
`2026-11-2x-be-star-rnaas-review.md`. Rank the issues **must fix / should fix / optional**. For each
one:
- quote the sentence;
- say what is wrong;
- propose the replacement text or the analysis needed.

If you think a claim needs a new analysis, write it as a step Claude Code can run. The data are
cached; the scripts are listed in the method notes below.

### Specific questions

1. **Headline and framing.** Is "within a star, variability rises with W1 disk excess" stated at
   the right strength? It is a population tendency: the effect is small next to the scatter, and
   48% of stars have the opposite slope. Is the title too strong?
2. **QLP vs SPOC.** 755 of the 1,084 stars are QLP-only. The cleaner SPOC/TESS-SPOC subset
   (329 stars) confirms B1, red noise and flicker, but not A1 (β 0.17, p 0.11). Does the note
   handle this honestly? Should the headline lead with B1 or flicker instead of A1?
3. **Dilution (k).** I report k = 0 slopes as the headline and call them conservative, given the
   measured k = 0.21. Is that the right choice and the right word? Is the sign argument correct?
   (Disk light adds optical flux, which lowers relative amplitudes. Correcting raises β by
   0.921·k.)
4. **The white-noise control rises too** (β 0.33). The SPOC formal flux error is flat, so I call
   it intrinsic. An alternative is a normalisation effect, which needs the optical to fade as W1
   brightens. The measured k says the optical brightens instead, which argues against it. Is the
   paragraph convincing? Is anything missing, such as scattered light, or a dependence of QLP
   photometry on epoch?
5. **The cross-star swap null** is the main defence against shared sector-era systematics. Is it
   explained clearly enough for a moderator?
6. **Frequency dependence.** This was post hoc. Is it fair to include it, with the caveat as
   written?
7. **Literature and novelty.** The note claims no priority explicitly. It links to
   Labadie-Bartz et al. 2022 and 2025, and reconciles with Niwano et al. 2024. Have I missed a
   paper that already did a within-star TESS × WISE (or other IR) disk test? Pay particular
   attention to 2025–2026 work. My search could not use ADS.
8. **Figure vs table.** RNAAS allows one figure or one table, not both. The figure shows B1 and
   flicker against disk excess. The alternative is the β-vs-frequency curve, which carries the
   frequency result. Which serves the note better?
9. **Phase terms.** I dropped the "dissipating sectors are less variable" claim because it is not
   robust to the interpolation. Building sectors are more variable, but the white-noise control
   shares that under linear interpolation. The discussion now has one sentence saying I do not
   claim a phase effect. Keep it, cut it, or strengthen it?
10. **Anything wrong with the physics.** For example: how Be-star disk growth maps to W1, the
    phase thresholds, the shell-star sign, or where Štefl frequencies and frequency groups sit.

## Context

### The project

Jon McBee is an independent researcher, with Claude doing the coding. Two RNAAS notes are already
with AAS (AAS82207R2 and AAS82225). This would be the third.

RNAAS limits:
- abstract ≤ 150 words;
- ≤ 1,350 counted words, counting everything after the abstract except the references;
- one figure or one table.

The current draft has a 149-word abstract and 1,113 counted words.

### The data

**Catalogue.** Jia et al. 2026 (arXiv:2611.02450), "An all-sky catalogue of Be transition stars
from WISE and Gaia". It lists about 3,300 Be stars whose W1 brightness changed state, with
China-VO W1/W2 light curves. I used the 3,237 stars not flagged as eclipsing binaries.

**TESS.** One product per star: QLP, TESS-SPOC, SPOC 120 s or GSFC-ELEANOR-LITE. All are cached
locally; none is committed.

**Per-sector metrics** (`src/metrics.py`). Each light curve is binned to 30 min and detrended
with a 5-d-knot spline. Then:
- the Lomb–Scargle amplitude spectrum at 0.25–24 c/d;
- A1, the amplitude of the strongest peak;
- B1–B3, the band rms in 0.25–0.8, 0.8–2 and 2–4 c/d;
- red_noise, the median residual amplitude at 0.5–3 c/d after prewhitening up to 10 peaks with
  SNR ≥ 4;
- flicker_rms, the rms of the residuals after 1-d smoothing;
- white_noise, the median amplitude at 15–24 c/d (the planned control).

**Disk excess** (`src/wise.py`).
- W1 visit medians come from China-VO plus later NEOWISE-R visits from IRSA.
- A FITPACK smoothing spline (s = 0.04) is fitted, as in Jia et al.
- excess = W1_quiet − W1_spline(sector mid-time), where W1_quiet is the faintest W1 level.
- A sector is used only if its bracketing visits are ≤ 250 d apart and the spline stays within
  0.1 mag of them. That guard was added after the spline was found to overshoot between sparse
  visits.
- Phase labels come from the spline slope: building < −0.02 mag/100 d, dissipating > +0.02, and
  quiescent where the excess is < 0.05 mag.

**Crowding cut.** CROWDSAP ≥ 0.9 for SPOC products; otherwise the TIC contamination ratio (or a
neighbour proxy) must be ≤ 0.1.

### The statistics (`src/stats4.py`, `src/s4_test.py`)

- **The fit:** within-star demeaned OLS of ln(metric) on excess, giving β per mag.
- **Uncertainties:** a bootstrap over stars (2,000 resamples), quoted as 68% intervals.
- **Within-star permutation null:** shuffle (excess, label) among each star's sectors.
- **Cross-star swap null:** evaluate a donor star's W1 spline at the target star's sector times.
  This keeps the sector-epoch structure and replaces only the disk history.
- 10,000 draws per null; p is floored at 1/10,001.
- **Synthetic test** (`tests/test_synthetic.py`): β recovery, null calibration and spline-guard
  checks all pass.

## Results behind the text

### Main sample: 1,084 stars, 4,108 sectors, k = 0

| Metric | β | 68% interval | ρ | p_perm | p_swap |
|---|---|---|---|---|---|
| A1 | 0.40 | 0.31–0.50 | 0.10 | < 1e-4 | < 1e-4 |
| B1 | 0.80 | 0.70–0.91 | 0.20 | < 1e-4 | < 1e-4 |
| B2 | 0.66 | 0.57–0.75 | 0.16 | < 1e-4 | < 1e-4 |
| B3 | 0.39 | 0.31–0.48 | 0.09 | < 1e-4 | < 1e-4 |
| red_noise | 0.74 | 0.65–0.84 | 0.19 | < 1e-4 | < 1e-4 |
| flicker_rms | 0.84 | 0.73–0.96 | 0.17 | < 1e-4 | < 1e-4 |
| white_noise (control) | 0.33 | 0.24–0.43 | 0.10 | < 1e-4 | < 1e-4 |

### Robustness (β for A1 / B1 / flicker)

| Variant | Stars | A1 | B1 | flicker | Notes |
|---|---|---|---|---|---|
| No crowding cut | 2,449 | 0.29 | 0.59 | 0.62 | all p < 1e-4 |
| Drop mixed-pipeline sectors | 1,011 | 0.44 | 0.82 | 0.85 | |
| Drop W1 < 8 stars | 1,016 | 0.49 | 0.85 | 0.93 | |
| QLP only | 755 | 0.52 | 0.86 | 0.90 | all p < 1e-3 |
| SPOC/TESS-SPOC only | 329 | **0.17 (p 0.11 / 0.15)** | 0.69 | 0.74 | B1 and flicker p < 1e-3 |
| Linear W1 interpolation | 1,086 | 0.44 | 0.87 | 0.89 | all p < 1e-4 |
| Without 4 shell candidates | 1,080 | 0.400 | 0.796 | 0.840 | |
| Without the Be X-ray binary | 1,083 | 0.41 | — | — | |

- **Jackknife** (drop one star): A1 stays within 0.38–0.46.
- **Noise floor:** rms Δln A1 between same-star sectors with |Δexcess| < 0.03 mag is 0.54. The
  brief's preset go/no-go threshold was 0.3. Jon relaxed it and treated this as intrinsic scatter.

### Controls

- **SPOC formal flux error:** β 0.09, p 0.56. QLP's error column tracks light-curve scatter, so
  it is not used as a control.
- **Ratios to white_noise:**

  | Metric | β | p (perm / swap) |
  |---|---|---|
  | B1 | 0.47 | < 1e-3 |
  | red_noise | 0.41 | < 1e-3 |
  | flicker | 0.51 | < 1e-3 |
  | A1 | 0.07 | 0.33 / 0.53 |
  | B3 | 0.06 | 0.37 / 0.61 |

### Frequency dependence (post hoc)

- 30 absolute bins:
  - below 1 c/d: β 0.71–0.92;
  - 1–2 c/d: 0.59–0.76;
  - 2–5 c/d: 0.39–0.51;
  - 5–24 c/d: flat at 0.31–0.37.
- Relative to each star's f1:
  - below 0.7 f1: 0.73–1.16;
  - at 0.9 f1: 0.53;
  - at 1.1 f1: 0.35 (the minimum);
  - above: 0.32–0.48.
- Every bin has p < 1e-3 under both nulls.

### Dilution (pre-submission check 3)

- **Sample:** 125 SPOC/TESS-SPOC stars, 564 sectors.
- **Optical magnitude:** Tmag = −2.5 log10(median SAP × CROWDSAP / FLFRCSAP) + 20.44.
- **Result:** pooled within-star dTmag/dexcess = −0.213 (95% −0.295 to −0.105), so k = 0.21.
  The optical brightens as W1 brightens.
- **Per-star slopes are mostly noise.**
- **Shell candidates** (optical fades as W1 brightens, p < 0.05): 4.
  - TIC 55295028 is strong: 34 sectors, slope +0.91 ± 0.07.
  - The other 3 have 3–4 sectors each.

### Phase terms (pre-submission checks 1 and 2)

Building and dissipating indicators were added to the within-star fit.

| | Spline (cross-star null) | Linear interpolation |
|---|---|---|
| Phase labels agree | — | 69% of sectors |
| Building term | survives for all six metrics (A1 +0.19) | stronger, but white_noise shows it too (+0.18) |
| Dissipating term | survives (A1 −0.09, p_swap 0.012) | mostly vanishes |

### Be X-ray binaries

- The cross-match was against Niwano et al.'s 5 sources plus a list of 33 known Galactic Be X-ray
  binaries from memory. SIMBAD was down, so the list is incomplete.
- 4 matched; 1 is in the main sample: RX J0440.9+4431 = TIC 410336237, 3 sectors, A1 slope −0.86.
- 48% of main-sample stars have negative per-star A1 slopes.
- If each star is negative with probability 0.48, P(≥ 4 of 5 negative) ≈ 0.16.

### Literature check (`results/full/LITERATURE_CHECK.md`)

- 8 web searches, plus forward citations of Labadie-Bartz 2022 and 2025 and of Niwano 2024 via
  OpenAlex and Semantic Scholar (87 works screened).
- 5 arXiv queries.
- No within-star, multi-year, population-scale TESS × IR-disk test was found.
- Closest work: NGC 7419 ZTF + NEOWISE classes (MNRAS 2025); PY Gem / HD 253659 (MNRAS 2026); three
  Be stars with TESS (Universe 2025); OGLE Be demographics (ApJ 2025).
- Not done: an ADS search (no token), and SIMBAD/VizieR (down on 2026-11-23).

## Open items Jon already knows about

- Dataset DOIs for GSFC-ELEANOR-LITE and NEOWISE-R. They are marked TODO in the source.
- The UAT keyword IDs (142, 235, 1625, 2109) are from memory and not yet checked.
- An ADS search, and a complete Be X-ray-binary cross-match.
- The current AAS policy on AI use needs rechecking before submission.
- The PDF will be compiled in Overleaf (aastex701). It has not been compiled yet.
- Every number traces to a file in the repo: `be-star-pulsations/rnaas/numbers.md`.

---

## Manuscript source (ms.tex, verbatim)

```latex
%% Research Note of the AAS (Be-star TESS variability vs WISE disk excess), AASTeX v7.0.1.
%% DRAFT 2026-11-23, not submitted. Provenance of every number: rnaas/numbers.md.
%% Code: https://github.com/coracle-us/duo-transits (be-star-pulsations/).
\documentclass{aastex701}

\begin{document}

\title{Within a Be Star, TESS Variability Rises with the WISE Disk Excess}

\author[orcid=0009-0000-7524-4470]{Jon McBee}
\affiliation{Independent researcher, USA}
\email{jon.mcbee@gmail.com}

\begin{abstract}
Be stars build and lose their decretion disks over months to years, and NEOWISE has tracked the
disk emission of thousands of them. I compare TESS light curves with the NEOWISE W1
light curves of 1,084 Be stars from the catalogue of \citet{Jia2026}, sector by sector and within
each star. When the W1 disk excess is larger, the star is more variable in the optical. The
amplitude of the strongest peak rises by 4\% per 0.1\,mag of excess, and the stochastic power
below 1\,c\,d$^{-1}$ rises by about 8\%. Both trends pass a null test that swaps disk histories
between stars ($p<10^{-4}$). The rise is largest below each star's dominant frequency and does not
single out the coherent pulsation. Within individual stars the trend is weak and about half go the
other way, so this is a population tendency, not a rule for each star.
\end{abstract}

\keywords{\uat{Be stars}{142} --- \uat{Circumstellar disks}{235} --- \uat{Stellar pulsations}{1625} --- \uat{Time domain astronomy}{2109}}

\section{Data and Method}

\citet{Jia2026} identified Be stars whose disks changed state during the WISE/NEOWISE mission
\citep{Wright2010,Mainzer2014}. Of their 3,237 stars not flagged as eclipsing binaries, 2,541 have
at least two TESS sectors within the NEOWISE baseline. For each star I took one product: QLP
\citep{Huang2020} for 2,020 stars, TESS-SPOC \citep{Caldwell2020} for 219, SPOC 120-s for 198 and
GSFC-ELEANOR-LITE for 104. Each sector was binned to 30\,min and detrended with a spline with 5-d
knots. From its amplitude spectrum (0.25--24\,c\,d$^{-1}$) I measured the amplitude $A_1$ of the
strongest peak; the rms in the bands 0.25--0.8, 0.8--2 and 2--4\,c\,d$^{-1}$ ($B_1$--$B_3$); the
median residual amplitude at 0.5--3\,c\,d$^{-1}$ after prewhitening (red noise); the rms of the
1-d-smoothed residuals (flicker); and the median amplitude at 15--24\,c\,d$^{-1}$ (white noise, a
control).

For W1 I combined the catalogue's visit magnitudes with later NEOWISE visits from IRSA and fitted
a smoothing spline as \citet{Jia2026} did. The disk excess of a sector is the spline at mid-sector,
measured from the star's faintest W1 level. A sector is used only if the visits on either side are
$\leq250$\,d apart and the spline stays within 0.1\,mag of them. The main sample keeps sectors with
CROWDSAP $\geq0.9$ (SPOC products) or a contamination ratio $\leq0.1$ (others), and stars with at
least two such sectors: 1,084 stars and 4,108 sectors, 755 of the stars from QLP.

For each metric $m$ I fitted $\ln m$ against disk excess with a per-star mean removed, so that
$\beta$ (per mag) uses only changes within each star. Uncertainties are 68\% bootstrap intervals
over stars. Two nulls, each of 10,000 draws, give $p$-values: a permutation of sectors within each
star, and a swap that gives each star another star's W1 spline, evaluated at its own sector times.
The swap keeps any trend shared by sectors of the same epoch, such as a change in QLP processing, so
it tests whether the star's own disk history matters. All analyses were checked on synthetic data
first.

\section{Results}

Every variability metric rises with disk excess (Figure~\ref{fig:beta}), with $p<10^{-4}$ under both
nulls. In the main sample $\beta=0.40\,[0.31,0.50]$ for $A_1$, $0.80\,[0.70,0.91]$ for $B_1$,
$0.66$ for $B_2$, $0.39$ for $B_3$, $0.74$ for red noise and $0.84\,[0.73,0.96]$ for flicker. The
swap nulls centre on zero (means between $-0.03$ and $+0.02$). The slopes hold without the crowding
cut (2,449 stars; $A_1$ 0.29, $B_1$ 0.59) and when dropping any one star, mixed-pipeline sectors, or
stars brighter than $W1=8$. They also hold with straight-line interpolation of W1 instead of the
spline, which raises every $\beta$ by 4--11\%.

Repeating the fit in 30 frequency bins shows where the rise sits. Below 1\,c\,d$^{-1}$
$\beta=0.71$--0.92; above 5\,c\,d$^{-1}$ it is flat at 0.31--0.37. Measured relative to each star's
dominant frequency $f_1$, $\beta$ is 0.73--1.16 below $0.7f_1$ and 0.35 at the dominant group
itself, the same as at high frequency. The coherent pulsation therefore follows the broadband level,
while the low-frequency stochastic signal rises about 2.5 times faster.

The white-noise control also rises ($\beta=0.33$). The formal flux error of the SPOC and TESS-SPOC
products does not ($\beta=0.09$, $p=0.56$), so this is not an instrumental noise floor. Relative to
the white noise, $B_1$, red noise and flicker still rise ($\beta=0.41$--0.51, $p<10^{-3}$), but $A_1$
does not ($0.07$, $p\geq0.33$).

Light from the disk dilutes relative amplitudes. A star whose $T$ magnitude changes by $-k$ per mag
of W1 excess has every $\beta$ raised by $0.921k$ once corrected. From SAP photometry of 125 SPOC and
TESS-SPOC stars I measure $k=0.21$ (95\% range 0.11--0.30), which raises $\beta$ to about 0.6 for $A_1$
and 1.0 for $B_1$. The uncorrected values above are therefore conservative. Four stars, one of them
(TIC~55295028) significantly, fade in the optical as W1 brightens, as shell stars should. Removing
them changes $\beta$ for $A_1$, $B_1$ and flicker by less than 0.01.

Two caveats limit what the data show. First, the effect is small against the scatter: at fixed disk
excess, $\ln A_1$ varies by 0.54 rms between sectors of the same star, the within-star rank
correlations are 0.09--0.20, and 48\% of stars have a negative $A_1$ slope of their own. Second,
the cleaner SPOC and TESS-SPOC photometry (329 stars) confirms the low-frequency trends ($B_1$
0.69, flicker 0.74, $p<10^{-3}$) but not $A_1$ ($0.17$, $p\geq0.11$). The frequency analysis was chosen
after seeing the first result, and I did not correct the $p$-values for the 7 metrics and several
samples; the main results sit at the $10^{-4}$ floor.

\begin{figure}[ht!]
\plotone{rnaas_fig1.pdf}
\caption{Within-star change in (a) $\ln B_1$ (0.25--0.8\,c\,d$^{-1}$ rms) and (b) $\ln$ flicker rms
against W1 disk excess, for 4,108 sectors of 1,084 Be stars. Both axes have each star's mean
removed. Symbols give the W1 disk phase from the spline slope. Squares are means in ten quantile
bins with standard errors; the line and band show the pooled $\beta$ and its 68\% bootstrap
interval. A few sectors lie outside the plotted range.\label{fig:beta}}
\end{figure}

\section{Discussion}

The rise sits where \citet{LabadieBartz2022,LabadieBartz2025} found extra power during individual
outbursts lasting days: on the low-frequency side of the frequency groups and in stochastic
variability. A plausible reading is that a star with a growing or dense disk has more ejection
episodes per sector. This note extends that picture to disk states lasting years and to a sample
of 1,084 stars. \citet{Niwano2024} found the main pulsation amplitude falling as the disk grew in
four of five Be X-ray binaries. One of them, RX~J0440.9+4431, is in my sample and also has a
negative slope. With 48\% of single stars going that way, four of five is consistent with the
positive population mean. Sectors where the W1 spline is brightening are more variable than the
slope alone predicts, but a sector's disk phase depends on how the W1 curve is interpolated, so I
do not claim a phase effect. Code and per-sector tables are at
\url{https://github.com/coracle-us/duo-transits}.

\begin{acknowledgments}
This work uses TESS data from MAST: SPOC
\dataset[doi:10.17909/t9-nmc8-f686]{https://doi.org/10.17909/t9-nmc8-f686}, TESS-SPOC
\dataset[doi:10.17909/t9-wpz1-8s54]{https://doi.org/10.17909/t9-wpz1-8s54}, QLP
\dataset[doi:10.17909/t9-r086-e880]{https://doi.org/10.17909/t9-r086-e880} and GSFC-ELEANOR-LITE
%% TODO verify the GSFC-ELEANOR-LITE HLSP DOI before submission
light curves; NEOWISE-R data from IRSA; and the WISE light curves of \citet{Jia2026} from China-VO.
%% TODO add the IRSA NEOWISE-R dataset DOI after verifying it
\end{acknowledgments}

\facilities{TESS, WISE, NEOWISE}
\software{lightkurve \citep{Lightkurve2018}, astroquery \citep{Ginsburg2019}, astropy, scipy}

\section*{AI disclosure statement}
I used Anthropic's Claude, an AI assistant, under my supervision to write and run the analysis code
(including synthetic-data tests and robustness checks) and to draft this note. I checked every number
against the saved outputs and every reference against arXiv or Crossref, and take full
responsibility for the content.

\begin{thebibliography}{}
\bibitem[Caldwell et al.(2020)]{Caldwell2020} Caldwell, D.~A., Tenenbaum, P., Twicken, J.~D., et al. 2020, RNAAS, 4, 201
\bibitem[Ginsburg et al.(2019)]{Ginsburg2019} Ginsburg, A., Sip{\H o}cz, B.~M., Brasseur, C.~E., et al. 2019, AJ, 157, 98
\bibitem[Huang et al.(2020)]{Huang2020} Huang, C.~X., Vanderburg, A., P\'al, A., et al. 2020, RNAAS, 4, 204
\bibitem[Jia et al.(2026)]{Jia2026} Jia, Q., Jian, M., Li, G., Cui, Q., \& Li, C. 2026, arXiv:2611.02450
\bibitem[Labadie-Bartz et al.(2022)]{LabadieBartz2022} Labadie-Bartz, J., Carciofi, A.~C., Henrique de Amorim, T., et al. 2022, AJ, 163, 226
\bibitem[Labadie-Bartz et al.(2025)]{LabadieBartz2025} Labadie-Bartz, J., Carciofi, A.~C., Rubio, A.~C., et al. 2025, A\&A, 699, A82
\bibitem[Lightkurve Collaboration et al.(2018)]{Lightkurve2018} Lightkurve Collaboration, Cardoso, J.~V.~d.~M., Hedges, C., et al. 2018, Lightkurve: Kepler and TESS time series analysis in Python, Astrophysics Source Code Library, ascl:1812.013
\bibitem[Mainzer et al.(2014)]{Mainzer2014} Mainzer, A., Bauer, J., Cutri, R.~M., et al. 2014, ApJ, 792, 30
\bibitem[Niwano et al.(2024)]{Niwano2024} Niwano, M., Fausnaugh, M.~M., Lau, R.~M., et al. 2024, MNRAS, 534, 2168
\bibitem[Wright et al.(2010)]{Wright2010} Wright, E.~L., Eisenhardt, P.~R.~M., Mainzer, A.~K., et al. 2010, AJ, 140, 1868
\end{thebibliography}

\end{document}
```
