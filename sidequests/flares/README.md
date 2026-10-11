# Side quest: flares on the 85 duotransit host stars

This is separate from the main duotransit period analysis. It uses only light curves already on
disk (`data/raw/lightcurves/`) and writes only inside `sidequests/flares/`.

## Bottom line

- **66 of 85 hosts** have SPOC 2-min data: 190 host-sectors, 4,252 d of good data and
  **3,849 d usable** after removing transit windows and segment edges. 13 of the 66 also have
  20-s data (40 sectors). The other **19 hosts have no 2-min/20-s data** and were skipped (list below).
- 4 raw candidates. After vetting: **1 accepted flare** (TIC 211409161, a K dwarf), 1 marginal
  event that is not counted, and 2 rejections: a Solar System object and a detrending artifact
  at a transit ingress.
- Accepted flare: TIC 211409161, S95, peak BTJD 3885.2053. Amplitude 3.4% (24σ), duration about
  28 min, ED = 16.8 ± 0.6 s, E_TESS ≈ 6.0e33 erg, E_bol ≈ 2.9e34 erg (9000 K blackbody). The
  injection tests show the ED is underestimated by about 27%, so the true E_TESS is closer to
  8e33 erg. The QLP 200-s FFI light curve shows the same event independently (+2.7%).
- No flares were found on the other 65 hosts. The median 90%-completeness limit is E_TESS ≈ 1.7e35 erg
  per star (range 5.9e33 to 2.2e36). Across the sample, no flares were seen above E_TESS = 1e35 erg in
  1,443 star-days from the 20 hosts that are complete at that energy. That gives a rate below 0.76 per
  star per year (95%). Above 10^35.5 erg the limit is below 0.37 per star per year (50 hosts, 2,927 star-days).
- These limits are much weaker than Kepler superflare statistics for solar-type stars.
  Treat this as a host-characterisation paragraph, not a stand-alone result.

## Data and quality handling

- Only `s<NNN>_SPOC_120s.csv.gz` is used for detection, along with `s<NNN>_SPOC_20s.csv.gz`, which is
  overplotted in the vetting figures. Every sector with 20-s data also has 2-min data, so detection
  runs on 2-min data throughout. This keeps the cadence, thresholds and completeness the same for every star. No file has
  sector > 200.
- The flux is PDCSAP, normalised per segment by the trend described below.
- Dropped SPOC bits (`BAD_BITS` in `scripts/flarelib.py`): 1, 2, 4, 8, 16, 32, 128, 2048, 4096,
  16384, 32768. Kept on purpose: 64/1024 (cosmic ray), 256, **512 (impulsive outlier)** and 8192.
  In practice PDCSAP is already NaN wherever the dropped bits are set. 112 good cadences carry bit
  512. 104 of them are single-cadence positive >3σ points, and none is part of a 3-cadence
  detection. **No candidate contains a 512-flagged cadence.**
- Momentum dumps are the times of all cadences with bit 32, including NaN-flux cadences. A
  candidate within 10 min of a dump is rejected automatically. None of the 4 was near a dump
  (nearest was 28 h).
- Transit mask: Tc1 and Tc2 ± 1.5 × Dur_h. Nothing else is masked.
- Segments are split at gaps > 0.5 d. Candidates peaking within 0.5 d of a segment edge are
  rejected, and those edge zones are also left out of the observing time.

## Method

1. **Detrending** (`flarelib.detrend`): per segment, a running median over 6 h with 3 iterations.
   Between iterations, points > 2.5σ above the trend are masked, plus 2 cadences before and
   5 after. The noise σ is a rolling 0.5-d robust (1.4826 × MAD) standard deviation of the
   residuals. A 2-h window is used instead when the residual scatter is more than 1.3× the
   point-to-point scatter and the shorter window reduces it. No host triggered this.
2. **Detection** (`flarelib.find_flares`): Chang et al. (2015) criteria as implemented in
   AltaiPony (Ilin et al. 2021). N1 = 3 means resid/σ > 3. N2 = 2 means (resid − flux_err)/σ > 2.
   N3 = 3 means at least 3 consecutive cadences, i.e. 6 min. Runs ≤ 3 cadences apart are merged.
   Each event is extended back while above 1σ (up to 5 cadences) and forward while above 0.5σ
   (up to 2 h), and overlapping events are merged.
   Automatic shape criterion (fast rise, slower decay): decay time ≥ rise time.
3. **Measurement**: peak time, amplitude (relative flux), duration, and
   ED = Σ resid × 120 s over the event (equivalent duration, e.g. Hawley et al. 2014).
   E_TESS = ED × L_TESS, where L_TESS = 4πR² · π∫B_λ(Teff) dλ over a **600–1000 nm top-hat**.
   This approximates the TESS response; the real response was not on disk. Catalog Teff and Rstar
   are used. E_bol = E_TESS × σT⁴/(π∫B_λ(T)dλ) for a T = 9000 K flare blackbody
   (Shibayama et al. 2013; Günther et al. 2020). The factor is 4.77.
4. **False-alarm check** (`false_alarm.py`): the same detector is run on the sign-flipped
   residuals. It found 8 negative events in 3,849 d, 2 of which pass the edge, dump and shape cuts.
   One is on TIC 211409161 (3.7σ), which is why C001 is called marginal.
5. **Injection-recovery** (`injrec.py`): all 66 hosts and all 190 sectors, with 20 trials × 5 flares
   per sector, for **19,000 injections**. The template is Davenport et al. (2014), integrated over the 2-min
   exposure. Amplitude is log-uniform in 0.5–50 × the sector σ. FWHM is log-uniform in 4–60 min.
   Injection times are random good, non-edge, out-of-transit times. The same detrending,
   detection and automatic cuts are applied. Recovery depends almost entirely on peak
   amplitude: about 30% at 3–5σ, about 85% at 5–8σ and over 95% above 8σ, nearly independent of FWHM (see
   `injection_recovery.png`).
   Recovered ED is biased low: median ratio to injected is 0.73 (68% range 0.5–0.9), because
   the tail is truncated at 0.5σ. ED50 and ED90 are the smallest ED bin (0.2 dex) above which the binned
   recovery stays at or above 50% and 90%. These depend on the assumed FWHM range. Three hosts
   (TIC 32179255, 289840544, 304339227) never stay above 90% within the injected range, so their
   E90 is blank. The pooled curve dips to 92% in its highest-ED bin, from one 50σ flare missed
   on TIC 396720998.
6. **Vetting** (`vet.py`, then a human decision in `manual_vetting.csv`): one plot per candidate
   in `vetting/`, showing PDCSAP + trend, residual/σ, SAP, SAP_BKG, centroids and quality flags,
   with dumps and edge lines marked and 20-s data overlaid where available. The script also
   records background and centroid excursion z-scores. As a cross-star check, it looks at every
   other host with 2-min data in the same sector **and camera** (camera from `tess-point`). A hit
   means ≥ 2 cadences above 3σ within ±10 min. This check is sparse: 0–3 comparison stars per
   candidate.

## Candidates (`candidates.csv`)

| id | TIC | S | peak BTJD | σ | ED (s) | decision | reason (short) |
|---|---|---|---|---|---|---|---|
| C000 | 211409161 | 95 | 3885.2053 | 24.3 | 16.8 | **accept** | classic impulsive rise and exponential decay; also in SAP and QLP; no bkg/centroid change (source within ~4″) |
| C001 | 211409161 | 95 | 3885.9400 | 4.0 | 4.7 | marginal | 3 points; same star shows negative events of similar size; not counted |
| C002 | 237605045 | 33 | 2220.2896 | 5.2 | 28.7 | reject | slow-rise, flat-topped 4.5-h bump after a background spike, with a centroid shift; cam 1 → asteroid |
| C003 | 394340183 | 93 | 3837.7596 | 4.3 | 5.4 | reject | detrending artifact at the ingress of an unmasked transit; fails shape test |

**Note for the main project (C003):** TIC 394340183 shows a full transit in S93. A trapezoid fit
gives mid-time BTJD **3837.926 ± 0.001**, depth about 10 ppt and duration about 8.1 h. The Tc1/Tc2 linear
ephemeris predicts 3837.652 for every alias, so this transit is **6.6 h late**. The ExoFOP third
transit (3113.67) is also about 0.29 d late relative to Tc1 + 2ΔT = 3113.38. The four epochs cannot
be strictly periodic, which points to TTVs or a problem with one of the epochs. This was a quick
check only; it should be verified in the main pipeline.

## Notable stars

- **TIC 211409161** (Teff 4660 K, 0.78 R☉, 1 sector, 18.9 d): the only flaring host.
  Its rate is 0.053 per day (1 flare; 95% upper limit 0.25 per day). The S95 light curve also has
  recurring ~3 ppt dips lasting 1–2 h, spaced about 0.7 d apart, which suggests an active or spotted star or a blend. Its earlier
  sectors (S1, S28, S68) exist only as FFI data and were not searched.
- **TIC 303317324** (late K dwarf, 4159 K, 40 pc, 5 sectors, 94.4 d): no flares. This is the most
  sensitive host, with E90 = 5.9e33 erg (E50 = 2.3e33). Rate < 0.032 per day (95%) above about 6e33 erg,
  so it is a quiet K dwarf.
- **TIC 396720998** (catalog Teff 50,000 K, 0.15 R☉, a hot subdwarf, Tmag 14): no flares in 185 d.
  The noise is high (13,000 ppm per 2 min), so E90 is about 8.8e35 erg. A flare here would point to a
  cool companion or a blend rather than to the subdwarf itself. The non-detection says little.
  The catalog Teff looks like a placeholder value, so the luminosity scaling is uncertain.

## Caveats

- The TESS band is a 600–1000 nm top-hat, not the real response. Expect L_TESS errors of order
  10–20% for FGK stars and more for the hot subdwarf. Catalog Teff and Rstar are used as given.
- The ED underestimate of about 27% (from the injection tests) is not corrected in `flares.csv`.
- Completeness depends on the injected FWHM range (4–60 min) and the template. Flares shorter
  than about 6 min cannot pass N3 = 3 at 2-min cadence.
- The cross-star check is sparse, and contamination is only checked through centroids. There was no
  Gaia neighbour search.
- FFI-only sectors and hosts (19 hosts with no 2-min data) are not searched.

## Files

| file | contents |
|---|---|
| `coverage.csv` | all 85 hosts: 2-min and 20-s sectors, good and usable days |
| `sectors.csv` | per host-sector: camera, noise, usable days, segments, dumps |
| `candidates_raw.csv` | raw detections with automatic flags |
| `candidates_vetmetrics.csv` | bkg/centroid z-scores, cross-star hits, distance to nearest dump |
| `manual_vetting.csv` | human decisions and reasons (input to `summarize.py`) |
| `candidates.csv` | all candidates with decision and reason |
| `flares.csv` | accepted flares with ED and energies |
| `per_star.csv` | per host: count, usable days, rate, Poisson 95% upper limit, L_TESS, analytic and injected ED/E limits |
| `sample_limits.csv` | sample-level rate limits above E_TESS thresholds |
| `false_alarm.csv` | negative-excursion events |
| `injection_recovery.csv`, `injrec_per_star.csv`, `injection_recovery.png` | completeness |
| `vetting/*.png` | candidate plots |

Hosts with no 2-min/20-s data: 9844069 38138512 67599025 71028120 101824521 140750416 145006304
153838604 176518126 207783865 215402824 251057075 256912435 275878706 287204963 332697924
437293313 439491923 442893494.

## Reproduce

```
cd sidequests/flares/scripts
PY=../../../.venv/Scripts/python.exe
$PY run_detection.py   # coverage.csv, sectors.csv, candidates_raw.csv, cache/
$PY vet.py             # vetting/*.png, candidates_vetmetrics.csv
$PY false_alarm.py     # false_alarm.csv
$PY injrec.py          # ~5 min; injection_recovery.*, injrec_per_star.csv (seeded)
$PY summarize.py       # candidates.csv, flares.csv, per_star.csv, sample_limits.csv
```
Uses only numpy, scipy, pandas, matplotlib and tess-point from the main venv. Nothing was installed.
