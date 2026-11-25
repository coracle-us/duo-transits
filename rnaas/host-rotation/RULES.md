# Host rotation retest: rules (written 2026-11-24, before any Santos number was computed)

Question: does the result of Binnenfeld, Sussholz, Faigler & Mazeh (arXiv 2611.08267; "B26" below) hold when the
rotation periods come from Santos et al. (2019, 2021) instead of McQuillan, Mazeh & Aigrain (2013, 2014; "MMA")?
The Santos result is reported whatever it shows. Every choice below was fixed before the Santos run; anything
changed afterwards is listed at the end of this file with its date and reason.

## Inputs (`fetch.py`; raw files in `data/`, not in git)

- **MMA:** J/ApJ/775/L11 table1 (KOIs; rows with Prot > 0) and J/ApJS/211/24 table1 (periods) and table2 (no
  period). B26 says the MMA sample combines the two. They share no KIC, so the union is used as is.
- **Santos:** J/ApJS/244/21 table3 (period) and table4 (no period); J/ApJS/255/17 table1 (period) and table2 (no
  period). They share no KIC. Primary: every published Prot. Sensitivity: drop stars with the CP/CB-candidate flag
  set (Fl1 = 1 in 2019, flag1 = 1 in 2021).
- **Gaia DR3** (B26 §2: GSP-Phot parameters, parallaxes, photometry, extinction): `teff_gspphot`, `mh_gspphot`,
  `ag_gspphot`, `ebpminrp_gspphot`, `parallax`, `phot_g_mean_mag`, `bp_rp`. Counterpart: the nearest Gaia DR3
  source with G < 17.5 within 1" of the KIC position, after moving Gaia positions back to J2000 with their proper
  motions; a KIC star with no such source is dropped. The Gaia sources come from a download of the whole Kepler
  field (`fetch.py gaia`), matched locally. (gaia-kepler.fun failed its TLS check, CDS XMatch was down, and the ESA
  archive timed out on joins with uploaded tables on 2026-11-24, so the match is our own; B26 matched by KIC, method
  not stated. This paragraph was revised on 2026-11-24 before any analysis had run, MMA or Santos.)
- **Hosts:** NASA Exoplanet Archive KOI cumulative table (accessed 2026-11-24; B26 used 2026-09-04). B26 says
  "stars associated with planets in the NASA Exoplanet Archive" without saying whether candidates count, so both
  are run: (C) confirmed only, koi_disposition = CONFIRMED; (C+K) CONFIRMED or CANDIDATE. Stars whose only KOIs are
  FALSE POSITIVE are non-hosts in both, as in B26 ("all remaining selected MMA stars are used as non-hosts"). The
  definition that reproduces B26's counts on MMA (step 2) is the primary one for Santos; the other is reported too.

## Cuts (B26 §2)

1. 4500 ≤ Teff_GSP-Phot ≤ 6000 K.
2. Main sequence: extinction-corrected CMD, M_G0 = G − A_G + 5 log10(ϖ/mas) − 10 against (BP−RP)0 = BP−RP −
   E(BP−RP). The main-sequence locus is a cubic in (BP−RP)0 fitted to the median M_G0 in 0.05-mag colour bins of
   all MMA-analysed Kepler targets (table1 + table2 + KOIs) that pass cut 1, iterated three times keeping stars
   within ±1 mag. A star is kept if |M_G0 − locus| ≤ 1 mag. Stars without parallax > 0 or extinction are dropped.
3. Distance d = 1000/ϖ ≤ 1300 pc (sensitivity: GSP-Phot distance).
4. [M/H]_GSP-Phot > −0.4.
5. Period ranges: fit on 1 < P ≤ 30 d; compare 1 < P ≤ 20 d and 20 < P ≤ 30 d.

## Statistic

- Expected period: ln P_exp = a0 + a1 x + a2 x² + a3 x³, x = (Teff − 5500 K)/1000 K, fitted as a median (L1,
  quantile 0.5) regression of ln P on x for all non-hosts with 1 < P ≤ 30 d in that catalogue's sample. B26 fitted
  "the median non-host rotation period"; median regression is our reading of that. Coefficient errors: standard
  deviation over 2,000 bootstrap resamples of the non-hosts.
- r = P_obs / P_exp. Offset = median r(hosts) / median r(non-hosts) − 1.
- CI: 2,000 bootstrap resamples, hosts and non-hosts resampled separately and the cubic refitted on each resampled
  non-host set; report the 16–84% and 2.5–97.5% intervals. Two-sided KS test of r (scipy.stats.ks_2samp).
- Santos: refit the cubic on the Santos non-hosts and also apply B26's published coefficients; report both.

## Step 2 pass criterion (MMA reproduction)

Host and non-host counts for 1 < P ≤ 20 d within 10% of 188 and 7,274, and the offset inside B26's interval,
13.7 (+2.2/−5.0)% (that is, 8.7–15.9%). If not: stop, report, and run nothing on Santos.

## Step 3 subsamples (Santos)

1. All Santos stars passing the cuts.
2. Stars with a period in both catalogues: the same test on each catalogue's periods, and P_Santos/P_MMA compared
   directly for hosts and non-hosts (medians, KS).
3. Stars with a Santos period and no MMA period.
Plus 20 < P ≤ 30 d for 1, and the detection-bias check: per 250-K Teff bin, the fraction of hosts and of non-hosts
that have a period, out of all stars each catalogue analysed (with and without a period) that pass cuts 1–4.

## Changes after the Santos run

(none yet)

## Step 2 outcome (2026-11-24): did not pass; Santos step not run

`results/mma.txt`, `results/mma.json`, `results/mma_cutflow.txt`. Confirmed-only hosts: 192 hosts (B26 188) but
8,992 non-hosts (B26 7,274, +24%); offset +16.7% (B26 interval 8.7-15.9%). Confirmed + candidates: 269 hosts
(+43%), offset +7.7%. Neither passes. Stopped here for Jon, as the brief asks.

## Step 2 diagnosis (2026-11-25): gate still not passed

`DIAGNOSIS.md`, `diagnose.py`, `results/diagnose.*`. Twelve single changes tested against all of B26's Kepler
checkpoints; none explains the +24% non-host excess consistently (a G < 15.5 limit passes the gate numerically but
breaks the host checkpoints and is not adopted). No rule above was changed. Santos step still held.
