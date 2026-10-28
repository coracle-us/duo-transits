# Local run of CLAUDE_CODE_BRIEF.md (2026-10-27, Jon's MacBook)

Env: `uv venv .venv` (Python 3.13) + `requirements.txt`. Nothing posted, submitted or signed in to.

## Step 1: Kepler alias scan: matches the cloud run
`python kepler_alias_scan.py 770 845` → `kepler_alias_scan_local.log`

    776.17-776.60  [(-4, 'hidden'), (-3, 'hidden')]
    790.50-790.66  [(-4, 'gap'), (-3, 'gap')]
    792.08-792.92  [(-4, 'gap'), (-3, 'hidden')]
    793.25-793.40  [(-4, 'gap'), (-3, 'hidden')]

`kepler_alias_survivors.csv` is byte-identical to the cloud copy.

## Step 2: TESS coverage: no TESS data in window 3
TESSCut (`tess_sectors.log`): sectors 14, 15, 41, 54, 55, 75, 81, 82. Nothing after S82.

tess-point 0.9.5 (pointings through Cycle 10) gives the same list (without S75, an edge case).
It then gives **S120, S121, S126, S127**: mid-times 2027-08-10, 2027-09-05, 2028-01-24
and 2028-02-20. Window 3 (2026 Oct 18-24) falls in **S109** (mid-time 2026-10-18), which
does not point at the star. **No TESS sector tests window 3.** The next TESS look is
Aug-Sep 2027.

## Step 2: ASAS-SN: blocked (server-side maintenance)
- `pyasassn` 0.6.4 (PyPI) is dead. Sky Patrol now returns Parquet and pyarrow dropped
  `deserialize`. The maintained client is the PyPI package **`skypatrol`** (0.6.21), which
  installs the same `pyasassn` module. `requirements.txt` is updated.
- With `skypatrol`, the light-curve download fails. The server reports zero block servers,
  and its banner (`ground_check_asassn_tess.log`) says: *"The light curve service will be
  down for several days while we run mainanence on our cluster."* (Kyle Hart, 13 Oct 2026).
  This is still the case on 27 Oct. Retry later; nothing to fix on our side.
- Bug fixed in `ground_check.py`: the 10″ cone holds **110 ASAS-SN ids**. Only
  60130747136 (0.23″ away, catalogue sources stellar_main/tic/aavsovsx) is the star; the
  rest are Swift-catalogue positions. The old code would have merged all 110 light curves
  into the same per-camera groups. The script now uses a 1″ cone and exits unless exactly
  one id comes back. Not tested end-to-end because the light-curve service is down.

## ZTF (optional): blocked (IRSA outage)
- Positional light-curve queries (`nph_light_curves?POS=...`) fail for every position,
  including IRSA's own documented example: `UsageFault: Internal error (invalid TAP query)`.
- Queries by object ID still work. However, the object-ID lookup (`ztf_objects_dr20`-`dr24`
  via IRSA TAP) fails with `ORA-12541: TNS:no listener`. So the star's ZTF ids can't be
  looked up. `ground_check_ztf.log` holds the resulting `KeyError: 'catflags'`; the
  response was an error VOTable, not CSV.
- Even when IRSA works, public ZTF releases lag by months, so ZTF could only test the
  2021/2023 epochs, not window 3.

## Step 3: AAVSO: Jon's step, not done
No `aavso_kic8462852.csv` yet. Download it while signed in to AAVSO (star "KIC 8462852",
JD 2459400 to today, all bands, CSV) and save it here. Then run:

    .venv/bin/python ground_check.py --aavso aavso_kic8462852.csv --asassn

Drop `--asassn` if Sky Patrol is still down. AAVSO is now the **only** source that can test
window 3.

## Bottom line for the channel
- Kepler survivors are reproduced exactly (776.17-776.60, 790.50-790.66, 792.08-792.92,
  793.25-793.40 d).
- No TESS data in window 3 (S109 misses the star). Next TESS coverage is S120/S121, Aug-Sep 2027.
- ASAS-SN is down for maintenance (since 13 Oct) and IRSA ZTF queries are broken today, so
  there is no per-window summary yet. It needs the AAVSO CSV (window 3) plus a rerun when
  ASAS-SN is back (2021/2023 epochs).

## Later on 2026-10-27 (round 2 session)
- `ground_check.py` statistic fixed (0968e1f). Points in the transit core were treated as
  independent, and AAVSO time-series observers share one nightly zero point across ~100 points
  a night. On synthetic nightly noise the old statistic gave 37% "excluded" and 26%
  "consistent" (false dips). Each group's core is now one mean, with the error taken from that
  group's nightly-mean scatter. The old numbers are still printed for comparison.
  `test_ground_check.py` checks this offline.
- ZTF: IRSA works again, but **ZTF has no source at the star** (past saturation). Within 10"
  there is only a G ~ 19-20 neighbour 8" away. ZTF can't test any window.
- ASAS-SN Sky Patrol: still down (same 13 Oct maintenance banner, no block servers).
- AAVSO: the Photometry (AID) search returned a server error ("An Error Occurred", AAVSO's side).
  The window 3 test is still waiting on `aavso_kic8462852.csv`.

## Rerun 2026-10-28 (brief steps 1 and 2; AAVSO export still down)
- `kepler_alias_scan.py 770 845` (`kepler_alias_scan_2026-10-28.log`): survivors 776.17-776.60,
  790.50-790.66, 792.08-792.92, 793.25-793.40 d, exactly the expected ranges.
- `ground_check.py --asassn --tess --ztf` (`ground_check_2026-10-28.log`): ASAS-SN Sky Patrol
  still in maintenance (banner of Oct 13; `asassn_status.py` says DOWN, no block servers), and
  the client dies with ZeroDivisionError. Rerun with `--tess --ztf`: ZTF returns no rows
  (no source at the star), so no photometry was loaded and every window is 100% untested
  (`epoch_coverage.csv`). The window plots were empty and are not kept. ASAS-SN nights and
  hours in window 3: unknown until Sky Patrol is back.
- TESSCut: sectors 14, 15, 41, 54, 55, 75, 81, 82. Nothing after S82, so none covers window 3
  (S109 misses the star, as above).
