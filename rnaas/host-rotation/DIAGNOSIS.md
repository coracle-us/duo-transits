# Why the MMA reproduction has too many non-hosts (2026-11-25)

Script: `diagnose.py` → `results/diagnose.txt`, `results/diagnose.json`. Each variant changes **one** thing in the
RULES.md baseline. It is scored against every Kepler-sample number B26 gives:

- the region before the period cut: 12,525 stars, 400 hosts;
- metal-rich stars beyond 1.3 kpc: 902 stars, 9 hosts;
- metal-poor stars within 1.3 kpc: 1,398 stars, 18 hosts;
- 1–20 d: 188 hosts, 7,274 non-hosts;
- median r: 0.909 for hosts, 0.800 for non-hosts;
- the offset;
- the KS test;
- the cubic coefficients;
- Fig. 2.

Confirmed-only hosts throughout. Nothing was run on Santos, and the authors were not contacted.

## Result: the reproduction still fails the step-2 gate, and the cause is not found

The baseline fails the gate:
- 192 hosts in 1–20 d (B26 188, +2%);
- 8,992 non-hosts (B26 7,274, +24%);
- offset +16.7% (gate 8.7–15.9%).

None of the cuts below fixes it in a way consistent with B26's other numbers.

## What the excess looks like

1. **It is in non-hosts only.** Our host counts match all three of B26's region checkpoints: 406/400, 10/9 and
   18/18. The non-host counts are +36%, +52% and +109%.
2. **The cubic still matches.** Our median cubic, fitted to the 1–30 d non-hosts, is within 1σ of B26's eq. 2.
   B26's non-hosts therefore follow the same period–temperature relation as ours; B26 simply has fewer of them.
3. **It does not depend on distance.** From Fig. 2, our stars per B26-implied star are ×1.25–1.48 in every distance
   bin from 0 to 1,800 pc. A brightness or parallax effect would grow with distance.
4. **It peaks at extreme metallicity.** By [M/H] bin, the excess is ×1.24–1.26 for −0.4 to 0.0. It is ×1.4–2.0 for
   0.0 to +0.8, where our non-hosts are cool, slow (median P 21–31 d) K dwarfs. It is ×1.7–3.1 below −0.4, where
   27–42% of our non-hosts have RUWE > 1.4 (unresolved binaries get low GSP-Phot [M/H]).
5. **It leans to long periods.** In B26, 60% of the region's non-hosts are at 1–20 d; in ours, 53%. We have
   about 2,100 non-hosts with P > 30 d.

## Variants tested, one at a time

| # | Change | 1–20 d hosts/non-hosts | Offset | Verdict |
|---|---|---|---|---|
| 0 | baseline | 192 / 8,992 | +16.7% | fail |
| 1 | KOI candidates out of the non-hosts (what B26 did for TARS) | 192 / 8,915 | +16.7% | removes only 1% |
| 2 | every non-confirmed KOI out of the non-hosts | 192 / 8,837 | +16.6% | removes only 2% |
| 3 | RUWE < 1.4 | 185 / 8,204 | +16.7% | metal-poor hosts fall from 18 to 8 (B26 18) |
| 4 | parallax/error > 5 or > 10 | 191 / 8,971; 190 / 8,927 | +17.1% | negligible |
| 5 | only one Gaia source within 2″ | 189 / 8,879 | +16.3% | negligible |
| 6 | MS locus fitted to all Kepler targets (Fig. 1 grey) | 192 / 8,992 | +16.7% | no change |
| 7 | 2013 KOI-table non-hosts out | 192 / 8,963 | +16.7% | negligible |
| 8 | logg_GSP-Phot ≥ 4.0 | 190 / 8,860 | +16.9% | negligible |
| 9 | MMA 2014 weight w > 0.40 (round 2) | 192 / 7,927 | +20.9% | region matches (12,409), but a0 falls to 2.72 (B26 2.846 ± 0.007) and the offset rises: ruled out |
| 10 | MMA 2014 Ref = 1 out (round 2) | 192 / 8,992 | +16.7% | no change |
| 11 | MMA 2014 n_Prot-flagged out (round 2) | 192 / 8,948 | +16.7% | negligible |
| 12 | G < 15.5 (round 2) | 170 / 7,588 | +14.8% | passes the gate by coincidence, see below |

Rows 1–8 were listed before any of them was run. Rows 9–12 were added after round 1 showed the long-period,
non-host-only pattern.

**Why the G < 15.5 pass doesn't count.** The cut removes hosts as fast as non-hosts:
- region hosts fall to 348 (B26 400);
- far metal-rich hosts fall to 5 (B26 9);
- the metal-poor count is still +66%.

It meets the two gate numbers while breaking three other checkpoints, and B26 states no magnitude limit. It is not
adopted.

**Also ruled out:**
- No KIC shares a Gaia source.
- The 2013 KOI table and 2014 table1 share no KIC.
- B26's starting sample cannot be Sibony et al. (2022)'s alone. That sample kept 32,049 of the ~34,800 MMA stars
  (Berger 2020 GOF < 0.99 cut; arXiv:2204.01421 §2.1), so it can remove at most ~8%. About 19% of the 1–20 d
  non-hosts have to go.

## Still unexplained, and what could test it

The missing step removes roughly a quarter of non-hosts, evenly in distance, and more at extreme GSP-Phot
metallicity and long period. It leaves the hosts and the period–temperature relation alone. B26 §2 lists no such
step.

Candidates that need columns we have not downloaded, one Gaia DR3 query on the matched source_ids:
- `libname_gspphot` (best-library choice);
- the GSP-Phot uncertainty bounds;
- `phot_variable_flag`;
- `non_single_star`;
- FLAME availability, since B26 says it tried FLAME ages.

The Berger 2020 GOF flag would show Sibony's share. Each would be tested one at a time against the same
checkpoints. Per RULES.md, the Santos step stays held until the gate passes for a stated reason.
