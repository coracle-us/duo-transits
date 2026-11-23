# DRAFT, NOT SENT: email to the Salinas et al. team

Status: draft for Jon to review. Nothing has been sent.

- **To:** H. Salinas, corresponding author of arXiv 2610.24882 (address on the paper).
- **Cc:** R. Brahm.
- **Before sending:**
  - Done 2026-11-23: no TESS data after S102 exists for any of the four stars, and none is on
    silicon in S103–S110 (`results/round3/independent/RESULTS.md`, "Predictions").
  - If S111 data are public by the time you send this, check the 275647095 phase-0.872 dip first.
  - Re-run `round3/s13_novelty.py` if more than a few days have passed.

Every number traces to a file:
- periods, depths and durations: `results/round3/oddeven_summary.csv`;
- next transits: `results/round3/predictions_2026-11-22.csv`;
- pixel and Gaia checks: `results/round3/PIXELS_2026-11-22.md` and `results/round3/independent/RESULTS.md`;
- novelty: `results/round3/NOVELTY_2026-11-22.md`.

---

**Subject:** Periods from later TESS sectors for four of your single-transit candidates (arXiv 2610.24882)

Dear Dr. Salinas,

I read your paper "Single-Transit Parameter Estimation with Calibrated Uncertainties" with interest. I
searched TESS Sectors 61–107 for second transits of your single-transit candidates, and for four of
them the transits now repeat on one period. As far as I can tell, these periods have not been reported:
- none of the four stars is a TOI or CTOI;
- ExoFOP has no follow-up files for them;
- I found nothing in the literature.

Since your group is following up its own candidates, I wanted you to have them first.

| TIC | T (mag) | P (d) | T0 (BTJD) | Depth (ppm) | T14 (h) | Transits used | Next transit (UTC) |
|---|---|---|---|---|---|---|---|
| 275647095 | 12.24 | 17.06402 ± 0.00003 | 2315.5125 ± 0.0018 | 9130 | 8.4 | 7 | 2026-12-05 18:11 ± 6 min |
| 241550976 | 12.17 | 21.19694 ± 0.00004 | 2339.4442 ± 0.0019 | 16350 | 5.1 | 5 | 2026-11-29 20:23 ± 6 min |
| 107113048 | 12.44 | 49.70152 ± 0.00017 | 2312.9429 ± 0.0033 | 8240 | 7.1 | 4 | 2026-12-25 21:44 ± 11 min |
| 230133616 | 12.20 | 23.75297 ± 0.00002 | 2396.1280 ± 0.0007 | 5880 | 9.3 | 31 | 2026-12-16 21:02 ± 3 min |

(T0 is your catalogued event, refitted. The depths are trapezoid fits to QLP and TESS-SPOC data.)

**What I checked for each star:**
- odd/even depths: all agree within 1.3σ (2.2σ in the second analysis for 107113048, which has only 2 + 2 transits);
- secondary eclipse at phase 0.5: none (but see TIC 275647095 below);
- transit shape and the stellar density it implies;
- TESScut difference images against Gaia DR3 neighbours;
- Gaia binarity indicators;
- the TESS EB catalogue and VSX: no matches within 30″.

A second analysis, with separate code written without access to the first, recovered all four periods
within 1.6σ.

**Notes per star:**

- **TIC 241550976** is the cleanest.
  - The dip is on target.
  - The Gaia radius (0.96 R☉) fits a circular orbit.
  - The companion is about 1.2 R_J.
- **TIC 107113048:** on target, and a circular orbit fits (about 0.9 R_J).
  - Only 4 transits so far. One of them (Sector 90) overlaps an eclipse of an uncatalogued eclipsing
    binary 97″ away (TIC 107116357).
  - That transit is on time to about 1 minute, but it is the only one that rules out 74.55 and 149.10 d.
  - TESS will not observe the 25 December transit (S112 does not cover the star), so a ground-based
    observation then would settle it.
  - Otherwise the first TESS test is Sector 114 (Feb–Mar 2027).
- **TIC 275647095:** on target, but the transit is about 1.4 times longer than a circular orbit allows
  for the Gaia radius (1.14 R☉). This suggests e ≳ 0.34, or a blend within ~2″.
  - There is a marginal dip at phase 0.872 (about 690 ± 150 ppm, 3–5σ depending on the data). The
    transiting object can't make it.
  - If it is real, the signal is more likely the primary eclipse of a blended eccentric eclipsing binary.
  - TESS Sectors 111–114 cover six of these dips and should decide it.
- **TIC 230133616:** this is not really a single.
  - It sits in the northern continuous viewing zone and transits in Sectors 16–86 (37 transits).
  - The host looks like a 2.45 R☉ subgiant with Gaia RUWE 3.8, so it probably has an unresolved
    companion.
  - The companion is 1.7–1.9 R_J, and the transit fits a circular orbit poorly.
  - There is no significant secondary eclipse at any phase we tested, which argues against an M-dwarf companion.
  - It is in the far north, so it may not suit southern follow-up.

I'd be glad to send the light curves, timing tables and code. They are in a public repository:
https://github.com/coracle-us/duo-transits, under `salinas-second-transits/`.

I may write the search up as a short Research Note of the AAS. If you are already working on any of
these stars, please let me know. I'd rather cite your work, or leave a star out, than step on it.

I am an independent researcher, and I did this work with the help of an AI assistant (Anthropic's
Claude), which wrote and ran most of the analysis code under my direction. I have checked the results
myself.

Best regards,
Jon McBee
