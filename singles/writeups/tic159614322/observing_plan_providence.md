# Observing plan: TIC 159614322 transit, night of 27–28 December 2026, from Providence, RI

Site: 41.82° N, 71.41° W.

Numbers come from `observing_plan.py`, with outputs `observing_plan.json` and `fig_providence.png`.
The ephemeris is from `ephemeris.py`: two full transits, with a planning σ equal to the formal
σ × 2.1.

**Bottom line (honest):**
- **27–28 December 2026:** a single 8–14″ telescope in Providence is unlikely to make a
  convincing detection of this transit.
  - The expected white-noise significance is about 2.6–3.8σ for the whole transit and about
    1.4–2.0σ for the ingress step. Real systematics will make it worse.
  - The star is only 17–20° up for the ingress and the pre-transit baseline.
  - The egress falls in morning twilight.
- **27–28 March 2027** is the better night from New England. See section 6.
- A December attempt is still worth doing as practice, or as part of a multi-site effort with
  western US or European observers.

## 1. The target

| | |
|---|---|
| Star | TIC 159614322 = Gaia DR3 1669046625055549312 |
| Coordinates (J2000) | RA 14h 40m 02.3s, Dec +65° 29′ 03″ (220.00977, +65.48404) |
| Where | In Draco, about 4° east of Thuban (α Dra) and 9° south of Kochab (β UMi). It never sets from Providence. |
| Brightness | V ≈ 11.8 (Gaia G 11.63). It needs a telescope. |
| Signal | 3.3 ppt (0.33%, about 3.6 mmag) dip, lasting about 9.7 h |
| Comparison stars | Six or more Gaia stars with G 11.0–12.3 and similar colour (BP−RP 0.68–0.94) lie 13–20′ away. A field of view of 40′ or more lets you use several of them. |

## 2. When: 27–28 December 2026 (both possible periods predict this transit)

| Event | UTC | EST (evening of 27 Dec) |
|---|---|---|
| Darkness begins (Sun −12°) | 22:42 | 5:42 pm |
| Ingress (start of transit) | 01:51 ± 45 min | 8:51 pm, window 8:06–9:36 pm |
| **Mid-transit** | **06:42 ± 45 min** | **1:42 am**, window 12:57–2:27 am |
| Darkness ends (Sun −12°) | 10:57 | 5:57 am |
| Egress (end of transit) | 11:33 ± 45 min | 6:33 am, window 5:48–7:18 am, in twilight |

The ± 45 min is the planning uncertainty: the formal ± 21 min, inflated because the third TESS
transit arrived about 13–26 minutes late.

## 3. Sky conditions through the night

| EST | Target altitude | Airmass | Moon altitude | State |
|---|---|---|---|---|
| 5:42 pm | 20° | 2.9 | below horizon | before transit |
| 7:42 pm | 17° | 3.4 | below horizon | before transit |
| 8:42 pm | 17° | 3.3 | below horizon | ingress window |
| 10:42 pm | 22° | 2.7 | 17° | in transit |
| 12:42 am | 30° | 2.0 | 38° | in transit |
| 2:42 am | 42° | 1.5 | 53° | in transit |
| 4:42 am | 54° | 1.24 | 54° | in transit |
| 5:42 am | 60° | 1.16 | 47° | in transit; twilight starts |
| 6:42 am | 64° | 1.11 | 39° | egress window; Sun −6° |

**Moon:** 76% illuminated and 72° from the target. It rises around 9 pm EST, so it brightens
the sky for the second half of the transit. Its distance from the target keeps the effect
manageable.

The star never transits the meridian during darkness, so a German equatorial mount needs **no
meridian flip** this night. That's one less systematic.

## 4. Expected precision and SNR

These are white-noise estimates. They include photon noise for V = 11.8, a moonlit suburban
sky, scintillation (Young's formula × 1.5) and four comparison stars. Exposures are 60 s.

| Telescope | σ per 60 s at 60° (X 1.2) | σ per 60 s at 18° (X 3) | σ per 10 min at X 1.2 / X 3 | Ingress-step SNR | Whole-transit SNR (vs 2.5 h of baseline at X ≈ 3) |
|---|---|---|---|---|---|
| 8″ | 3.4 ppt | 12.9 ppt | 1.1 / 4.3 ppt | 1.4 | 2.6 |
| 11″ | 2.6 ppt | 10.4 ppt | 0.9 / 3.5 ppt | 1.7 | 3.2 |
| 14″ | 2.2 ppt | 8.8 ppt | 0.7 / 2.9 ppt | 2.0 | 3.8 |

**Why it's hard:**
- **Scintillation at low altitude.** At 17–20°, atmospheric twinkling dominates, and the only
  pre-transit baseline is taken there.
- **Airmass trend.** Airmass falls from 3.4 to 1.1 across the night. Even with well-matched
  comparison stars, differential extinction can leave trends of 1–3 ppt, and a smooth trend
  looks very much like a 9.7-h, 3.3-ppt box. The sharp ingress step is the only feature a trend
  can't fake, and it happens at the worst altitude.
- **Egress in twilight.** The egress falls in morning twilight (Sun −11° to +1°). Photometry
  of a V = 11.8 star is still possible in nautical twilight, but the sky brightens fast.

## 5. If you try it anyway

1. **Rehearse.** Practise on a known hot Jupiter with a 5–10 ppt transit at a similar
   magnitude. The Exoplanet Watch target list has suitable ones. Doing this before December
   shows what precision your setup actually reaches.
2. **Telescope and camera:**
   - an 8–14″ telescope with a cooled CMOS camera;
   - a field of view of 40′ or more;
   - autoguiding, so the star stays on the same pixels all night;
   - slight defocus, keeping peak counts at about 50% of saturation;
   - 30–60 s exposures, with a clear or luminance filter (or Sloan r′).
3. **Calibration:** take flats, darks and bias frames, and leave focus fixed all night.
4. **When to record:** start as soon as it's dark (5:45 pm EST) and keep going until the sky is
   too bright, around 7 am. That's about 13 hours.
5. **Analysis:** use AstroImageJ, which can fit airmass and other detrending terms together with
   the transit model, or EXOTIC.
6. **Report:** give the time of any ingress you see, and the depth. Send it to the people we
   emailed.

## 6. The better night from New England: 27–28 March 2027 (89.965-d period only)

| Event | UTC | EDT (evening of 27 Mar) |
|---|---|---|
| Darkness begins | 00:21 | 8:21 pm |
| Ingress | 01:00 ± 47 min | 9:00 pm |
| Mid-transit | 05:51 ± 47 min | 1:51 am |
| Darkness ends | 09:21 | 5:21 am |
| Egress | 10:42 ± 47 min | 6:42 am, after sunrise |

The star is at 36–67° through the transit. Expected ingress-step SNR is 3.2 (8″), 4.1 (11″) and
4.8 (14″), white noise only. There is very little dark pre-ingress baseline, and the egress is
lost to daylight.

Mid-transit is near the meridian, so a German equatorial mount will flip mid-event. Treat the
data before and after the flip as separate segments.

**What each outcome means:**
- **A dip** says the period is 89.965 d.
- **Flat data**, if the December transit has been confirmed by someone, says 179.930 d.

A combined effort makes this much stronger: several amateurs, or a 0.4–1 m telescope run by a
follow-up group.
