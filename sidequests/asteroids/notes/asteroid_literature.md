# Asteroid rotation periods from TESS: literature check (2026-10-11)

Compiled from a web search (ADS/arXiv/journal pages) by a sub-agent, with the key items
re-checked by hand on arXiv (marked ✓). LCDB coverage numbers in §5 are our own
(`scripts/lcdb_coverage.py` → `data/lcdb_coverage.csv`). Unverified items say so.

## 1. Bottom line

- **TESS asteroid photometry is a crowded field.** Year 1 (S1–13) has been analysed three times
  (Pál et al. 2020 TSSYS-DR1; McNeill et al. 2023; Vavilov & Carry 2025). The Konkoly group
  (Pál, Kiss, Takács, Szabó, ...) holds a processed "Early Data Release 2" through Sector 48
  (Takács et al. 2025 ✓ for the paper; the EDR2 statistics, ~750k targets and ~1.76M light curves,
  are quoted from the paper body by the sub-agent) and is publishing population papers from it
  (Trojan fast rotators, Hilda fast rotators). The full EDR2 release is promised but not public.
- NASA Ames' open `tess-asteroids` pipeline (Tuson, Martínez-Palomera, Hedges; AAS 245 abstract)
  makes moving-target light curves easy for anyone; its stated goal is the rotation-period
  distribution.
- Survey inversion keeps absorbing TESS: Ďurech & Hanuš 2026 (arXiv:2610.06082 ✓, posted
  5 Oct 2026) combine nine surveys including TESS DR1 into ~30,000 models, including
  **247 Jupiter Trojans and 157 Hildas**. Gaia DR4 (2 Dec 2026) will add 66 months of SSO data.
  Rubin is forecast to deliver ~15,500 Hilda periods (Fleming et al. 2025).
- **The real TESS edge is slow rotators (P ≈ 24–300 h)**: no diurnal aliasing, 27-d continuous
  baseline. For Trojans and Hildas, the gap is at **H ≳ 12 (V ≳ 18 at opposition)**, which is at
  or beyond TESS's useful FFI limit. The bright end is already done (LCDB) or being done (Konkoly).
- Our only uncontested data: sectors after S48 (Konkoly EDR2 ends there), in particular the
  200-s sectors S56+ and the ecliptic-heavy sectors S70–72, S91–92.

## 2. TESS papers

| Paper | Data | Result |
|---|---|---|
| Pál, Molnár & Kiss 2018, PASP 130, 114503 (arXiv:1809.07403) | prospects | camera 1 at ecliptic latitude 18° in Y1–2 |
| **Pál et al. 2020, ApJS 247, 26** (arXiv:2001.05822), TSSYS-DR1 | S1–13, camera 1, 30-min FFIs | 9,912 light curves; scan to 24 c/d; claims long periods underestimated by ground surveys; data at archive.konkoly.hu/pub/tssys/dr1 |
| McNeill et al. 2020, ApJS (arXiv:1911.01495) | pilot | 300 periods; 43 with P > 100 h; predicted ~10x more P > 100 h objects |
| **McNeill et al. 2023, AJ 166, 152** (VizieR J/AJ/166/152) | Cycle 1, own photometry | 37,965 objects; 3,492 "highly reliable" periods |
| Szabó et al. 2022, A&A (arXiv:2202.04372) | DR1 | 2,859 members of 16 families |
| **Vavilov & Carry 2025, A&A 693, A66** (arXiv:2501.07189) | DR1 re-analysis | 4,521 periods; extrema-count alias test; deficit at 50–100 h |
| **Takács et al. 2025, PASP** (arXiv:2503.13332 ✓) "EDR2 I" | S1–48 processed | spin/shape for 44 MBAs; full EDR2 release "in a separate paper" |
| Kiss et al. 2025, A&A 694, L17 (arXiv:2501.18440) | EDR2 | 3 fast-rotating Trojans (2.9–4.8 h) |
| Takács et al. 2025, ApJL 986, L33 (arXiv:2506.00144) | S42–46 | 3 fast Hildas; their TESS survey contains ~400 Hildas |
| DPS 2023 abstract (BAAS 55(8), 2023n8i113p05) | TESS Trojans | >500 Trojan light curves, amplitude limits for ~3,000 (page 403'd; authors unverified) |
| Gowanlock et al. 2024, AJ 168, 181 (arXiv:2408.04096) | ZTF+TESS | 222 objects; ~85% recovered; 48-h ZTF alias |
| Tuson et al., `tess-asteroids` (2025AAS...24510315T) | S1–3 demo | open moving-target pipeline; poster: 402 bright asteroids |

Not found: any "TSSYS-DR2" catalogue paper, any Pál 2023/2025 first-author catalogue, any public
TESS Trojan or Hilda period catalogue.

## 3. Gaia

- Ďurech & Hanuš 2023, A&A 675, A24 (arXiv:2305.10798): DR3 SSO photometry, ~8,600 spin states;
  the spin–size gap.
- Cellino et al. 2024, A&A 687, A277: 22,815 triaxial-ellipsoid inversions.
- MacLennan et al. 2026, A&A 707, A131: >8,000 from DR3 alone (unverified by hand).
- Zhou et al. 2025 (arXiv:2502.19631): most slow rotators are probably tumbling, so a slow-rotator
  paper now needs a non-principal-axis check.
- Gaia DR4: 2 Dec 2026, with `sso_source` for 436,551 objects.

## 4. Ground surveys, K2 and others

- ATLAS: Ďurech et al. 2020, A&A 643, A59 (~2,750 models); Ďurech et al. 2022, Front. (>5,000 periods).
- ATLAS+ZTF super-slow rotators: Erasmus et al. 2021, MNRAS 506, 3872 (39 with P > 1000 h).
- **ZTF Trojans: Li et al. 2025, Icarus (arXiv:2504.15530 ✓)**: 216 Trojan periods, 80 new,
  4.6–447.8 h. They state that ~240 reliable Trojan periods were known in 2024.
- Pan-STARRS1 (ApJS 241, 6, 2019): 876 periods. iPTF: Chang et al. 2015/2016 (~2,500).
- K2: Szabó et al. 2017 (56 L4 Trojans; ~20% with P > 50 h); Kalup et al. 2021, ApJS 254, 7
  (101 Trojans; 25% with P ≥ 30 h); Szabó et al. 2020 (125 Hildas; 18% with P > 100 h);
  Molnár et al. 2018 (608 MBAs).
- FOSSIL/Subaru: 53 small Trojans, 17 Hildas.
- Rubin First Look (Greenstreet et al. 2026): 76 periods, short baselines.
  Euclid 2026: 889 periods. KMTNet 2026: 96 periods.

## 5. LCDB, and how underserved Trojans and Hildas are

- LCDB (Warner, Harris & Pravec 2009, Icarus 202, 134): **the current public file is still the
  2023 Oct 1 release** (`LCLIST_PUB_CURRENT.zip` downloaded 2026-10-11 contains 2023-Oct files).
  It therefore lacks Li et al. 2025, Vavilov & Carry 2025, Ďurech & Hanuš 2026 and most TESS work.
- U codes: 1 = possibly wrong; 2 = ambiguous or partial coverage (may be off by ~30%); 3 = secure;
  +/− refine these; 0 = proven wrong (Butkiewicz-Bąk et al. 2017, MNRAS 470, 1314).
- Our cross-match of numbered asteroids (JPL SBDB, 2026-10-11) against LCDB 2023-Oct
  (`data/lcdb_coverage.csv`):

| population (numbered) | N | LCDB U ≥ 2 | fraction | U ≥ 3− |
|---|---|---|---|---|
| Jupiter Trojans (SBDB class TJN) | 10,633 | 376 | 3.5% | 135 |
| Hildas (a 3.7–4.2, e < 0.3, i < 20°) | 3,821 | 184 | 4.8% | 76 |
| outer MB (a 2.82–3.70) | 292,903 | 7,539 | 2.6% | 1,061 |
| middle MB | 316,118 | 9,199 | 2.9% | 1,259 |
| inner MB | 240,584 | 8,864 | 3.7% | 1,854 |

  By H, Trojans: H < 10, 99% have U ≥ 2; H 10–11, 71%; H 11–12, 31%; **H 12–13, 7%
  (75/1,031)**; H 13–14, 1%. Hildas: H < 12, ≥ 90%; **H 12–13, 40%**; H 13–14, 11%.
  Adding ZTF (80 new Trojans), K2 (~100 Trojans, ~125 Hildas), Ďurech & Hanuš 2026 (247 Trojans,
  157 Hildas) and Konkoly's unpublished TESS work, the bright end (H < 11 for Trojans, H < 12 for
  Hildas) is essentially complete. The open territory is Trojans at H ≈ 11.5–13.5, i.e.
  **V ≈ 18–20 at TESS geometry** (Δ ≈ 4.3–5 au).
- Slow rotators in LCDB (U ≥ 2): P > 24 h is 22–32% of entries in every population and P > 100 h
  is 6–9% (`data/lcdb_slow_rotators.csv`). These are known to be biased low by ground-survey
  windows (Pál 2020; McNeill 2023).
- ALCDEF (alcdef.org): 24,816 objects; simple text upload (S-ALCDEF).
- Minor Planet Bulletin: still publishing (Vol. 53, 2026); quarterly; ADS-indexed; no page charges.
  The natural outlet for a handful of new amateur-style periods, with an ALCDEF deposit.

## 6. Which population for the pilot, and why

Jupiter Trojans crossing a post-S48, 200-s sector, which is outside Konkoly EDR2:

1. Low LCDB coverage (3.5%) and a high slow-rotator / binary fraction, where 27 d of continuous
   coverage beats ground surveys.
2. S71/S72 have the most Trojans on silicon in S56–S107: 916 and 1,132 with H ≤ 13.5 and
   ≥ 10 d (`data/sector_counts_trojan_H13.5.csv`).
3. The open question is whether TESS reaches the V ≈ 18–19 objects that lack periods. The pilot
   is designed to measure exactly that.

Hildas give far fewer objects per sector (≤ 78 with H ≤ 13.5 in the best sector, S72;
`data/sector_counts_hilda_H13.5.csv`).

## References

Pál+2018 arXiv:1809.07403 · Pál+2020 2020ApJS..247...26P arXiv:2001.05822 · McNeill+2020
arXiv:1911.01495 · McNeill+2023 2023AJ....166..152M · Szabó+2022 arXiv:2202.04372 · Vavilov & Carry
2025 arXiv:2501.07189 · Takács+2025 arXiv:2503.13332 · Kiss+2025 arXiv:2501.18440 · Takács+2025
arXiv:2506.00144 · Gowanlock+2024 arXiv:2408.04096 · Tuson+ 2025AAS...24510315T · Ďurech & Hanuš 2023
arXiv:2305.10798 · Cellino+2024 A&A 687, A277 · Zhou+2025 arXiv:2502.19631 · Ďurech & Hanuš 2026
arXiv:2610.06082 · Ďurech+2020 arXiv:2010.01820 · Erasmus+2021 arXiv:2106.16066 · Li+2025
arXiv:2504.15530 · Szabó+2017 arXiv:1609.02760 · Szabó+2020 arXiv:2001.06656 · Kalup+2021
arXiv:2102.09447 · Molnár+2018 arXiv:1706.06056 · Nesvorný+2020 arXiv:2003.12634 · Fleming+2025
arXiv:2512.00624 · Warner, Harris & Pravec 2009 Icarus 202, 134 · Butkiewicz-Bąk+2017 arXiv:1707.04427 ·
LCDB https://minplanobs.org/mpinfo/php/lcdb.php · ALCDEF https://alcdef.org · MPB https://mpbulletin.org
