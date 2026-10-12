# Literature check: the three NCVZ timing signals that B26 does not list

Checked 2026-10-12 for TIC 353854078, 441804888 and 233193944. These three were flagged or noted in the northern
CVZ false-alarm set (`data/ncvz_fa_vetting.csv`).

## What was searched

- **arXiv LaTeX sources** of Mitnyan et al. 2024 (arXiv:2402.01486) and Borkovits et al. 2026 (arXiv:2604.04004).
  - The TIC numbers were searched both as plain digits and split by separators (`353\,854\,078`, superscripts, etc.).
  - **None of the three appears in either paper.**
- **SIMBAD:** names and full reference lists.
- **Tokovinin's Multiple Star Catalog** (MSC; VizieR J/ApJS/235/6).
- **The original period-study papers**, where these exist.

## Caveat on "not in B26"

B26 did **not** re-time all ~3,500 NCVZ EBs on the S14–86 baseline. They returned only to the **351 EBs** that Mitnyan+2024
had flagged for non-linear ETVs in S14–60 (B26 §2: "we returned to the previously selected narrower sample of the
above-mentioned 351 targets").

That list of 351 is not published. Mitnyan+2024 tabulates only the systems it solved. So "not in B26" means one of two things:
- the EB was never selected from S14–60, which is plausible for P3 ≳ 1,300 d;
- it was selected and then rejected.

We cannot tell which from the papers. Mitnyan's input list was the GSFC machine-learning EB list (Powell+2021), restricted
to ≥ 8 Year-2 sectors. All three stars sit at ecliptic latitude 82–87°, so they easily meet that condition.

## Per system

| TIC | Identity | Prior work | Verdict |
|---|---|---|---|
| 353854078 | **SX Dra** (HD 238811), Algol, P = 5.17 d, semi-detached, δ Sct primary | Century of ground-based minima. **Soydugan & Kaçar 2013, AJ 145, 87:** parabola (mass transfer) plus a cyclic term, P ≈ 82–88 yr, A ≈ 0.137 d. Their LTTE reading needs m3,min = 3.4 M☉, which they reject in favour of magnetic activity. Also Soydugan et al. 2002, Ap&SS 282, 399. | **Known period variable.** The ~1,330 d, ~60 s TESS "orbit" is most likely the local curvature of the known ~85-yr variation and/or activity, and is not a new third body. B26 would reasonably exclude it. |
| 233193944 | **HD 160861** (HIP 86187), F5, P = 2.506 d | **MSC 17364+6820** is a quadruple: Aa1,Aa2 = the 2.506 d SB1 (Carquillat 1976), which is our EB; Aa,Ab = a **visual orbit, P = 24.5 yr, e = 0.9, a = 0.10″** (grade 9, poor); C at 5.5″; B is optical (GJ 687). Gaia DR3 RUWE 5.9. Also Tokovinin et al. 2006 (tertiary companions) and Kervella+2019 / Brandt 2021 (proper-motion anomaly). | **Known tertiary.** At plx ≈ 10 mas, a ≈ 10 AU. With masses ~1.9 + 1.0 M☉ (MSC), the EB's LTTE semi-amplitude is of order 10–20 min over the 24.5-yr orbit. The ±5 min TESS curvature fits a section of that orbit. **TESS ETVs could improve the grade-9 visual orbit**, but this is not a discovery. |
| 441804888 | **TYC 4436-905-1**, P = 2.738 d, Tmag 11.3 | 7 SIMBAD refs, all catalogues: Prša+2022 (TESS EBs), Gaia DR3 variability and EB catalogue, APOGEE/ATL. **No timing or multiplicity work.** | **Unreported.** The signal is weak: a smooth common O−C whose late upturn depends on the S73–86 season, and it fails the v2 jackknife. **Candidate only.** |

## Net effect on the README claim

"3/49 NCVZ EBs show real long-period ETVs not in B26" stands as a statement about the ETVs. But only **1 of the 3 is
unreported in the literature**, and it is the weakest of the three.
- One is a famous period-variable Algol.
- One has a known 24.5-yr visual tertiary.

So the false-alarm test does not show that B26 is incomplete for *new* triples. It also tightens the procedure for the
south: **every candidate needs a SIMBAD / MSC / WDS / ground-based period-study check, not just a check against the ETV papers.**

## The same check on the four southern candidates (done 2026-10-12)

| TIC | Identity | MSC (J/ApJS/235/6) | WDS (B/wds) | SIMBAD refs | Verdict |
|---|---|---|---|---|---|
| 32000625 | TYC 9155-22-1 | none within 20″ | none | 8: RAVE, APOGEE, Prša+2022, Justesen & Albrecht 2021 EB catalogue, Shi+2022 pulsators in EA EBs | No prior timing or multiplicity work. |
| 279323859 | UCAC4 155-008982 | none | none | 5: TESS-HERMES, GALAH, Prša+2022, Gaia DR3 variability and EB catalogue | No prior work. |
| 278861826 | TYC 8544-1497-1 | none | none | 5: catalogues only | No prior work. |
| 294327312 | HD 56910 (HIP 35095), Am | none | none | 27, including Makarov & Kaplan 2005, Frankowski+2007, Kervella+2019 and Brandt 2021 (HGCA): Hipparcos/Gaia **astrometric acceleration** | No published outer orbit. The astrometric acceleration supports the ETV companion; it is not a prior detection of its orbit. |

All four still look unreported. Remaining gap: no ADS full-text search, and no check of the Kreiner O−C atlas
(not relevant for these faint, recently discovered EBs, but it matters for bright, named ones like SX Dra).
