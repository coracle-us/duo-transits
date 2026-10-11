"""
Plan step 2 (triage): assign a status to every candidate and write data/catalog/status.csv.

Status values
  resolved        a single period has been claimed with evidence (see `tier` for how strong)
  false-positive  dispositioned as not a planet (TFOPWG, published FP catalog, or an ExoFOP
                  note documenting an eclipsing/spectroscopic binary)
  unclear         conflicting claims, a period claim with no stated evidence, or binary hints
                  that we could not confirm
  open            no period resolution or disposition found (may still have partial follow-up)

tier (for resolved): peer-reviewed | catalog-paper | thesis | exofop-community | toi-catalog
flags: concerns that do not change the status (e.g. likely-FP statistics, Gaia orbits)

Assignments are curated by hand from: data/catalog/crossmatch.csv (archive/TOI/CTOI),
data/catalog/exofop_planet_params.csv + data/raw/exofop_targets/ (ExoFOP uploads/notes),
data/catalog/simbad_refs.csv, and the papers/theses cited in `source`. Every period listed
was checked to be an integer alias of dT (see crossmatch.csv / exofop_planet_params.csv).
Literature checked through 2026-10-11.

Usage:  python scripts/build_status.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "catalog"

EXOFOP = "https://exofop.ipac.caltech.edu/tess/target.php?id={}"
OBSNOTES = "https://exofop.ipac.caltech.edu/tess/edit_obsnotes.php?id={}"
THESIS = ("Hawthorn 2024, PhD thesis (Univ. Warwick), Ch. 6 "
          "https://warwick.ac.uk/fac/sci/physics/research/astro/theses/faithhawthorn.phd.pdf")
CATALOG_PAPER = "Hawthorn et al. 2024, MNRAS 528, 1841 (arXiv:2310.17268)"

# TIC: (status, tier, period_d, n_alias, evidence, source, flags)
S = {
    # ---- resolved: peer-reviewed / refereed or validation papers
    124029677: ("resolved", "peer-reviewed", 20.33003, 36, "TOI-5153 b, confirmed (RVs)",
                "Ulmer-Moll et al. 2022, A&A 666, A46 (2022A&A...666A..46U)", ""),
    157698565: ("resolved", "peer-reviewed", 61.6277, 12, "TOI-2589 b, confirmed",
                "Brahm et al. 2023, AJ 165, 227 (2023AJ....165..227B)", ""),
    269333648: ("resolved", "peer-reviewed", 64.5949, 11, "TOI-2529 b, confirmed",
                "Jones et al. 2024, A&A 683, A192 (2024A&A...683A.192J)", ""),
    303317324: ("resolved", "peer-reviewed", 29.566602, 25,
                "TOI-4310 b, statistically validated; period from extra transits (S69, QLP S96)",
                "Collier et al. 2026, arXiv:2608.21560; ExoFOP TOI comment", ""),
    333736132: ("resolved", "peer-reviewed", 25.241192, 29, "NGTS-35 b = TOI-6669 b, confirmed",
                "Kendall et al. 2026, MNRAS 547 (2026MNRAS.547f2189K, arXiv:2512.07716)", ""),
    404518509: ("resolved", "peer-reviewed", 25.1292, 28, "HD 21520 b = TOI-4320, confirmed",
                "Nies et al. 2024, MNRAS 534, 3744 (2024MNRAS.534.3744N)", ""),
    466206508: ("resolved", "peer-reviewed", 75.12375, 5, "TOI-5542 b, confirmed",
                "Grieves et al. 2022, A&A 668, A29 (2022A&A...668A..29G)", ""),
    77437543: ("resolved", "peer-reviewed", 60.33, 12, "TOI-2490 b, eccentric brown dwarf "
               "(not in pscomppars)", "Henderson et al. 2024, MNRAS 533, 2823 (2024MNRAS.533.2823H)",
               "brown dwarf, not a planet"),
    118339710: ("resolved", "peer-reviewed", 80.4, 9, "TOI-6695 b; both duotransit events are "
                "planet b (c is RV-only). Archive 80.389 d vs static dT/9=80.433 d (TTVs, near-resonant pair)",
                "Eberhardt et al. 2025, AJ 169, 298 (2025AJ....169..298E)", "TTVs"),
    349091983: ("resolved", "catalog-paper", 549.065, 1, "only one alias allowed by Cycle 1+3 "
                "coverage (CVZ)", CATALOG_PAPER + " Sec. 4.2 / Table A3 (N_alias=1)", ""),
    # ---- resolved: PhD thesis (NGTS photometry), not peer reviewed
    224279805: ("resolved", "thesis", 22.9840, 32, "two NGTS ingresses confirmed period; "
                "TRICERATOPS 'likely planet' (FPP 0.049); ExoFOP upload (mason 2023) agrees: 22.984184 d",
                THESIS + " Sec. 6.4.4; " + EXOFOP.format(224279805), ""),
    142278054: ("resolved", "thesis", 26.4937, 27, "NGTS ingress confirmed period; 'likely planet' "
                "(FPP 0.42); extra transit in S63", THESIS + " Sec. 6.4.3", ""),
    439491923: ("resolved", "thesis", 27.3608, 26, "three NGTS partial transits; 'likely planet' "
                "(FPP 0.24)", THESIS + " Sec. 6.4.7", ""),
    39167176: ("resolved", "thesis", 21.7883, 34, "one NGTS ingress; TRICERATOPS FPP 0.96, NGTS depth "
               "shallower than TESS (chromatic?)", THESIS + " Sec. 6.4.2", "likely FP (thesis)"),
    339399841: ("resolved", "thesis", 27.2487, 14, "four NGTS partial transits; TRICERATOPS FPP 0.70, "
                "most likely NEB on TIC 339399847", THESIS + " Sec. 6.4.6", "likely FP (thesis)"),
    # ---- resolved: TFOP / ExoFOP community analyses (unrefereed)
    323295479: ("resolved", "exofop-community", 185.62807, 4, "TOI-1861: EXOFASTv2 alias analysis; "
                "ground detection 2459736 + non-detections exclude all n except 4",
                OBSNOTES.format(323295479) + " (lhill, 2026-07-14); Lofvenholm 2026 MSc thesis, Chalmers "
                "https://odr.chalmers.se/items/50878040-03ee-4d18-a2f4-64202738fd83",
                "TOI table still lists 742.5 d"),
    81089255: ("resolved", "exofop-community", 46.8905, 16, "TOI-6480: transits in S7,35,61,62",
               EXOFOP.format(81089255) + " (upload: mason, 2023-06-01)", ""),
    188620407: ("resolved", "exofop-community", 43.448716, 25, "TOI-7869: transits in S2,42,69 and "
                "absence elsewhere; TOI table adopts this period", EXOFOP.format(188620407) +
                " (upload: mason, 2023-11-17); ExoFOP TOI comment", ""),
    242241304: ("resolved", "exofop-community", 26.7517, 27, "transits in S11,38,65,102",
                EXOFOP.format(242241304) + " (upload: slee, 2026-06-09)", ""),
    265465927: ("resolved", "exofop-community", 31.1287617, 23, "transits in S1,27,67,94,101,102 "
                "(last two consecutive); LCO ground detections 2025", EXOFOP.format(265465927) +
                " (upload: mason, 2026-05-29; TFOP notes files)", ""),
    289840544: ("resolved", "exofop-community", 25.9308, 29, "transits in S1,28,68,92,95",
                EXOFOP.format(289840544) + " (upload: slee, 2026-04-20)", ""),
    304339227: ("resolved", "exofop-community", 51.0497543, 14, "TOI-4629: transits in "
                "S13,39,66,101,103", EXOFOP.format(304339227) + " (upload: mason, 2026-06-14)",
                "TOI table still lists 714.7 d"),
    457649900: ("resolved", "exofop-community", 26.74387, 28, "TOI-4958: transits in S11,38,65,101",
                EXOFOP.format(457649900) + " (upload: vincent, 2026-08-15)", "TOI table still lists 748.8 d"),
    # ---- resolved: TOI catalog period from QLP faint-star multi-sector search
    193318850: ("resolved", "toi-catalog", 15.717856, 47, "TOI-7762.01 found in QLP faint-star search "
                "at an alias period", "ExoFOP TOI table (data/raw/exofop_toi.csv)", "unverified"),
    306249066: ("resolved", "toi-catalog", 26.162212, 28, "TOI-6596.01 found in QLP faint-star search; "
                "agent-reported on-time ground transits (ASTEP 2024-03-29, 2025-02-06)",
                "ExoFOP TOI table; " + EXOFOP.format(306249066) + " (TFOP files; not opened by us)", "unverified"),
    412386707: ("resolved", "toi-catalog", 39.295951, 19, "TOI-7619.01 found in QLP faint-star search",
                "ExoFOP TOI table (data/raw/exofop_toi.csv)", "unverified"),
    # ---- false positives
    70561926: ("false-positive", "", 21.49257, 35, "KELT FP catalog: SB1, RV amp 6.6 km/s, P=21.4926 d "
               "= dT/35, depth 14.5 mmag (duo depth 16 ppt); Gaia DR3 SB1 21.5 d",
               "Collins et al. 2018, AJ 156, 234 (VizieR J/AJ/156/234, KS22C028401)", ""),
    157119927: ("false-positive", "", None, None, "Gaia DR3 SB1 P=21.865+/-0.003 d (K~30 km/s; near but "
                "~3 sigma from dT/32=21.856 d); ExoFOP note: "
                "'The TESS detection shows this system is an eclipsing binary'",
                OBSNOTES.format(157119927) + " (tarrants, 2024-03-10); Gaia DR3 NSS", ""),
    161169240: ("false-positive", "", None, None, "TOI-6667.01 retired by TFOPWG as APC; spectroscopic "
                "disposition SB1 (PFS: 8 km/s shift over 25 d)", OBSNOTES.format(161169240) +
                " (mlund, 2024-08-16)", ""),
    # ---- unclear
    13072758: ("unclear", "", None, None, "CONFLICT: thesis gives 23.7940 d (n=31; one NGTS ingress + "
               "WINE full transit) but ExoFOP upload gives 49.1746 d (n=15; QLP transits S5,32,98, "
               "'no transit later in S98 rules out P/2'). Not mutually consistent.",
               THESIS + " Sec. 6.4.1; " + EXOFOP.format(13072758) + " (upload: mason, 2026-05-28)",
               "good pipeline test case"),
    25194908: ("unclear", "", 337.45806, 2, "CTOI .02 upload claims the shorter of the 2 allowed aliases; "
               "no supporting evidence stated", EXOFOP.format(25194908) + " (upload: sedgley, 2025-12-15)", "CVZ"),
    294097549: ("unclear", "", 199.53498, 3, "CTOI .02 upload claims the shortest of the 3 allowed aliases; "
                "no supporting evidence stated", EXOFOP.format(294097549) + " (upload: sedgley, 2025-12-15)", ""),
    71028120: ("unclear", "", 19.58205, 38, "CTOI .02 upload at an alias period; note only says "
               "'Impact parameter and R2 high in this fit'", EXOFOP.format(71028120) +
               " (upload: sedgley, 2025-12-15)", ""),
    2669681: ("unclear", "", 45.975, 8, "K2 notes: 'possible eccentric binary', primary depth ~4.66%, "
              "P=45.975 d = dT/8 (TESS duo depth ~0.8%)", OBSNOTES.format(2669681) +
              " (ExoFOP-K2 notes for EPIC 213689594, 2016-2017)", "possible EB"),
    32179255: ("unclear", "", None, None, "Gaia DR3 SB1 orbit P=72.13+/-0.11 d, K~27 km/s, close to dT/10=72.26 d; "
               "Hawthorn+24 lists Gaia NSS=2, RV amp 49 km/s", "Gaia DR3 nss_two_body_orbit (our inference)",
               "possible binary"),
    100776118: ("unclear", "", None, None, "Gaia DR3 AstroSpectroSB1 P=251.3+/-0.8 d vs dT/3=249.0 d (~3 sigma); "
                "RUWE 2.37; Hawthorn+24 Gaia NSS=3", "Gaia DR3 nss_two_body_orbit (our inference)", "possible binary"),
}

# Open candidates that have relevant notes but no resolution
OPEN_NOTES = {
    140215502: ("TFOP SG1 ground campaign Dec 2024 on the 146.777 d (n=5) alias; mixed/unsettled detections",
                EXOFOP.format(140215502) + " (TFOP files)"),
    218977148: ("TRAPPIST-South 2025-10-25 tested P=31.94 d (n=23): possible early egress; later LCO runs, "
                "no results posted", EXOFOP.format(218977148) + " (TFOP files)"),
    394340183: ("third transit at Tc1+2dT (BTJD 3113.67) reported; consistent with every alias, does not narrow",
                EXOFOP.format(394340183) + " (upload: essack, 2025-05-28)"),
    296737508: ("TOI-6267: 3 transits ~734 d apart (S8,35,62); does not narrow", "ExoFOP TOI comment"),
    313671132: ("TOI-6559: 3 transits in S11,38,65 at 733 d spacing; does not narrow. exoplanet.eu "
                "lists it as confirmed at 733 d citing arXiv:2512.00967, which only reuses the TOI value",
                "ExoFOP TOI comment"),
    396720998: ("TOI-709.02 (duotransit) still PC. Separate 32 d signal TOI-709.01 retired as APC/BEB "
                "(chromatic). Hot subdwarf host with K companion", OBSNOTES.format(396720998)),
    156716001: ("TOI-7602.01 (QLP, 10.30 d, 0.8 ppt, 3.8 h) is NOT an alias of dT (dT/P=71.8) and differs in "
                "depth/duration: a different signal; duotransit itself unresolved", "ExoFOP TOI table"),
    207783865: ("Gaia DR3 SB1 7.545 d (K~4.8 km/s), not an alias of dT", "Gaia DR3 NSS"),
    275878706: ("Hawthorn+24 lists Gaia NSS=1", CATALOG_PAPER),
}


def main():
    cand = pd.read_csv(CAT / "candidates.csv")
    rows = []
    for _, c in cand.iterrows():
        t = int(c.TIC)
        if t in S:
            st, tier, P, n, ev, src, fl = S[t]
        else:
            st, tier, P, n, fl = "open", "", None, None, ""
            ev, src = OPEN_NOTES.get(t, ("no period resolution or disposition found",
                                         "pscomppars, ExoFOP TOI/CTOI + target page, SIMBAD, "
                                         "ADS/arXiv searches (2026-10-11)"))
        if P is not None and n is not None:
            # sanity: claimed periods must be integer aliases of dT (TTVs allowed for TOI-6695)
            assert abs(c.dT / P - n) < 0.01, (t, c.dT / P, n)
        rows.append(dict(TIC=t, status=st, evidence=ev, source=src, tier=tier,
                         period_claim_d=P, n_alias=n, flags=fl, N_alias_paper=int(c.N_alias),
                         Pmin_paper=round(c.Pmin, 4), Pmax_paper=round(c.Pmax, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(CAT / "status.csv", index=False)
    print(out.status.value_counts().to_string())
    print()
    print(out[out.status == "resolved"].tier.value_counts().to_string())


if __name__ == "__main__":
    main()
