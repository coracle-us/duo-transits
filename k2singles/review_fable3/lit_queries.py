"""Fable review 3: SIMBAD bibliography, Gaia DR3 variability tables and VizieR period catalogues
for the eight table rows with 'none' in the Prior-period column (plus HD 220967 by name)."""
import warnings, sys
warnings.filterwarnings("ignore")
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

d = pd.read_csv("data/lead_status.csv")
stars = {246194159: "HD 220967", 248811085: None, 246089124: None, 246696483: None, 211390677: None,
         247967714: None, 212012030: None, 211503363: None}
for epic, name in stars.items():
    r = d[d.EPIC == epic].iloc[0]
    sid = int(r.gaia_source_id)
    print(f"\n===== EPIC {epic} TIC {r.TIC} Gaia DR3 {sid} G={r.phot_g_mean_mag:.2f} RUWE={r.ruwe:.2f}")
    # SIMBAD: main id, sptype, bibcodes
    try:
        s = Simbad()
        s.add_votable_fields("sp_type", "otype", "ids")
        t = s.query_object(f"EPIC {epic}")
        if t is None or len(t) == 0:
            t = s.query_object(f"Gaia DR3 {sid}")
        if t is not None and len(t):
            row = t[0]
            cols = {c.lower(): c for c in t.colnames}
            print("  SIMBAD:", row[cols.get("main_id")], "| sptype:", row[cols.get("sp_type")] if cols.get("sp_type") else "", "| otype:", row[cols.get("otype")] if cols.get("otype") else "")
            bib = s.query_bibobj(row[cols.get("main_id")])
            if bib is not None and len(bib):
                bc = [c for c in bib.colnames if "bibcode" in c.lower()]
                print("  SIMBAD refs:", len(bib), list(bib[bc[0]])[:30] if bc else bib.colnames)
            else:
                print("  SIMBAD refs: 0")
        else:
            print("  SIMBAD: not found")
    except Exception as ex:
        print("  SIMBAD error", ex)
    # Gaia DR3 variability tables
    try:
        for tab in ("gaiadr3.vari_eclipsing_binary", "gaiadr3.vari_classifier_result", "gaiadr3.vari_summary", "gaiadr3.nss_two_body_orbit", "gaiadr3.nss_non_linear_spectro", "gaiadr3.nss_acceleration_astro"):
            q = Gaia.launch_job(f"SELECT TOP 5 * FROM {tab} WHERE source_id = {sid}").get_results()
            if len(q):
                keep = [c for c in q.colnames if any(k in c for k in ("period", "frequency", "class", "solution", "eccentricity", "best"))]
                print(f"  {tab}: {len(q)} row(s):", {c: q[c][0] for c in keep})
            else:
                print(f"  {tab}: none")
    except Exception as ex:
        print("  Gaia error", ex)
    # VizieR: period-bearing catalogues within 3 arcsec
    try:
        c = SkyCoord.from_name(f"Gaia DR3 {sid}")
        v = Vizier(columns=["**"], row_limit=5)
        res = v.query_region(c, radius=3 * u.arcsec, catalog=["B/vsx", "J/ApJS/258/16", "J/A+A/691/A242", "J/ApJS/244/11", "J/ApJS/244/43", "J/ApJS/281/52", "II/366", "J/AJ/170/287", "J/ApJS/272/1", "J/MNRAS/519/5271", "J/A+A/674/A22", "J/other/RAA/21.292", "J/ApJS/249/22", "J/ApJS/282/55", "J/ApJS/280/26", "J/A+A/579/A19", "J/MNRAS/452/3561", "J/MNRAS/495/1252", "J/ApJS/256/14", "J/ApJS/258/26", "J/ApJS/266/40"])
        for k in res.keys():
            t = res[k]
            pc = [cn for cn in t.colnames if any(s_ in cn.lower() for s_ in ("per", "period", "type", "class", "epic"))]
            print(f"  VizieR {k}: {len(t)} row(s)", {cn: str(t[cn][0]) for cn in pc[:10]})
    except Exception as ex:
        print("  VizieR error", ex)
