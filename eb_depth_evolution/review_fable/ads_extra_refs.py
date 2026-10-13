"""
ADS check of references the draft does NOT cite but arguably should (switching-eclipse systems, EDV catalogues,
pipelines named in the text). Token from ADS_DEV_KEY or the Windows user environment (never printed).
Output: eb_depth_evolution/review_fable/ads_extra_refs.txt
"""
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).with_name("ads_extra_refs.txt")
QUERIES = {
    "SS Lac third body (Torres 2001)": 'author:"^Torres, G" year:2001 title:"SS Lacertae"',
    "V907 Sco revisited": 'title:"V907 Sco" year:2000-2026',
    "HS Hya vanishing eclipses (Zasche & Paschke 2012)": 'title:"HS Hya" year:2012',
    "HS Hya TESS (Davenport 2021)": 'author:"^Davenport" title:"HS Hya"',
    "Borkovits 2016 Kepler triples": 'author:"^Borkovits" year:2016 bibstem:MNRAS title:"triple"',
    "Borkovits 2026 northern CVZ EDVs": 'author:"^Borkovits" year:2025-2026 abs:"eclipse depth"',
    "Baum 2025 KELT+TESS EDVs": 'author:"^Baum" year:2025 abs:"eclipse depth variation"',
    "Kostov 2025 TESS Ten-Thousand": 'author:"^Kostov" year:2025 title:"ten thousand"',
    "Quadruples changing inclination 2025": 'year:2025 abs:"changing inclination" abs:"quadruple"',
    "Juryšek 2018 dynamical EDV (V907 Sco / others)": 'author:"^Jurysek" year:2018 bibstem:A&A',
    "KIC 2856960 Marsh 2014": 'author:"^Marsh, T" year:2014 title:"KIC 2856960"',
    "Butters 2010 WASP public archive": 'author:"^Butters" year:2010 bibstem:A&A',
    "Luger 2016 EVEREST": 'author:"^Luger" year:2016 title:"EVEREST"',
    "Vanderburg & Johnson 2014 K2SFF": 'author:"^Vanderburg" year:2014 bibstem:PASP',
    "Shappee 2014 ASAS-SN": 'author:"^Shappee" year:2014 bibstem:ApJ title:"ASAS-SN"',
    "Soderhjelm 1975 nodal precession": 'author:"^Soderhjelm" year:1975',
    "AAS journals AI policy / RNAAS guidelines (sanity)": 'bibstem:RNAAS year:2026 title:"Research Notes"',
    "Eclipse depth variation TESS 2022-2026": 'abs:"eclipse depth variations" year:2022-2026 property:refereed',
}


def token():
    t = os.environ.get("ADS_DEV_KEY")
    if t:
        return t
    import winreg
    with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as k:
        return winreg.QueryValueEx(k, "ADS_DEV_KEY")[0]


def main():
    tok = token()
    lines = []
    for key, q in QUERIES.items():
        url = ("https://api.adsabs.harvard.edu/v1/search/query?" +
               urllib.parse.urlencode({"q": q, "fl": "bibcode,author,title,page,volume,year,pub", "rows": 8, "sort": "citation_count desc"}))
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                docs = json.load(r)["response"]["docs"]
        except Exception as ex:  # noqa: BLE001
            lines.append(f"## {key}: ERROR {ex}")
            continue
        lines.append(f"## {key}: {len(docs)} match(es) for {q}")
        for d in docs:
            au = "; ".join((d.get("author") or [])[:3])
            lines.append(f"  {d['bibcode']}  {d.get('pub', '')} v{d.get('volume', '?')} p.{(d.get('page') or ['?'])[0]}  {au}  | {(d.get('title') or [''])[0][:110]}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines).encode("ascii", "replace").decode())


if __name__ == "__main__":
    main()
