"""
Verify the bibliography of eb_depth_evolution/rnaas/rnaas.tex against ADS (token from ADS_DEV_KEY or the Windows
user environment; never printed). For each reference: query by first author, year, journal and volume or page, and
print the matching bibcode, authors (first 3) and title.
Output: eb_depth_evolution/rnaas/refs_check.txt
"""
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "rnaas" / "refs_check.txt"
REFS = {
    "Torres2000": 'author:"^Torres, G" year:2000 bibstem:AJ volume:119',
    "Lacy1999": 'author:"^Lacy" year:1999 bibstem:AJ volume:117',
    "Pollacco2006": 'author:"^Pollacco" year:2006 bibstem:PASP volume:118',
    "Pepper2007": 'author:"^Pepper, J" year:2007 bibstem:PASP volume:119',
    "Grindlay2009": 'author:"^Grindlay" year:2009 bibstem:ASPC volume:410',
    "LaCourse2015": 'author:"^LaCourse" year:2015 bibstem:MNRAS volume:452',
    "Shang2022": 'author:"^Shang" year:2022 bibstem:ApJS volume:259',
    "IJspeert2024": 'author:"^IJspeert" year:2024 bibstem:A&A volume:691',
    "Howell2014": 'author:"^Howell, S" year:2014 bibstem:PASP volume:126',
    "Caldwell2020": 'author:"^Caldwell, D" year:2020 bibstem:RNAAS volume:4',
    "Huang2020": 'author:"^Huang, C" year:2020 bibstem:RNAAS volume:4 page:204',
    "Kreidberg2015": 'author:"^Kreidberg" year:2015 bibstem:PASP volume:127',
    "Ginsburg2019": 'author:"^Ginsburg" year:2019 bibstem:AJ volume:157',
    "Gaia2023": 'year:2023 bibstem:A&A volume:674 page:A1',
    "Zasche2023": 'author:"^Zasche" year:2023 bibstem:AJ volume:165',
    "Davenport2021": 'author:"^Davenport" year:2021 bibstem:AJ volume:162',
    "Luger2016": 'author:"^Luger" year:2016 bibstem:AJ volume:152 page:100',
    "Vanderburg2014": 'author:"^Vanderburg" year:2014 bibstem:PASP volume:126 page:948',
    "Butters2010": 'author:"^Butters" year:2010 bibstem:A&A volume:520',
    "Borkovits2016": 'author:"^Borkovits" year:2016 bibstem:MNRAS volume:455',
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
    for key, q in REFS.items():
        url = ("https://api.adsabs.harvard.edu/v1/search/query?" +
               urllib.parse.urlencode({"q": q, "fl": "bibcode,author,title,page,volume", "rows": 5}))
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
        with urllib.request.urlopen(req, timeout=60) as r:
            docs = json.load(r)["response"]["docs"]
        lines.append(f"## {key}: {len(docs)} match(es) for {q}")
        for d in docs:
            au = "; ".join((d.get("author") or [])[:3])
            lines.append(f"  {d['bibcode']}  p.{(d.get('page') or ['?'])[0]}  {au}  | {(d.get('title') or [''])[0][:90]}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines).encode("ascii", "replace").decode())


if __name__ == "__main__":
    main()
