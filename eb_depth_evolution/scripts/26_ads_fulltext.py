"""
ADS full-text novelty search for EPIC 202126877 / TIC 86491361.
The token comes from the ADS_DEV_KEY environment variable, or the Windows user-level variable in the registry
(HKCU\\Environment). It is never printed or written anywhere.
Output: eb_depth_evolution/ads_fulltext_results.txt
Usage:  python eb_depth_evolution/scripts/26_ads_fulltext.py
"""
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "ads_fulltext_results.txt"

QUERIES = [
    'full:"202126877"',
    'full:"86491361"',
    'full:"1883-1721-1"',
    'full:"3383995253529149056"',
    'full:"EPIC 202126877"',
    'full:"J063155.2+252855"',
    'object:"EPIC 202126877"',
    # topical: eclipse depth variations / eclipses appearing in TESS vs K2
    'abs:("eclipse depth variation" OR "disappearing eclipses" OR "appearing eclipses") year:2020-2026',
]


def token():
    t = os.environ.get("ADS_DEV_KEY")
    if t:
        return t
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as k:
            return winreg.QueryValueEx(k, "ADS_DEV_KEY")[0]
    except OSError:
        return None


def search(q, tok, rows=50):
    url = ("https://api.adsabs.harvard.edu/v1/search/query?" +
           urllib.parse.urlencode({"q": q, "fl": "bibcode,title,first_author,year", "rows": rows, "sort": "date desc"}))
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["response"]


def main():
    tok = token()
    if not tok:
        raise SystemExit("ADS_DEV_KEY not found")
    lines = []
    for q in QUERIES:
        try:
            r = search(q, tok)
        except Exception as e:  # noqa: BLE001
            lines.append(f"## {q}\n  ERROR {type(e).__name__}: {str(e)[:120]}")
            continue
        lines.append(f"## {q}  ->  {r['numFound']} hits")
        for d in r["docs"][:50]:
            lines.append(f"  {d.get('bibcode')}  {d.get('first_author', '')}  {(d.get('title') or [''])[0][:110]}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines).encode("ascii", "replace").decode())


if __name__ == "__main__":
    main()
