"""
ADS API novelty check for the southern ETV candidates (and positive controls).

For each star, query the ADS full-text index (`full:`, which covers abstract + body + ack where ADS has
full text) with every catalogue identifier we know. (`object:` queries are not accepted by the plain search
API; the SIMBAD reference lists were checked directly instead, see notes/ncvz_unreported_check.md.)
Highlights from the body/abstract are saved, so each hit can be judged without opening the paper.

Positive controls are stars known to be tabulated in Mitnyan+2024 / Borkovits+2026 / Marcadon & Prsa 2024;
they measure how well the full-text index finds identifiers that sit in tables.

Token: ADS_DEV_KEY from the environment, else the Windows user environment (HKCU\\Environment).
It is never printed or written out.

Output: notes/ads_api_check.json
Usage:  python sidequests/eb_timing/scripts/ads_check.py
"""
import json
import os
import time
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "notes" / "ads_api_check.json"
API = "https://api.adsabs.harvard.edu/v1/search/query"

STARS = {
    "TIC 32000625": ["32000625", "TYC 9155-22-1", "TYC 9155-0022-1", "4642053155197039104",
                     "2MASS J03435710-7128445", "J03435710-7128445"],
    "TIC 279323859": ["279323859", "UCAC4 155-008982", "5480636311333638272",
                      "2MASS J06565890-5909400", "J06565890-5909400"],
    "TIC 278861826": ["278861826", "TYC 8544-1497-1", "GSC 08544-01497", "5484778171633161088",
                      "2MASS J06423805-5632306", "J06423805-5632306"],
    "TIC 294327312": ["294327312", "HD 56910", "HD56910", "HIP 35095", "HIP35095", "CD-58 1696",
                      "CPD-58 885", "5485783679313899136", "2MASS J07152908-5832334", "J07152908-5832334"],
}
# Positive controls: (star, identifiers, papers that are known to list it)
CONTROLS = {
    "TIC 424461577": (["424461577"], ["Mitnyan+2024 (A&A)", "Borkovits+2026 (Universe)", "Marcadon & Prsa 2024 (ApJ)"]),
    "TIC 229771234": (["229771234"], ["Mitnyan+2024", "Borkovits+2026", "Marcadon & Prsa 2024"]),
    "TIC 167692429": (["167692429"], ["Borkovits+2020 (MNRAS)"]),
    "TIC 220397947": (["220397947"], ["Borkovits+2020 (MNRAS)"]),
}


def token():
    t = os.environ.get("ADS_DEV_KEY") or os.environ.get("ADS_API_TOKEN")
    if t:
        return t
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as k:
            return winreg.QueryValueEx(k, "ADS_DEV_KEY")[0]
    except OSError:
        raise SystemExit("no ADS_DEV_KEY found")


def query(q, tok):
    params = {"q": q, "fl": "id,bibcode,title,year,doctype,pub", "rows": 200, "sort": "date desc",
              "hl": "true", "hl.fl": "body,abstract,ack", "hl.snippets": 3, "hl.fragsize": 200}
    for attempt in range(3):
        r = requests.get(API, params=params, headers={"Authorization": f"Bearer {tok}"}, timeout=60)
        if r.status_code == 200:
            break
        time.sleep(5)
    r.raise_for_status()
    j = r.json()
    hl = j.get("highlighting", {})
    docs = []
    for d in j["response"]["docs"]:
        h = hl.get(d.get("id", ""), {}) if hl else {}
        docs.append({"bibcode": d["bibcode"], "year": d.get("year"), "pub": d.get("pub"),
                     "doctype": d.get("doctype"), "title": (d.get("title") or [""])[0],
                     "highlights": [s for v in h.values() for s in v]})
    return j["response"]["numFound"], docs


def run(name, ids, tok):
    full = " OR ".join(f'full:"{i}"' for i in ids)
    res = {}
    for label, q in [("full", full)]:
        n, docs = query(q, tok)
        res[label] = {"query": q, "n": n, "docs": docs}
        print(f"  {label:6s} {n:3d}  {q[:110]}")
        for d in docs:
            print(f"      {d['bibcode']}  {d['title'][:90]}")
    return res


def main():
    tok = token()
    out = {"date": time.strftime("%Y-%m-%d"), "stars": {}, "controls": {}}
    print("== candidates")
    for name, ids in STARS.items():
        print(name)
        out["stars"][name] = run(name, ids, tok)
    print("== positive controls")
    for name, (ids, known) in CONTROLS.items():
        print(name, "known in:", "; ".join(known))
        out["controls"][name] = run(name, ids, tok) | {"known_in": known}
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("->", OUT)


if __name__ == "__main__":
    main()
