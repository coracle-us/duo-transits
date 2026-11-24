"""
Check the bibliography of salinas-second-transits/rnaas/rnaas.tex against ADS, and Salinas et al. 2026 against
arXiv too (for a journal reference). ADS token from ADS_DEV_KEY or ~/.ads/dev_key (never printed).
For each reference: query by first author, year, journal and volume or page; print bibcode, first 3 authors, title.
Output: salinas-second-transits/rnaas/checks/refs_check.txt
Usage:  python salinas-second-transits/rnaas/checks/verify_refs.py
"""
import json
import os
import re
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / "refs_check.txt"
REFS = {
    "Brasseur2019": 'bibstem:ascl year:2019 "astrocut"',
    "Caldwell2020": 'author:"^Caldwell, D" year:2020 bibstem:RNAAS volume:4 page:201',
    "Gaia2023": 'year:2023 bibstem:A&A volume:674 page:A1',
    "Ginsburg2019": 'author:"^Ginsburg" year:2019 bibstem:AJ volume:157 page:98',
    "Hippke2019": 'author:"^Hippke" year:2019 bibstem:AJ volume:158 page:143',
    "Huang2020": 'author:"^Huang, C" year:2020 bibstem:RNAAS volume:4 page:204',
    "Lightkurve2018": 'bibstem:ascl year:2018 title:"Lightkurve"',
    "Prsa2022": 'author:"^Prsa" year:2022 bibstem:ApJS volume:258 page:16',
    "Salinas2026": 'identifier:"arXiv:2610.24882"',
    "Seager2003": 'author:"^Seager" year:2003 bibstem:ApJ volume:585 page:1038',
    "Watson2006": 'author:"^Watson, C" year:2006 bibstem:SASS volume:25 page:47',
}


def token():
    t = os.environ.get("ADS_DEV_KEY") or os.environ.get("ADS_API_TOKEN")
    return t or (Path.home() / ".ads" / "dev_key").read_text().strip()


def ads(q, tok):
    url = ("https://api.adsabs.harvard.edu/v1/search/query?" +
           urllib.parse.urlencode({"q": q, "fl": "bibcode,author,title,page,volume,pub,doi,identifier", "rows": 5}))
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["response"]["docs"]


def arxiv(aid):
    with urllib.request.urlopen(f"https://export.arxiv.org/api/query?id_list={aid}", timeout=60) as r:
        x = r.read().decode()
    get = lambda tag: (re.findall(rf"<{tag}[^>]*>(.*?)</{tag}>", x, re.S) or [""])
    return dict(title=" ".join(get("title")[-1].split()), updated=get("updated")[-1],
                journal_ref=get("arxiv:journal_ref")[0], doi=get("arxiv:doi")[0],
                authors=[a.strip() for a in re.findall(r"<name>(.*?)</name>", x)][:4])


def main():
    tok = token()
    lines = []
    for key, q in REFS.items():
        docs = ads(q, tok)
        lines.append(f"## {key}: {len(docs)} match(es) for {q}")
        for d in docs:
            au = "; ".join((d.get("author") or [])[:3])
            lines.append(f"  {d['bibcode']}  p.{(d.get('page') or ['?'])[0]}  {au}  | {(d.get('title') or [''])[0][:90]}"
                         f"  | pub: {d.get('pub')}  doi: {(d.get('doi') or [''])[0]}")
    a = arxiv("2610.24882")
    lines.append(f"## arXiv 2610.24882 (export.arxiv.org): {a}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT.read_text())


if __name__ == "__main__":
    main()
