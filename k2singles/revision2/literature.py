"""
AAS82225 revision, literature search (editor's point 1; PLAN.md) for all 163 LaCourse & Jacobs (2018) stars.

Per star:
  identifiers  every SIMBAD identifier of the object (TAP on the Harvard mirror, simbad.harvard.edu; the CDS host
               returned "Page unavailable" on 2026-11-15), found from "EPIC <n>", else "TIC <n>"; plus EPIC and TIC
               from the catalogue in any case
  ADS          one full-text query per star: full:("<id1>" OR "<id2>" ...) over all identifiers, plus the bare
               EPIC number (9 digits) and bare Gaia source_id (19 digits), which tables often use; with body
               highlights (snippets). If the query finds more than 150 papers, an identifier is too generic as a
               phrase; each is counted alone and those with more than 150 hits are dropped (listed in per_star.csv). Token from ADS_DEV_KEY or ~/.ads/dev_key (never printed).
  SIMBAD       the object's bibliography (has_ref)
Every paper is one row: EPIC, bibcode, year, title, first author, found_by (ads / simbad / both), ADS snippets, and
a keyword flag `period_words` (snippet mentions period, ephemeris, orbit or eclipse timing) for the manual review.
Output: k2singles/revision2/data/literature/hits.csv, identifiers.csv, per_star.csv. Resumable (per-star json cache).
Usage:  python k2singles/revision2/literature.py
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "data" / "literature"
CACHE = OUT / "cache"
SIMBAD_TAP = "https://simbad.harvard.edu/simbad/sim-tap/sync"
ADS = "https://api.adsabs.harvard.edu/v1/search/query"
AMBIGUOUS = 150
PERIOD_RE = re.compile(r"period|ephemer|orbit|P\s*=|\bdays?\b|eclips|transit", re.I)


def token():
    t = os.environ.get("ADS_DEV_KEY")
    if not t:
        f = Path.home() / ".ads" / "dev_key"
        t = f.read_text().strip() if f.exists() else None
    if not t:
        raise SystemExit("no ADS token (ADS_DEV_KEY or ~/.ads/dev_key)")
    return t


def tap(q):
    url = SIMBAD_TAP + "?" + urllib.parse.urlencode(dict(request="doQuery", lang="adql", format="json", query=q))
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.loads(r.read())
            return [row for row in d["data"]]
        except Exception:  # noqa: BLE001
            time.sleep(5 * (k + 1))
    raise RuntimeError("SIMBAD TAP failed")


def simbad(name):
    n = name.replace("'", "''")
    ids = [r[0] for r in tap(f"SELECT id2.id FROM ident AS id1 JOIN ident AS id2 ON id1.oidref = id2.oidref WHERE id1.id = '{n}'")]
    refs = [r[0] for r in tap(f"SELECT r.bibcode FROM ident AS i JOIN has_ref AS h ON h.oidref = i.oidref "
                               f"JOIN ref AS r ON r.oidbib = h.oidbibref WHERE i.id = '{n}'")] if ids else []
    return ids, refs


def ads(q, tok):
    params = dict(q=q, fl="id,bibcode,title,year,first_author,doctype", rows=200, hl="true", **{"hl.fl": "body,abstract"},
                  **{"hl.snippets": 3})
    url = ADS + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                d = json.loads(r.read())
            docs = d["response"]["docs"]
            hl = d.get("highlighting", {})
            return docs, hl, int(d["response"]["numFound"])
        except Exception as e:  # noqa: BLE001
            if "400" in str(e):
                raise RuntimeError(f"ADS rejected the query: {q[:300]}") from e
            if "429" in str(e):
                time.sleep(60)
            time.sleep(5 * (k + 1))
    raise RuntimeError("ADS query failed")


def ads_meta(bibcodes, tok):
    if not bibcodes:
        return {}
    out = {}
    for i in range(0, len(bibcodes), 50):
        chunk = bibcodes[i:i + 50]
        q = "bibcode:(" + " OR ".join(f'"{b}"' for b in chunk) + ")"
        docs, _, _ = ads(q, tok)
        out.update({d["bibcode"]: d for d in docs})
    return out


def ads_phrase(name):
    """SIMBAD name -> ADS-safe phrase: drop SIMBAD's object-type prefixes (V*, **, SV*, *, EM*, NAME), brackets,
    asterisks and quotes, which the ADS parser rejects (HTTP 400), and collapse spaces."""
    n = re.sub(r"^(V\*|\*\*|SV\*|\*|EM\*|NAME)\s+", "", name.strip())
    n = re.sub(r"[\[\]\*\"():{}^~\\]", " ", n)
    return re.sub(r"\s+", " ", n).strip()


def clean_snip(s):
    return re.sub(r"\s+", " ", re.sub(r"</?em>", "**", s)).strip()


def star(row, tok):
    epic, tic = int(row.EPIC), int(row.TIC)
    f = CACHE / f"{epic}.json"
    if f.exists():
        return json.loads(f.read_text())
    ids, refs = simbad(f"EPIC {epic}")
    if not ids:
        ids, refs = simbad(f"TIC {tic}")
    ids = sorted(set(re.sub(r"\s+", " ", i).strip() for i in ids) | {f"EPIC {epic}", f"TIC {tic}"})
    phrases = set(ids)
    phrases.add(str(epic))
    for i in ids:
        m = re.fullmatch(r"Gaia DR[23] (\d{15,20})", i)
        if m:
            phrases.add(m.group(1))
    phrases = sorted({c for c in (ads_phrase(p) for p in phrases) if len(c) >= 6})
    q = "full:(" + " OR ".join('"' + p.replace('"', "") + '"' for p in phrases) + ")"
    try:
        docs, hl, nfound = ads(q, tok)
    except RuntimeError as e:
        if "rejected" in str(e):
            raise
        # the combined query timed out (HTTP 504 on 2026-11-15 for EPIC 201207683): one query per identifier, merged
        docs, hl, seen = [], {}, set()
        for p in phrases:
            d1, h1, _ = ads(f'full:"{p}"', tok)
            for d in d1:
                if d["bibcode"] not in seen:
                    seen.add(d["bibcode"])
                    docs.append(d)
            hl.update(h1)
        nfound = len(docs)
        q = q + "  [split into one query per identifier after a timeout]"
    dropped = []
    if nfound > AMBIGUOUS:
        # some identifier is too generic as a phrase (e.g. "K2-3" tokenises to K2 3): count each one alone and drop
        # those with more than AMBIGUOUS hits, then query again
        for p in phrases:
            _, _, n1 = ads(f'full:"{p}"', tok)
            if n1 > AMBIGUOUS:
                dropped.append(f"{p} ({n1})")
        phrases = [p for p in phrases if not any(dd.startswith(p + " (") for dd in dropped)]
        q = "full:(" + " OR ".join('"' + p.replace('"', "") + '"' for p in phrases) + ")"
        docs, hl, nfound = ads(q, tok)
    hits = {}
    for d in docs:
        # ADS keys highlighting by its internal document id, not the bibcode (keyed by bibcode, every snippet was empty)
        sn = [clean_snip(x) for v in hl.get(str(d.get("id")), {}).values() for x in v]
        hits[d["bibcode"]] = dict(bibcode=d["bibcode"], year=d.get("year"), title=(d.get("title") or [""])[0],
                                  first_author=d.get("first_author"), doctype=d.get("doctype"), found_by="ads",
                                  snippets=" || ".join(sn))
    meta = ads_meta([b for b in refs if b not in hits], tok)
    for b in refs:
        if b in hits:
            hits[b]["found_by"] = "both"
        else:
            d = meta.get(b, {})
            hits[b] = dict(bibcode=b, year=d.get("year"), title=(d.get("title") or [""])[0], first_author=d.get("first_author"),
                           doctype=d.get("doctype"), found_by="simbad", snippets="")
    res = dict(EPIC=epic, TIC=tic, identifiers=ids, ads_query=q, ads_numfound=nfound, dropped_identifiers=dropped,
               n_simbad_refs=len(refs),
               hits=list(hits.values()))
    f.write_text(json.dumps(res, indent=1))
    return res


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    tok = token()
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC")
    rows, idrows, per = [], [], []
    for k, r in enumerate(d.itertuples(), 1):
        res = star(r, tok)
        for h in res["hits"]:
            h = dict(h)
            h["period_words"] = bool(PERIOD_RE.search(h["snippets"] or ""))
            rows.append(dict(EPIC=res["EPIC"], TIC=res["TIC"], **h))
        idrows.append(dict(EPIC=res["EPIC"], TIC=res["TIC"], identifiers=" | ".join(res["identifiers"])))
        per.append(dict(EPIC=res["EPIC"], TIC=res["TIC"], n_identifiers=len(res["identifiers"]), ads_numfound=res["ads_numfound"],
                        dropped_identifiers="; ".join(res.get("dropped_identifiers", [])),
                        n_simbad_refs=res["n_simbad_refs"], n_papers=len(res["hits"])))
        print(f"[{k}/{len(d)}] EPIC {res['EPIC']}: {len(res['identifiers'])} ids, ADS {res['ads_numfound']}, "
              f"SIMBAD {res['n_simbad_refs']}, papers {len(res['hits'])}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "hits.csv", index=False)
    pd.DataFrame(idrows).to_csv(OUT / "identifiers.csv", index=False)
    pd.DataFrame(per).to_csv(OUT / "per_star.csv", index=False)


if __name__ == "__main__":
    main()
