"""
AAS82225 revision, Search 1 (PLAN.md): blind single-event search of every K2 campaign and TESS sector of all 163
LaCourse & Jacobs (2018) stars. No catalogued depth, duration or epoch enters detection.

Detection code and rules are those of singles/revision_search/search_any.py (imported, not copied):
free depth, durations 1-48 h, red-noise SNR >= 8 and the other rules in its docstring, identical search on 2 - flux
as the control.
  TESS  best product per sector searched (SPOC 2-min > TESS-SPOC > QLP > OWN), each detection re-measured in the
        other products (search_any.annotate)
  K2    EVEREST, K2SFF and mission PDCSAP each searched separately; detections merged across products and each
        re-measured in the other two ("confirmed": SNR >= 4 and depth 0.5-2x, as search_any)
The catalogued event is used only afterwards: blind_hit = a blind event within max(0.1 d, D_fit/2) of the
catalogued epoch; and, for reference, a free measurement at the catalogued time and duration in each K2 product.

Usage:  python k2singles/revision2/search_blind.py [--epics ...] [--workers N]
Output: k2singles/revision2/data/per_star/<EPIC>.json
"""
import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "singles" / "revision_search"))
import search_any as sa  # noqa: E402

DATA = HERE / "data"
OUT = DATA / "per_star"
K2CACHE = ROOT / "k2singles" / "cache"
HLSP = DATA / "hlsp"
DURS_H = [1, 2, 4, 8, 16, 24, 48]
K2_PRODUCTS = ["EVEREST", "K2SFF", "K2PDC"]
NO_DUMPS = np.array([])


def k2_products(epic):
    """prods[1000 + campaign][product] = dict(t, f, seg), as search_any.load_products."""
    prods = {}

    def add(camp, key, t, f):
        o = np.argsort(t)
        t, f = np.asarray(t, float)[o], np.asarray(f, float)[o]
        ok = np.isfinite(t) & np.isfinite(f)
        t, f = t[ok], f[ok]
        if len(t) < 50:
            return
        sec = 1000 + camp
        seg = np.concatenate([[0], np.cumsum(np.diff(t) > sa.SEG_GAP)])
        prods.setdefault(sec, {})[key] = dict(t=t, f=f / np.nanmedian(f), seg=sec * 1000 + seg)

    for fn in sorted((K2CACHE / str(epic)).glob("k2_c*.npz")):
        z = np.load(fn)
        add(int(fn.stem.split("_c")[1]), "K2PDC", z["t"], z["f"])
    for fn in sorted((HLSP / str(epic)).glob("*_c*.npz")):
        key, c = fn.stem.split("_c")
        z = np.load(fn)
        add(int(c), key, z["t"], z["f"])
    return prods


def lc_of(prods, key):
    parts = [(p[key]["t"], p[key]["f"], p[key]["seg"], np.full(len(p[key]["t"]), s))
             for s, p in sorted(prods.items()) if key in p]
    if not parts:
        return None, {}
    t = np.concatenate([x[0] for x in parts])
    o = np.argsort(t)
    lc = dict(t=t[o], f=np.concatenate([x[1] for x in parts])[o], seg=np.concatenate([x[2] for x in parts])[o],
              sec=np.concatenate([x[3] for x in parts])[o])
    return lc, {s: key for s, p in prods.items() if key in p}


def k2_search(prods, durs, sign):
    """Search each K2 product separately, merge across products, refine, cross-check."""
    allev, why = [], {}
    segs = {}
    for key in K2_PRODUCTS:
        lc, used = lc_of(prods, key)
        if lc is None:
            continue
        segs[key] = sa.segments(lc)
        evs, w = sa.search(lc, used, NO_DUMPS, durs, sign)
        for k, v in w.items():
            why[f"{key}:{k}"] = why.get(f"{key}:{k}", 0) + v
        allev += evs
    merged = sa.merge(allev)
    for ev in merged:
        ev["found_in"] = sorted({e["product"] for e in allev if abs(e["T"] - ev["T"]) < (e["D"] + ev["D"]) / 2})
        sa.refine(ev, segs[ev["product"]], sign)
        alt = {}
        for key in K2_PRODUCTS:
            if key == ev["product"] or key not in prods.get(ev["sector"], {}):
                continue
            p = prods[ev["sector"]][key]
            j = np.searchsorted(p["t"], ev["T"])
            if j <= 0 or j >= len(p["t"]):
                alt[key] = dict(status="no-coverage")
                continue
            m = p["seg"] == p["seg"][j]
            t, f = p["t"][m], (p["f"][m] if sign > 0 else 2.0 - p["f"][m])
            r = sa.test_free(t, f, ev["T"], ev["D"], scan=ev["D"] / 8)
            if r is None or r["cover"] < 0.5:
                alt[key] = dict(status="no-coverage")
                continue
            conf = r["snr"] >= sa.CONF_SNR and sa.CONF_LO * ev["depth"] <= r["depth"] <= sa.CONF_HI * ev["depth"]
            alt[key] = dict(status="confirmed" if conf else "not-confirmed", depth_ppm=round(r["depth"] * 1e6, 1),
                            snr=round(r["snr"], 1))
        ev["alt"] = alt
        conf = [k for k, v in alt.items() if v.get("status") == "confirmed"]
        ev["n_reductions"] = 1 + len(conf)
        ev["n_alt_tested"] = sum(v.get("status") in ("confirmed", "not-confirmed") for v in alt.values())
        ev["independent_confirm"] = len(conf) > 0
        ev["mission"] = "K2"
    return merged, why


def tess_search(prods, dumps, durs, sign):
    lc, used = sa.primary(prods)
    if lc is None:
        return [], {}, {}
    evs, why = sa.search(lc, used, dumps, durs, sign)
    evs = sa.annotate(evs, prods, sign, sa.segments(lc))
    for ev in evs:
        ev["mission"] = "TESS"
        ev["found_in"] = [ev["product"]]
    return evs, why, used


def catalogued_free(prods_k2, Tev, Dcat):
    """For reference only: free depth at the catalogued time and duration in each K2 product."""
    out = {}
    for sec, p in prods_k2.items():
        for key, d in p.items():
            if not (d["t"][0] < Tev < d["t"][-1]):
                continue
            j = np.searchsorted(d["t"], Tev)
            m = d["seg"] == d["seg"][min(j, len(d["t"]) - 1)]
            r = sa.test_free(d["t"][m], d["f"][m], Tev, Dcat, scan=Dcat / 4)
            if r is not None:
                ok, reason = sa.passes(r, NO_DUMPS)
                out[key] = dict(T=r["T"], depth_ppm=round(r["depth"] * 1e6, 1), snr=round(r["snr"], 1),
                                cover=round(r["cover"], 2), passes=ok, reason=reason)
    return out


def clean(evs):
    for e in evs:
        e.pop("_tin", None)
        e.pop("_rin", None)
    return evs


def run_star(spec):
    epic, tic = int(spec["EPIC"]), int(spec["TIC"])
    t0 = time.time()
    durs = [d / 24 for d in DURS_H]
    prods_t, dumps = sa.load_products(tic)
    prods_k = k2_products(epic)
    rep = dict(spec, tess_sectors=sorted(prods_t), k2_campaigns=sorted(s - 1000 for s in prods_k),
               k2_products={str(s - 1000): sorted(p) for s, p in prods_k.items()})
    ev_k, why_k = k2_search(prods_k, durs, +1.0)
    ct_k, why_kc = k2_search(prods_k, durs, -1.0)
    ev_t, why_t, used = tess_search(prods_t, dumps, durs, +1.0)
    ct_t, why_tc, _ = tess_search(prods_t, dumps, durs, -1.0)
    rep["tess_products_used"] = {str(k): v for k, v in used.items()}
    rep["events"] = clean(sorted(ev_k + ev_t, key=lambda e: e["T"]))
    rep["control_events"] = clean(sorted(ct_k + ct_t, key=lambda e: e["T"]))
    rep["near_miss"] = dict(k2=why_k, tess=why_t, k2_control=why_kc, tess_control=why_tc)
    events = []
    for i, (T, D) in enumerate(zip(spec["epochs"], spec["durs_h"])):
        hits = [e for e in rep["events"] if abs(e["T_fit"] - T) < max(0.1, e["D_fit"] / 2, D / 48)]
        events.append(dict(event_idx=i, epoch=T, blind_hit=bool(hits),
                           blind_T=hits[0]["T_fit"] if hits else None,
                           blind_depth_ppm=round(hits[0]["depth_fit"] * 1e6, 1) if hits else None,
                           blind_D_h=round(hits[0]["D_fit"] * 24, 2) if hits else None,
                           blind_products=hits[0]["found_in"] if hits else [],
                           catalogued_free=catalogued_free(prods_k, T, D / 24)))
    rep["catalogued_events"] = events
    rep["runtime_s"] = round(time.time() - t0, 1)
    rep["status"] = "ok"
    (OUT / f"{epic}.json").write_text(json.dumps(rep, indent=1, default=float))
    return rep


def run_safe(spec):
    try:
        return run_star(spec)
    except Exception:  # noqa: BLE001
        rep = dict(spec, status="error", error=traceback.format_exc()[-1500:])
        (OUT / f"{int(spec['EPIC'])}.json").write_text(json.dumps(rep, indent=1, default=float))
        return rep


def specs():
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    out = []
    for epic, g in d.groupby("EPIC", sort=False):
        out.append(dict(EPIC=int(epic), TIC=int(g.TIC.iloc[0]), epochs=[float(x) for x in g.epoch_btjd],
                        durs_h=[float(x) for x in g.dur_h], depths_ppm=[float(x) for x in g.depth_ppm],
                        campaign=int(g.campaign.iloc[0]), comment=str(g.comment.iloc[0])))
    return out


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--epics", nargs="*", type=int)
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--force", action="store_true")
    args = a.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(DATA / "fetch_status.csv")
    fetched = set(st[st.status == "done"].EPIC)
    todo = [s for s in specs() if s["EPIC"] in fetched and (not args.epics or s["EPIC"] in args.epics)
            and (args.force or not (OUT / f"{s['EPIC']}.json").exists())]
    print(f"{len(todo)} stars to search, {args.workers} workers", flush=True)
    with ProcessPoolExecutor(args.workers, max_tasks_per_child=1) as ex:
        futs = [ex.submit(run_safe, s) for s in todo]
        for k, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            print(f"[{k}/{len(todo)}] EPIC {r['EPIC']} {r['status']} K2 {r.get('k2_campaigns')} "
                  f"TESS {len(r.get('tess_sectors', []))} sectors | events {len(r.get('events', []))} "
                  f"control {len(r.get('control_events', []))} | catalogued hit "
                  f"{[e['blind_hit'] for e in r.get('catalogued_events', [])]} | {r.get('runtime_s')} s "
                  f"{r.get('error', '')[-300:]}", flush=True)


if __name__ == "__main__":
    main()
