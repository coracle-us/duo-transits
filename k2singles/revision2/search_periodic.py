"""
AAS82225 revision, stage 2 (PLAN.md): significance from the inverted control, then two blind periodic searches.

  thresholds  for K2 and TESS separately: S* = smallest SNR at which eligible CONTROL events number <= 0.05 per
              searched star (eligible = n_reductions >= 2, or n_alt_tested == 0). Written to thresholds.json.
  significant an event that is eligible, has SNR >= S* of its mission and is not common-mode (within 0.5 d of a
              first-stage event on >= 2 other stars in the same campaign/sector)
  linking     per star, significant events grouped into families (depth within x2, duration within 1/2-2x of the
              family's first, highest-SNR, member). For every pair in a family, P = dT/n >= 1 d; a period is
              supported if every family event lies on the ephemeris within max(0.1 d, D/4); it is excluded if a
              clean predicted epoch (coverage >= 0.8, not within D of a segment edge) has box depth
              < 0.5 x family depth - 3 sigma. Box depths come from a grid (step D/8) over each star's primary
              light curves (TESS best product per sector; K2 EVEREST, else K2SFF, else PDCSAP).
  BLS         astropy BoxLeastSquares per data block (data separated by gaps > 60 d), each block divided by a
              running median (window 3 d) per segment; periods 1 d to min(100 d, block span / 2), durations
              1, 2, 4, 8 and 16 h; SDE = (peak - median) / std of the power. Threshold: the 99th percentile of the
              maximum SDE over all blocks of the INVERTED light curves. Up to three signals per block, masking
              each (3 durations wide, at most half the period) before the next. (Per-block, not on the combined light curve as first
              planned: a combined 12-yr baseline needs ~4e8 trial periods; see PLAN.md.)
Usage:  python k2singles/revision2/search_periodic.py [--workers N]
Output: k2singles/revision2/data/thresholds.json, periodic/<EPIC>.json, periodic_summary.csv
"""
import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "singles" / "revision_search"))
import search_any as sa  # noqa: E402
import search_blind as sb  # noqa: E402

DATA = HERE / "data"
PER = DATA / "per_star"
OUTP = DATA / "periodic"
FA_PER_STAR = 0.05
CM_WINDOW = 0.5
P_MIN = 1.0
BLS_DURS_H = [1, 2, 4, 8, 16]
BLOCK_GAP = 60.0
BLS_PMAX = 100.0
K2_RANK = ["EVEREST", "K2SFF", "K2PDC"]


# ------------------------------------------------------------------ stage 2a: significance
def eligible(e):
    return e.get("n_reductions", 1) >= 2 or e.get("n_alt_tested", 0) == 0


def thresholds(reps):
    out = {}
    for mission in ("K2", "TESS"):
        n_stars = sum(1 for r in reps if any(e["mission"] == mission for e in r["events"] + r["control_events"])
                      or (mission == "K2" and r.get("k2_campaigns")) or (mission == "TESS" and r.get("tess_sectors")))
        snr = np.sort([e["snr"] for r in reps for e in r["control_events"] if e["mission"] == mission and eligible(e)])[::-1]
        allowed = int(np.floor(FA_PER_STAR * n_stars))
        # S* is the smallest SNR with at most `allowed` control events at or above it
        s_star = float(snr[allowed] * (1 + 1e-9)) if len(snr) > allowed else 8.0
        s_star = max(s_star, 8.0)
        out[mission] = dict(S_star=round(s_star, 2), n_stars=n_stars, n_control_eligible=int(len(snr)),
                            allowed_control_events=allowed,
                            control_snr_top=[round(x, 1) for x in snr[:allowed + 5]])
    return out


CM_P = 1e-3


def common_mode(reps, key="events"):
    """(EPIC, T) of common-mode events (PLAN.md, 2026-11-15 replacement rule): k = number of OTHER stars with an
    event within CM_WINDOW in the same sector/campaign; expected lambda = sum over other stars of
    n_events * 2 CM_WINDOW / span (span = the sector's event time range); flagged if k >= 2 and
    P(K >= k | lambda) < CM_P."""
    from scipy.stats import poisson
    rows = [(r["EPIC"], e["sector"], e["T"]) for r in reps for e in r[key]]
    df = pd.DataFrame(rows, columns=["EPIC", "sector", "T"])
    flag = set()
    for _, g in df.groupby("sector"):
        span = max(g["T"].max() - g["T"].min(), 1.0)
        n_by = g.groupby("EPIC").size()
        for epic, T in zip(g.EPIC, g["T"]):
            k = g[(g.EPIC != epic) & (np.abs(g["T"] - T) < CM_WINDOW)].EPIC.nunique()
            if k < 2:
                continue
            lam = float(n_by.drop(epic, errors="ignore").sum()) * 2 * CM_WINDOW / span
            if poisson.sf(k - 1, lam) < CM_P:
                flag.add((epic, round(T, 5)))
    return flag


# ------------------------------------------------------------------ light curves for epoch tests and BLS
def primary_lcs(rep):
    tic, epic = int(rep["TIC"]), int(rep["EPIC"])
    prods_t, _ = sa.load_products(tic)
    lc_t, _ = sa.primary(prods_t)
    prods_k = sb.k2_products(epic)
    parts = []
    if lc_t is not None:
        parts.append(lc_t)
    for sec, p in sorted(prods_k.items()):
        key = next((k for k in K2_RANK if k in p), None)
        if key:
            d = p[key]
            parts.append(dict(t=d["t"], f=d["f"], seg=d["seg"], sec=np.full(len(d["t"]), sec)))
    if not parts:
        return None
    t = np.concatenate([p["t"] for p in parts])
    o = np.argsort(t)
    return dict(t=t[o], f=np.concatenate([p["f"] for p in parts])[o], seg=np.concatenate([p["seg"] for p in parts])[o])


def box_grid(lc, D):
    """Box depth, sigma and coverage on a D/8 grid in every segment (local linear baseline, as search_any.fast_pass,
    sigma from the 1.4826*MAD point scatter times a red-noise factor at timescale D)."""
    Ts, dep, sig, cov = [], [], [], []
    for u in np.unique(lc["seg"]):
        m = lc["seg"] == u
        t, f = lc["t"][m], lc["f"][m]
        if len(t) < 20 or t[-1] - t[0] < 3 * D:
            continue
        cad = float(np.median(np.diff(t)))
        T = np.arange(t[0] + D, t[-1] - D + 1e-9, D / 8)
        if not len(T):
            continue
        B = max(0.75, 2 * D)
        cf, ct = np.concatenate([[0], np.cumsum(f)]), np.concatenate([[0], np.cumsum(t)])

        def box(a, b):
            i, j = np.searchsorted(t, a), np.searchsorted(t, b)
            n = j - i
            with np.errstate(invalid="ignore", divide="ignore"):
                return n, (cf[j] - cf[i]) / n, (ct[j] - ct[i]) / n

        n_in, m_in, _ = box(T - sa.FLAT * D / 2, T + sa.FLAT * D / 2)
        mh = D / 2 + 0.02
        nL, mL, tL = box(T - mh - B, T - mh)
        nR, mR, tR = box(T + mh, T + mh + B)
        sp = 1.4826 * np.median(np.abs(np.diff(f) - np.median(np.diff(f)))) / np.sqrt(2)
        r = f / np.median(f) - 1
        beta = sa.red_noise_beta(r - pd.Series(r).rolling(max(3, int(2 * B / cad)), center=True, min_periods=3).median().to_numpy(),
                                 max(1, int(round(sa.FLAT * D / cad))))
        beta = 1.0 if not np.isfinite(beta) else beta
        with np.errstate(invalid="ignore", divide="ignore"):
            base = np.where((nL >= 3) & (nR >= 3), mL + (mR - mL) * (T - tL) / (tR - tL),
                            np.where(nL >= 3, mL, mR))
            d = base - m_in
            s = beta * sp * np.sqrt(1 / n_in + 0.25 / np.maximum(nL, 1) + 0.25 / np.maximum(nR, 1))
        c = np.minimum(1.0, n_in * cad / (sa.FLAT * D))
        ok = (n_in >= 2) & ((nL >= 3) | (nR >= 3))
        Ts.append(T[ok]); dep.append(d[ok]); sig.append(s[ok]); cov.append(c[ok])
    if not Ts:
        return np.array([]), np.array([]), np.array([]), np.array([])
    T = np.concatenate(Ts)
    o = np.argsort(T)
    return T[o], np.concatenate(dep)[o], np.concatenate(sig)[o], np.concatenate(cov)[o]


def link(sig_events, lc):
    """Families and periods (stage 2a). Returns a list of family dicts."""
    evs = sorted(sig_events, key=lambda e: -e["snr"])
    fams = []
    for e in evs:
        for fam in fams:
            h = fam["members"][0]
            if 0.5 <= e["depth_fit"] / h["depth_fit"] <= 2 and 0.5 <= e["D_fit"] / h["D_fit"] <= 2:
                fam["members"].append(e)
                break
        else:
            fams.append(dict(members=[e]))
    out = []
    for fam in fams:
        mem = sorted(fam["members"], key=lambda e: e["T_fit"])
        Tm = np.array([e["T_fit"] for e in mem])
        D = float(np.median([e["D_fit"] for e in mem]))
        d0 = float(np.median([e["depth_fit"] for e in mem]))
        d0_err = d0 / float(np.median([e["snr"] for e in mem]))
        res = dict(n_members=len(mem), T=[round(x, 4) for x in Tm], D_h=round(D * 24, 2), depth_ppm=round(d0 * 1e6, 1),
                   missions=[e["mission"] for e in mem], periods=[])
        if len(mem) >= 2 and lc is not None:
            Tg, dg, sg, cg = box_grid(lc, D)
            tol = max(0.1, D / 4)
            seen = set()
            for i in range(len(Tm)):
                for j in range(i + 1, len(Tm)):
                    dT = Tm[j] - Tm[i]
                    if dT < P_MIN:
                        continue
                    for n in range(1, int(dT / P_MIN) + 1):
                        P = dT / n
                        ph = (Tm - Tm[i]) / P
                        if np.any(np.abs(ph - np.round(ph)) * P > tol):
                            continue
                        key = round(P, 3)
                        if key in seen:
                            continue
                        seen.add(key)
                        # predicted epochs over the data
                        excl, n_clean = False, 0
                        if len(Tg):
                            k0, k1 = np.ceil((Tg[0] - Tm[i]) / P), np.floor((Tg[-1] - Tm[i]) / P)
                            Tk = Tm[i] + np.arange(k0, k1 + 1) * P
                            Tk = Tk[np.min(np.abs(Tk[:, None] - Tm[None, :]), axis=1) > tol] if len(Tk) else Tk
                            idx = np.clip(np.searchsorted(Tg, Tk), 1, len(Tg) - 1)
                            idx = np.where(np.abs(Tg[idx - 1] - Tk) < np.abs(Tg[idx] - Tk), idx - 1, idx)
                            near = np.abs(Tg[idx] - Tk) <= D / 8
                            clean = near & (cg[idx] >= 0.8)
                            n_clean = int(clean.sum())
                            tot = np.hypot(sg[idx], d0_err)
                            excl = bool(np.any(clean & (dg[idx] < 0.5 * d0 - 3 * tot)))
                        if not excl:
                            res["periods"].append(dict(P=round(P, 5), n_support=len(mem), n_clean_flat_or_dip=n_clean))
            res["n_allowed"] = len(res["periods"])
        out.append(res)
    return out


# ------------------------------------------------------------------ stage 2b: BLS
def blocks(lc):
    t = lc["t"]
    cut = np.concatenate([[0], np.where(np.diff(t) > BLOCK_GAP)[0] + 1, [len(t)]])
    return [(cut[i], cut[i + 1]) for i in range(len(cut) - 1)]


def bls_block(t, f, seg, sign=1.0, n_sig=3):
    from astropy.timeseries import BoxLeastSquares
    from wotan import flatten
    ff = f if sign > 0 else 2.0 - f
    det = np.empty_like(ff)
    for u in np.unique(seg):
        m = seg == u
        if m.sum() < 10:
            det[m] = np.nan
            continue
        det[m] = flatten(t[m], ff[m], method="median", window_length=3.0)
    ok = np.isfinite(det)
    t, y = t[ok], det[ok]
    if len(t) > 1 and np.median(np.diff(t)) < 10 / 1440:  # bin to 10 min for speed
        b = np.floor((t - t[0]) / (10 / 1440)).astype(int)
        cnt = np.bincount(b)
        keep = cnt > 0
        t = (np.bincount(b, t) / np.maximum(cnt, 1))[keep]
        y = (np.bincount(b, y) / np.maximum(cnt, 1))[keep]
    span = t[-1] - t[0] if len(t) else 0
    pmax = min(BLS_PMAX, span / 2)
    sigs = []
    if len(t) < 200 or pmax <= P_MIN * 1.5:
        return sigs, span
    durs = np.array(BLS_DURS_H) / 24
    for _ in range(n_sig):
        bls = BoxLeastSquares(t, y)
        periods = bls.autoperiod(durs, minimum_period=P_MIN, maximum_period=pmax, frequency_factor=2.0)
        pw = bls.power(periods, durs, objective="snr")
        p = np.asarray(pw.power)
        sde = (p.max() - np.median(p)) / np.std(p)
        k = int(np.argmax(p))
        st = dict(P=float(pw.period[k]), T0=float(pw.transit_time[k]), D_h=float(pw.duration[k] * 24),
                  depth_ppm=float(pw.depth[k] * 1e6), sde=float(sde), n_trial=int(len(periods)))
        sigs.append(st)
        mask = bls.transit_mask(t, pw.period[k], min(3 * pw.duration[k], 0.5 * pw.period[k]), pw.transit_time[k])
        t, y = t[~mask], y[~mask]
        if len(t) < 200:
            break
    return sigs, span


def run_star(rep, thr, cm):
    epic = int(rep["EPIC"])
    lc = primary_lcs(rep)
    sig = []
    for e in rep["events"]:
        s_star = thr[e["mission"]]["S_star"]
        e["eligible"] = eligible(e)
        e["common_mode"] = (epic, round(e["T"], 5)) in cm
        e["significant"] = bool(e["eligible"] and e["snr"] >= s_star and not e["common_mode"])
        if e["significant"]:
            sig.append(e)
    out = dict(EPIC=epic, TIC=int(rep["TIC"]), n_events=len(rep["events"]), n_significant=len(sig),
               significant=[dict(T=round(e["T_fit"], 4), mission=e["mission"], sector=e["sector"], product=e["product"],
                                 depth_ppm=round(e["depth_fit"] * 1e6, 1), D_h=round(e["D_fit"] * 24, 2),
                                 snr=round(e["snr"], 1), found_in=e.get("found_in"), n_reductions=e.get("n_reductions"))
                            for e in sig],
               families=link(sig, lc))
    bl, blc = [], []
    if lc is not None:
        for a, b in blocks(lc):
            s, span = bls_block(lc["t"][a:b], lc["f"][a:b], lc["seg"][a:b], +1.0)
            c, _ = bls_block(lc["t"][a:b], lc["f"][a:b], lc["seg"][a:b], -1.0, n_sig=1)
            blk = dict(t0=round(float(lc["t"][a]), 2), t1=round(float(lc["t"][b - 1]), 2), span=round(span, 1))
            bl.append(dict(blk, signals=s))
            blc.append(dict(blk, signals=c))
    out["bls"], out["bls_control"] = bl, blc
    (OUTP / f"{epic}.json").write_text(json.dumps(out, indent=1, default=float))
    return out


def run_safe(rep, thr, cm):
    try:
        return run_star(rep, thr, cm)
    except Exception:  # noqa: BLE001
        return dict(EPIC=int(rep["EPIC"]), error=traceback.format_exc()[-1500:])


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--epics", nargs="*", type=int)
    args = a.parse_args()
    OUTP.mkdir(parents=True, exist_ok=True)
    reps = [json.loads(f.read_text()) for f in sorted(PER.glob("*.json"))]
    reps = [r for r in reps if r.get("status") == "ok"]
    thr = thresholds(reps)
    cm = common_mode(reps)
    (DATA / "thresholds.json").write_text(json.dumps(dict(thr, n_common_mode_events=len(cm)), indent=1))
    print(json.dumps({k: {kk: v[kk] for kk in ("S_star", "n_stars", "n_control_eligible")} for k, v in thr.items()}))
    todo = [r for r in reps if not args.epics or r["EPIC"] in args.epics]
    res = []
    with ProcessPoolExecutor(args.workers, max_tasks_per_child=1) as ex:
        futs = [ex.submit(run_safe, r, thr, cm) for r in todo]
        for k, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            res.append(r)
            print(f"[{k}/{len(todo)}] EPIC {r['EPIC']} sig {r.get('n_significant')} "
                  f"families {[(f['n_members'], f.get('n_allowed')) for f in r.get('families', [])]} "
                  f"{r.get('error', '')[-300:]}", flush=True)
    # BLS threshold from the inverted control, then flag signals
    ctl = [s["sde"] for r in res for b in r.get("bls_control", []) for s in b["signals"][:1]]
    sde_thr = float(np.percentile(ctl, 99)) if ctl else np.inf
    rows = []
    for r in res:
        if "error" in r:
            rows.append(dict(EPIC=r["EPIC"], error=True))
            continue
        f = OUTP / f"{r['EPIC']}.json"
        j = json.loads(f.read_text())
        j["bls_sde_threshold"] = sde_thr
        nb = 0
        for b in j["bls"]:
            for s in b["signals"]:
                s["detected"] = bool(s["sde"] >= sde_thr)
                nb += s["detected"]
        f.write_text(json.dumps(j, indent=1, default=float))
        best = max((p for fam in j["families"] for p in fam["periods"]), key=lambda p: p["n_support"], default=None)
        rows.append(dict(EPIC=j["EPIC"], TIC=j["TIC"], n_events=j["n_events"], n_significant=j["n_significant"],
                         n_families=len(j["families"]),
                         max_family=max((fam["n_members"] for fam in j["families"]), default=0),
                         linked_families=sum(1 for fam in j["families"] if fam["n_members"] >= 2 and fam.get("n_allowed")),
                         periodic_3plus=sum(1 for fam in j["families"] if fam["n_members"] >= 3 and fam.get("n_allowed")),
                         n_bls_detected=nb))
    pd.DataFrame(rows).to_csv(DATA / "periodic_summary.csv", index=False)
    t = json.loads((DATA / "thresholds.json").read_text())
    t["BLS"] = dict(sde_threshold=sde_thr, n_control_blocks=len(ctl))
    (DATA / "thresholds.json").write_text(json.dumps(t, indent=1))
    print(f"BLS SDE threshold (99th pct of inverted max): {sde_thr:.2f} from {len(ctl)} blocks")


if __name__ == "__main__":
    main()
