"""Period constraints from detected events + 'covered and flat' epochs.

Covered & flat (per product, per trial centre t0, from s4_search.scan):
  * >= 50 % of the expected in-event cadences (central 70 % of the duration) present,
  * out-of-event baseline on at least one side (polynomial order drops to 0/1 if one-sided),
  * measured depth < 0.5 * D_eff and the expected depth D_eff excluded at > 3 sigma
    (sigma includes a red-noise beta factor), with D_eff = catalogue depth for K2 and
    0.5 x catalogue depth for TESS (allowance for dilution in 21" pixels),
  * product is one of K2 PDCSAP / EVEREST / K2SFF, SPOC, TESS-SPOC, QLP, TGLC, eleanor,
    own TESScut photometry (K2 raw SAP, TARS, T16, CDIPS are not used for exclusion),
  * no product shows a detection (SNR > 5, depth > 0.3 D_eff) within +-D of that t0.
A trial period P is excluded if any T_ref + k P (k != 0) lands on such an epoch.
Excluded P-sets are built exactly by interval arithmetic: for every merged excluded time
interval [a, b] and every integer k, P in [(a - T)/k, (b - T)/k] is excluded."""
import sys, os, pickle, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import lib
from targets import T

SCAN = os.path.join(lib.ROOT, "data", "scan")


def merge(iv, tol=0.0):
    if len(iv) == 0: return np.zeros((0, 2))
    iv = np.asarray(iv, float); iv = iv[np.argsort(iv[:, 0])]
    out = [iv[0].copy()]
    for a, b in iv[1:]:
        if a <= out[-1][1] + tol:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append(np.array([a, b]))
    return np.array(out)


def time_masks(scan, dur, kinds_ok=None):
    """Excluded-epoch intervals and detection-protected intervals.
    For long events (> 0.5 d) the spline/roll-detrended products (QLP KSPSAP/DET, K2SFF,
    EVEREST) are not trusted to exclude, since they partially remove long dips.
    Protection uses detections with a two-sided baseline from non-raw products only."""
    step = dur / 8
    if kinds_ok is None:
        kinds_ok = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "qlpk", "tglc", "own", "eleanor"}
        if dur > 0.5:
            kinds_ok = {"K2pdc", "spoc", "tspoc", "qlp", "tglc", "own", "eleanor"}
    prot_kinds = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "qlpk", "tglc", "own", "eleanor"}
    ex, det = [], []
    for x in scan:
        t0 = x["t0"]
        if x["kind"] in prot_kinds:
            for i in np.where(x["det"] & (x["nb"] >= 3) & (x["na"] >= 3))[0]:
                det.append((t0[i] - dur, t0[i] + dur))
        if x["kind"] not in kinds_ok: continue
        for i in np.where(x["excl"])[0]:
            ex.append((t0[i] - step / 2, t0[i] + step / 2))
    ex = merge(ex, tol=1e-6); det = merge(det)
    # subtract detection-protected intervals from excluded intervals
    out = []
    for a, b in ex:
        segs = [(a, b)]
        for c, d in det:
            new = []
            for (u, v) in segs:
                if d <= u or c >= v: new.append((u, v)); continue
                if c > u: new.append((u, c))
                if d < v: new.append((d, v))
            segs = new
        out += segs
    return merge(out, tol=1e-6), det


def excluded_P(ex, Tref, Pmin=0.3, Pmax=8000.0):
    ivs = []
    for a, b in ex:
        if a <= Tref <= b:
            continue  # the reference epoch itself (should not happen)
        if a > Tref:
            lo, hi = a - Tref, b - Tref
        else:
            lo, hi = Tref - b, Tref - a
        kmax = int(np.floor(hi / Pmin))
        if kmax < 1: continue
        k = np.arange(1, kmax + 1)
        pl, ph = lo / k, hi / k
        m = (ph >= Pmin) & (pl <= Pmax)
        ivs.append(np.c_[pl[m], ph[m]])
    if not ivs: return np.zeros((0, 2))
    return merge(np.vstack(ivs))


def allowed_from_excluded(exP, Pmin=0.3, Pmax=8000.0):
    al = []; cur = Pmin
    for a, b in exP:
        if a > cur: al.append((cur, min(a, Pmax)))
        cur = max(cur, b)
        if cur >= Pmax: break
    if cur < Pmax: al.append((cur, Pmax))
    return np.array(al)


def is_excluded(P, exP):
    i = np.searchsorted(exP[:, 0], P, side="right") - 1
    return i >= 0 and exP[i, 1] >= P


def load(e):
    return pickle.load(open(os.path.join(SCAN, f"{e}.pkl"), "rb"))


def analyse(e, events, Pmin=0.3, Pmax=8000.0, verbose=True):
    """events: list of (t0, err) of confirmed events of the SAME kind (e.g. all primaries)."""
    P = load(e); r = P["row"]; dur = r["dur"] / 24
    ex, det = time_masks(P["scan"], dur)
    Tref = events[0][0]
    exP = excluded_P(ex, Tref, Pmin, Pmax)
    out = dict(epic=e, n_events=len(events))
    if len(events) == 1:
        al = allowed_from_excluded(exP, Pmin, Pmax)
        out["allowed"] = al
        if verbose:
            print(f"EPIC {e}: single event at {Tref:.4f}. Allowed P windows (first 15):")
            for a, b in al[:15]:
                print(f"    {a:9.3f} - {b:9.3f}")
            print(f"    ... total {len(al)} windows; smallest allowed P = {al[0][0] if len(al) else np.nan:.3f}")
        return out
    # multiple events: candidate periods = gcd-like set
    ts = np.array([t for t, _ in events]); es = np.array([s for _, s in events])
    o = np.argsort(ts); ts, es = ts[o], es[o]
    span = ts[-1] - ts[0]
    cands = []
    nmax = int(span / Pmin)
    for n in range(1, min(nmax, 200000) + 1):
        Pn = span / n
        if Pn < Pmin: break
        # all events must fall on integer cycles
        cyc = (ts - ts[0]) / Pn
        resid = (cyc - np.round(cyc)) * Pn
        tol = 3 * np.sqrt(es ** 2 + es[0] ** 2) + 0.02 * dur * 0  # timing tolerance
        if np.all(np.abs(resid) <= np.maximum(tol, 0.25 * dur)):
            # refine P by linear fit
            N = np.round(cyc)
            if len(ts) > 2:
                A = np.c_[np.ones_like(N), N]
                w = 1 / es ** 2
                cov = np.linalg.inv(A.T @ (A * w[:, None]))
                sol = cov @ (A.T @ (w * ts))
                Pfit, ePfit = sol[1], np.sqrt(cov[1, 1])
                chi2 = np.sum(w * (ts - A @ sol) ** 2)
            else:
                Pfit, ePfit, chi2 = Pn, np.sqrt(es[0] ** 2 + es[-1] ** 2) / n, 0.0
            exd = is_excluded(Pfit, excluded_P(ex, ts[0], Pfit * 0.999, Pfit * 1.001)) if True else False
            cands.append(dict(n=n, P=Pfit, eP=ePfit, chi2=chi2, excluded=bool(exd)))
    out["cands"] = cands
    if verbose:
        ok = [c for c in cands if not c["excluded"]]
        print(f"EPIC {e}: {len(ts)} events, span {span:.4f} d; {len(cands)} commensurate periods >= {Pmin}, {len(ok)} not excluded:")
        for c in ok[:40]:
            print(f"    n={c['n']:5d}  P={c['P']:.5f} +- {c['eP']:.5f}  chi2={c['chi2']:.1f}")
    return out


if __name__ == "__main__":
    pass


def why_excluded(e, T0, P, show=True):
    """Return the epochs (and products) that exclude period P (using the same rules)."""
    D = load(e); r = D["row"]; dur = r["dur"] / 24
    ex, det = time_masks(D["scan"], dur)
    tmin = min(s["t"][0] for s in D["segs"]); tmax = max(s["t"][-1] for s in D["segs"])
    ks = np.arange(int(np.floor((tmin - T0) / P)) - 1, int(np.ceil((tmax - T0) / P)) + 2)
    hits = []
    for k in ks:
        tc = T0 + k * P
        i = np.searchsorted(ex[:, 0], tc, side="right") - 1
        if i >= 0 and ex[i, 1] >= tc:
            # which products flag it
            prods = []
            for x in D["scan"]:
                j = np.argmin(np.abs(x["t0"] - tc)) if len(x["t0"]) else None
                if j is not None and abs(x["t0"][j] - tc) <= dur / 8 and x["excl"][j]:
                    prods.append(f"{x['src']}({x['depth'][j]*1e6:.0f}+-{x['err'][j]*1e6:.0f})")
            hits.append((k, tc, prods))
    if show:
        print(f"  P={P:.5f}: {len(hits)} excluding epochs" + "".join(f"\n     E={k} t={tc:.3f}: " + "; ".join(p[:3]) for k, tc, p in hits[:4]))
    return hits
