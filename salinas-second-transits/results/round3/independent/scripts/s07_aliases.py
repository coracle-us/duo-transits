"""Alias scan: which periods are consistent with the detected transits AND the rest of the data?

Candidates: P_N = (T_last - T_first) / N for N = 1..Nmax (P_N >= Pmin), where T_first/T_last
are the earliest/latest detected transit mid-times. A candidate is "event-consistent" if every
detected transit falls within `tol` of an integer epoch.  For event-consistent candidates,
every other predicted epoch that lands in data is tested: the transit depth scale s (shape
fixed from s05, Tc fixed at the prediction) is fitted on the *raw* normalised flux with a
quadratic baseline in a +-1.5 duration window, so the test does not depend on detrending.
Expected s = 1. A predicted epoch with good coverage (in-transit >= 0.5) and s inconsistent with
1 at > 5 sigma rules the candidate out.

Usage: s07_aliases.py TIC  (reads out/times_<tic>.csv for detected events: the timing-source row of each epoch;
and out/shape_<tic>.json)
"""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import TARGETS, OUT, load_lc_files, load_tesscut_phot, batman_model

RANK = {"TESS-SPOC": 0, "QLP": 1, "TESScut": 2}


def shape_model(sh, tc, t, cad):
    inc = np.degrees(np.arccos(np.clip(sh["b"] / sh["a"], 0, 1)))
    return batman_model(t, tc, sh["P0"], sh["rp"], sh["a"], inc, sh["u"], exptime=cad)


def depth_at(lcs, sh, tc, dur):
    """Best-covered measurement of depth scale at tc. Returns (s, err, cov, sector, source)."""
    best = None
    for lc in lcs:
        t = lc["t"]
        if t.min() > tc + dur or t.max() < tc - dur:
            continue
        w = np.abs(t - tc) < 1.5 * dur
        cad = lc["exptime"]
        inn = np.abs(t - tc) < 0.5 * dur
        cov = inn.sum() * cad / dur
        oot = w & ~(np.abs(t - tc) < 0.6 * dur)
        if cov < 0.3 or oot.sum() < 10:
            continue
        tt, ff = t[w], lc["f"][w]
        m0 = shape_model(sh, tc, tt, cad)
        dt = tt - tc
        A = np.vstack([np.ones_like(dt), dt, dt**2, -(1 - m0)]).T
        x, *_ = np.linalg.lstsq(A, ff, rcond=None)
        res = ff - A @ x
        sig = np.std(res)
        C = np.linalg.inv(A.T @ A) * sig**2
        # beta for red noise at ~dur/4
        nb = max(2, int(round(dur / 4 / cad)))
        k = (len(res) // nb) * nb
        beta = 1.0
        if k >= 3 * nb:
            beta = max(1.0, np.std(res[:k].reshape(-1, nb).mean(1)) / (sig / np.sqrt(nb)))
        err = np.sqrt(C[3, 3]) * beta
        cand = (x[3], err, cov, lc["sector"], lc["source"])
        key = (-min(cov, 1.0), RANK[lc["source"]])
        if best is None or key < best[0]:
            best = (key, cand)
    return None if best is None else best[1]


def run(tic, Pmin=1.0, tol_min=45.0, verbose=True):
    tgt = TARGETS[tic]
    dur = tgt["dur"]
    sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
    # detected events = the timing-source mid-times (one per epoch, chosen in s05)
    rows = [r for r in csv.DictReader(open(os.path.join(OUT, f"times_{tic}.csv")))
            if r["tc"] and r["timing_source"] == "1"]
    ev = {}
    for r in rows:  # one time per epoch, prefer timing source / best coverage
        e = int(r["epoch"])
        key = (r["timing_source"] != "1", -float(r["cov_in"]))
        if e not in ev or key < ev[e][0]:
            ev[e] = (key, float(r["tc"]))
    times = np.array(sorted(v[1] for v in ev.values()))
    lcs = [lc for lc in load_lc_files(tic) + load_tesscut_phot(tic) if len(lc["t"]) > 100]
    tmin = min(lc["t"].min() for lc in lcs); tmax = max(lc["t"].max() for lc in lcs)
    span = times[-1] - times[0]
    Nmax = int(span / Pmin)
    tol = tol_min / 1440
    cache = {}
    results = []
    for N in range(1, Nmax + 1):
        P = span / N
        ph = (times - times[0]) / P
        if np.max(np.abs(ph - np.round(ph))) * P > tol:
            continue
        # refine with least squares on the events
        E = np.round(ph)
        A = np.vstack([np.ones_like(E), E]).T
        T0, P = np.linalg.lstsq(A, times, rcond=None)[0]
        n = np.arange(np.ceil((tmin - T0) / P), np.floor((tmax - T0) / P) + 1)
        tests = []
        for k in n:
            tc = T0 + k * P
            if np.min(np.abs(times - tc)) < 0.5 * dur:
                continue  # a detected event
            key = round(tc, 3)
            if key not in cache:
                cache[key] = depth_at(lcs, sh, tc, dur)
            m = cache[key]
            if m is not None:
                tests.append((int(k), tc, *m))
        z = [(s - 1) / e for (_, _, s, e, cov, _, _) in tests if cov >= 0.5]
        worst = min(z) if z else np.nan
        ruled = bool(bool(z) and worst < -5)
        results.append(dict(N=N, P=float(P), T0=float(T0), ntest=len(z), worst_z=float(worst), ruled_out=ruled, tests=tests))
    if verbose:
        print(f"TIC {tic}: {len(times)} detected events, span {span:.3f} d; "
              f"{len(results)} event-consistent periods >= {Pmin} d (tol {tol_min} min)")
        for r in results:
            print(f"  N={r['N']:5d} P={r['P']:11.5f}  extra covered epochs={r['ntest']:3d}  "
                  f"worst z={r['worst_z']:+7.1f}  {'RULED OUT' if r['ruled_out'] else 'ALLOWED'}")
    return results, times


def detail(tic, P, T0):
    """Print depth tests at every predicted epoch in data for a given ephemeris."""
    tgt = TARGETS[tic]
    dur = tgt["dur"]
    sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
    lcs = [lc for lc in load_lc_files(tic) + load_tesscut_phot(tic) if len(lc["t"]) > 100]
    tmin = min(lc["t"].min() for lc in lcs); tmax = max(lc["t"].max() for lc in lcs)
    print(f"Predicted epochs for P={P:.6f} T0={T0:.5f}: depth scale s (1 = transit, 0 = flat)")
    for k in np.arange(np.ceil((tmin - T0) / P), np.floor((tmax - T0) / P) + 1):
        tc = T0 + k * P
        m = depth_at(lcs, sh, tc, dur)
        if m is not None:
            s, e, cov, sec, src = m
            print(f"   E={int(k):5d} Tc={tc:10.4f} S{sec:<3d} {src:<9} cov={cov:.2f}  s={s:5.2f} ± {e:.2f}")


if __name__ == "__main__":
    tic = int(sys.argv[1])
    res, times = run(tic)
    json.dump([{k: v for k, v in r.items() if k != "tests"} for r in res],
              open(os.path.join(OUT, f"aliases_{tic}.json"), "w"), indent=1)
    allowed = [r for r in res if not r["ruled_out"]]
    for r in allowed:
        detail(tic, r["P"], r["T0"])
    # the strongest ruled-out shorter aliases: show evidence
    for r in res:
        if r["ruled_out"] and r["N"] <= 4 * (allowed[0]["N"] if allowed else 1):
            print(f"\nEvidence against P={r['P']:.5f} (N={r['N']}):")
            for (k, tc, s, e, cov, sec, src) in r["tests"]:
                if cov >= 0.5 and (s - 1) / e < -5:
                    print(f"   E={k} Tc={tc:.4f} S{sec} {src} cov={cov:.2f} s={s:.2f}±{e:.2f}")
