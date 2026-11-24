"""Independent re-implementation of the duotransit alias table (driver check).

Written from the text of
  catalogue/prereg_catalogue.yaml  (what_is_computed.step_4_alias_table.hawthorn_duotransits)
  calibration/mono/round3/FROZEN_RULE3.md  ("Alias table (T1)")
  calibration/mono/round2/FROZEN_RULE.md   (definitions)
only. It uses the shared loading / detrending / reference / event-search functions
(monocal.load_star, monocal.prepare, st2core helpers) and nothing of catrun, rule2, rule3 or
st2core.alias_table. Choices made where the text left one are in AMBIGUITIES.md.

Usage: python ind_duo.py [n_processes<=3]
"""
import math
import socket
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

socket.setdefaulttimeout(120)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "calibration" / "mono"))
sys.path.insert(0, str(ROOT / "salinas-second-transits"))
import monocal as mc  # noqa: E402
import st2core as core  # noqa: E402

CACHE = Path("/Users/mcbee/Documents/GitHub/coracle/duo-transits/salinas-second-transits/work/lc_cache")
DEFS = HERE.parent / "duo_check_defs.csv"

SIG_T_MIN = 15.0                                    # error of each transit time, minutes
SIG_DT_D = math.hypot(SIG_T_MIN, SIG_T_MIN) / 1440  # days
P_FLOOR = 1.0
N_MAX = 3000
H_CAP = 16
CLEAN_FRAC, CLEAN_EDGE = 0.8, 0.5
COVERED_FRAC = 0.5


def boxes(T, F, c, dur, sig, edges, depth0, derr):
    """depth, frac, err, thr, clean, covered of the boxes at centres c."""
    c = np.atleast_1d(np.asarray(c, float))
    depth, frac, _ = core.box_stats(T, F, c, dur)
    err = sig / np.sqrt(np.maximum(frac, 0.001))
    thr = 0.5 * depth0 - 3.0 * np.hypot(err, derr)
    clean = (frac >= CLEAN_FRAC) & (core.edge_distance(edges, c) >= CLEAN_EDGE)
    covered = frac >= COVERED_FRAC
    return depth, frac, err, thr, clean, covered


def alias_rows(T, F, t_a, t_b, dur, depth0, derr, measured, w0, pmin):
    t1, t2 = sorted([float(t_a), float(t_b)])
    dT = t2 - t1
    sig = core.duration_scale_noise(T, F, dur)
    edges = core.gap_edges(T)
    tmin, tmax = float(T.min()), float(T.max())
    step = dur / 8.0
    out = []
    n = 1
    while n <= N_MAX and dT / n >= P_FLOOR:
        P = dT / n
        k = np.arange(int(math.floor((tmin - t1) / P)) - 1, int(math.ceil((tmax - t1) / P)) + 2)
        k = k[(k != 0) & (k != n)]
        ep = t1 + k * P
        far = np.minimum(np.abs(k), np.abs(k - n))
        h = w0 + np.ceil(2.0 * far * (SIG_DT_D / n) / step - 1e-9).astype(int)
        constrained = h <= H_CAP
        depth, frac, err, thr, clean, covered = boxes(T, F, ep, dur, sig, edges, depth0, derr)
        with np.errstate(invalid="ignore"):
            cand = constrained & clean & (depth >= -3.0 * err) & (depth < thr) & bool(measured)
            supp = constrained & clean & (np.abs(depth - depth0) < 3.0 * np.hypot(err, derr))
        excluded = False
        if cand.any():
            ci = np.where(cand)[0]
            j = np.arange(-H_CAP, H_CAP + 1)
            cc = ep[ci][:, None] + j[None, :] * step
            inwin = np.abs(j)[None, :] <= h[ci][:, None]
            d2, _, _, thr2, _, cov2 = boxes(T, F, cc.ravel(), dur, sig, edges, depth0, derr)
            with np.errstate(invalid="ignore"):
                veto = (cov2 & (d2 >= thr2)).reshape(cc.shape) & inwin
            excluded = bool((~veto.any(axis=1)).any())
        if excluded:
            status = "excluded"
        elif supp.any():
            status = "supported"
        elif np.isfinite(pmin) and P < pmin:
            status = "allowed_below_pmin"
        else:
            status = "allowed"
        out.append((n, P, status))
        n += 1
    return out


def second_depth_error(prep, t2_cat):
    """derr of the second transit: pipeline event, else refit at t2_cat, else the reference's."""
    dur = prep["dur_ref"]
    near = [e for e in prep["events"]
            if "error" not in e and e.get("snr") is not None and abs(e["tc"] - t2_cat) <= dur / 2]
    if near:
        ev = max(near, key=lambda e: e["snr"])
        return float(ev["depth_err"]), "event"
    for s in sorted(prep["det"]):
        t, f = prep["det"][s]
        if t.min() <= t2_cat <= t.max():
            r = core.refine_event(t, f, t2_cat, dur)
            if r and np.isfinite(r["depth_err"]) and r["depth_err"] > 0:
                return float(r["depth_err"]), "refit"
            break
    return float(prep["derr_ref"]), "reference"


def do_row(row):
    rid = row["row_id"]
    last = None
    for attempt in range(3):
        try:
            sectors, errors = mc.load_star(int(row["tic"]), CACHE)
            if not sectors:
                raise RuntimeError(f"no sectors loaded ({errors[:3]})")
            break
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            time.sleep(5)
    else:
        return dict(row_id=rid, status=f"load failed: {last}"), []
    try:
        prep = mc.prepare(sectors, float(row["t_cat"]), float(row["dur_cat_d"]), float(row["depth_cat"]),
                          float(row["depth_cat_err"]))
        if prep is None:
            return dict(row_id=rid, status="reference not in data", load_errors="; ".join(errors)), []
        T, F = prep["T"], prep["F"]
        dur, depth0 = prep["dur_ref"], prep["depth_ref"]
        measured = bool(prep["out"]["ref_measured"])
        derr, src = second_depth_error(prep, float(row["t2_cat"]))
        multi = str(row["multi"]).strip().lower() in ("true", "1")
        w0 = 8 if multi else 4
        mstar = float(row["mstar"]) if pd.notna(row["mstar"]) else np.nan
        rstar = float(row["rstar"]) if pd.notna(row["rstar"]) else np.nan
        pmin = core.pmin_from_duration(dur, mstar, rstar, depth0)
        tab = alias_rows(T, F, float(row["t_cat"]), float(row["t2_cat"]), dur, depth0, derr, measured, w0, pmin)
        info = dict(row_id=rid, ref_measured=measured, ref_snr=prep["out"]["ref_snr"], dur_d=dur,
                    depth0=depth0, derr_used=derr, derr_from=src, n_alias=len(tab),
                    n_excluded=sum(s == "excluded" for _, _, s in tab), status="ok",
                    n_sectors=len(sectors), load_errors="; ".join(errors), pmin_d=pmin, w0=w0)
        return info, [dict(row_id=rid, n=n, period_d=P, status=s) for n, P, s in tab]
    except Exception:
        return dict(row_id=rid, status="error: " + traceback.format_exc()[-400:]), []


def main():
    nproc = min(3, int(sys.argv[1])) if len(sys.argv) > 1 else 3
    defs = pd.read_csv(DEFS)
    rows = [r.to_dict() for _, r in defs.iterrows()]
    if nproc > 1:
        from multiprocessing import Pool
        with Pool(nproc) as p:
            res = p.map(do_row, rows, chunksize=1)
    else:
        res = [do_row(r) for r in rows]
    infos = pd.DataFrame([a for a, _ in res])
    tabs = pd.DataFrame([x for _, b in res for x in b], columns=["row_id", "n", "period_d", "status"])
    cols = ["row_id", "ref_measured", "ref_snr", "dur_d", "depth0", "derr_used", "derr_from", "n_alias",
            "n_excluded"]
    infos = infos[cols + [c for c in infos.columns if c not in cols]]
    infos.to_csv(HERE / "ind_duo_rows.csv", index=False)
    tabs.to_csv(HERE / "ind_duo_tables.csv", index=False)
    print(infos.to_string())
    print(tabs.groupby(["row_id", "status"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
