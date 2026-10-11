"""
Single-transit project, step 3: per-target analysis of a single-transit candidate.

For each target (row of singles/data/pilot_targets.csv or a validation list):
  1. fetch TESS light curves (reuses scripts/fetch_data.py; shared store data/raw/lightcurves)
  2. confirm the catalogued event in our light curves (scripts/alias_pipeline.test_window)
  3. repeat search: test a transit window of the catalogued depth/duration on a grid of trial
     times (step D/4) over ALL data; robust detections (SNR >= 8, depth 0.5-2x, coverage >= 0.8,
     baseline both sides; same cuts as the validated alias pipeline) are repeat candidates
  4. period constraints: for trial periods 5-1500 d (log grid), predicted epochs that fall in a
     FLAT window (pipeline 'excluded' verdict) rule the period out; report the excluded fraction
     of log-period space and the allowed ranges
  5. difference image at the event (and at any repeat) with scripts/diff_image.py, source offset
     measured from the target's out-of-transit centroid (validated ~1-3 arcsec floor)

Usage:  python singles/scripts/pilot.py [targets.csv] [out_subdir] [--no-diff]
        (defaults: singles/data/pilot_targets.csv, singles/results/pilot)
Output: <out>/<TIC>.json, <out>/<TIC>_windows.png, <out>/<TIC>_diffimg_s*.png, <out>/summary.csv
"""
import json
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402
import diff_image as di  # noqa: E402
import fetch_data as fd  # noqa: E402

DETECT_SNR, DMIN, DMAX = 8.0, 0.5, 2.0
P_GRID = np.exp(np.linspace(np.log(5), np.log(1500), 6000))
EPOCH_ERR = 0.01  # d, used to scan the catalogued event time


def robust(r, delta):
    if r.get("verdict") != "supported":
        return False
    side_min = ap.DETECT_SIDE * ap.BASE_W / r["cadence"]
    return (r["snr"] >= DETECT_SNR and DMIN * delta <= r["depth"] <= DMAX * delta
            and r["cover"] >= ap.DETECT_COVER and r["n_oot_left"] >= side_min
            and r["n_oot_right"] >= side_min)


def grid_windows(lc, D, delta, epoch):
    t, f, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
    step = D / 4
    T_all = np.arange(t.min(), t.max(), step)
    j = np.searchsorted(t, T_all)
    j0, j1 = np.clip(j - 1, 0, len(t) - 1), np.clip(j, 0, len(t) - 1)
    near = np.minimum(np.abs(t[j0] - T_all), np.abs(t[j1] - T_all)) < D / 2
    rows = []
    for T in T_all[near]:
        r = ap.test_window(t, f, T, 0.0, D, delta, sg)
        r["sector"] = int(lc.sector.iloc[min(np.searchsorted(t, T), len(t) - 1)])
        r["product"] = lc["product"].iloc[min(np.searchsorted(t, T), len(t) - 1)]
        rows.append(r)
    g = pd.DataFrame(rows)
    g["near_event"] = (g.T_pred - epoch).abs() < D
    return g


def period_constraints(g, epoch, D):
    """Excluded/allowed trial periods from flat windows on the grid."""
    Tg = g.T_pred.to_numpy()
    flat = (g.verdict == "excluded").to_numpy()
    conclusive = g.verdict.isin(["excluded", "supported"]).to_numpy()
    tmin, tmax = Tg.min(), Tg.max()
    status = []
    for P in P_GRID:
        ks = np.arange(np.ceil((tmin - epoch) / P), np.floor((tmax - epoch) / P) + 1)
        ks = ks[ks != 0]
        ex = False
        for k in ks:
            i = np.searchsorted(Tg, epoch + k * P)
            for ii in (i - 1, i):
                if 0 <= ii < len(Tg) and abs(Tg[ii] - (epoch + k * P)) < D / 8 and flat[ii]:
                    ex = True
        status.append(ex)
    status = np.array(status)
    allowed = P_GRID[~status]
    # summarise allowed ranges
    ranges, start = [], None
    for i, ok in enumerate(~status):
        if ok and start is None:
            start = P_GRID[i]
        if (not ok or i == len(P_GRID) - 1) and start is not None:
            end = P_GRID[i - 1] if not ok else P_GRID[i]
            ranges.append((round(start, 2), round(end, 2)))
            start = None
    return dict(frac_logP_excluded=float(status.mean()),
                min_allowed_P=float(allowed.min()) if len(allowed) else None,
                n_allowed_ranges=len(ranges),
                allowed_ranges=ranges[:15] + ([("...", len(ranges))] if len(ranges) > 15 else []),
                frac_conclusive_windows=float(conclusive.mean()))


def run_target(row, out, do_diff=True):
    tic = int(row.TIC)
    D, delta, epoch = row.dur_h / 24, row.depth_ppm * 1e-6, float(row.epoch_btjd)
    if not (fd.LC_DIR / str(tic) / "search_results.csv").exists():
        fd._fetch_one(0, tic)
    lc = ap.load_lightcurves(tic)
    rep = dict(TIC=tic, CTOI=row.CTOI, Tmag=row.Tmag, depth_ppm=row.depth_ppm, dur_h=row.dur_h,
               epoch_btjd=epoch, n_sectors=int(lc.sector.nunique()) if len(lc) else 0)
    if len(lc) == 0:
        rep["status"] = "no light curves"
        return rep
    t, f, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
    ev = ap.test_window(t, f, epoch, EPOCH_ERR, D, delta, sg)
    rep["event_verdict"] = ev.get("verdict")
    rep["event_depth_ppm"] = round(ev.get("depth", np.nan) * 1e6, 1) if "depth" in ev else None
    rep["event_snr"] = round(ev.get("snr", np.nan), 1) if "snr" in ev else None
    ev_T = epoch + (ev.get("offset", 0.0) or 0.0)
    g = grid_windows(lc, D, delta, ev_T)
    reps = []
    for _, r in g[~g.near_event].sort_values("snr", ascending=False).iterrows():
        if robust(r, delta) and all(abs(r.T_pred - x["T"]) > D for x in reps):
            reps.append(dict(T=round(r.T_pred, 4), sector=int(r.sector), product=r["product"],
                             depth_ppm=round(r.depth * 1e6), snr=round(r.snr, 1)))
    rep["repeats"] = sorted(reps, key=lambda x: x["T"])
    rep["n_repeats"] = len(reps)
    rep.update(period_constraints(g, ev_T, D))
    cands = periods_from_repeats(g, ev_T, rep["repeats"], D)
    rep["candidate_periods"] = cands
    rep["n_candidate_periods"] = len(cands)
    # difference images: event sector + each repeat sector
    if not do_diff:
        plot_windows(tic, lc, ev_T, reps, D, delta, out)
        return rep
    try:
        gaia = di.gaia_field(row.RA, row.Dec)
        tgt = gaia.iloc[0]
        diffs = []
        ev_secs = sorted(set(lc.sector[np.abs(lc.time - ev_T) < 1]))
        targets = [(ev_T, s, "event") for s in ev_secs] + [(x["T"], x["sector"], "repeat") for x in reps[:3]]
        for T, s, lab in targets:
            r = di.analyse(tic, 1e6, T, D, int(s), tgt.ra, tgt.dec, gaia, out_dir=out)
            if r and r.get("n_transits"):
                diffs.append(dict(kind=lab, sector=int(s), offset_arcsec=round(r["offset_from_target_arcsec"], 1),
                                  offset_err_arcsec=round(r["offset_err_px"] * 21, 1),
                                  offset_sigma=round(r["offset_sigma"], 1), snr_peak=round(r["snr_peak"], 1),
                                  nearest_gaia_G=round(r["nearest_gaia_to_source"]["G"], 2),
                                  nearest_gaia_sep=round(r["nearest_gaia_to_source"]["sep_from_target_arcsec"], 1)))
            else:
                diffs.append(dict(kind=lab, sector=int(s), note="no usable difference image"))
        rep["diff_images"] = diffs
    except Exception as e:  # noqa: BLE001
        rep["diff_images"] = [dict(error=str(e)[:200])]
    plot_windows(tic, lc, ev_T, reps, D, delta, out)
    return rep


def plot_windows(tic, lc, ev_T, reps, D, delta, out):
    wins = [("event", ev_T)] + [(f"repeat s{x['sector']}", x["T"]) for x in reps[:5]]
    fig, axs = plt.subplots(1, len(wins), figsize=(3.6 * len(wins), 3), squeeze=False)
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    for ax, (lab, T) in zip(axs[0], wins):
        m = np.abs(t - T) < max(1.5 * D, 0.6)
        if m.sum():
            tt, ff = t[m], f[m]
            oot = np.abs(tt - T) > D / 2 + 0.05
            if oot.sum() > 3:
                ff = ff / np.polyval(np.polyfit(tt[oot] - T, ff[oot], 1), tt - T)
            ax.plot((tt - T) * 24, (ff - 1) * 1e3, ".", ms=2, color="0.5")
        ax.plot([-D * 12, -D * 12, D * 12, D * 12], [0, -delta * 1e3, -delta * 1e3, 0], "r-", lw=1)
        ax.set_title(f"TIC {tic} {lab} T={T:.2f}", fontsize=8)
        ax.set_ylim(-3 * delta * 1e3, 2 * delta * 1e3)
    fig.tight_layout()
    fig.savefig(out / f"{tic}_windows.png", dpi=80)
    plt.close(fig)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sys.argv = [sys.argv[0]] + args + [a for a in sys.argv[1:] if a.startswith("--")]
    targets = Path(args[0]) if len(args) > 0 else ROOT / "singles" / "data" / "pilot_targets.csv"
    out = ROOT / "singles" / "results" / (args[1] if len(args) > 1 else "pilot")
    out.mkdir(parents=True, exist_ok=True)
    do_diff = "--no-diff" not in sys.argv
    rows = []
    for _, row in pd.read_csv(targets).iterrows():
        try:
            rep = run_target(row, out, do_diff)
        except Exception:  # noqa: BLE001
            rep = dict(TIC=int(row.TIC), status="error", error=traceback.format_exc()[-400:])
        (out / f"{int(row.TIC)}.json").write_text(json.dumps(rep, indent=1, default=float))
        rows.append({k: (json.dumps(v, default=float) if isinstance(v, (list, dict)) else v)
                     for k, v in rep.items()})
        print(f"TIC {rep['TIC']}: event {rep.get('event_verdict')} snr {rep.get('event_snr')} | "
              f"repeats {rep.get('n_repeats')} | periods {rep.get('candidate_periods', [])[:6]} | "
              f"logP excluded {rep.get('frac_logP_excluded')} | "
              f"diff {rep.get('diff_images')}", flush=True)
    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)



def excluded_by_flat(g, epoch, P, D):
    """True if any predicted epoch (k != 0) of period P falls in a flat (excluded) window."""
    Tg = g.T_pred.to_numpy()
    flat = (g.verdict == "excluded").to_numpy()
    ks = np.arange(np.ceil((Tg.min() - epoch) / P), np.floor((Tg.max() - epoch) / P) + 1)
    for k in ks[ks != 0]:
        i = np.searchsorted(Tg, epoch + k * P)
        for ii in (i - 1, i):
            if 0 <= ii < len(Tg) and abs(Tg[ii] - (epoch + k * P)) < D / 8 and flat[ii]:
                return True
    return False


def periods_from_repeats(g, epoch, reps, D, p_min=5.0):
    """Candidate periods implied by detected repeats. Trial periods P = (T_rep - epoch)/n from
    the repeat closest to the event; for each, assign cycle numbers to the event and all
    repeats, fit a linear ephemeris (absorbing the trial period's error, since repeat times
    are only known to the grid step D/4), and keep it if every residual is < D/3, the cycle
    numbers are distinct, and no predicted epoch falls in a flat window. Returns [(n, P_fit)]."""
    if not reps:
        return []
    ts = np.array(sorted([x["T"] for x in reps], key=lambda T: abs(T - epoch)))
    times = np.concatenate([[epoch], ts])
    dT = abs(ts[0] - epoch)
    out = []
    for n in range(1, int(dT / p_min) + 1):
        P0 = dT / n
        k = np.round((times - epoch) / P0)
        if len(np.unique(k)) < len(k):
            continue
        A = np.column_stack([np.ones_like(k), k])
        (T0, P), *_ = np.linalg.lstsq(A, times, rcond=None)
        if np.max(np.abs(times - (T0 + k * P))) >= D / 3:
            continue
        if not excluded_by_flat(g, T0, P, D):
            out.append((n, round(float(P), 5)))
    return out

if __name__ == "__main__":
    main()
