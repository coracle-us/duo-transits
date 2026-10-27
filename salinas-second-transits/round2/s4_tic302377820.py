#!/usr/bin/env python3
"""Round 2, section 4: TIC 302377820 (candidates a and b) and CTOI 302377820.01.

1. CTOI ephemeris residuals at candidate a, candidate b, the S68 dip and the S95 dips.
2. Alias sets for a and b rebuilt from the S68 transit (the brief's choice;
   the clean-event rule picks the S95 dip at 3890.2, which has higher SNR).
3. Fold of all sectors on the CTOI period, with and without the single events.
4. Plots of each dip with the CTOI predicted epochs marked.

Run after the search and alias steps (uses work/lc/302377820_s*.npz).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2core as core  # noqa: E402

TIC = 302377820
WORK, OUT = HERE / "work", HERE / "results" / "round2"


def ctoi():
    c = pd.read_csv(WORK / "exofop_ctoi.csv", comment="#", low_memory=False)
    r = c[c["TIC ID"] == TIC].iloc[0]
    return dict(ctoi=r["CTOI"], E=r["Transit Epoch (BJD)"] - 2457000, sE=r["Transit Epoch (BJD) err"],
                P=r["Period (days)"], sP=r["Period (days) Error"], depth=r["Depth ppm"] * 1e-6,
                dur=r["Duration (hrs)"] / 24, user=r["User"], notes=r["Notes"])


def load_all():
    T, F, S = [], [], []
    for p in sorted((WORK / "lc").glob(f"{TIC}_s*.npz")):
        if "tesscut" in p.name:
            continue
        d = np.load(p)
        T.append(d["t"]); F.append(d["f"]); S.append(np.full(d["t"].size, int(p.stem.split("_s")[1])))
    return np.concatenate(T), np.concatenate(F), np.concatenate(S)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c = ctoi()
    cands = pd.read_csv(HERE / "candidates.csv")
    cands = cands[cands.tic == TIC].set_index("planet")
    ev = pd.read_csv(WORK / "events_all.csv")
    ev = ev[(ev.tic == TIC) & (ev.passes == True)]
    t, f, sec = load_all()

    # 1. ephemeris residuals
    marks = [("candidate a (paper)", cands.loc["a", "t0"]), ("candidate b (paper)", cands.loc["b", "t0"])]
    for t2, g in ev.groupby(ev.t2.round(1)):
        marks.append((f"S{int(g.sector.iloc[0])} {g.author.iloc[0]} dip", float(g.t2.mean())))
    rows = []
    for lab, tm in marks:
        n = (tm - c["E"]) / c["P"]
        k = int(np.round(n))
        pred = c["E"] + k * c["P"]
        sig = np.hypot(c["sE"], k * c["sP"])
        rows.append(dict(event=lab, t_btjd=tm, ctoi_cycle=k, ctoi_pred_btjd=pred,
                         resid_h=(tm - pred) * 24, ctoi_sigma_h=sig * 24, resid_sigma=(tm - pred) / sig,
                         within_half_ctoi_duration=abs(tm - pred) < c["dur"] / 2))
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tic302377820_ctoi_residuals.csv", index=False)
    print(f"CTOI {c['ctoi']}: P={c['P']} +/- {c['sP']}, E={c['E']:.5f} +/- {c['sE']}, depth {c['depth']*1e6:.0f} ppm, "
          f"duration {c['dur']*24:.2f} h ({c['user']})")
    print(res.to_string(index=False))
    b_s68 = ev[(ev.cand_id.str.startswith(f"{TIC}b")) & (ev.sector == 68)].t2.iloc[0]
    n = (b_s68 - cands.loc["b", "t0"]) / c["P"]
    print(f"S68 dip minus candidate b = {b_s68 - cands.loc['b', 't0']:.4f} d = {n:.3f} CTOI periods; "
          f"P from b and S68 at n={round(n)}: {(b_s68 - cands.loc['b', 't0']) / round(n):.6f} d")

    # 2. aliases from the S68 transit
    summ = []
    for pl in ("a", "b"):
        e68 = ev[(ev.cand_id.str.startswith(f"{TIC}{pl}")) & (ev.sector == 68)].iloc[0]
        dur, depth, derr = e68.dur_ref_h / 24, e68.depth_ref_ppm * 1e-6, e68.depth_err_ppm * 1e-6
        cov = pd.read_csv(WORK / "coverage.csv").set_index("tic").loc[TIC]
        pmin = core.pmin_from_duration(dur, cov.mstar, cov.rstar, depth)
        at = pd.DataFrame(core.alias_table(e68.t1, e68.t2, dur, depth, derr, t, f, p_floor=1.0, pmin_density=pmin))
        at.to_csv(OUT / f"tic302377820{pl}_aliases_from_S68.csv", index=False)
        ok = at[at.status.isin(["allowed", "supported"])]
        near = at.iloc[(at.period_d - c["P"]).abs().argmin()]
        summ.append(dict(candidate=pl, t1=e68.t1, t2_S68=e68.t2, dT_d=e68.t2 - e68.t1, pmin_circ_d=pmin,
                         n_allowed=int((at.status == "allowed").sum()), n_supported=int((at.status == "supported").sum()),
                         n_excluded=int((at.status == "excluded").sum()),
                         n_below_pmin=int((at.status == "allowed_below_pmin").sum()),
                         longest_allowed_d=ok.period_d.max(), shortest_allowed_d=ok.period_d.min(),
                         supported=";".join(f"{p:.4f}" for p in at.loc[at.status == "supported", "period_d"]),
                         alias_nearest_ctoi_P=near.period_d, its_status=near.status, its_n=int(near.n)))
    sm = pd.DataFrame(summ)
    sm.to_csv(OUT / "tic302377820_aliases_from_S68_summary.csv", index=False)
    print(sm.T.to_string())

    # 3. fold on the CTOI ephemeris
    singles = [cands.loc["a", "t0"], cands.loc["b", "t0"]] + sorted(ev.t2.round(2).unique())
    near_single = np.zeros(t.size, bool)
    for ts in singles:
        near_single |= np.abs(t - ts) < 0.5
    ph = ((t - c["E"]) / c["P"] + 0.5) % 1 - 0.5
    fold = []
    for lab, m in (("all data", np.ones(t.size, bool)), ("single events masked (+/-0.5 d)", ~near_single)):
        inn = m & (np.abs(ph * c["P"]) < c["dur"] / 2 * 0.8)
        out = m & (np.abs(ph * c["P"]) > c["dur"]) & (np.abs(ph * c["P"]) < 3 * c["dur"])
        d = np.median(f[out]) - np.mean(f[inn])
        err = core.robust_std(f[out]) / np.sqrt(inn.sum())
        fold.append(dict(sample=lab, n_in=int(inn.sum()), depth_ppm=d * 1e6, err_ppm=err * 1e6))
    fd = pd.DataFrame(fold)
    fd.to_csv(OUT / "tic302377820_ctoi_fold.csv", index=False)
    print(fd.to_string(index=False))
    plot(t, f, sec, c, marks, ph, near_single)


def plot(t, f, sec, c, marks, ph, near_single):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = marks
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(3, 4)
    for i, (lab, tm) in enumerate(show):
        a = fig.add_subplot(gs[i // 4, i % 4])
        w = 0.6
        m = np.abs(t - tm) < w
        x, y = (t[m] - tm) * 24, (f[m] - 1) * 1e6
        a.plot(x, y, ".", ms=2, color="0.6")
        if m.sum():
            e = np.arange(-w, w + 1 / 48, 1 / 48) * 24
            idx = np.digitize(x, e)
            a.plot([x[idx == j].mean() for j in np.unique(idx)], [y[idx == j].mean() for j in np.unique(idx)],
                   "o", ms=3, color="C0")
        k0 = int(np.round((tm - c["E"]) / c["P"]))
        for k in (k0 - 1, k0, k0 + 1):
            tp = c["E"] + k * c["P"]
            if abs(tp - tm) < w:
                s = np.hypot(c["sE"], k * c["sP"])
                a.axvline((tp - tm) * 24, color="C3")
                a.axvspan((tp - s - tm) * 24, (tp + s - tm) * 24, color="C3", alpha=0.15)
        a.axhline(-c["depth"] * 1e6, color="C3", ls=":", lw=0.8)
        a.set_title(f"{lab}\nT={tm:.3f}", fontsize=8)
        a.set_xlabel("hours")
    a = fig.add_subplot(gs[2, 2:])
    for lab, m, col in (("all", np.ones(t.size, bool), "0.7"), ("singles masked", ~near_single, "C0")):
        x = ph[m] * c["P"] * 24
        sel = np.abs(x) < 15
        e = np.arange(-15, 15.25, 0.5)
        idx = np.digitize(x[sel], e)
        a.plot([x[sel][idx == j].mean() for j in np.unique(idx)],
               [(f[m][sel][idx == j].mean() - 1) * 1e6 for j in np.unique(idx)], "o-", ms=3, color=col, label=lab)
    a.axhline(-c["depth"] * 1e6, color="C3", ls=":", lw=0.8, label="CTOI depth")
    a.legend(fontsize=8)
    a.set_title(f"fold on CTOI P={c['P']:.6f} d (30-min bins)", fontsize=9)
    a.set_xlabel("hours from CTOI epoch")
    fig.suptitle(f"TIC {TIC}: dips with CTOI {c['ctoi']} predicted epochs (red line, shaded +/-1 sigma)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "tic302377820_ctoi_check.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
