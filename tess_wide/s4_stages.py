"""
tess_wide step 4 (PLAN.md sections 5.1-5.3 and 7): detections -> single-event candidates, for the real search,
the inverted control and the injections.

  5.1 common-mode veto  per camera/CCD, event centres in 0.5-d bins between the CCD's first and last event;
                        a bin is flagged if P(N >= n | lambda) < 1e-6, lambda = mean count per bin, recomputed
                        once without the flagged bins. Real and control each use their own flags; injections use
                        the real run's flags (PLAN 7).
  5.2 one event/star    exactly one event left -> single; two with depth and duration within 0.5-2x of each other
                        -> repeating; otherwise multi-event.
  5.3 size              Rp = R* sqrt(depth (1 + contratio)) > 2 R_J -> EB-size.
Injections: recovered at detection = an event within max(0.1 d, D_inj/2) of T_inj; end to end also needs the star
to pass 5.1-5.3 with that event (the CNN stage is added by s5_cnn.py).

Output (outside git): tess_wide/data/s<NNN>_stage_{real,ctrl,inj}.csv, s<NNN>_cm_bins.csv, s<NNN>_stages.json
Usage:  python tess_wide/s4_stages.py --sector 104
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
BIN, P_FLAG = 0.5, 1e-6
RATIO = (0.5, 2.0)
RP_MAX_RJ, RSUN_RJ = 2.0, 9.731


def flag_bins(ev):
    rows = []
    for (cam, ccd), g in ev.groupby(["camera", "ccd"]):
        t0, t1 = np.floor(g["T"].min() / BIN) * BIN, g["T"].max()
        edges = np.arange(t0, t1 + BIN, BIN)
        n, _ = np.histogram(g["T"], edges)
        lam = n.mean()
        flagged = poisson.sf(n - 1, lam) < P_FLAG
        if (~flagged).any():
            lam = n[~flagged].mean()
            flagged = poisson.sf(n - 1, lam) < P_FLAG
        for i in np.flatnonzero(flagged):
            rows.append(dict(camera=cam, ccd=ccd, t_lo=edges[i], t_hi=edges[i + 1], n=int(n[i]), lam=lam))
    return pd.DataFrame(rows, columns=["camera", "ccd", "t_lo", "t_hi", "n", "lam"])


def vetoed(ev, bins):
    v = np.zeros(len(ev), bool)
    for b in bins.itertuples():
        v |= ((ev.camera == b.camera) & (ev.ccd == b.ccd) & (ev["T"] >= b.t_lo) & (ev["T"] < b.t_hi)).to_numpy()
    return v


def per_star(ev):
    """Label each surviving event single / repeating / multi by its star's surviving events."""
    lab = pd.Series("", index=ev.index)
    for tic, g in ev.groupby("tic"):
        if len(g) == 1:
            lab[g.index] = "single"
        elif len(g) == 2:
            a, b = g.iloc[0], g.iloc[1]
            rd, rD = a.depth / b.depth, a.D / b.D
            same = RATIO[0] <= rd <= RATIO[1] and RATIO[0] <= rD <= RATIO[1]
            lab[g.index] = "repeating" if same else "multi"
        else:
            lab[g.index] = "multi"
    return lab


def stage(ev, bins, stars):
    ev = ev.copy()
    ev["cm_veto"] = vetoed(ev, bins) if len(ev) else []
    keep = ev[~ev.cm_veto]
    ev["star_class"] = ""
    ev.loc[keep.index, "star_class"] = per_star(keep)
    st = stars.set_index("id")
    R = ev.tic.map(st.rad)
    c = ev.tic.map(st.contratio).fillna(0.0)
    ev["rp_rj"] = R * np.sqrt(ev.depth.clip(lower=0) * (1 + c)) * RSUN_RJ
    ev["size_ok"] = ev.rp_rj <= RP_MAX_RJ
    ev["survivor"] = (~ev.cm_veto) & (ev.star_class == "single") & ev.size_ok
    return ev


def counts(ev, n_stars):
    return dict(stars_searched=n_stars, detections=len(ev), stars_with_detection=int(ev.tic.nunique()),
                after_common_mode=int((~ev.cm_veto).sum()),
                star_class=ev.loc[~ev.cm_veto, "star_class"].value_counts().to_dict(),
                single_events=int((ev.star_class == "single").sum()),
                after_size=int(ev.survivor.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    s = ap.parse_args().sector
    stars = pd.read_csv(DATA / f"s{s:04d}_stars.csv")
    n_searched = sum(1 for p in (DATA / f"lc_s{s:04d}").glob("*.npz") if ".tmp" not in p.name)
    out = {}
    real = pd.read_csv(DATA / f"s{s:04d}_events_real.csv")
    ctrl = pd.read_csv(DATA / f"s{s:04d}_events_ctrl.csv")
    b_real, b_ctrl = flag_bins(real), flag_bins(ctrl)
    pd.concat([b_real.assign(kind="real"), b_ctrl.assign(kind="ctrl")]).to_csv(DATA / f"s{s:04d}_cm_bins.csv",
                                                                                 index=False)
    out["cm_bins"] = dict(real=len(b_real), ctrl=len(b_ctrl))
    for kind, ev, b in (("real", real, b_real), ("ctrl", ctrl, b_ctrl)):
        e = stage(ev, b, stars)
        e.to_csv(DATA / f"s{s:04d}_stage_{kind}.csv", index=False)
        out[kind] = counts(e, n_searched)
    fi = DATA / f"s{s:04d}_events_inj.csv"
    if fi.exists():
        inj = pd.read_csv(fi)
        tr = pd.read_csv(DATA / f"s{s:04d}_inj_truth.csv")
        e = stage(inj, b_real, stars) if len(inj) else inj
        tr["rec_detect"], tr["rec_stage"], tr["ev_index"] = False, False, -1
        for i, r in tr.iterrows():
            g = e[e.tic == r.tic]
            m = g[np.abs(g["T"] - r["T"]) <= max(0.1, r.dur_h / 48)]
            if len(m):
                j = m.snr.idxmax()
                tr.loc[i, ["rec_detect", "ev_index"]] = True, j
                tr.loc[i, "rec_stage"] = bool(e.loc[j, "survivor"])
        st = stars.set_index("id")
        tr["rp_true_rj"] = tr.tic.map(st.rad) * np.sqrt(tr.depth_true) * RSUN_RJ
        e.to_csv(DATA / f"s{s:04d}_stage_inj.csv", index=True, index_label="ev_index")
        tr.to_csv(DATA / f"s{s:04d}_inj_truth.csv", index=False)
        out["inj"] = dict(injected=len(tr), rec_detect=int(tr.rec_detect.sum()), rec_stage=int(tr.rec_stage.sum()),
                          rp_gt_2=int((tr.rp_true_rj > RP_MAX_RJ).sum()))
    (DATA / f"s{s:04d}_stages.json").write_text(json.dumps(out, indent=1, default=int))
    print(json.dumps(out, indent=1, default=int))


if __name__ == "__main__":
    main()
