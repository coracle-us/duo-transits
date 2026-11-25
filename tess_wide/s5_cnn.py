"""
tess_wide step 5 (PLAN.md section 6): the Salinas et al. 2026 CNN (helemysm/singletransiters) as a re-measurement
vetter. Runs in that repository's own environment; its code is imported from the local clone, never copied here.

  .../salinas-second-transits/.venv-cnn/bin/python tess_wide/s5_cnn.py --sector 104 \
      --cnn-repo .../salinas-second-transits/singletransiters

Input per event: the star's QLP SAP_FLUX where QUALITY == 0 (finite), as a lightkurve LightCurve, through their
preprocess_lc (30-min bins) and _run_inference_on_lc at the box time. Injected events get the same injection
(inject.py) first. Events: real and control survivors of s4 (survivor == True) and the s4 injection events that
survived (rec_stage).

Pass rule (PLAN 6, all four): depth/sigma_depth >= 5; CNN depth within 0.5-2x the box depth; CNN duration within
0.5-2x the box duration; |CNN mid-time - box time| <= max(1 h, D_box/2).
Out of range (not passed or failed): box depth > 15,000 ppm, or box duration outside 1.5-15 h.

Output (outside git): tess_wide/data/s<NNN>_cnn_{real,ctrl,inj}.csv, s<NNN>_cnn.json
"""
import argparse
import contextlib
import functools
import io
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE))
import inject  # noqa: E402

SIG_MIN, RATIO, DT_MIN_H = 5.0, (0.5, 2.0), 1.0
DEPTH_MAX_PPM, DUR_RANGE_H = 15000.0, (1.5, 15.0)


def verdict(r, box_depth_ppm, box_D_h, box_T):
    if box_depth_ppm > DEPTH_MAX_PPM or not (DUR_RANGE_H[0] <= box_D_h <= DUR_RANGE_H[1]):
        return "out_of_range"
    ok = (r["sigma_depth"] > 0 and r["depth_ppm"] / r["sigma_depth"] >= SIG_MIN
          and RATIO[0] <= r["depth_ppm"] / box_depth_ppm <= RATIO[1]
          and RATIO[0] <= r["duration_h"] / box_D_h <= RATIO[1]
          and abs(r["t_mid_pred_btjd"] - box_T) * 24 <= max(DT_MIN_H, box_D_h / 2))
    return "pass" if ok else "fail"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    ap.add_argument("--cnn-repo", required=True)
    ap.add_argument("--kinds", default="real,ctrl,inj")
    a = ap.parse_args()
    s = a.sector
    repo = Path(a.cnn_repo).resolve()
    sys.path.insert(0, str(repo))
    os.chdir(repo)  # their config uses repo-relative paths
    import lightkurve as lk
    from computeparams import param_inference_core as pic
    pic.load_bundle = functools.lru_cache(maxsize=2)(pic.load_bundle)  # load the weights once, not per event
    summ = {}
    for kind in a.kinds.split(","):
        st = pd.read_csv(DATA / f"s{s:04d}_stage_{kind}.csv")
        if kind == "inj":
            tr = pd.read_csv(DATA / f"s{s:04d}_inj_truth.csv")
            idx = tr.loc[tr.rec_stage, "ev_index"].astype(int)
            ev = st.set_index("ev_index").loc[idx].reset_index()
            injs = tr.set_index("tic")
        else:
            ev = st[st.survivor].reset_index(drop=True)
        rows, t0 = [], time.time()
        for e in ev.itertuples():
            z = np.load(DATA / f"lc_s{s:04d}" / f"{int(e.tic)}.npz")
            t, f, q = z["time"], z["sap"].astype(float), z["quality"]
            if kind == "inj":
                f = f * inject.model(t, injs.loc[int(e.tic)].to_dict())
            m = np.isfinite(t) & np.isfinite(f) & (q == 0)
            if kind == "ctrl":
                f = 2 * np.nanmedian(f[m]) - f  # the control searched 2 - flux
            row = dict(tic=int(e.tic), T=e.T, D_h=e.D * 24, depth_ppm_box=e.depth * 1e6)
            try:
                lc = pic.preprocess_lc(lk.LightCurve(time=t[m], flux=f[m] / np.nanmedian(f[m])))
                with contextlib.redirect_stdout(io.StringIO()):
                    r = pic._run_inference_on_lc(lc, float(e.T), plot=False)
                row.update({k: r[k] for k in ("depth_ppm", "sigma_depth", "duration_h", "sigma_duration",
                                              "t_mid_pred_btjd")})
                row["cnn"] = verdict(r, row["depth_ppm_box"], row["D_h"], e.T)
            except Exception as ex:  # noqa: BLE001
                row.update(cnn="error", error=str(ex)[:200])
            rows.append(row)
        out = pd.DataFrame(rows)
        out.to_csv(DATA / f"s{s:04d}_cnn_{kind}.csv", index=False)
        summ[kind] = dict(events=len(out), seconds=round(time.time() - t0, 1),
                          verdicts=out.cnn.value_counts().to_dict() if len(out) else {})
        print(kind, summ[kind], flush=True)
    (DATA / f"s{s:04d}_cnn.json").write_text(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
