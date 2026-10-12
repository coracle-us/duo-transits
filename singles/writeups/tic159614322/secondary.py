"""
Secondary-eclipse check for TIC 159614322 at phase 0.5 of the two surviving ephemerides
(P = 89.965 d and 179.929 d, ephemeris.json). Box of the primary duration at each phase-0.5
epoch that falls in TESS data (any product, best-ranked per sector), local linear baseline;
weighted mean depth over epochs.
Output: secondary.json
Run:    python singles/writeups/tic159614322/secondary.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402

TIC = 159614322
D = 9.66 / 24


def main():
    eph = json.load(open(HERE / "ephemeris.json"))["ephemerides"]
    lc = ap.load_lightcurves(TIC)
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    out = {}
    for key in ("n8", "n4"):
        e = eph[key]
        rows = []
        for k in range(-30, 30):
            T = e["T0"] + (k + 0.5) * e["P"]
            m = np.abs(t - T) < D / 2 + 0.6
            if m.sum() < 20:
                continue
            tt, ff = t[m], f[m]
            inn = np.abs(tt - T) < D / 2
            oot = ~inn
            if inn.sum() < 0.6 * D / np.median(np.diff(tt)) or (tt[oot] < T).sum() < 5 or (tt[oot] > T).sum() < 5:
                continue
            p = np.polyfit(tt[oot] - T, ff[oot], 1)
            r = ff / np.polyval(p, tt - T) - 1
            d = -r[inn].mean()
            s = r[oot].std() / np.sqrt(inn.sum())
            rows.append(dict(T=float(T), sector=int(lc.sector.to_numpy()[m][0]), depth_ppm=float(d * 1e6), err_ppm=float(s * 1e6)))
        if rows:
            w = np.array([1 / r["err_ppm"] ** 2 for r in rows])
            dm = float(np.sum(w * [r["depth_ppm"] for r in rows]) / w.sum())
            out[key] = dict(P=e["P"], epochs=rows, mean_depth_ppm=dm, mean_err_ppm=float(1 / np.sqrt(w.sum())))
        else:
            out[key] = dict(P=e["P"], epochs=[], note="no phase-0.5 epoch in TESS data")
    (HERE / "secondary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
