"""Stage 3: pulsation and flicker metrics per star-sector."""
import warnings
import numpy as np, pandas as pd
from astropy.io import fits
from common import TESS_DIR, SPEC_DIR, PILOT
from metrics import prepare, sector_metrics

warnings.filterwarnings("ignore")
MIN_BINS = 400          # ~8.3 d of 30-min bins
CROWD_MIN, CONTRATIO_MAX = 0.9, 0.1


def load(path, col):
    with fits.open(TESS_DIR / path, memmap=False) as h:
        d = h[1].data
        return (np.asarray(d["TIME"], float), np.asarray(d[col], float),
                np.asarray(d["QUALITY"], int))


def main():
    sec = pd.read_csv(PILOT / "sector_wise.csv")
    pr = pd.read_csv(PILOT / "crowding_proxy.csv").set_index("tic")
    rows, peaks = [], []
    for i, s in sec.iterrows():
        r = s.to_dict()
        t, f, q = load(s.file, s.flux_col)
        r["n_points"] = int(np.sum(np.isfinite(f) & (q == 0)))
        tb, yb = prepare(t, f, q)
        if len(tb) >= MIN_BINS:
            m, pk, spec = sector_metrics(tb, yb)
            r.update(m)
            np.save(SPEC_DIR / f"{s.tic}_{s.sector}.npy", spec.astype(np.float32))
            for j, p in enumerate(pk):
                peaks.append(dict(tic=s.tic, sector=s.sector, rank=j + 1, **p))
        else:
            r["n_bins"] = len(tb)
        rows.append(r)
        if i % 25 == 0:
            print(i, len(sec), flush=True)
    out = pd.DataFrame(rows)
    # crowding: CROWDSAP for SPOC / TESS-SPOC; TIC contratio, else neighbour proxy, otherwise
    spoc = out["product"].isin(["SPOC120", "TESS-SPOC"])
    cr = out.tic.map(pr.contratio_proxy)
    out["crowd_source"] = np.where(spoc, "CROWDSAP",
                          np.where(out.contratio.notna(), "tic_contratio", "tic_proxy"))
    ctr = out.contratio.where(out.contratio.notna(), cr)
    out["contratio_used"] = np.where(spoc, np.nan, ctr)
    out["crowd"] = np.where(spoc, out.crowdsap, 1.0 - ctr)
    out["crowd_ok"] = np.where(spoc, out.crowdsap >= CROWD_MIN, ctr <= CONTRATIO_MAX)
    out["metrics_ok"] = out.n_bins >= MIN_BINS
    out["usable"] = out.wise_usable & out.metrics_ok
    out.to_csv(PILOT / "sector_metrics.csv", index=False)
    pd.DataFrame(peaks).to_csv(PILOT / "peaks.csv", index=False)
    print("sectors", len(out), "metrics_ok", int(out.metrics_ok.sum()), "usable",
          int(out.usable.sum()), "usable & crowd_ok", int((out.usable & out.crowd_ok).sum()))


if __name__ == "__main__":
    main()
