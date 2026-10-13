"""
Overall recovery of the catalogued events for the 306 classified stars (editor Q6), split by the
product the pipeline used at the epoch and by own-photometry use; plus the recovered depth ratio.

Inputs: final_classification_2026-10-11_rnaas.csv (frozen), per-star JSONs in singles/results/<run>/,
        the light-curve files on disk (to identify the product alias_pipeline.load_lightcurves picks
        for the epoch sector; PRODUCT_RANK + OWN, as in pilot.py).
Output: singles/revision/recovery_by_product.csv (per star) and a printed summary.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REV = Path(__file__).resolve().parent
sys.path.insert(0, str(REV))
import revlib as rl  # noqa: E402

ROOT = rl.ROOT


def main():
    cl = pd.read_csv(ROOT / "singles" / "results" / "final_classification_2026-10-11_rnaas.csv")
    rows = []
    for _, r in cl.iterrows():
        tic = int(r.TIC)
        js = ROOT / "singles" / "results" / r.run / f"{tic}.json"
        d = json.loads(js.read_text())
        ep = d.get("event_T", d["epoch_btjd"])
        lc = rl.ap.load_lightcurves(tic, products=rl.ap.PRODUCT_RANK + ["OWN"])
        prod, sec = "none", None
        if len(lc):
            j = np.argmin(np.abs(lc.time.to_numpy() - ep))
            if abs(lc.time.iloc[j] - ep) < 1.0:
                prod, sec = lc["product"].iloc[j], int(lc.sector.iloc[j])
        rows.append(dict(TIC=tic, run=r.run, cls=r.cls, event_verdict=d.get("event_verdict"),
                         recovered=r.cls != "event-not-recovered", event_product=prod, event_sector=sec,
                         own_used_any=bool(d.get("own_photometry_sectors")),
                         depth_ratio=(d["event_depth_ppm"] / d["depth_ppm"]) if d.get("event_depth_ppm") is not None else np.nan,
                         event_snr=d.get("event_snr")))
    df = pd.DataFrame(rows)
    df.to_csv(REV / "recovery_by_product.csv", index=False)
    print(pd.crosstab(df.event_product, df.recovered, margins=True))
    print(pd.crosstab(df.own_used_any, df.recovered, margins=True))
    print(df.event_verdict.value_counts())
    rec = df[df.recovered]
    print("recovered: depth ratio median %.2f, 16-84%% %.2f-%.2f; SNR median %.1f" % (
        rec.depth_ratio.median(), rec.depth_ratio.quantile(0.16), rec.depth_ratio.quantile(0.84), rec.event_snr.median()))
    print("event products among recovered:", rec.event_product.value_counts().to_dict())


if __name__ == "__main__":
    main()
