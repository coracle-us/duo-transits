"""
EPIC 202126877 (TIC 86491361): eclipse depths over time, like for like, and K2 aperture checks
for the RNAAS draft (review_fable3 item 4).

Inputs:
- the blind analysis's per-event, per-product fits (k2singles/independent_table/data/final_details.json;
  data/ is gitignored, so the extracted rows are saved here);
- the K2 C0 light-curve header (k2singles/cache/202126877/dl/ktwo202126877-c00_llc.fits).
Output: k2singles/rnaas/epic202126877_depths.csv (one row per event and product) and printed
summary lines (own TESScut r=1.5 px series; secondary/primary ratios; CROWDSAP, FLFRCSAP).
Usage:  python k2singles/scripts/epic202126877_depths.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]


def main():
    d = json.loads((ROOT / "k2singles" / "independent_table" / "data" / "final_details.json").read_text())["202126877"]
    rows = []
    for ev in d["events"]:
        for prod, (t0, depth_ppm) in ev["per_prod"].items():
            rows.append(dict(kind=ev["label"], t0_btjd=round(t0, 5), product=prod, depth_ppm=int(depth_ppm)))
    df = pd.DataFrame(rows).sort_values(["kind", "t0_btjd", "product"])
    df.to_csv(ROOT / "k2singles" / "rnaas" / "epic202126877_depths.csv", index=False)
    own = df[df["product"].str.contains("own-TESScut")]
    for kind in ("P", "S"):
        sub = own[own.kind == kind]
        print(kind, "own TESScut r=1.5:", [(r.t0_btjd, round(r.depth_ppm / 1e4, 1)) for r in sub.itertuples()])
    for era, (lo, hi) in {"2021": (2400, 2600), "2023": (3200, 3300)}.items():
        p = own[(own.kind == "P") & own.t0_btjd.between(lo, hi)].depth_ppm.mean()
        s = own[(own.kind == "S") & own.t0_btjd.between(lo, hi)].depth_ppm.mean()
        print(f"{era}: S/P = {s / p:.2f} (own TESScut)")
    k2 = df[df["product"].str.startswith("K2")]
    print("K2 C0 depths by product (%):", [(r.product, round(r.depth_ppm / 1e4, 2)) for r in k2.itertuples()])
    with fits.open(ROOT / "k2singles" / "cache" / "202126877" / "dl" / "ktwo202126877-c00_llc.fits") as h:
        hd = h[1].header
        ap = h[2].data
        n_pix = int(np.sum((ap & 2) > 0))
        print(f"K2 C0 aperture: {n_pix} pixels; CROWDSAP {hd.get('CROWDSAP')}; FLFRCSAP {hd.get('FLFRCSAP')}; "
              f"Kp {h[0].header.get('KEPMAG')}")


if __name__ == "__main__":
    main()
