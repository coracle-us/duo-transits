"""Work estimate: target-sectors for SNR>=7 candidates, and how many have any pipeline light curve
(SPOC, TESS-SPOC, QLP) vs need our own TESScut photometry. Output: data/work_estimate.txt"""
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parents[1]
f = pd.read_csv(HERE / "data" / "feasibility.csv")
inv = pd.read_csv(HERE / "data" / "tess_lc_inventory.csv")
pipe = inv[inv.author.isin(["SPOC", "TESS-SPOC", "QLP"])].dropna(subset=["sector"])
have = set(zip(pipe.TIC, pipe.sector.astype(int)))
L = []
for lab, m in [("planet-like singles SNR>=7", f.planet_like & f.snr_ok & (f.n_kepler_transits == 1)),
               ("all singles SNR>=7 (incl. EB-like)", f.snr_ok & (f.n_kepler_transits == 1)),
               ("all candidates SNR>=7", f.snr_ok)]:
    s = f[m].drop_duplicates("TIC")
    ts = [(r.TIC, int(x)) for r in s.itertuples() for x in str(r.tess_sectors_observed).split()]
    n_pipe = sum(t in have for t in ts)
    L.append(f"{lab}: {len(s)} stars, {len(ts)} target-sectors; with a SPOC/TESS-SPOC/QLP LC: {n_pipe}; "
             f"own TESScut needed: {len(ts) - n_pipe}; stars with no pipeline LC at all: "
             f"{sum(all((r.TIC, int(x)) not in have for x in str(r.tess_sectors_observed).split()) for r in s.itertuples())}")
txt = "\n".join(L)
(HERE / "data" / "work_estimate.txt").write_text(txt + "\n", encoding="utf-8")
print(txt)
