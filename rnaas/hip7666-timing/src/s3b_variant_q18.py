"""Sensitivity variant (logged in RULES.md change log): S18 eclipse times from QLP instead of
TESS-SPOC, S58 and S85 unchanged. Writes results/eclipse_times_mixq18.csv for s3_models.py."""
import pandas as pd
import hip
m = pd.read_csv(hip.OUTDIR / "eclipse_times_main.csv")
q = pd.read_csv(hip.OUTDIR / "eclipse_times_qlp.csv")
pd.concat([q[q.sector == 18], m[m.sector != 18]]).to_csv(hip.OUTDIR / "eclipse_times_mixq18.csv", index=False)
