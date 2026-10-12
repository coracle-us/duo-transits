"""Light-curve loading helpers (independent re-analysis of TIC 159614322)."""
import glob, os, re
import numpy as np
import pandas as pd

RAW = r"C:\git\duo-transits\data\raw\lightcurves\159614322"
OUT = os.path.dirname(os.path.abspath(__file__))


def list_files():
    out = []
    for f in sorted(glob.glob(os.path.join(RAW, "s*_*.csv.gz"))):
        m = re.match(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", os.path.basename(f))
        out.append((int(m.group(1)), m.group(2), int(m.group(3)), f))
    return out


def load(sector, product, which="flux", strict_quality=True):
    """Return dataframe with time, f (normalised), e, q, plus raw columns.
    which: 'flux' (PDCSAP for SPOC/TESS-SPOC, KSPSAP-like detrended for QLP 'flux'),
           'sap' for raw SAP flux.
    """
    fs = [x for x in list_files() if x[0] == sector and x[1] == product]
    if not fs:
        return None
    df = pd.read_csv(fs[0][3])
    col = {"flux": "flux", "sap": "sap_flux"}.get(which, which)
    df = df[np.isfinite(df["time"]) & np.isfinite(df[col])].copy()
    if strict_quality:
        df = df[df["quality"] == 0]
    med = np.nanmedian(df[col])
    df["f"] = df[col] / med
    if "flux_err" in df and which == "flux":
        df["e"] = df["flux_err"] / med
    else:
        df["e"] = np.nan
    df["cad"] = fs[0][2]
    return df.reset_index(drop=True)
