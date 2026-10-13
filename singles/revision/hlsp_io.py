"""Download a MAST light-curve product by dataURI into singles/revision/data/lc/<TIC>/ and read
it into a uniform table (time in BTJD, several flux columns, quality). Used for products the
main fetch never stored (TARS, TASOC, CDIPS, T16) and for any sector file missing locally."""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

REV = Path(__file__).resolve().parent
DATA = REV / "data" / "lc"

# flux-like columns per HLSP, in order of preference (first = primary)
CANDIDATE_FLUX = ["PDCSAP_FLUX", "SAP_FLUX", "KSPSAP_FLUX", "FLUX_CORR", "FLUX_RAW", "FLUX",
                  "CAL_PSF_FLUX", "CAL_APER_FLUX", "PSF_FLUX", "APER_FLUX", "CORR_FLUX", "RAW_FLUX",
                  "PCA_FLUX", "AP_FLUX", "FLUX_AP", "PSF_FLUX_COR", "AP2_FLUX_COR", "FCOR",
                  "IRM1", "IRM2", "IRM3", "TFA1", "TFA2", "TFA3", "PCA1", "PCA2", "PCA3"]
QUAL = ["QUALITY", "QUALITY_FLAG", "FLAG", "DQUALITY", "IRQ1"]


def download(uri, tic):
    from astroquery.mast import Observations
    d = DATA / str(tic)
    d.mkdir(parents=True, exist_ok=True)
    p = d / uri.split("/")[-1]
    if not p.exists():
        Observations.download_file(uri, local_path=str(p), cache=False)
    return p


def read_any(path):
    """Return (DataFrame with 'time', flux columns, 'quality'), list of flux column names."""
    with fits.open(path) as h:
        hdu = next(x for x in h[1:] if isinstance(x, fits.BinTableHDU))
        d = hdu.data
        names = [n.upper() for n in d.columns.names]
        col = {n.upper(): n for n in d.columns.names}
        tname = next((n for n in ("TIME", "TMID_BJD", "BTJD", "BJD") if n in names), None)
        if tname is None:
            raise ValueError(f"no time column in {path.name}: {names}")
        t = np.asarray(d[col[tname]], float)
        if np.nanmedian(t) > 2.4e6:
            t = t - 2457000.0
        out = {"time": t}
        fl = []
        for n in CANDIDATE_FLUX:
            if n in names:
                v = np.asarray(d[col[n]], float)
                if v.ndim != 1:
                    continue
                if n.startswith(("IRM", "TFA", "PCA")) and np.nanmedian(v) > 0 and np.nanmedian(v) < 30:
                    v = 10 ** (-0.4 * (v - np.nanmedian(v)))  # magnitudes -> relative flux
                out[n.lower()] = v
                fl.append(n.lower())
        qn = next((n for n in QUAL if n in names), None)
        q = np.asarray(d[col[qn]]) if qn else np.zeros(len(t))
        if q.dtype.kind in "SUO":  # CDIPS IRQ flags are strings
            q = np.array([0 if str(x).strip() in ("", "G", "0") else 1 for x in q])
        out["quality"] = np.nan_to_num(np.asarray(q, float)).astype(int)
    return pd.DataFrame(out), fl
