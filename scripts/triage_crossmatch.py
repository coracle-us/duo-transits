"""
Triage step 1 (automated): cross-match the 85 candidates against
  - NASA Exoplanet Archive pscomppars (confirmed planets, by tic_id)
  - ExoFOP TOI table (dispositions, periods, comments)
  - ExoFOP CTOI table (periods, notes)
For each catalog entry, test whether its period is consistent with an alias of the duotransit
(P = dT / n for integer n, within tolerance) and whether its epoch lines up with Tc1, so we know
the catalog entry refers to the *same* signal rather than another planet in the system.

Output: data/catalog/crossmatch.csv (one row per catalog entry matched to a candidate TIC).
This is evidence for status.csv, not the status itself (that is assigned by hand with sources).

Usage:  python scripts/triage_crossmatch.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
CAT = ROOT / "data" / "catalog"


def alias_check(P, T0, cand, ptol=0.002):
    """Is period P an alias dT/n of this candidate, and does epoch T0 (BTJD) phase onto Tc1?
    Returns (n_nearest, frac_period_err, phase_offset_hours)."""
    if not np.isfinite(P) or P <= 0:
        return (np.nan, np.nan, np.nan)
    n = cand.dT / P
    n_int = max(1, int(round(n)))
    frac = abs(cand.dT / n_int - P) / P
    ph = np.nan
    if np.isfinite(T0):
        k = round((cand.Tc1 - T0) / P)
        ph = (cand.Tc1 - (T0 + k * P)) * 24
    return (n_int, frac, ph)


def main():
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC")
    rows = []

    ps = pd.read_csv(RAW / "nasa_pscomppars.csv")
    ps["tic"] = pd.to_numeric(ps["tic_id"].str.replace("TIC ", ""), errors="coerce")
    for _, r in ps[ps.tic.isin(cand.index)].iterrows():
        c = cand.loc[int(r.tic)]
        T0 = r.pl_tranmid - 2457000 if np.isfinite(r.pl_tranmid) else np.nan
        n, frac, ph = alias_check(r.pl_orbper, T0, c)
        rows.append(dict(TIC=int(r.tic), source="pscomppars", name=r.pl_name, disposition="CONFIRMED",
                         period=r.pl_orbper, epoch_btjd=T0, n_alias=n, frac_P_err=frac,
                         Tc1_phase_off_h=ph, ref=r.pl_orbper_reflink, note=r.disc_refname))

    toi = pd.read_csv(RAW / "exofop_toi.csv")
    for _, r in toi[toi["TIC ID"].isin(cand.index)].iterrows():
        c = cand.loc[int(r["TIC ID"])]
        T0 = r["Epoch (BJD)"] - 2457000 if np.isfinite(r["Epoch (BJD)"]) else np.nan
        n, frac, ph = alias_check(r["Period (days)"], T0, c)
        rows.append(dict(TIC=int(r["TIC ID"]), source="TOI", name=f"TOI-{r['TOI']}",
                         disposition=f"TFOPWG={r['TFOPWG Disposition']} TESS={r['TESS Disposition']}",
                         period=r["Period (days)"], epoch_btjd=T0, n_alias=n, frac_P_err=frac,
                         Tc1_phase_off_h=ph, ref=r.get("Planet Name"), note=r.get("Comments")))

    ctoi = pd.read_csv(RAW / "exofop_ctoi.csv")
    for _, r in ctoi[ctoi["TIC ID"].isin(cand.index)].iterrows():
        c = cand.loc[int(r["TIC ID"])]
        ep = r["Transit Epoch (BJD)"]
        T0 = ep - 2457000 if np.isfinite(ep) else np.nan
        n, frac, ph = alias_check(r["Period (days)"], T0, c)
        rows.append(dict(TIC=int(r["TIC ID"]), source="CTOI", name=f"CTOI-{r['CTOI']}",
                         disposition=f"TFOPWG={r['TFOPWG Disposition']} user={r['User Disposition']} "
                                     f"promoted={r['Promoted to TOI']}",
                         period=r["Period (days)"], epoch_btjd=T0, n_alias=n, frac_P_err=frac,
                         Tc1_phase_off_h=ph, ref=r.get("Paper"),
                         note=f"{r.get('Group')}|{r.get('Tag')}|{r.get('Notes')}"))

    out = pd.DataFrame(rows).sort_values(["TIC", "source"])
    # Same signal if the period is an alias (<0.2% off) and the epoch phases onto Tc1 within 12 h
    out["same_signal"] = (out.frac_P_err < 0.002) & (out.Tc1_phase_off_h.abs() < 12)
    out.to_csv(CAT / "crossmatch.csv", index=False)
    print(f"{len(out)} catalog entries for {out.TIC.nunique()} of {len(cand)} candidates "
          f"-> {CAT / 'crossmatch.csv'}")
    pd.set_option("display.width", 250, "display.max_colwidth", 70)
    print(out.drop(columns=["epoch_btjd"]).to_string(index=False))


if __name__ == "__main__":
    main()
