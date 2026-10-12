"""
Collect the per-system fit results (results/oc/<TIC>_fit.json), Gaia/TIC flags
(data/gaia_tic_flags.csv) and the by-eye classification (data/manual_classification.csv, if
present) into one table.

Minimum tertiary mass: f(m3) from the LTTE fit; m12 = 2 x TIC mass (twin assumption; the TIC
mass is a single-star estimate from the blended Teff), or 2.0 Msun when TIC has no mass.
This is only meaningful when the ETV is LTTE-dominated (P3 >> P1 and pri/sec in phase).

Per-product timing offsets: weighted mean of the LTTE-model residuals per product
(SPOC_120 / TESS-SPOC / QLP), as a check on cadence-dependent bias.

Usage:  python sidequests/eb_timing/scripts/summarize.py LIST_CSV OUT_CSV
  LIST_CSV: a csv with a TIC column (e.g. data/pilot_list.csv)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import etvlib as el  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


def product_offsets(tic):
    f = HERE / "results" / "oc" / f"{tic}_oc.csv"
    oc = pd.read_csv(f)
    if "resid_ltte" not in oc:
        return {}
    oc = oc[oc.flag.isna() | (oc.flag == "")].dropna(subset=["resid_ltte"])
    out = {}
    for p, g in oc.groupby("product"):
        w = 1 / g.sig_s ** 2
        mu = np.sum(w * g.resid_ltte) / w.sum()
        # error from scatter (red noise makes the formal error too small)
        err = g.resid_ltte.std() / np.sqrt(max(g.sector.nunique(), 1))
        out[p] = (mu * 86400, err * 86400, len(g))
    return out


def main():
    lst = pd.read_csv(HERE / sys.argv[1], comment="#")
    gaia = pd.read_csv(HERE / "data" / "gaia_tic_flags.csv").set_index("tic")
    man_f = HERE / "data" / "manual_classification.csv"
    man = pd.read_csv(man_f).set_index("tic") if man_f.exists() else pd.DataFrame()
    rows = []
    for tic in lst.TIC.astype(int):
        f = HERE / "results" / "oc" / f"{tic}_fit.json"
        if not f.exists():
            rows.append(dict(tic=tic, note="no fit"))
            continue
        r = json.loads(f.read_text())
        L = r["ltte"]
        g = gaia.loc[tic] if tic in gaia.index else pd.Series(dtype=float)
        m1 = g.get("tic_mass", np.nan)
        m12 = 2 * m1 if np.isfinite(m1) else 2.0
        po = product_offsets(tic)
        row = dict(tic=tic, P=round(r["P"], 7), n_pri=r["n_pri"], n_sec=r["n_sec"],
                   n_sectors=r["n_sectors"], span_d=round(r["span_d"]),
                   sig_pri_s=round(r["median_sig_pri_s"], 1),
                   chi2r_quad=round(r["chi2r_quad"], 2), chi2r_best=round(r["chi2r_best"], 2),
                   dBICc_quad_vs_lin=round(r["dBICc_quad_vs_lin"], 1),
                   dBICc_aps_vs_quad=round(r["dBICc_aps_vs_quad"], 1),
                   dBICc_ltte_vs_quad=round(r["dBICc_ltte_vs_quad"], 1),
                   ltte_base=r["ltte_base"], P3=round(L["P3"], 1), P3_err=round(L["P3_err"], 1),
                   P3_over_span=round(L["P3_over_span"], 2), A_ltte_s=round(L["A_ltte_s"], 1),
                   asini_err_s=round(L["asini_err_s"], 1), e3=round(L["e3"], 2),
                   rho_pri_sec=round(r["rho_pri_sec"], 2) if r["rho_pri_sec"] == r["rho_pri_sec"] else None,
                   f_m3=float(f"{L['f_m3']:.3g}"), m12_assumed=round(m12, 2),
                   m3_min=round(el.min_m3(L["f_m3"], m12), 2),
                   oconnell_sd_ppt=round(r["oconnell_sd"] * 1e3, 2) if r["oconnell_sd"] else None,
                   ruwe=round(g.get("ruwe", np.nan), 2), nss=g.get("non_single_star", np.nan),
                   nss_type=g.get("nss_nss_solution_type", ""), rv_var=g.get("rv_variable", ""),
                   prod_offsets_s="; ".join(f"{k}:{v[0]:+.1f}±{v[1]:.1f}(n={v[2]})" for k, v in po.items()))
        if tic in man.index:
            row.update(man.loc[tic].to_dict())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / sys.argv[2], index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(df.drop(columns=["prod_offsets_s"]).to_string())
        print(df[["tic", "prod_offsets_s"]].to_string())


if __name__ == "__main__":
    main()
