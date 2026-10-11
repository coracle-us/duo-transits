"""Merge detection, automatic flags, vetting metrics and the human vetting log.

Inputs : candidates_raw.csv, candidates_vetmetrics.csv, manual_vetting.csv,
         sectors.csv, coverage.csv, injrec_per_star.csv
Outputs: candidates.csv (every candidate, decision + reason)
         flares.csv     (accepted flares with ED and energies)
         per_star.csv   (count, observing time, rate, limits per host with 2-min data)
"""
import numpy as np
import pandas as pd

import flarelib as F

# Poisson one-sided 95% upper limits on the mean for n = 0..5 observed events
POISSON_UL95 = {0: 3.00, 1: 4.74, 2: 6.30, 3: 7.75, 4: 9.15, 5: 10.51}


def auto_reason(r):
    rs = []
    if r.near_dump:
        rs.append("within 10 min of momentum dump")
    if r.edge:
        rs.append("within 0.5 d of segment edge")
    if not r.shape_ok:
        rs.append("fails fast-rise/slow-decay")
    if isinstance(r.compare_hits, str) and r.compare_hits.strip():
        rs.append(f"simultaneous excursion in TIC {r.compare_hits}")
    return "; ".join(rs) if rs else "passes automatic checks"


def main():
    cat = F.load_catalog()
    raw = pd.read_csv(F.SQ / "candidates_raw.csv")
    vm = pd.read_csv(F.SQ / "candidates_vetmetrics.csv", dtype={"compare_hits": str})
    man = pd.read_csv(F.SQ / "manual_vetting.csv")
    c = raw.merge(vm, on="cand_id", how="left").merge(man, on="cand_id", how="left")
    c["auto_result"] = c.apply(auto_reason, axis=1)
    c["decision"] = c["decision"].fillna("unvetted")
    c = c.drop(columns=["i0", "i1", "ipk"])
    c.to_csv(F.SQ / "candidates.csv", index=False)

    bc = F.bol_correction(9000.0)
    acc = c[c.decision == "accept"].copy()
    acc["Teff"] = acc.tic.map(cat.Teff)
    acc["Rstar"] = acc.tic.map(cat.Rstar)
    acc["L_tess_erg_s"] = [F.l_tess(t, r) for t, r in zip(acc.Teff, acc.Rstar)]
    acc["E_tess_erg"] = acc.ed_s * acc.L_tess_erg_s
    acc["E_tess_err_erg"] = acc.ed_err_s * acc.L_tess_erg_s
    acc["E_bol_erg_9000K"] = acc.E_tess_erg * bc
    cols = ["cand_id", "tic", "sector", "camera", "t_start", "t_peak", "t_stop", "amp", "amp_sig", "dur_min",
            "rise_min", "decay_min", "ed_s", "ed_err_s", "Teff", "Rstar", "L_tess_erg_s", "E_tess_erg",
            "E_tess_err_erg", "E_bol_erg_9000K", "n_flag512", "n_flag_cr", "reason"]
    acc[cols].to_csv(F.SQ / "flares.csv", index=False)

    sec = pd.read_csv(F.SQ / "sectors.csv")
    cov = pd.read_csv(F.SQ / "coverage.csv")
    ir = pd.read_csv(F.SQ / "injrec_per_star.csv").set_index("tic")
    rows = []
    for _, r in cov[cov.has_2min].iterrows():
        tic = r.TIC
        s = sec[sec.tic == tic]
        n = int((acc.tic == tic).sum())
        nm = int(((c.tic == tic) & (c.decision == "marginal")).sum())
        T = float(s.obs_days.sum())
        L = F.l_tess(float(r.Teff), float(r.Rstar))
        sig = float(np.median(s.sigma_ppm)) * 1e-6
        ed_min = F.N1 * sig * F.N3 * 120.0  # 3 cadences at exactly N1 sigma
        rows.append(dict(
            TIC=tic, Teff=r.Teff, Rstar=r.Rstar, Tmag=r.Tmag,
            sectors_2min=r.sectors_2min, n_sectors=len(s), obs_days=round(T, 2),
            sigma_ppm_median=round(sig * 1e6, 0),
            n_flares=n, n_marginal=nm,
            rate_per_day=round(n / T, 4),
            rate_ul95_per_day=round(POISSON_UL95[n] / T, 4),
            L_tess_erg_s=L,
            ED_min_analytic_s=round(ed_min, 3), E_min_analytic_erg=ed_min * L,
            ED50_s=ir.loc[tic, "ed50_s"] if tic in ir.index else np.nan,
            ED90_s=ir.loc[tic, "ed90_s"] if tic in ir.index else np.nan,
            E50_erg=ir.loc[tic, "E50_erg"] if tic in ir.index else np.nan,
            E90_erg=ir.loc[tic, "E90_erg"] if tic in ir.index else np.nan,
            flare_E_tess_erg=";".join(f"{e:.2e}" for e in acc[acc.tic == tic].E_tess_erg),
        ))
    per = pd.DataFrame(rows)
    per.to_csv(F.SQ / "per_star.csv", index=False)
    pd.set_option("display.width", 250)
    print(c[["cand_id", "tic", "sector", "amp_sig", "ed_s", "auto_result", "decision"]].to_string())
    print(acc[["cand_id", "tic", "ed_s", "E_tess_erg", "E_bol_erg_9000K"]].to_string())
    print(per.drop(columns=["sectors_2min"]).describe().T.to_string())
    print("total obs days:", per.obs_days.sum())

    # sample-level limit: stars complete (>=90%) above E_th contribute their time
    srows = []
    for lth in (33.5, 34.0, 34.5, 35.0, 35.5, 36.0):
        eth = 10 ** lth
        ok = per[per.E90_erg <= eth]
        nfl = int(sum((acc.tic.isin(ok.TIC)) & (acc.E_tess_erg >= eth)))
        T = float(ok.obs_days.sum())
        srows.append(dict(log10_E_tess_threshold=lth, n_hosts_complete=len(ok), star_days=round(T, 1),
                          n_flares_above=nfl,
                          rate_ul95_per_star_day=POISSON_UL95[nfl] / T if T > 0 else np.nan,
                          rate_ul95_per_star_year=365.25 * POISSON_UL95[nfl] / T if T > 0 else np.nan))
    sl = pd.DataFrame(srows)
    sl.to_csv(F.SQ / "sample_limits.csv", index=False)
    print(sl.to_string())


if __name__ == "__main__":
    main()
