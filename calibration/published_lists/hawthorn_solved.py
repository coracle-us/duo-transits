"""
Do the published Hawthorn et al. (2024) alias limits keep the periods that were later solved?

Phase 3, first question of docs/STRATEGY_living_alias_catalogue.md: our own former epoch rule excluded
the true period for most confirmed planets. Does a published, MonoTools-based alias list do the same?

Exploratory, not pre-registered: the solved periods were already in data/catalog/status.csv and had
been looked at before this script was written. Every choice is stated here.

Truth  : data/catalog/status.csv, rows with status == 'resolved' (26) and a period_claim_d.
Paper  : data/catalog/candidates.csv (Tables A2, A3): Tc1, Tc2, Pmin, Pmax, N_alias, P_marg and its errors.
For each solved candidate:
  n_true      = round(dT / P_true); alias_resid = |dT / n_true - P_true|
  n_max_paper = round(dT / Pmin)
  kept        = n_true <= n_max_paper              (the true period is not below the paper's Pmin)
  contiguous  = N_alias in (n_max_paper, n_max_paper + 1): the paper's list then has no gaps between
                Pmin and Pmax, so 'kept' means the true alias is in the list. Otherwise the full
                online list is needed to say so.
  in_pmarg    = P_marg - errm <= P_true <= P_marg + errp
  rank        = n_max_paper - n_true + 1           (1 = the shortest allowed period)
Period priors: for weights w_n proportional to P_n**(-a) over n = 1..n_max_paper, the probability of the
true alias and the mid-PIT (weight on shorter periods plus half the true alias's weight). a = 8/3 is the
exponent I recall MonoTools using (check Osborn 2022 before citing); 5/3 and 2/3 are priors A and B of
arXiv 2611.03019; 0 is flat.

Outputs: hawthorn_solved.csv, hawthorn_solved_summary.txt (this folder)
Usage  : python calibration/published_lists/hawthorn_solved.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXPONENTS = {"8/3": 8 / 3, "5/3": 5 / 3, "2/3": 2 / 3, "flat": 0.0}


def alias_weights(dt, n_max, a):
    n = np.arange(1, n_max + 1)
    w = (dt / n) ** (-a)
    return n, w / w.sum()


def analyse(row):
    dt = row.Tc2 - row.Tc1
    n_true = int(round(dt / row.period_claim_d))
    n_max = int(round(dt / row.Pmin))
    out = dict(TIC=int(row.TIC), tier=row.tier, P_true=row.period_claim_d, dT=dt, n_true=n_true,
               alias_resid_d=abs(dt / n_true - row.period_claim_d), Pmin_paper=row.Pmin,
               n_max_paper=n_max, N_alias_paper=int(row.N_alias), kept=n_true <= n_max,
               contiguous=int(row.N_alias) in (n_max, n_max + 1), rank=n_max - n_true + 1,
               P_marg=row.P_marg,
               in_pmarg=bool(row.P_marg - row.P_marg_errm <= row.period_claim_d
                             <= row.P_marg + row.P_marg_errp))
    for name, a in EXPONENTS.items():
        n, w = alias_weights(dt, n_max, a)
        if n_true <= n_max:
            out[f"p_true_{name}"] = w[n_true - 1]
            out[f"pit_{name}"] = w[n > n_true].sum() + 0.5 * w[n_true - 1]   # shorter periods = larger n
        else:
            out[f"p_true_{name}"] = 0.0
            out[f"pit_{name}"] = np.nan
    return out


def main():
    s = pd.read_csv(ROOT / "data/catalog/status.csv")
    c = pd.read_csv(ROOT / "data/catalog/candidates.csv")
    m = s[(s.status == "resolved") & s.period_claim_d.notna()].merge(c, on="TIC")
    r = pd.DataFrame([analyse(x) for x in m.itertuples()])
    r.to_csv(HERE / "hawthorn_solved.csv", index=False)
    lines = [
        f"solved Hawthorn duotransits: {len(r)}",
        f"claimed period is an alias of dT to within 0.05 d: {(r.alias_resid_d < 0.05).sum()} "
        f"(largest residual {r.alias_resid_d.max():.3f} d, TIC {int(r.loc[r.alias_resid_d.idxmax(), 'TIC'])})",
        f"true period not below the paper's Pmin (kept): {int(r.kept.sum())} of {len(r)}",
        f"paper's list contiguous (N_alias = n_max or n_max + 1): {int(r.contiguous.sum())} of {len(r)}",
        f"kept and contiguous, so the true alias is in the published list: "
        f"{int((r.kept & r.contiguous).sum())} of {len(r)}",
        f"not contiguous (full online list needed): TIC {r.TIC[~r.contiguous].tolist()}",
        f"true period inside P_marg's quoted errors: {int(r.in_pmarg.sum())} of {len(r)}",
        f"true period is the shortest allowed (rank 1): {int((r['rank'] == 1).sum())}; "
        f"rank <= 3: {int((r['rank'] <= 3).sum())}; median rank {r['rank'].median():.0f} "
        f"of a median {r.n_max_paper.median():.0f} aliases",
        "period prior exponent a (w ~ P^-a): sum of ln p(true alias), mean mid-PIT, KS p against uniform",
    ]
    for name in EXPONENTS:
        k = r[r.kept]
        ks = stats.kstest(k[f"pit_{name}"], "uniform")
        lines.append(f"  a = {name:>4}: sum ln p = {np.log(k[f'p_true_{name}']).sum():7.2f}; "
                     f"mean PIT = {k[f'pit_{name}'].mean():.2f}; KS p = {ks.pvalue:.3f}")
    (HERE / "hawthorn_solved_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(r[["TIC", "tier", "P_true", "n_true", "n_max_paper", "N_alias_paper", "kept", "contiguous",
             "rank", "P_marg", "in_pmarg"]].to_string(index=False))


if __name__ == "__main__":
    main()
