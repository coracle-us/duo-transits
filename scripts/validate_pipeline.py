"""
Plan step 4: validate the alias pipeline (scripts/alias_pipeline.py) before trusting it.

A. Known periods. For every candidate in data/catalog/status.csv with status 'resolved' and
   a claimed alias n, run the pipeline on ALL TESS data and check:
     PASS      true alias not excluded (verdict 'supported' or 'allowed')
     FAIL      true alias 'excluded' or 'conflict'
   and report how many aliases survive and whether the true one is uniquely recovered.
   Tiers are kept separate: peer-reviewed / catalog-paper / thesis periods are the ground
   truth; exofop-community and toi-catalog periods are claims being tested at the same time.

B. Paper comparison. Run all 85 candidates on paper-era data only (Sectors <= 39, plus each
   candidate's own Tc sectors) and compare the surviving aliases with Hawthorn+24 Table A3.
   Note: in Table A3, P_min = dT/(N_alias - 1) for most rows, i.e. P_min and N_alias are
   off by one relative to each other. Our survivor sets match n = 1..N_alias, so we compare
   against N_alias read as the count of allowed aliases.

Outputs: results/validation_known.csv, results/validation_paper.csv (+ per-TIC files under
results/aliases/all and results/aliases/maxs39).

Usage:  python scripts/validate_pipeline.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

import alias_pipeline as ap

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "catalog"
RES = ROOT / "results"


def run(c, max_sector, label):
    lc = ap.load_lightcurves(int(c.TIC), max_sector=None if max_sector is None
                             else max(max_sector, int(c.S2)))
    ep, al = ap.run_candidate(c, lc)
    out = RES / "aliases" / label
    out.mkdir(parents=True, exist_ok=True)
    ep.to_csv(out / f"{c.TIC}_epochs.csv", index=False)
    al.to_csv(out / f"{c.TIC}_aliases.csv", index=False)
    return ep, al


def part_a(cand, status):
    truth = status[(status.status == "resolved") & status.n_alias.notna()]
    rows = []
    for _, s in truth.iterrows():
        c = cand.loc[s.TIC]
        ep, al = run(c, None, "all")
        n_true = int(s.n_alias)
        v = al.set_index("n").verdict
        survivors = al[al.verdict.isin(["supported", "allowed"])].n.tolist()
        known = ep[ep.known & (ep.n == 1)].verdict.tolist() if len(ep) else []
        tv = v.get(n_true, "not-tested")
        rows.append(dict(
            TIC=s.TIC, tier=s.tier, period=s.period_claim_d, n_true=n_true, verdict_true=tv,
            result="PASS" if tv in ("supported", "allowed") else "FAIL",
            n_survivors=len(survivors), unique=(survivors == [n_true]),
            survivors=" ".join(map(str, survivors[:30])) + (" ..." if len(survivors) > 30 else ""),
            true_supported_at=al.set_index("n").supported_at.get(n_true, ""),
            true_excluded_by=al.set_index("n").excluded_by.get(n_true, ""),
            true_misses_corroborated=al.set_index("n").misses_corroborated.get(n_true, ""),
            uncorroborated_dips=al.uncorroborated_dips.iloc[0] if len(al) else "",
            known_transits=" ".join(known), flags=s.flags))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "validation_known.csv", index=False)
    return out


def part_b(cand):
    rows = []
    for _, c in cand.iterrows():
        ep, al = run(c, 39, "maxs39")
        surv = al[al.verdict.isin(["supported", "allowed"])]
        ours = set(surv.n)
        paper = set(range(1, int(c.N_alias) + 1))
        known = ep[ep.known & (ep.n == 1)].verdict.tolist() if len(ep) else []
        rows.append(dict(
            TIC=c.TIC, N_alias_paper=int(c.N_alias), our_n_survivors=len(ours),
            dN=len(ours) - int(c.N_alias), same_set=(ours == paper),
            ours_not_paper=" ".join(map(str, sorted(ours - paper))),
            paper_not_ours=" ".join(map(str, sorted(paper - ours))),
            our_Pmin=surv.P.min() if len(surv) else np.nan, Pmin_paper=c.Pmin,
            depth_ppt=np.nanmean([c.Depth1_ppt, c.Depth2_ppt]), b=c.b,
            known_transits=" ".join(known)))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "validation_paper.csv", index=False)
    return out


def main():
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC", drop=False)
    status = pd.read_csv(CAT / "status.csv")
    pd.set_option("display.width", 250, "display.max_colwidth", 60, "display.max_rows", 200)

    a = part_a(cand, status)
    print("=== A. known periods (all TESS data) ===")
    print(a[["TIC", "tier", "period", "n_true", "verdict_true", "result", "n_survivors", "unique",
             "known_transits"]].to_string(index=False))
    print(a.groupby("tier").result.value_counts().unstack(fill_value=0))

    b = part_b(cand)
    print("\n=== B. paper-era data (S<=39) vs Hawthorn+24 Table A3 ===")
    # The paper publishes only the count; its sets need not be contiguous, so compare counts.
    print("our count - N_alias:", b.dN.value_counts().sort_index().to_dict())
    print(f"count equal: {(b.dN == 0).sum()}, within +1: {(b.dN.between(0, 1)).sum()}, "
          f"within +3: {(b.dN.between(0, 3)).sum()}, fewer than paper: {(b.dN < 0).sum()} (of {len(b)})")
    print(b[(b.dN > 3) | (b.dN < 0)][["TIC", "N_alias_paper", "our_n_survivors", "depth_ppt", "b",
                                      "known_transits"]].to_string(index=False))
    print("known transits recovered as 'supported':",
          b.known_transits.str.count("supported").sum(), "of", 2 * len(b))


if __name__ == "__main__":
    main()
