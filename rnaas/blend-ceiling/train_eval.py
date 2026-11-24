"""
Steps 4-5 of the blend-ceiling test: cross-fitted gradient boosting (RULES.md amendment A3/A4), the AUCs with
bootstrap CIs (section 4), the confounder and sanity checks, and the figure.

Input:  data/features.csv (make_windows.py), data/blends.csv (tag_blends.py), data/sample.csv, data/eb_sample.csv
Output: results/window_scores.csv, results/object_scores.csv, results/auc.txt, fig1.pdf, fig1.png
Usage:  .venv/bin/python -I train_eval.py [--cnn]
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir; our modules live here

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve

import lcwin

HERE = os.path.dirname(os.path.abspath(__file__))
D, R = os.path.join(HERE, "data"), os.path.join(HERE, "results")
SEED, NFOLD, NBOOT = 20261124, 5, 2000
TRAIN_KINDS = ["pos", "eb", "tce"]  # FP/FA never trained on (amendment A2)


def folds(tics):
    u = np.array(sorted(set(tics)))
    perm = np.random.default_rng(SEED).permutation(u)
    return {int(t): i % NFOLD for i, t in enumerate(perm)}


def balanced_weights(y):
    w = np.ones(len(y))
    for c in (0, 1):
        w[y == c] = len(y) / (2 * max((y == c).sum(), 1))
    return w


def crossfit_gbm(F, cols):
    s = np.full(len(F), np.nan)
    for k in range(NFOLD):
        tr = (F.fold != k) & F.kind.isin(TRAIN_KINDS)
        te = F.fold == k
        y = F.label[tr].values
        if len(set(y)) < 2:
            print(f"fold {k}: one class only in training; not scored (A6)")
            continue
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=SEED)
        m.fit(F.loc[tr, cols].values, y, sample_weight=balanced_weights(y))
        s[te.values] = m.predict_proba(F.loc[te, cols].values)[:, 1]
    return s


def crossfit_cnn(F, X):
    import torch
    import torch.nn as nn
    torch.manual_seed(SEED)
    so = F.sig_out.values[:, None]
    xin = np.stack([np.broadcast_to(lcwin.T_REL_H / (lcwin.WIN_H / 2), X.shape),
                    np.clip(X / np.where(so > 0, so, 1), -50, 50)], axis=1).astype(np.float32)
    s = np.full(len(F), np.nan)
    for k in range(NFOLD):
        torch.manual_seed(SEED + k)
        net = nn.Sequential(nn.Conv1d(2, 16, 7, padding=3), nn.ReLU(), nn.MaxPool1d(2),
                            nn.Conv1d(16, 32, 5, padding=2), nn.ReLU(), nn.MaxPool1d(2),
                            nn.Conv1d(32, 64, 5, padding=2), nn.ReLU(), nn.AdaptiveAvgPool1d(1),
                            nn.Flatten(), nn.Linear(64, 1))
        tr = ((F.fold != k) & F.kind.isin(TRAIN_KINDS)).values
        te = (F.fold == k).values
        xt, yt = torch.tensor(xin[tr]), torch.tensor(F.label.values[tr], dtype=torch.float32)
        wt = torch.tensor(balanced_weights(F.label.values[tr]), dtype=torch.float32)
        opt = torch.optim.Adam(net.parameters(), lr=1e-3)
        g = torch.Generator().manual_seed(SEED + k)
        for ep in range(30):
            net.train()
            for idx in torch.randperm(len(xt), generator=g).split(64):
                opt.zero_grad()
                loss = nn.functional.binary_cross_entropy_with_logits(net(xt[idx]).squeeze(1), yt[idx], weight=wt[idx])
                loss.backward()
                opt.step()
        net.eval()
        with torch.no_grad():
            s[te] = torch.sigmoid(net(torch.tensor(xin[te])).squeeze(1)).numpy()
    return s


def auc_ci(pos, neg, rng):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) < 2 or len(neg) < 2:
        return np.nan, np.nan, np.nan
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    a = roc_auc_score(y, np.r_[pos, neg])
    b = [roc_auc_score(y, np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]) for _ in range(NBOOT)]
    return a, np.percentile(b, 2.5), np.percentile(b, 97.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnn", action="store_true")
    args = ap.parse_args()
    os.makedirs(R, exist_ok=True)
    F = pd.read_csv(os.path.join(D, "features.csv"))
    X = np.load(os.path.join(D, "windows_flux.npz"))["x"]
    fmap = folds(F.tic)
    F["fold"] = F.tic.map(fmap)
    F[lcwin.FEATURES] = F[lcwin.FEATURES].astype(float)
    F["s_gbm"] = crossfit_gbm(F, lcwin.FEATURES)
    F["s_shape"] = crossfit_gbm(F, [c for c in lcwin.FEATURES if c not in lcwin.DEPTH_SCALE])
    models = ["s_gbm", "s_shape"]
    if args.cnn:
        F["s_cnn"] = crossfit_cnn(F, X)
        models.append("s_cnn")
    F.to_csv(os.path.join(R, "window_scores.csv"), index=False)

    # per object: max (main), mean, and first window by time (RULES section 2)
    F = F.sort_values(["obj_id", "t_c"])
    agg = {f"{m}_{a}": (m, f) for m in models for a, f in (("max", "max"), ("mean", "mean"), ("first", "first"))}
    O = F.groupby(["obj_id", "kind", "tic"], as_index=False).agg(n_win=("t_c", "size"), **agg)
    s = pd.read_csv(os.path.join(D, "sample.csv"))
    s["obj_id"] = "TOI-" + s.toi.astype(str)
    bl = pd.read_csv(os.path.join(D, "blends.csv"))
    bl["obj_id"] = "TOI-" + bl.toi.astype(str)
    O = O.merge(s[["obj_id", "disposition", "depth_ppm", "tmag", "period_d"]], on="obj_id", how="left")
    O = O.merge(bl[["obj_id", "blend_C", "blend_G", "blend_CG"]], on="obj_id", how="left")
    O.to_csv(os.path.join(R, "object_scores.csv"), index=False)

    rng = np.random.default_rng(SEED)
    pos, fp = O[O.kind == "pos"], O[O.kind == "fpfa"]
    eb, tce, san = O[O.kind == "eb"], O[O.kind == "tce"], O[O.kind == "sanity"]
    out = [f"Scored objects: {len(pos)} CP/KP, {len(fp)} FP/FA ({int(fp.blend_C.sum())} blends by C, "
           f"{int(fp.blend_CG.sum())} by C or G), {len(eb)} EBs, {len(tce)} TCE-like, {len(san)} sanity windows",
           f"Unscored (fold with one class): {int(O.s_gbm_max.isna().sum())}", ""]
    comps = [("all FP/FA", fp), ("non-blend FP/FA (C)", fp[~fp.blend_C.astype(bool)]),
             ("blend FP/FA (C)", fp[fp.blend_C.astype(bool)]),
             ("non-blend FP/FA (C or G)", fp[~fp.blend_CG.astype(bool)]),
             ("blend FP/FA (C or G)", fp[fp.blend_CG.astype(bool)]),
             ("catalogued EBs (Prsa)", eb), ("TCE-like windows", tce)]
    res = {}
    for m in models:
        for agg_ in ("max", "mean", "first"):
            col = f"{m}_{agg_}"
            out.append(f"== {m} ({agg_} over windows)")
            for name, neg in comps:
                a, lo, hi = auc_ci(pos[col].dropna(), neg[col].dropna(), rng)
                res[(col, name)] = (a, lo, hi, len(neg))
                out.append(f"  CP/KP ({pos[col].notna().sum()}) vs {name} ({neg[col].notna().sum()}): "
                           f"AUC {a:.3f} [{lo:.3f}, {hi:.3f}]")
    out.append("\n== depth bins (table depth), s_gbm_max, CP/KP vs all FP/FA")
    edges = [0, 1000, 3000, 10000, np.inf]
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        p = pos[(pos.depth_ppm >= lo_) & (pos.depth_ppm < hi_)].s_gbm_max.dropna()
        n = fp[(fp.depth_ppm >= lo_) & (fp.depth_ppm < hi_)].s_gbm_max.dropna()
        a, l, h = auc_ci(p, n, rng)
        flag = "" if len(p) >= 10 and len(n) >= 10 else "  (fewer than 10 of a class: not interpreted)"
        out.append(f"  {lo_:>6.0f}-{hi_:<6.0f} ppm: {len(p)} pos, {len(n)} neg, AUC {a:.3f} [{l:.3f}, {h:.3f}]{flag}")
    out.append("\n== depth and Tmag of the scored CP/KP and FP/FA")
    for c in ("depth_ppm", "tmag"):
        a_, b_ = pos[c].dropna(), fp[c].dropna()
        ks = ks_2samp(a_, b_)
        out.append(f"  {c}: CP/KP median {a_.median():.1f} (IQR {a_.quantile(.25):.1f}-{a_.quantile(.75):.1f}); "
                   f"FP/FA median {b_.median():.1f} (IQR {b_.quantile(.25):.1f}-{b_.quantile(.75):.1f}); "
                   f"KS D = {ks.statistic:.3f}, p = {ks.pvalue:.2g}")
    out.append("\n== sanity windows (out of transit, positives' light curves), s_gbm window scores")
    q = lambda v: f"median {np.median(v):.3f}, 10-90% {np.percentile(v, 10):.3f}-{np.percentile(v, 90):.3f}"
    out.append(f"  sanity ({len(san)}): {q(san.s_gbm_max.dropna())}")
    out.append(f"  FP/FA objects: {q(fp.s_gbm_max.dropna())}")
    out.append(f"  TCE-like: {q(tce.s_gbm_max.dropna())}")
    out.append(f"  CP/KP objects: {q(pos.s_gbm_max.dropna())}")
    out.append("\nSalinas et al. overlap: their CNN was trained on synthetic light curves only and its 120 real "
               "validation TOIs are not listed in the paper, so the overlap could not be checked (RULES section 4).")
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(R, "auc.txt"), "w").write(txt + "\n")
    make_fig(pos, fp, res)


def make_fig(pos, fp, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from tag_blends import wilson
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.15]))
    y = lambda p, n: (np.r_[np.ones(len(p)), np.zeros(len(n))], np.r_[p, n])
    sets = [("all FP/FA", fp, "k", "-"), ("non-blend (C)", fp[~fp.blend_C.astype(bool)], "C0", "--"),
            ("blend (C)", fp[fp.blend_C.astype(bool)], "C3", ":")]
    for name, n, c, ls in sets:
        if len(n) < 2:
            continue
        yy, ss = y(pos.s_gbm_max.dropna(), n.s_gbm_max.dropna())
        fpr, tpr, _ = roc_curve(yy, ss)
        a1.plot(fpr, tpr, color=c, ls=ls, lw=1.4, label=f"{name}, AUC {roc_auc_score(yy, ss):.2f}")
    a1.plot([0, 1], [0, 1], color="0.6", lw=0.8)
    a1.set(xlabel="False-positive rate", ylabel="True-positive rate", xlim=(0, 1), ylim=(0, 1))
    a1.legend(fontsize=7, loc="lower right", frameon=False)
    names = ["all FP/FA", "non-blend FP/FA (C)", "blend FP/FA (C)", "catalogued EBs (Prsa)"]
    lab = ["all FP/FA", "non-blend (C)", "blend (C)", "Prša EBs"]
    for i, nm in enumerate(names):
        a, lo, hi, _ = res[("s_gbm_max", nm)]
        a2.errorbar(i, a, yerr=[[a - lo], [hi - a]], fmt="o", color="k", ms=4, capsize=3)
    for col, c, nm in (("blend_C", "C3", "ceiling, comments"), ("blend_CG", "C1", "ceiling, comments or Gaia")):
        f = fp[col].astype(bool).mean()
        _, flo, fhi = wilson(fp[col].astype(bool).sum(), len(fp))
        a2.axhline(1 - f / 2, color=c, lw=1.2, label=f"{nm}: {1 - f / 2:.2f}")
        a2.axhspan(1 - fhi / 2, 1 - flo / 2, color=c, alpha=0.12, lw=0)
    a2.axhline(0.5, color="0.6", lw=0.8, ls="--", label="chance")
    a2.set_xticks(range(len(lab)), lab, fontsize=8)
    a2.set(ylabel="AUC vs CP/KP (95% CI)", ylim=(0.3, 1.02), xlim=(-0.5, len(lab) - 0.5))
    a2.legend(fontsize=7, loc="lower left", frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"fig1.{ext}"), dpi=200)


if __name__ == "__main__":
    main()
