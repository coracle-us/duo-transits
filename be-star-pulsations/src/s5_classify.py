"""Stage 5: unsupervised sorting of amplitude-spectrum shapes (PCA + Gaussian mixture)."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from common import PILOT, FIGURES, SPEC_DIR, SEED, figpath
from metrics import LOG_GRID


def features(sm):
    X = []
    for _, s in sm.iterrows():
        a = np.load(SPEC_DIR / f"{s.tic}_{s.sector}.npy").astype(float)
        a = np.log10(np.clip(a / a.max(), 1e-4, None))
        a = np.convolve(np.r_[a[0], a, a[-1]], np.ones(3) / 3, mode="valid")
        X.append(a)
    return np.array(X)


def chi2(tab):
    e = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
    m = e > 0
    return float(np.sum((tab[m] - e[m]) ** 2 / e[m]))


def table(cl, lab, K, labels):
    t = np.zeros((K, len(labels)))
    for c, l in zip(cl, lab):
        t[c, labels.index(l)] += 1
    return t


def main():
    sm = pd.read_csv(PILOT / "sector_metrics.csv")
    sm = sm[sm.metrics_ok.astype(bool) & sm.wise_usable.astype(bool)].reset_index(drop=True)
    X = features(sm)
    Z = PCA(n_components=10, random_state=SEED).fit_transform(X)
    fits = {k: GaussianMixture(k, covariance_type="full", n_init=5, random_state=SEED).fit(Z)
            for k in range(3, 7)}
    bic = {k: g.bic(Z) for k, g in fits.items()}
    K = min(bic, key=bic.get)
    cl = fits[K].predict(Z)
    sm["cluster"] = cl
    sm[["tic", "sector", "cluster", "phase_label", "disk_excess", "A1", "f1", "n_peaks",
        "red_noise"]].to_csv(PILOT / "clusters.csv", index=False)
    print("BIC", {k: round(v, 1) for k, v in bic.items()}, "-> K =", K)
    print("cluster sizes", np.bincount(cl).tolist())
    for c in range(K):
        g = sm[sm.cluster == c]
        print(f"cluster {c}: n={len(g)}, median n_peaks {g.n_peaks.median():.0f}, "
              f"median snr1 {sm.loc[g.index, 'snr1'].median():.1f}, median f1 {g.f1.median():.2f}")
    # within-star changes and association with phase label
    multi = sm.groupby("tic").filter(lambda g: len(g) >= 2)
    changes = multi.groupby("tic").cluster.nunique()
    print(f"stars with >=2 usable sectors: {len(changes)}; cluster changes within star: "
          f"{int((changes > 1).sum())}")
    labels = sorted(sm.phase_label.unique())
    lab = multi.phase_label.values
    clm = multi.cluster.values
    obs = chi2(table(clm, lab, K, labels))
    rng = np.random.default_rng(SEED)
    import stats4
    from wise import DiskCurve
    from common import WISE_DIR
    lc = stats4.code(lab)
    mt = multi.tic.values
    ustars, gid = np.unique(mt, return_inverse=True)
    P = np.argsort(gid, kind="stable")

    def within_perm():
        Q = np.lexsort((rng.random(len(gid)), gid))
        out = np.empty_like(lc)
        out[P] = lc[Q]
        return out

    def chi2c(cl_, lc_, use_=None):
        m = (lc_ >= 0) if use_ is None else (use_ & (lc_ >= 0))
        t = np.bincount(cl_[m] * 4 + lc_[m], minlength=4 * K).reshape(K, 4).astype(float)
        return chi2(t[:, t.sum(0) > 0])

    null = np.array([chi2c(clm, within_perm()) for _ in range(10000)])
    p = (1 + np.sum(null >= obs)) / (1 + len(null))
    # Pair statistic (purely within-star): over pairs of sectors of the same star,
    # P(cluster differs | phase differs) - P(cluster differs | phase same).
    pi, pj = [], []
    for g in range(len(ustars)):
        r = np.where(gid == g)[0]
        a_, b_ = np.triu_indices(len(r), 1)
        pi.append(r[a_]); pj.append(r[b_])
    pi, pj = np.concatenate(pi), np.concatenate(pj)
    dcl = clm[pi] != clm[pj]

    def pair_stat(lc_, use_):
        v = use_[pi] & use_[pj]
        dl = lc_[pi] != lc_[pj]
        a_, b_ = v & dl, v & ~dl
        if not a_.any() or not b_.any():
            return np.nan
        return dcl[a_].mean() - dcl[b_].mean()

    allu = np.ones(len(lc), bool)
    obs_pair = pair_stat(lc, allu)
    pn = np.array([pair_stat(within_perm(), allu) for _ in range(2000)])
    # cross-star null: labels from another star's WISE curve at this star's sector times
    # (keeps sector-era systematics, e.g. the QLP 1800/600/200 s cadence eras)
    tics = sorted(sm.tic.unique())
    disks = []
    for t in tics:
        v = pd.read_csv(WISE_DIR / "visits" / f"{t}.csv")
        disks.append(DiskCurve(v.mjd.values, v.w1mpro_ep.values))
    own = np.array([tics.index(t) for t in ustars])

    def ev(d, tt):
        o = d.at(tt)
        return o["disk_excess"], o["phase_label"], o["wise_usable"]

    _, Lm, Um = stats4.build_swap_table(multi.mid_mjd.values, disks, ev)
    rows_ = np.arange(len(mt))
    sn = np.full(2000, np.nan)
    for i in range(len(sn)):
        donors = rng.integers(0, len(disks) - 1, size=len(ustars))
        donors += donors >= own
        dj = donors[gid]
        sn[i] = pair_stat(Lm[rows_, dj], Um[rows_, dj])

    def p2(o, n):
        n = n[np.isfinite(n)]; c = n.mean()
        return (1 + np.sum(np.abs(n - c) >= abs(o - c))) / (1 + len(n))
    p_pair_perm, p_pair_swap = p2(obs_pair, pn), p2(obs_pair, sn)
    print(f"pair statistic {obs_pair:.3f}: within-star perm p = {p_pair_perm:.4f} "
          f"(null mean {pn.mean():.3f}), cross-star p = {p_pair_swap:.4f} "
          f"(null mean {np.nanmean(sn):.3f})")
    tab = pd.DataFrame(table(clm, lab, K, labels).astype(int), columns=labels)
    tab.index.name = "cluster"
    tab.to_csv(PILOT / "cluster_vs_phase.csv")
    print(tab.to_string())
    print(f"within-star permutation chi2 {obs:.2f}, p = {p:.4f}")
    with open(PILOT / "classifier_summary.txt", "w") as fh:
        fh.write(f"K={K} BIC={bic}\nsizes={np.bincount(cl).tolist()}\n"
                 f"stars_multi={len(changes)} changes={int((changes > 1).sum())}\n"
                 f"chi2={obs:.3f} p_perm={p:.4f}\n"
                 f"pair_stat={obs_pair:.4f} p_perm={p_pair_perm:.4f} p_swap={p_pair_swap:.4f}\n"
                 f"{tab.to_string()}\n")
    # figure
    fig, axs = plt.subplots(1, K, figsize=(3.2 * K, 3.4), sharey=True)
    for c, ax in enumerate(np.atleast_1d(axs)):
        idx = np.where(cl == c)[0]
        for j in rng.choice(idx, size=min(5, len(idx)), replace=False):
            ax.plot(LOG_GRID, X[j], lw=0.6, alpha=0.6)
        ax.plot(LOG_GRID, X[idx].mean(0), "k-", lw=2)
        ax.set_xscale("log")
        ax.set_title(f"cluster {c} (n={len(idx)})", fontsize=9)
        ax.set_xlabel("frequency (c/d)")
    np.atleast_1d(axs)[0].set_ylabel("log10 A / A_max")
    fig.tight_layout()
    fig.savefig(figpath("clusters"), dpi=110)


if __name__ == "__main__":
    main()
