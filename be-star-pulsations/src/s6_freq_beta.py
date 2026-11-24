"""Stage 6 (follow-up, post hoc): which frequencies grow with the disk?

Within-star slope beta(f) of ln(amplitude) against disk_excess, frequency by frequency,
from the saved per-sector amplitude spectra (600-point log grid, 0.25-24 c/d, rms of A in
each log bin). The grid is merged into N_BIN coarser log bins (rms). Two axes:

  absolute:  bins in f (c/d);
  relative:  bins in f / f1_star, where f1_star is the star's median f1 over its usable
             sectors. Bin 1 is then the star's own dominant pulsation group, so this
             separates the coherent peak from the continuum around it.

Each bin gets the stage-4 test (within-star permutation and cross-star disk-swap nulls,
bootstrap over stars). Variants: main (crowd-cut), no crowding cut, SPOC/TESS-SPOC only.
Outputs: results/<run>/freq_beta.csv, figures/<run>_freq_beta.png.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import PILOT, SPEC_DIR, SEED, figpath
from metrics import LOG_GRID
import s4_test as T, stats4

N_BIN = 30
N_NULL, N_BOOT = 10000, 2000
VARIANTS = ["main (crowd-cut)", "no crowding cut", "SPOC/TESS-SPOC stars only"]


def load_spectra(sm):
    S = np.full((len(sm), len(LOG_GRID)), np.nan)
    for i, (tic, sec) in enumerate(zip(sm.tic, sm.sector)):
        p = SPEC_DIR / f"{tic}_{sec}.npy"
        if p.exists():
            S[i] = np.load(p).astype(float)
    return S


def rebin(S, lf, edges):
    """rms of S over grid points whose log10 frequency falls in each edge bin (per row)."""
    out = np.full((S.shape[0], len(edges) - 1), np.nan)
    for b in range(len(edges) - 1):
        m = (lf >= edges[b]) & (lf < edges[b + 1])
        if np.ndim(m) == 1:
            if m.any():
                out[:, b] = np.sqrt(np.nanmean(S[:, m] ** 2, axis=1))
        else:
            m = m & np.isfinite(S)
            n = m.sum(1)
            s = np.sum(np.where(m, S, 0.0) ** 2, axis=1)
            ok = n >= 3
            out[ok, b] = np.sqrt(s[ok] / n[ok])
    return out


def main():
    sm, sw, disks, dos = T.load()
    S = load_spectra(sm)
    print("spectra loaded:", int(np.isfinite(S[:, 0]).sum()), "of", len(sm), flush=True)
    sid, t, x = sm.tic.values, sm.mid_mjd.values, sm.disk_excess.values
    lab, wu = sm.phase_label.values.astype(object), sm.wise_usable.values.astype(bool)
    lg = np.log10(LOG_GRID)
    # absolute bins
    e_abs = np.linspace(lg[0], lg[-1] + 1e-9, N_BIN + 1)
    Y_abs = rebin(S, lg, e_abs)
    # relative bins: log10(f / f1_star), f1_star = median f1 over the star's usable sectors
    use = sm.usable.astype(bool).values
    f1s = sm[use].groupby("tic").f1.median()
    lrel = lg[None, :] - np.log10(sm.tic.map(f1s).values)[:, None]
    e_rel = np.linspace(-1.2, 1.2, 25)
    Y_rel = rebin(S, lrel, e_rel)
    table = stats4.build_swap_table(t, disks, T.evaluate)
    vs = T.variants(sm)
    rows = []
    for axis, Y, edges in (("absolute", Y_abs, e_abs), ("relative", Y_rel, e_rel)):
        for vname in VARIANTS:
            sel = vs[vname]
            for b in range(Y.shape[1]):
                y = np.where(sel, np.log(Y[:, b]), np.nan)
                r = stats4.run_test(sid, x, y, lab, wu, t, disks, dos, T.evaluate,
                                    n_null=N_NULL, n_boot=N_BOOT, seed=SEED,
                                    swap_table=table, keep_nulls=True)
                sn = r.get("_swap_null", np.array([np.nan]))
                sn = sn[np.isfinite(sn)]
                rows.append(dict(axis=axis, variant=vname, bin=b,
                                 lo=10 ** edges[b], hi=10 ** edges[b + 1],
                                 centre=10 ** (0.5 * (edges[b] + edges[b + 1])),
                                 n_stars=r.get("n_stars"), n_sectors=r.get("n_sectors"),
                                 beta=r.get("beta"), beta_lo=r.get("beta_lo"),
                                 beta_hi=r.get("beta_hi"), p_perm=r.get("p_perm"),
                                 p_swap=r.get("p_swap"),
                                 swap_q025=np.quantile(sn, 0.025) if len(sn) else np.nan,
                                 swap_q975=np.quantile(sn, 0.975) if len(sn) else np.nan))
                print(axis, vname, b, round(rows[-1]["centre"], 3), rows[-1]["n_stars"],
                      None if rows[-1]["beta"] is None else round(rows[-1]["beta"], 3),
                      rows[-1]["p_swap"], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(PILOT / "freq_beta.csv", index=False)
    figure(out)


def figure(out):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    cols = {"main (crowd-cut)": "C0", "no crowding cut": "0.5", "SPOC/TESS-SPOC stars only": "C3"}
    for ax, axis in zip(axs, ("absolute", "relative")):
        for v, c in cols.items():
            d = out[(out.axis == axis) & (out.variant == v)]
            ls = "--" if v == "no crowding cut" else "-"
            ax.plot(d.centre, d.beta, ls, c=c, lw=1.6, label=v)
            if v != "no crowding cut":
                ax.fill_between(d.centre, d.beta_lo, d.beta_hi, color=c, alpha=0.2, lw=0)
        d = out[(out.axis == axis) & (out.variant == "main (crowd-cut)")]
        ax.fill_between(d.centre, d.swap_q025, d.swap_q975, color="k", alpha=0.12, lw=0,
                        label="cross-star null 95% (main)")
        ax.axhline(0, c="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_ylabel("beta: d ln A / d disk_excess (per mag)")
    axs[0].set_xlabel("frequency (c/d)")
    axs[0].set_title("absolute frequency", fontsize=10)
    axs[1].set_xlabel("f / f1 of the star (median over its sectors)")
    axs[1].axvline(1, c="k", ls=":", lw=0.8)
    axs[1].set_title("relative to the star's dominant frequency", fontsize=10)
    axs[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(figpath("freq_beta"), dpi=120)


if __name__ == "__main__":
    main()
