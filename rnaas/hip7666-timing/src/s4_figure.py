"""Figure 1 (fig1.pdf/png) and the S58 check figure (check_S58.png).

Top row: eclipse O-C against our linear ephemeris (M0, the BIC choice) for S18, S58, S85, with
pulsation time delays (common tau per orbit, linear term per sector removed) on the same time axis.
Bottom: the 95% LTTE semi-amplitude limit vs P3, with the minimum tertiary mass on the right axis.
"""
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import hip

oc = pd.read_csv(hip.OUTDIR / "oc_main.csv")
pdl = pd.read_csv(hip.OUTDIR / "pulsation_delays.csv")
lim = pd.read_csv(hip.OUTDIR / "ltte_limit_main.csv")
limq = pd.read_csv(hip.OUTDIR / "ltte_limit_mixq18.csv")
pdl = pdl[pdl["mode"] == "all"]

fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1], hspace=0.38, wspace=0.08)
ylim = 70
for j, s in enumerate((18, 58, 85)):
    ax = fig.add_subplot(gs[0, j])
    for ty, mk, lab in (("p", "o", "primary eclipse"), ("s", "s", "secondary eclipse")):
        g = oc[(oc.sector == s) & (oc.type == ty)]
        ax.errorbar(g.t_btjd, g.oc_own_s, g.err_s, fmt=mk, ms=3.5, mfc="C0" if ty == "p" else "white",
                    color="C0", lw=0.8, label=lab)
    g = pdl[pdl.sector == s]
    ax.errorbar(g.t_mid, g.tau_detr, g.err_s, fmt="D", ms=3, color="C3", lw=0.8, alpha=0.8,
                label="pulsation delay (per orbit)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylim(-ylim, ylim)
    lo, hi = oc[oc.sector == s].t_btjd.min(), oc[oc.sector == s].t_btjd.max()
    ax.set_xlim(lo - 1.5, hi + 1.5)
    ax.set_title(f"Sector {s}", fontsize=9)
    ax.tick_params(labelsize=7)
    if j == 0:
        ax.set_ylabel("O−C or delay (s)", fontsize=8)
        ax.legend(fontsize=6, loc="lower left", frameon=False)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("BTJD", fontsize=8)

ax = fig.add_subplot(gs[1, :])
ax.loglog(lim.P3_d, lim.A95_s, "C0-", lw=1.5, label="excluded above (all main products)")
ax.loglog(limq.P3_d, limq.A95_s, "C0--", lw=1, label="with S18 from QLP")
ax.fill_between(lim.P3_d, lim.A95_s, 3e3, color="C0", alpha=0.12)
ax.axvline(oc.t_btjd.max() - oc.t_btjd.min(), color="0.5", ls=":", lw=1)
ax.text((oc.t_btjd.max() - oc.t_btjd.min()) * 0.93, 1200, "TESS baseline", ha="right", fontsize=7, color="0.4")
ax.set_xlabel("Tertiary period $P_3$ (d)", fontsize=8)
ax.set_ylabel("LTTE semi-amplitude (s)", fontsize=8)
ax.set_ylim(2, 3e3)
ax.tick_params(labelsize=7)
for p, a, m in lim[["P3_d", "A95_s", "m3min_msun"]].values[::20]:
    ax.annotate(f"{m:.2f}", (p, a), textcoords="offset points", xytext=(0, 5), fontsize=6, color="C0",
                ha="center")
ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.text(0.99, 0.04, "labels: minimum tertiary mass (M$_\\odot$) at the limit", transform=ax.transAxes,
        ha="right", fontsize=6.5, color="C0")
fig.savefig(hip.HERE / "fig1.pdf", bbox_inches="tight")
fig.savefig(hip.HERE / "fig1.png", dpi=200, bbox_inches="tight")

# ---- S58 check figure (brief step 2)
P = pickle.load(open(hip.OUTDIR / "prepared_main.pkl", "rb"))
TPL = pickle.load(open(hip.OUTDIR / "template_S58.pkl", "rb"))
d = P[58]
t = d["t"]
shape = TPL(hip.phase(t))
fig, ax = plt.subplots(4, 1, figsize=(8, 9))
ax[0].plot(t, d["f"], "k.", ms=0.6); ax[0].plot(t, d["trend"], "r-", lw=0.8)
ax[0].set_title("S58 SPOC 120 s PDCSAP, with the trend (red)", fontsize=9)
ph = hip.phase(t)
ax[1].plot(ph, d["y"], "k.", ms=0.6); ax[1].plot(np.sort(ph), TPL(np.sort(ph)), "r-", lw=0.8)
ax[1].set_title("pulsations subtracted, folded on Lee Eq. 1, with the binary template (red)", fontsize=9)
z = (t > 2890) & (t < 2892.5)
r = d["f"] / d["trend"] - shape
ax[2].plot(t[z], r[z] * 1e3, "k.", ms=1.5); ax[2].plot(t[z], d["puls"][z] * 1e3, "r-", lw=0.8)
ax[2].set_title("binary removed (ppt), one orbit, with the 14-frequency pulsation model (red)", fontsize=9)
res = d["y"] - shape
ax[3].plot(ph, res * 1e3, "k.", ms=0.6)
ax[3].set_title(f"residual after both (ppt); out-of-eclipse rms {np.std(res[d['oot']])*1e3:.2f} ppt", fontsize=9)
ax[3].set_ylim(-5, 5)
for a in ax:
    a.tick_params(labelsize=7)
plt.tight_layout()
fig.savefig(hip.HERE / "check_S58.png", dpi=130)
print("wrote fig1.pdf, fig1.png, check_S58.png")
