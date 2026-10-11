"""Step 5: injection-recovery on the real 2-min light curves of every host.

Per host-sector: N_TRIALS realizations, each injecting N_PER synthetic flares
(Davenport et al. 2014 template, integrated over the 2-min exposure) at random
good, non-edge, out-of-transit times separated by >= 0.5 d and >= 0.3 d from
real candidates. Amplitude log-uniform in [0.5, 50] x sector noise sigma;
FWHM log-uniform in [4, 60] min. The SAME detrend + detection + automatic flags
(edge, momentum dump, rise/decay shape) are applied. An injection counts as
recovered if an accepted event's extent (+/-2 cadences) contains its peak.

Writes injection_recovery.csv (one row per injected flare),
injrec_per_star.csv (ED/energy at 50% and 90% completeness) and
injection_recovery.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import flarelib as F

N_TRIALS, N_PER = 20, 5
AMP_RANGE = (0.5, 50.0)      # in units of sector median sigma
FWHM_MIN = (4.0, 60.0)       # minutes
SEED = 20261011
ED_BINS = np.arange(-1.0, 3.61, 0.2)  # log10 ED [s]


def completeness_level(ed, rec, level):
    """Smallest ED bin centre above which the binned recovery fraction stays >= level."""
    lg = np.log10(ed)
    idx = np.digitize(lg, ED_BINS) - 1
    cen = 0.5 * (ED_BINS[1:] + ED_BINS[:-1])
    frac = np.array([rec[idx == k].mean() if np.sum(idx == k) >= 5 else np.nan for k in range(len(cen))])
    ok = np.isfinite(frac)
    if not ok.any():
        return np.nan
    good = frac >= level
    # last bin that fails, then the next finite bin
    for k in range(len(cen) - 1, -1, -1):
        if ok[k] and not good[k]:
            nxt = [j for j in range(k + 1, len(cen)) if ok[j]]
            return 10 ** cen[nxt[0]] if nxt else np.nan
    return 10 ** cen[np.where(ok)[0][0]]


def main():
    rng = np.random.default_rng(SEED)
    cat = F.load_catalog()
    sectors = pd.read_csv(F.SQ / "sectors.csv")
    real = pd.read_csv(F.SQ / "candidates_raw.csv")
    ed_unit = F.davenport_ed_per_amp_fwhm()
    rows = []
    for n_s, s in enumerate(sectors.itertuples()):
        base = F.load_sector(int(s.tic), int(s.sector), cat)
        base_d = F.detrend(base)
        sig = float(np.nanmedian(base_d["sigma"][base["good"]]))
        t = base["t"]
        allowed = base["good"] & ~base["edge"]
        for tp in real[(real.tic == s.tic) & (real.sector == s.sector)].t_peak:
            allowed &= np.abs(t - tp) > 0.3
        tpool = t[allowed]
        if len(tpool) < 100:
            continue
        for trial in range(N_TRIALS):
            # pick N_PER peak times separated by >= 0.5 d
            peaks = []
            for _ in range(200):
                if len(peaks) == N_PER:
                    break
                tp = rng.choice(tpool) + rng.uniform(-0.5, 0.5) * F.CADENCE_D
                if all(abs(tp - p) > 0.5 for p in peaks):
                    peaks.append(tp)
            amps = sig * 10 ** rng.uniform(*np.log10(AMP_RANGE), len(peaks))
            fwhm = 10 ** rng.uniform(*np.log10(FWHM_MIN), len(peaks)) / 1440
            model = np.zeros(len(t))
            fin = np.isfinite(t)
            for tp, a, w in zip(peaks, amps, fwhm):
                near = fin & (t > tp - 2 * w - 0.01) & (t < tp + 40 * w + 0.01)
                model[near] += F.binned_template(t[near], tp, w, a)
            lc = dict(base)
            lc["f"] = base["f"] * (1 + model)
            lc = F.detrend(lc)
            fl = F.find_flares(lc)
            if len(fl):
                fl = F.auto_flags(fl, lc)
                acc = fl[(~fl.near_dump) & (~fl.edge) & fl.shape_ok]
            else:
                acc = fl
            for tp, a, w in zip(peaks, amps, fwhm):
                hit = None
                if len(acc):
                    m = (acc.t_start - 2 * F.CADENCE_D <= tp) & (acc.t_stop + 2 * F.CADENCE_D >= tp)
                    if m.any():
                        hit = acc[m].iloc[0]
                rows.append(dict(tic=s.tic, sector=s.sector, trial=trial, t_peak=tp, amp=a, amp_sig=a / sig,
                                 fwhm_min=w * 1440, ed_inj_s=ed_unit * a * w * 86400,
                                 recovered=hit is not None,
                                 ed_rec_s=float(hit.ed_s) if hit is not None else np.nan,
                                 amp_rec=float(hit.amp) if hit is not None else np.nan))
        if n_s % 20 == 0:
            print(f"{n_s}/{len(sectors)} done")
    ir = pd.DataFrame(rows)
    ir.to_csv(F.SQ / "injection_recovery.csv", index=False)

    # ------------------------------------------------ per-star completeness
    per = []
    for tic, g in ir.groupby("tic"):
        row = cat.loc[tic]
        L = F.l_tess(float(row.Teff), float(row.Rstar))
        ed50 = completeness_level(g.ed_inj_s.values, g.recovered.values, 0.5)
        ed90 = completeness_level(g.ed_inj_s.values, g.recovered.values, 0.9)
        per.append(dict(tic=tic, n_injected=len(g), frac_recovered=round(g.recovered.mean(), 3),
                        ed50_s=ed50, ed90_s=ed90, L_tess_erg_s=L,
                        E50_erg=ed50 * L, E90_erg=ed90 * L))
    per = pd.DataFrame(per)
    per.to_csv(F.SQ / "injrec_per_star.csv", index=False)
    print(per.describe().T)

    # ----------------------------------------------------------------- plot
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    lg = np.log10(ir.ed_inj_s)
    cen = 0.5 * (ED_BINS[1:] + ED_BINS[:-1])
    idx = np.digitize(lg, ED_BINS) - 1
    frac = [ir.recovered[idx == k].mean() if np.sum(idx == k) else np.nan for k in range(len(cen))]
    a = ax[0]
    for tic, g in ir.groupby("tic"):
        gi = np.digitize(np.log10(g.ed_inj_s), ED_BINS) - 1
        fr = [g.recovered[gi == k].mean() if np.sum(gi == k) >= 5 else np.nan for k in range(len(cen))]
        a.plot(10 ** cen, fr, color="grey", lw=0.5, alpha=0.5)
    a.plot(10 ** cen, frac, "o-", color="C3", lw=2, label="all hosts pooled")
    a.set_xscale("log")
    a.axhline(0.9, ls="--", color="k", lw=0.8)
    a.axhline(0.5, ls=":", color="k", lw=0.8)
    a.set_xlabel("injected ED [s]")
    a.set_ylabel("recovery fraction")
    a.set_title("Recovery vs ED (grey: individual hosts)")
    a.legend()

    a = ax[1]
    ir["E_inj"] = ir.ed_inj_s * ir.tic.map(per.set_index("tic").L_tess_erg_s)
    eb = np.arange(30, 37.01, 0.25)
    ec = 0.5 * (eb[1:] + eb[:-1])
    ei = np.digitize(np.log10(ir.E_inj), eb) - 1
    ef = [ir.recovered[ei == k].mean() if np.sum(ei == k) >= 10 else np.nan for k in range(len(ec))]
    a.plot(10 ** ec, ef, "o-", color="C0")
    a.set_xscale("log")
    a.axhline(0.9, ls="--", color="k", lw=0.8)
    a.set_xlabel("injected energy E_TESS [erg]")
    a.set_ylabel("recovery fraction (pooled)")
    a.set_title("Recovery vs TESS-band energy (pooled)")

    a = ax[2]
    xb = np.log10(np.geomspace(0.5, 50, 11))
    yb = np.log10(np.geomspace(4, 60, 9))
    H = np.full((len(xb) - 1, len(yb) - 1), np.nan)
    xi = np.digitize(np.log10(ir.amp_sig), xb) - 1
    yi = np.digitize(np.log10(ir.fwhm_min), yb) - 1
    for i in range(len(xb) - 1):
        for j in range(len(yb) - 1):
            m = (xi == i) & (yi == j)
            if m.sum():
                H[i, j] = ir.recovered[m].mean()
    im = a.pcolormesh(10 ** xb, 10 ** yb, H.T, vmin=0, vmax=1, cmap="viridis")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("injected amplitude / sigma")
    a.set_ylabel("injected FWHM [min]")
    a.set_title("Recovery vs amplitude and FWHM")
    fig.colorbar(im, ax=a, label="recovery fraction")
    fig.tight_layout()
    fig.savefig(F.SQ / "injection_recovery.png", dpi=110)


if __name__ == "__main__":
    main()
