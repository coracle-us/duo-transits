"""
Plan step 6: vet a candidate period P = dT/n found by the alias pipeline.

Checks (all written to results/vetting/<TIC>_n<n>.{txt,png,json}):
  1. Joint transit fit (batman, quadratic LD fixed) to every epoch on the alias grid that has
     data, P and T0 free (starting at dT/n, Tc1). Windows are detrended with the same
     segment-aware polynomial baseline as the pipeline (transit masked). Delta-chi2 / BIC for
     transit vs flat, for the new-data epochs alone (sectors not used by Hawthorn+24).
  2. Per-epoch depth scale s_k (model shape fixed, amplitude free): consistency chi2 of the
     s_k about their mean; the depths must agree with the two discovery transits.
  3. Odd/even: weighted mean s over odd vs even epochs (k counted from Tc1 on THIS alias).
  4. Secondary eclipse: box of the transit duration at phase 0.5, all epochs.
  5. Independent reductions: per-epoch scale in every product available for that sector
     (SPOC 2-min PDCSAP and SAP, TESS-SPOC PDCSAP and SAP, QLP SAP).
  6. Centroids (SPOC/TESS-SPOC mom_centr1/2): in- minus out-of-transit shift per epoch,
     combined significance.
  7. Momentum dumps (quality bit 32) and data-segment edges within +/- D of each epoch.
  8. Gaia DR3 neighbours within 2.5 arcmin: depth the neighbour would need to have to produce
     the signal if it were the source (> 100% = impossible).

Usage:  python scripts/vet_candidate.py TIC n [TIC n ...]
"""
import json
import re
import sys
from pathlib import Path

import batman
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import alias_pipeline as ap

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "catalog"
LC = ROOT / "data" / "raw" / "lightcurves"
OUT = ROOT / "results" / "vetting"
PAPER_SECTORS = set(range(1, 40))
LD = (0.32, 0.22)            # quadratic LD, typical F/G star in the TESS band (fixed)
GOOD_QMASK = 1 | 2 | 4 | 8 | 32 | 128 | 1024 | 2048 | 4096 | 8192 | 32768  # bad-data bits
WIN = 0.6                     # d of data either side of each transit window


# ------------------------------------------------------------------ data per product
def load_product(path):
    """One light-curve file -> dict of arrays (good cadences), with segment ids."""
    df = pd.read_csv(path)
    m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", path.name)
    sec, author, exp = int(m.group(1)), m.group(2), int(m.group(3))
    key = f"SPOC_{exp}" if author == "SPOC" else author
    dumps = df.time[(df.quality & 32) > 0].to_numpy()
    good = np.isfinite(df.time) & np.isfinite(df.flux) & (df.quality == 0)
    df = df[good].sort_values("time")
    seg = np.concatenate([[0], np.cumsum(np.diff(df.time.values) > ap.SEG_GAP)])
    if key == "QLP":
        st = pd.Series(df.time.values).groupby(seg).transform("min").values
        keep = df.time.values - st >= ap.QLP_TRIM
        df, seg = df[keep], seg[keep]
    out = dict(sector=sec, product=key, exptime=exp, time=df.time.to_numpy(),
               flux=(df.flux / np.nanmedian(df.flux)).to_numpy(), seg=seg, dumps=dumps)
    for col in ("sap_flux", "mom_centr1", "mom_centr2"):
        if col in df and np.isfinite(df[col]).sum() > 10:
            v = df[col].to_numpy()
            out[col] = v / np.nanmedian(v) if col == "sap_flux" else v
    seg_t = pd.Series(df.time.values).groupby(seg).agg(["min", "max"])
    out["seg_edges"] = seg_t.to_numpy()
    return out


def all_products(tic):
    prods = []
    for f in sorted((LC / str(tic)).glob("s*.csv.gz")):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if not m or int(m.group(1)) > 200 or m.group(2) in ("TGLC", "GSFC-ELEANOR-LITE"):
            continue
        if m.group(2) == "SPOC" and m.group(3) == "20":
            continue  # 2-min covers the same sectors
        try:
            p = load_product(f)
        except Exception as e:   # noqa: BLE001
            print(f"  skip {f.name}: {e}")
            continue
        if len(p["time"]) > 50:
            prods.append(p)
    return prods


def window(p, T, D, col="flux"):
    """Detrended data around T from product p: polynomial + per-segment offsets fitted to
    out-of-transit points. Returns t, f (normalised), sigma, or None."""
    t = p["time"]
    m = np.abs(t - T) < D / 2 + WIN
    if m.sum() < 10:
        return None
    tt, ff, sg = t[m], p[col][m], p["seg"][m]
    ok = np.isfinite(ff)
    tt, ff, sg = tt[ok], ff[ok], sg[ok]
    oot = np.abs(tt - T) > D / 2 + 0.03
    if oot.sum() < 8:
        return None
    keep = np.isin(sg, np.unique(sg[oot]))
    tt, ff, sg, oot = tt[keep], ff[keep], sg[keep], oot[keep]
    nl, nr = np.sum(oot & (tt < T)), np.sum(oot & (tt > T))
    deg = 2 if min(nl, nr) >= 10 else (1 if min(nl, nr) >= 3 else 0)
    segs = np.unique(sg)
    X = np.column_stack([(tt - T) ** k for k in range(1, deg + 1)] + [(sg == u) * 1.0 for u in segs])
    c, *_ = np.linalg.lstsq(X[oot], ff[oot], rcond=None)
    fn = ff / (X @ c)
    sig = 1.4826 * np.median(np.abs(fn[oot] - np.median(fn[oot])))
    return tt, fn, max(sig, 1e-6), nl, nr


# ------------------------------------------------------------------ transit model
def model(t, T0, P, rp, b, aR, exptime):
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a = T0, P, rp, aR
    pm.inc = np.degrees(np.arccos(np.clip(b / aR, 0, 1)))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, list(LD), "quadratic"
    ss = max(1, int(round(exptime / 120)))
    m = batman.TransitModel(pm, t, supersample_factor=ss, exp_time=exptime / 86400.0)
    return m.light_curve(pm)


def fit_joint(data, T0, P, D, delta, b0):
    """data: list of (t, f, sig, exptime). Fit T0, P, rp, b, aR."""
    rp0 = np.sqrt(delta)
    aR0 = max(2.0, P / (np.pi * D) * np.sqrt((1 + rp0) ** 2 - min(b0, 0.9) ** 2))

    def resid(x):
        r = []
        for t, f, s, ex in data:
            r.append((f - model(t, x[0], x[1], x[2], x[3], x[4], ex)) / s)
        return np.concatenate(r)

    x0 = [T0, P, rp0, min(b0, 0.9), aR0]
    lo = [T0 - 0.1, P * 0.999, 0.005, 0.0, 1.5]
    hi = [T0 + 0.1, P * 1.001, 0.5, 1.3, max(300.0, 3 * aR0)]  # long trial periods need large a/R*
    r = least_squares(resid, x0, bounds=(lo, hi), x_scale=[0.003, 1e-4, 0.005, 0.1, 2])
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, np.sum(r.fun ** 2) / (len(r.fun) - 5))
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(5, np.nan)
    return r.x, err, np.sum(r.fun ** 2), len(r.fun)


def scale_fit(t, f, s, mdl):
    """Best amplitude multiplier for a fixed transit shape: f = 1 - k (1 - mdl)."""
    dm = 1 - mdl
    w = 1 / s ** 2
    den = np.sum(w * dm * dm)
    if den <= 0:
        return np.nan, np.nan
    k = np.sum(w * dm * (1 - f)) / den
    return k, 1 / np.sqrt(den)


def wmean(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if not ok.any():
        return np.nan, np.nan
    w = 1 / e[ok] ** 2
    return np.sum(w * v[ok]) / np.sum(w), 1 / np.sqrt(np.sum(w))


# ------------------------------------------------------------------ neighbours
def gaia_neighbours(ra, dec, depth, radius_arcsec=150):
    try:
        from astroquery.gaia import Gaia
        Gaia.ROW_LIMIT = 200
        q = f"""SELECT source_id, ra, dec, phot_g_mean_mag,
                DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep
                FROM gaiadr3.gaia_source
                WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{radius_arcsec/3600}))
                ORDER BY phot_g_mean_mag"""
        tab = Gaia.launch_job(q).get_results().to_pandas()
    except Exception as e:   # noqa: BLE001
        return None, f"Gaia query failed: {e}"
    if tab.empty:
        return tab, "no Gaia sources"
    tab = tab.sort_values("sep")
    tgt = tab.iloc[0]
    tab["dG"] = tab.phot_g_mean_mag - tgt.phot_g_mean_mag
    # depth the neighbour must have if it alone causes the dip (all flux in aperture)
    fr = 10 ** (-0.4 * tab.dG)
    tab["required_depth"] = depth * (1 + 1 / fr)
    return tab, f"target G={tgt.phot_g_mean_mag:.2f} (sep {tgt.sep:.1f}\")"


# ------------------------------------------------------------------ main vetting
def vet(tic, n, cand):
    c = cand.loc[tic]
    dT = c.Tc2 - c.Tc1
    P0, D = dT / n, c.Dur_h / 24
    delta = 1e-3 * np.nanmean([c.Depth1_ppt, c.Depth2_ppt])
    prods = all_products(tic)
    best_rank = {}
    for p in prods:
        r = ap.PRODUCT_RANK.index(p["product"]) if p["product"] in ap.PRODUCT_RANK else 9
        if p["sector"] not in best_rank or r < best_rank[p["sector"]][0]:
            best_rank[p["sector"]] = (r, p)
    best = {s: v[1] for s, v in best_rank.items()}

    # epochs with data (best product per sector)
    tmin = min(p["time"].min() for p in prods)
    tmax = max(p["time"].max() for p in prods)
    ks = range(int(np.ceil((tmin - c.Tc1) / P0)) - 1, int(np.floor((tmax - c.Tc1) / P0)) + 2)
    epochs = []
    for k in ks:
        T = c.Tc1 + k * P0
        for s, p in best.items():
            if p["time"].min() - 1 < T < p["time"].max() + 1:
                w = window(p, T, D)
                if w is not None and np.sum(np.abs(w[0] - T) < D / 2) >= 3:
                    epochs.append(dict(k=k, T=T, sector=s, p=p, w=w))
    fitdata = [(e["w"][0], e["w"][1], e["w"][2], e["p"]["exptime"]) for e in epochs]
    x, xe, chi2, npts = fit_joint(fitdata, c.Tc1, P0, D, delta, c.b if np.isfinite(c.b) else 0.5)
    T0, P, rp, b, aR = x
    rep = {"TIC": int(tic), "n": int(n), "P_alias": P0, "fit": dict(zip(
        ["T0", "P", "rp", "b", "aR"], map(float, x))), "fit_err": dict(zip(
        ["T0", "P", "rp", "b", "aR"], map(float, xe)))}
    # duration T14 from fit
    inc = np.arccos(min(b / aR, 1))
    arg = np.sqrt(max((1 + rp) ** 2 - b ** 2, 0)) / (aR * np.sin(inc))
    rep["fit"]["T14_h"] = float(P / np.pi * np.arcsin(min(arg, 1)) * 24)
    rep["fit"]["depth_ppt_model"] = float((1 - model(np.array([T0]), T0, P, rp, b, aR, 120)[0]) * 1e3)

    # per-epoch scales in best product; delta-chi2 for new data
    rows = []
    dchi_new, npts_new = 0.0, 0
    for e in epochs:
        t, f, s, nl, nr = e["w"]
        mdl = model(t, T0, P, rp, b, aR, e["p"]["exptime"])
        k_, ke = scale_fit(t, f, s, mdl)
        chi_flat = np.sum(((f - 1) / s) ** 2)
        chi_tr = np.sum(((f - mdl) / s) ** 2)
        new = e["sector"] not in PAPER_SECTORS and e["k"] not in (0, n)
        if new:
            dchi_new += chi_flat - chi_tr
            npts_new += len(t)
        Tk = T0 + e["k"] * P
        near_dump = bool(np.any(np.abs(e["p"]["dumps"] - Tk) < D))
        edges = e["p"]["seg_edges"]
        near_edge = bool(np.any((np.abs(edges[:, 0] - Tk) < D) | (np.abs(edges[:, 1] - Tk) < D)))
        cen = {}
        for cc in ("mom_centr1", "mom_centr2"):
            if cc in e["p"]:
                tm = np.abs(e["p"]["time"] - Tk)
                v = e["p"][cc]
                inn, out = (tm < 0.35 * D), (tm > D / 2 + 0.05) & (tm < D / 2 + WIN)
                if inn.sum() > 3 and out.sum() > 10:
                    d = np.nanmean(v[inn]) - np.nanmean(v[out])
                    se = np.nanstd(v[out]) * np.sqrt(1 / inn.sum() + 1 / out.sum())
                    cen[cc] = (float(d), float(se))
        # other products for the same epoch (independent reductions, SAP)
        alt = {}
        for p in prods:
            if p["sector"] != e["sector"]:
                continue
            for col in ("flux", "sap_flux"):
                if col not in p or (p["product"] == "QLP" and col == "sap_flux"):
                    continue  # QLP 'flux' is already SAP
                w = window(p, Tk, D, col=col)
                if w is None or np.sum(np.abs(w[0] - Tk) < D / 2) < 3:
                    continue
                kk, kke = scale_fit(w[0], w[1], w[2], model(w[0], T0, P, rp, b, aR, p["exptime"]))
                lab = p["product"] + (":SAP" if col == "sap_flux" or p["product"] == "QLP" else ":PDC")
                alt[lab] = (float(kk), float(kke))
        rows.append(dict(k=e["k"], T=Tk, sector=e["sector"], product=e["p"]["product"],
                         new=new, scale=k_, scale_err=ke, odd=(e["k"] % 2 == 1),
                         n_oot_l=nl, n_oot_r=nr, near_dump=near_dump, near_edge=near_edge,
                         centroid=json.dumps(cen), alt_products=json.dumps(alt)))
    ep = pd.DataFrame(rows)
    rep["n_epochs"] = len(ep)
    rep["n_new_epochs"] = int(ep.new.sum())
    rep["new_data_delta_chi2"] = float(dchi_new)
    rep["new_data_delta_BIC"] = float(dchi_new - 5 * np.log(max(npts_new, 1)))
    m, me = wmean(ep.scale, ep.scale_err)
    ok = ep.scale_err > 0
    rep["depth_scale_mean"] = [m, me]
    rep["depth_consistency_chi2"] = float(np.sum(((ep.scale[ok] - m) / ep.scale_err[ok]) ** 2))
    rep["depth_consistency_dof"] = int(ok.sum() - 1)
    disc = ep[ep.k.isin([0, n])]
    rep["discovery_scale"] = list(wmean(disc.scale, disc.scale_err))
    newe = ep[ep.new]
    rep["new_scale"] = list(wmean(newe.scale, newe.scale_err))
    o, oe = wmean(ep.scale[ep.odd], ep.scale_err[ep.odd])
    ev, eve = wmean(ep.scale[~ep.odd], ep.scale_err[~ep.odd])
    rep["odd_even"] = {"odd": [o, oe], "even": [ev, eve],
                       "diff_sigma": float(abs(o - ev) / np.hypot(oe, eve)) if np.isfinite(o + ev) else None}

    # secondary eclipse at phase 0.5: box depth from best products
    sd, se = [], []
    for k in ks:
        Ts = T0 + (k + 0.5) * P
        for s, p in best.items():
            if p["time"].min() < Ts < p["time"].max():
                w = window(p, Ts, D)
                if w is not None:
                    inn = np.abs(w[0] - Ts) < 0.35 * D
                    if inn.sum() >= 3:
                        sd.append(1 - np.mean(w[1][inn]))
                        se.append(w[2] / np.sqrt(inn.sum()))
    rep["secondary_depth_ppt"] = [x * 1e3 for x in wmean(sd, se)] + [len(sd)]

    # centroids combined
    for cc in ("mom_centr1", "mom_centr2"):
        v = [json.loads(z).get(cc) for z in ep.centroid]
        v = [z for z in v if z]
        if v:
            mm, mme = wmean([z[0] for z in v], [z[1] for z in v])
            rep[f"{cc}_shift_pix"] = [mm, mme, float(abs(mm) / mme) if mme else None]
    # alt products: mean scale per product label
    alts = {}
    for z in ep.alt_products:
        for lab, (kk, kke) in json.loads(z).items():
            alts.setdefault(lab, []).append((kk, kke))
    rep["scale_by_product"] = {lab: list(wmean([a for a, _ in v], [b_ for _, b_ in v])) + [len(v)]
                               for lab, v in alts.items()}
    rep["epochs_near_dump"] = int(ep.near_dump.sum())
    rep["epochs_near_edge"] = int(ep.near_edge.sum())
    tab, note = gaia_neighbours(c.RA, c.Dec, delta)
    rep["gaia_note"] = note
    if tab is not None and len(tab):
        nb = tab.iloc[1:]
        rep["gaia_neighbours_within_150as"] = int(len(nb))
        risky = nb[nb.required_depth < 0.5]
        rep["gaia_neighbours_could_host_signal"] = risky[["source_id", "sep", "phot_g_mean_mag",
                                                          "dG", "required_depth"]].round(4).to_dict("records")

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / f"{tic}_n{n}"
    ep.to_csv(f"{base}_epochs.csv", index=False)
    with open(f"{base}.json", "w") as fh:
        json.dump(rep, fh, indent=1, default=float)
    plot(tic, n, c, epochs, ep, x, D, base)
    return rep, ep


def plot(tic, n, c, epochs, ep, x, D, base):
    T0, P, rp, b, aR = x
    nE = len(epochs)
    ncol = 5
    nrow = 1 + int(np.ceil(nE / ncol))
    fig = plt.figure(figsize=(3.2 * ncol, 2.6 * nrow))
    ax0 = fig.add_subplot(nrow, 2, 1)
    tt = np.linspace(-1.2 * D, 1.2 * D, 400)
    for e in epochs:
        t, f = e["w"][0], e["w"][1]
        ph = t - (T0 + e["k"] * P)
        col = "C1" if (e["sector"] not in PAPER_SECTORS and e["k"] not in (0, n)) else "C0"
        ax0.plot(ph * 24, (f - 1) * 1e3, ".", ms=1.5, color=col, alpha=0.5)
    ax0.plot(tt * 24, (model(tt + T0, T0, P, rp, b, aR, 120) - 1) * 1e3, "k-")
    ax0.set_xlim(-1.2 * D * 24, 1.2 * D * 24)
    ax0.set_title(f"TIC {tic} n={n} P={P:.5f} d  folded (blue: discovery/paper-era, orange: new)", fontsize=8)
    ax1 = fig.add_subplot(nrow, 2, 2)
    ax1.errorbar(ep["T"], ep.scale, ep.scale_err, fmt="o", ms=3)
    ax1.axhline(1, color="k", lw=0.5)
    ax1.set_title("per-epoch depth scale (1 = joint model)", fontsize=8)
    for i, e in enumerate(epochs):
        ax = fig.add_subplot(nrow, ncol, ncol + i + 1)
        t, f = e["w"][0], e["w"][1]
        Tk = T0 + e["k"] * P
        ax.plot((t - Tk) * 24, (f - 1) * 1e3, ".", ms=2, color="0.5")
        tb = np.arange(t.min(), t.max(), 1 / 48)
        idx = np.digitize(t, tb)
        bx = [t[idx == j].mean() for j in np.unique(idx)]
        by = [f[idx == j].mean() for j in np.unique(idx)]
        ax.plot((np.array(bx) - Tk) * 24, (np.array(by) - 1) * 1e3, "o", ms=2.5, color="C0")
        ax.plot((t - Tk) * 24, (model(t, T0, P, rp, b, aR, e["p"]["exptime"]) - 1) * 1e3, "r-", lw=0.8)
        r = ep.iloc[i]
        ax.set_title(f"k={e['k']} s{e['sector']} {e['p']['product']} s={r.scale:.2f}+/-{r.scale_err:.2f}"
                     + (" DUMP" if r.near_dump else "") + (" EDGE" if r.near_edge else ""), fontsize=7)
        ax.set_ylim(-3 * c.Depth1_ppt, 2 * c.Depth1_ppt)
        ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(f"{base}.png", dpi=85)
    plt.close(fig)


def main():
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC", drop=False)
    args = list(map(int, sys.argv[1:]))
    for tic, n in zip(args[::2], args[1::2]):
        rep, ep = vet(tic, n, cand)
        print(json.dumps(rep, indent=1, default=float))
        print(ep[["k", "T", "sector", "product", "new", "scale", "scale_err", "near_dump",
                  "near_edge", "alt_products"]].to_string(index=False))


if __name__ == "__main__":
    main()
