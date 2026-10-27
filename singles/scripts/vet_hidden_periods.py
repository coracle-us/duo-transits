"""
Hand vetting of the repeats-inconsistent stars flagged by recheck_periods_robust.py as having a
period that explains all but one repeat (singles/results/period_recheck/README.md, 2026-10-12):
140215502, 162031311, 76140400.

Per star:
  1. scripts/vet_candidate.vet on each (event, repeat) pair at n = 1. This gives per-epoch depth scale,
     depths in every product (SPOC PDC/SAP, TESS-SPOC PDC/SAP, QLP), centroid shifts, momentum
     dumps and segment edges, and the Gaia neighbours. The singles inputs replace the Hawthorn
     catalogue row: Tc1, Tc2 = the earlier and later of event_T and the repeat, the catalogue depth and duration, b unknown.
  2. vet_candidate.vet on the specific period hypotheses listed in HYPOTHESES.
  3. Cross-star coincidence: every other singles star (all runs) with a repeat within 0.25 d of
     this star's repeat. Shared timing points to a spacecraft or scattered-light systematic.
  4. Our own TESScut aperture photometry (singles/scripts/tesscut_photometry.py, cached cutouts) for
     the sectors in OWN_SECTORS: an independent reduction of every event, written under OUT/own/.
     For the events in KEEP_STRAYLIGHT a second version ('OWNSL') keeps stray-light-flagged cadences.
  5. Per-event trapezoid fit (Tc, depth, T14, ingress fraction) in every product, own photometry
     included, so each event's depth and duration can be compared without a shared model.
  6. Quality flags at every event (FFI cutout, SPOC, TESS-SPOC, QLP); 2048 = stray light.
  7. Wide-window plots (+-1.5 d, every product, dumps marked): <TIC>_events_wide.png.
  8. Predicted transits of the HYPOTHESES fits for the next 15 months: <TIC>_<label>_predictions.csv.

Usage:  python singles/scripts/vet_hidden_periods.py [--lc-dir PATH]
Output: singles/results/vetting/hidden_periods/<TIC>_<label>.{json,png,_epochs.csv}, summary.csv,
        coincidences.csv, events.csv
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import alias_pipeline as ap  # noqa: E402
import recheck_periods_robust as rr  # noqa: E402
import vet_candidate as vc  # noqa: E402

RES = ROOT / "singles" / "results"
OUT = RES / "vetting" / "hidden_periods"
STARS = [140215502, 162031311, 76140400]
# (label, Tc1, Tc2, n): P = (Tc2 - Tc1) / n
HYPOTHESES = {
    140215502: [("P297", 2194.4489, 3681.0676, 5)],   # event + S87 + S107 (all-but-one set, n = 1 of 297.33/n)
    76140400: [("P14.74", 2398.0984, 2412.8404, 1)],  # the two S40 dips as one periodic signal
}
COINC_D = 0.25
OWN_SECTORS = {140215502: [5, 32, 87], 162031311: [16, 56, 83], 76140400: [26, 40]}


def cand_frame(j, ctoi):
    return dict(TIC=j["TIC"], Tc1=j["event_T"], Tc2=None, Depth1_ppt=j["depth_ppm"] / 1e3,
                Depth2_ppt=j["depth_ppm"] / 1e3, Dur_h=j["dur_h"], b=np.nan,
                RA=float(ctoi.RA), Dec=float(ctoi.Dec))


def run_vet(row, label, n):
    cand = pd.DataFrame([row]).set_index("TIC", drop=False)
    rep, ep = vc.vet(int(row["TIC"]), n, cand)
    # vet_candidate names files <TIC>_n<n>; rename to the hypothesis label
    for suf in (".json", ".png", "_epochs.csv"):
        src = OUT / f"{row['TIC']}_n{n}{suf}"
        dst = OUT / f"{row['TIC']}_{label}{suf}"
        if dst.exists():
            dst.unlink()
        src.rename(dst)
    return rep, ep


def own_photometry(tic, ra, dec, lc_dir):
    """tesscut_photometry.photometry on cached cutouts, written to OUT/own/ (not the shared cache)."""
    import tesscut_photometry as tp
    cut = lc_dir.parent / "tesscut" / str(tic)
    tp.di.get_cutout = lambda t, r, d, sector: next(cut.rglob(f"*s{sector:04d}*.fits"), None)
    tp.ROOT = OUT / "own"
    (tp.ROOT / "data" / "raw" / "lightcurves" / str(tic)).mkdir(parents=True, exist_ok=True)
    return [tp.photometry(tic, s, ra, dec) for s in OWN_SECTORS.get(tic, [])]


STRAYLIGHT = 2048
# (TIC, sector): events whose FFI cadences all carry the stray-light flag, so own photometry with the
# default quality == 0 cut has no data there
KEEP_STRAYLIGHT = [(162031311, 16), (76140400, 26), (76140400, 40)]


def own_photometry_straylight(tic, sector, ra, dec, lc_dir):
    """tesscut_photometry.photometry with stray-light-flagged cadences kept (all other flags cut).
    Same aperture and background recipe; its per-frame background (faintest 40% of pixels) is
    independent of SPOC's, so a background-subtraction artefact need not reproduce. Written as
    product 'OWNSL' under OUT/own/."""
    from astropy.io import fits
    from astropy.wcs import WCS
    path = next((lc_dir.parent / "tesscut" / str(tic)).rglob(f"*s{sector:04d}*.fits"))
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, flux, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"].astype(int)
    tx, ty = wcs.all_world2pix(ra, dec, 0)
    yy, xx = np.mgrid[0:flux.shape[1], 0:flux.shape[2]]
    med = np.nanmedian(flux, axis=0)
    bkg = np.nanmedian(flux[:, med <= np.nanpercentile(med, 40)], axis=1)
    aperture = (np.hypot(xx - tx, yy - ty) <= 1.6) & (med > np.nanpercentile(med, 60))
    f = np.nansum(flux[:, aperture] - bkg[:, None], axis=1)
    good = np.isfinite(t) & np.isfinite(f) & ((q & ~STRAYLIGHT) == 0) & (f > 0)
    cad = int(round(np.nanmedian(np.diff(t[good])) * 86400))
    out = OUT / "own" / "data" / "raw" / "lightcurves" / str(tic) / f"s{sector:03d}_OWNSL_{cad}s.csv.gz"
    pd.DataFrame({"time": t[good], "flux": f[good], "flux_err": np.sqrt(np.abs(f[good])), "quality": 0}).to_csv(
        out, index=False, compression="gzip")
    return out


def trapezoid(t, Tc, depth, T14, f):
    """Box with linear ingress/egress of duration f*T14 (0 < f <= 0.5)."""
    x = np.abs(t - Tc)
    half, ing = T14 / 2, max(f * T14, 1e-4)
    return 1 - depth * np.clip((half - x) / ing, 0, 1)


def fit_event(w, T, D):
    from scipy.optimize import least_squares
    t, f, s = w[0], w[1], w[2]
    x0 = [T, max(1 - np.median(f[np.abs(t - T) < D / 4]), 1e-4), D, 0.15]
    lo, hi = [T - 0.15, 0, 0.3 * D, 0.02], [T + 0.15, 0.05, 2 * D, 0.5]
    r = least_squares(lambda x: (f - trapezoid(t, *x)) / s, np.clip(x0, lo, hi), bounds=(lo, hi))
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * max(1.0, np.sum(r.fun ** 2) / (len(t) - 4))
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    return r.x, err, int(np.sum(np.abs(t - r.x[0]) < r.x[2] / 2))


def event_shapes(tic, j, own_files):
    D = j["dur_h"] / 24
    prods = vc.all_products(tic) + [vc.load_product(Path(f)) for f in own_files]
    rows = []
    for lab, T in [("event", j["event_T"])] + [(f"repeat S{x['sector']}", x["T"]) for x in j["repeats"]]:
        for p in prods:
            if not p["time"].min() - 0.5 < T < p["time"].max() + 0.5:
                continue
            for col in ("flux", "sap_flux"):
                if col not in p or (p["product"] in ("QLP", "OWN", "OWNSL") and col == "sap_flux"):
                    continue
                w = vc.window(p, T, D, col=col)
                if w is None or np.sum(np.abs(w[0] - T) < D / 2) < 5:
                    continue
                x, e, n_in = fit_event(w, T, D)
                rows.append(dict(TIC=tic, label=lab, T_grid=T, sector=p["sector"],
                                 product=p["product"] + (":SAP" if col == "sap_flux" else ""),
                                 Tc=round(x[0], 4), Tc_err=round(e[0], 4), depth_ppt=round(x[1] * 1e3, 2),
                                 depth_err=round(e[1] * 1e3, 2), T14_h=round(x[2] * 24, 2),
                                 T14_err=round(e[2] * 24, 2), ingress_frac=round(x[3], 2), n_in=n_in))
    return rows


def quality_bits(q):
    return " ".join(f"{1 << b}:{int(np.sum((q & (1 << b)) > 0))}" for b in range(17) if np.any(q & (1 << b)))


def event_flags(tic, j, lc_dir):
    from astropy.io import fits
    D = j["dur_h"] / 24
    rows = []
    for lab, T in [("event", j["event_T"])] + [(f"repeat S{x['sector']}", x["T"]) for x in j["repeats"]]:
        for f in (lc_dir.parent / "tesscut" / str(tic)).rglob("*.fits"):
            with fits.open(f) as h:
                t, q = h[1].data["TIME"], h[1].data["QUALITY"].astype(int)
            m = np.abs(t - T) < D / 2
            if m.any():
                rows.append(dict(TIC=tic, label=lab, T=T, source="FFI cutout", n_in=int(m.sum()),
                                 n_q0=int((q[m] == 0).sum()), bits=quality_bits(q[m])))
        for f in sorted((lc_dir / str(tic)).glob("s*.csv.gz")):
            if not any(k in f.name for k in ("_SPOC_120", "_TESS-SPOC_", "_QLP_")):
                continue
            d = pd.read_csv(f, usecols=["time", "quality"])
            m = (np.abs(d.time - T) < D / 2).to_numpy()
            if m.any():
                q = d.quality.to_numpy()[m].astype(int)
                rows.append(dict(TIC=tic, label=lab, T=T, source=f.name, n_in=int(m.sum()),
                                 n_q0=int((q == 0).sum()), bits=quality_bits(q)))
    return rows


def plot_wide(tic, j, own_files):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = j["dur_h"] / 24
    prods = vc.all_products(tic) + [vc.load_product(Path(f)) for f in own_files]
    evs = [("event", j["event_T"])] + [(f"repeat S{x['sector']}", x["T"]) for x in j["repeats"]]
    fig, axs = plt.subplots(1, len(evs), figsize=(4.2 * len(evs), 3.4), squeeze=False)
    for ax, (lab, T) in zip(axs[0], evs):
        off = 0
        for p in prods:
            m = np.abs(p["time"] - T) < 1.5
            if m.sum() < 10:
                continue
            tt, ff = p["time"][m], p["flux"][m]
            tb = np.floor((tt - T) * 24)    # 1-h bins
            b = pd.DataFrame({"t": tt, "f": ff}).groupby(tb).mean()
            ax.plot((b.t - T) * 24, (b.f / np.median(ff) - 1) * 1e3 - off, ".-", ms=3, lw=0.5,
                    label=f"{p['product']} S{p['sector']}")
            for dmp in p["dumps"][np.abs(p["dumps"] - T) < 1.5]:
                ax.axvline((dmp - T) * 24, color="r", lw=0.5, ls=":")
            off += 3 * j["depth_ppm"] / 1e3
        ax.axvspan(-D * 12, D * 12, color="0.9")
        ax.set_title(f"TIC {tic} {lab} T={T:.2f}", fontsize=8)
        ax.set_xlabel("hours from T", fontsize=7)
        ax.legend(fontsize=5)
        ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(OUT / f"{tic}_events_wide.png", dpi=90)
    plt.close(fig)


def predictions(tic, label, divisors, t_from, t_to):
    """Transit times (UTC) of P_fit/m for each surviving divisor m, from the hypothesis fit."""
    from astropy.time import Time
    r = json.loads((OUT / f"{tic}_{label}.json").read_text())
    T0, P = r["fit"]["T0"], r["fit"]["P"]
    eT0, eP = r["fit_err"]["T0"], r["fit_err"]["P"]
    rows = []
    for m in divisors:
        Pm = P / m
        for k in range(int(np.ceil((t_from - T0) / Pm)), int(np.floor((t_to - T0) / Pm)) + 1):
            T = T0 + k * Pm
            err = np.hypot(eT0, k * eP / m)
            rows.append(dict(TIC=tic, alias_P=round(Pm, 4), m=m, BTJD=round(T, 4), err_h=round(err * 24, 2),
                             utc=Time(T + 2457000, format="jd", scale="tdb").utc.iso[:16]))
    df = pd.DataFrame(rows).sort_values("BTJD")
    df["n_aliases_predicting"] = df.groupby(df.BTJD.round(1)).alias_P.transform("count")
    df.to_csv(OUT / f"{tic}_{label}_predictions.csv", index=False)
    return df


def all_repeats():
    rows = []
    for run in rr.RUNS:
        for f in (RES / run).glob("*.json"):
            j = json.loads(f.read_text())
            for x in j.get("repeats") or []:
                rows.append(dict(TIC=int(j["TIC"]), run=run, T=x["T"], sector=x["sector"], product=x["product"],
                                 snr=x["snr"]))
    return pd.DataFrame(rows).drop_duplicates(["TIC", "T"])


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--lc-dir", type=Path, default=None)
    lc_dir = a.parse_args().lc_dir or rr.default_lc_dir()
    vc.LC = ap.LC = lc_dir
    vc.OUT = OUT
    vc.PAPER_SECTORS = set()   # Hawthorn-era flag is meaningless here; 'new' = not one of the two input events
    OUT.mkdir(parents=True, exist_ok=True)
    ctoi = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv").drop_duplicates("TIC").set_index("TIC")
    reps_all = all_repeats()
    summary, coinc, shapes, flags = [], [], [], []
    for tic in STARS:
        j = json.loads((RES / "full" / f"{tic}.json").read_text())
        ra, dec = float(ctoi.loc[tic].RA), float(ctoi.loc[tic].Dec)
        own = own_photometry(tic, ra, dec, lc_dir)
        own += [own_photometry_straylight(tic, s, ra, dec, lc_dir) for t_, s in KEEP_STRAYLIGHT if t_ == tic]
        shapes += event_shapes(tic, j, own)
        flags += event_flags(tic, j, lc_dir)
        plot_wide(tic, j, own)
        base = cand_frame(j, ctoi.loc[tic])
        jobs = [(f"pair_S{x['sector']}_{x['T']:.1f}", dict(base, Tc1=min(j["event_T"], x["T"]),
                                                          Tc2=max(j["event_T"], x["T"])), 1)
                for x in j["repeats"]]
        jobs += [(lab, dict(base, Tc1=t1, Tc2=t2), n) for lab, t1, t2, n in HYPOTHESES.get(tic, [])]
        for lab, row, n in jobs:
            rep, ep = run_vet(row, lab, n)
            print(f"\n=== TIC {tic} {lab}: P = {rep['fit']['P']:.4f} d, depth {rep['fit']['depth_ppt_model']:.2f} ppt, "
                  f"T14 {rep['fit']['T14_h']:.1f} h, b {rep['fit']['b']:.2f}")
            print(f"depth chi2 {rep['depth_consistency_chi2']:.1f}/{rep['depth_consistency_dof']}, "
                  f"secondary {np.round(rep['secondary_depth_ppt'], 3).tolist()}, "
                  f"centroids {[rep.get(c) for c in ('mom_centr1_shift_pix', 'mom_centr2_shift_pix')]}")
            print(f"by product: { {k: np.round(v, 2).tolist() for k, v in rep['scale_by_product'].items()} }")
            print(f"Gaia: {rep.get('gaia_note')}; could host: {rep.get('gaia_neighbours_could_host_signal')}")
            print(ep[["k", "T", "sector", "product", "scale", "scale_err", "near_dump", "near_edge",
                      "centroid", "alt_products"]].to_string(index=False))
            summary.append(dict(TIC=tic, label=lab, P=rep["fit"]["P"], depth_ppt=rep["fit"]["depth_ppt_model"],
                                T14_h=rep["fit"]["T14_h"], b=rep["fit"]["b"],
                                depth_chi2=rep["depth_consistency_chi2"], dof=rep["depth_consistency_dof"],
                                epochs=" ".join(f"{r.sector}:{r.scale:.2f}+-{r.scale_err:.2f}" for r in ep.itertuples())))
        for x in j["repeats"]:
            m = reps_all[(reps_all.TIC != tic) & ((reps_all["T"] - x["T"]).abs() < COINC_D)]
            coinc.append(dict(TIC=tic, T=x["T"], sector=x["sector"], product=x["product"], n_other_stars=len(m),
                              others=" ".join(f"{r.TIC}({r.T:.2f},{r.product})" for r in m.itertuples())))
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    ev = pd.DataFrame(shapes)
    ev.to_csv(OUT / "events.csv", index=False)
    print("\nper-event trapezoid fits (every product):")
    print(ev.drop(columns=["T_grid"]).to_string(index=False))
    fl = pd.DataFrame(flags)
    fl.to_csv(OUT / "quality_flags.csv", index=False)
    print("\nquality flags at each event (2048 = stray light, 32 = momentum dump):")
    print(fl.to_string(index=False))
    # 140215502, hypothesis P297: surviving divisors m of 297.32 d (recheck: 297.3, 148.7, 99.1, 74.3, 49.6 d)
    pr = predictions(140215502, "P297", [1, 2, 3, 4, 6], 4325.0, 4790.0)
    print("\n140215502 P297 predicted transits (next 15 months):")
    print(pr.to_string(index=False))
    c = pd.DataFrame(coinc)
    c.to_csv(OUT / "coincidences.csv", index=False)
    print("\ncross-star coincidences (other stars with a repeat within 0.25 d):")
    print(c.to_string(index=False))


if __name__ == "__main__":
    main()
