"""
Leads triage, step 2: cheap checks on every signal of the 47 lead stars of the AAS82225 revision's blind search
(leads.py). Rules written 2026-11-23 before this script was first run; any later change is logged in README.md.

Data: the same light curves the search used (search_any.load_products for TESS, search_blind.k2_products for K2),
plus, for centroids, the SAP flux and flux-weighted centroids in the K2 mission llc files and in the SPOC 2-min /
TESS-SPOC light curves already in the cache. No new files are downloaded.

Per epoch and product: search_any.test_free (local quadratic baseline, red-noise sigma) at the signal's duration D,
centre scanned +/- min(D/4, 0.1 d) for event-based signals and +/- D/4 for BLS predicted epochs.
"Confirmed" in a product = SNR >= 4, depth 0.5-2x the signal depth, coverage >= 0.5 (search_any's rule).

Checks
  R  reductions.  Event signals (family, post-hoc): pass if >= 2 epochs are confirmed in >= 2 products; fail if
     fewer than 2 epochs could reach that even counting epochs covered by only one product; else untestable.
     BLS: per product, inverse-variance stack of the predicted epochs; confirmed if stacked SNR >= 5 and depth
     0.5-2x. Pass if confirmed in >= 2 products; fail if >= 2 products were tested and < 2 confirm; else untestable.
     BLS dominance: share of the primary-product stacked chi^2 (sum (d/sigma)^2) carried by the single strongest
     epoch. >= 0.6 means one event carries the peak ("one-event BLS"); that signal is not periodic evidence.
  C  own control.  control-level if the star's own inverted-flux statistic >= 0.8 x the signal's statistic:
     post-hoc support vs max_support_control (link_calibrated.json); BLS SDE vs the star's maximum control SDE
     (bls_control, any block); family: n_members vs the number of the star's eligible control events at SNR >= the
     family's weakest member.
  X  centroid.  Per epoch, regress SAP flux and both moment centroids in +/- max(0.75 d, 2D) on a time quadratic
     (K2: plus POS_CORR1/2 to second order), out-of-event cadences only; in-event (|t-T| < 0.35 D) mean residuals
     give the SAP depth d and the centroid shift dc. Source offset from the out-of-event centroid = -dc (1-d)/d,
     sigma from the out-of-event scatter (red-noise inflated). Epochs with SAP d/sigma < 3 are skipped. Weighted mean
     per mission, in arcsec (K2 3.98"/px, TESS 21"/px).
     off-target if |offset| >= 3 sigma and >= 0.5 px; on-target if not off-target and sigma <= 0.5 px;
     else uninformative.  Added 2026-11-24 after the first full run: |offset| >= 5 px is "excursion" (a source that
     far off cannot shift a moment centroid in proportion to the depth; it means a pointing or systematic jump),
     which is neither off- nor on-target.  This is a light-curve centroid test, not a pixel difference image.
  E  eclipsing binary.  Any of: LaCourse & Jacobs comment says "secondary eclipse"; odd-even depth difference
     >= 3 sigma; box depth at phase 0.5 >= 3 sigma (circular orbit, single-period signals only); sqrt(depth) R* >
     2 R_Jup (R* from the K2 catalogue; box depth); VSX type starting with E; Gaia DR3 classifier ECL; TESS EB
     catalogue (Prsa+2022) entry.
  Flags (no verdict by themselves): Gaia DR3 NSS orbit or RUWE > 1.4 (binary host); VSX / Gaia variable that is not
     eclipsing (if the signal period is within 1% of the variable's period, x1/2 or x2, the signal is "systematic").

Verdict, first that applies: systematic (R fail, or period = host variability) > control-level (C) > EB (E) >
off-target (X off in any mission) > on-target candidate (R pass and X on-target in at least one mission, no binary
flag, not a one-event BLS) > unclear.
Per star: the best verdict of its signals, in the order on-target candidate > unclear > off-target > EB >
control-level > systematic.

Output: k2singles/leads_triage/signals.csv (one row per signal), triage.csv (one row per lead star)
Usage:  python k2singles/leads_triage/triage.py
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "k2singles" / "revision2"))
sys.path.insert(0, str(ROOT / "singles" / "revision_search"))
import search_any as sa  # noqa: E402
import search_blind as sb  # noqa: E402
import search_periodic as sp  # noqa: E402
from leads import DATA, REV, lead_epics, signals  # noqa: E402

K2_RANK = sp.K2_RANK
PIX = {"K2": 3.98, "TESS": 21.0}
RJ_RSUN = 0.10045
K2DL = ROOT / "k2singles" / "cache"
LCDIR = ROOT / "data" / "raw" / "lightcurves"


# ---------------------------------------------------------------- data
def load(epic, tic):
    pt, _ = sa.load_products(tic)
    pk = sb.k2_products(int(epic))
    prods = dict(pt)
    prods.update(pk)
    return prods


def where(prods, T):
    """(sector, {product: (t, f)}) of the segment containing T in each product."""
    for sec, avail in prods.items():
        out = {}
        for key, p in avail.items():
            t = p["t"]
            if not (t[0] - 0.5 <= T <= t[-1] + 0.5):
                continue
            j = int(np.clip(np.searchsorted(t, T), 0, len(t) - 1))
            m = p["seg"] == p["seg"][j]
            out[key] = (t[m], p["f"][m])
        if out:
            return sec, out
    return None, {}


def rank_key(sec, keys):
    order = K2_RANK if sec >= 1000 else sa.RANK
    return next((k for k in order if k in keys), None)


def measure(prods, T, D, scan):
    sec, segs = where(prods, T)
    res = {}
    for key, (t, f) in segs.items():
        r = sa.test_free(t, f, T, D, scan)
        if r is not None and r["cover"] >= 0.5:
            res[key] = dict(depth=r["depth"], sig=r["beta"] * r["sigma_white"] / np.sqrt(r["n_in"]), snr=r["snr"],
                            T=r["T"], cover=r["cover"])
    if len(res) >= 2:  # product spread: K2 detrendings disagree on depth by up to ~40%
        ds = [m["depth"] for m in res.values()]
        for m in res.values():
            m["sig_sys"] = float(np.hypot(m["sig"], 0.5 * (max(ds) - min(ds))))
    for m in res.values():
        m.setdefault("sig_sys", m["sig"])
    return sec, res


def confirmed(m, depth):
    return m["snr"] >= sa.CONF_SNR and sa.CONF_LO * depth <= m["depth"] <= sa.CONF_HI * depth


# ---------------------------------------------------------------- R
def check_R(sig, prods):
    D, dep = sig["D"], sig["depth"]
    scan = D / 4 if sig["kind"] == "BLS" else min(D / 4, 0.1)
    ep = []
    for T in sig["T_events"]:
        sec, res = measure(prods, T, D, scan)
        if res:
            ep.append(dict(T=T, sec=sec, res=res))
    out = dict(n_epochs_covered=len(ep))
    if sig["kind"] != "BLS":
        n_multi = sum(sum(confirmed(m, dep) for m in e["res"].values()) >= 2 for e in ep)
        n_one = sum(len(e["res"]) == 1 for e in ep)
        n_contra = sum(len(e["res"]) >= 2 and sum(confirmed(m, dep) for m in e["res"].values()) == 1 for e in ep)
        n_none = sum(sum(confirmed(m, dep) for m in e["res"].values()) == 0 for e in ep)
        R = "pass" if n_multi >= 2 else ("fail" if n_multi + n_one < 2 else "untestable")
        prods_seen = sorted({k for e in ep for k in e["res"]})
        out.update(R=R, R_detail=f"{n_multi} epochs confirmed in >=2 products, {n_contra} in only 1 of >=2, "
                                 f"{n_one} covered by 1 product, {n_none} confirmed nowhere; products {' '.join(prods_seen)}")
        out["dominance"] = np.nan
    else:
        stack = {}
        for e in ep:
            for k, m in e["res"].items():
                fam = "K2" if e["sec"] >= 1000 else "TESS"
                stack.setdefault(f"{fam}:{k}", []).append(m)
        conf, tested, txt = 0, 0, []
        for k, ms in sorted(stack.items()):
            w = np.array([1 / m["sig"] ** 2 for m in ms])
            d = np.array([m["depth"] for m in ms])
            ds, ss = float(np.sum(w * d) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))
            ok = ds / ss >= 5 and 0.5 * dep <= ds <= 2 * dep
            tested += 1
            conf += ok
            txt.append(f"{k} {ds * 1e6:.0f}+-{ss * 1e6:.0f}ppm n={len(ms)}{' ok' if ok else ''}")
        R = "pass" if conf >= 2 else ("fail" if tested >= 2 else "untestable")
        # dominance in the primary product of each epoch
        chi = []
        for e in ep:
            k = rank_key(e["sec"], e["res"])
            m = e["res"][k]
            chi.append(max(m["depth"], 0) / m["sig"] if m["sig"] > 0 else 0)
        chi2 = np.array(chi) ** 2
        out["dominance"] = float(chi2.max() / chi2.sum()) if chi2.sum() > 0 else np.nan
        out["n_epochs_snr3"] = int(np.sum(np.array(chi) >= 3))
        out.update(R=R, R_detail="; ".join(txt))
    out["_epochs"] = ep
    return out


# ---------------------------------------------------------------- odd/even, secondary
def check_E_lc(sig, prods, ep):
    P = sig["P"]
    out = dict(oddeven_sigma=np.nan, oddeven_detail="", sec_depth_ppm=np.nan, sec_sigma=np.nan)
    if P is None or not ep:
        out["oddeven_detail"] = "no single period"
        return out
    T0 = sig["T0"]
    od, ev = [], []
    for e in ep:
        k = rank_key(e["sec"], e["res"])
        m = e["res"][k]
        n = int(round((e["T"] - T0) / P))
        (od if n % 2 else ev).append(m)

    def comb(ms):
        w = np.array([1 / m["sig_sys"] ** 2 for m in ms])
        return float(np.sum(w * np.array([m["depth"] for m in ms])) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))

    if od and ev:
        (a, sa_), (b, sb_) = comb(od), comb(ev)
        out["oddeven_sigma"] = abs(a - b) / np.hypot(sa_, sb_)
        out["oddeven_detail"] = f"odd {a * 1e6:.0f}+-{sa_ * 1e6:.0f} (n={len(od)}), even {b * 1e6:.0f}+-{sb_ * 1e6:.0f} (n={len(ev)})"
    else:
        out["oddeven_detail"] = f"only {'odd' if od else 'even'} epochs covered"
    # phase 0.5 over all data
    tmin = min(p["t"][0] for a in prods.values() for p in a.values())
    tmax = max(p["t"][-1] for a in prods.values() for p in a.values())
    ns = np.arange(np.floor((tmin - T0) / P) - 1, np.ceil((tmax - T0) / P) + 1)
    ms = []
    for n in ns:
        T = T0 + (n + 0.5) * P
        sec, res = measure(prods, T, sig["D"], 0.0)
        if res:
            ms.append(res[rank_key(sec, res)])
    if ms:
        d, s = comb(ms)
        out.update(sec_depth_ppm=d * 1e6, sec_sigma=d / s, sec_n=len(ms))
    return out


# ---------------------------------------------------------------- X centroid
_cent_cache = {}


def cent_data(epic, tic, sec):
    key = (epic, sec)
    if key in _cent_cache:
        return _cent_cache[key]
    out = None
    if sec >= 1000:
        from astropy.io import fits
        fs = sorted((K2DL / str(epic) / "dl").glob(f"ktwo{epic}-c{sec - 1000:02d}_llc.fits"))
        if fs:
            with fits.open(fs[0]) as h:
                d = h[1].data
                ok = (d["SAP_QUALITY"] == 0) & np.isfinite(d["TIME"]) & np.isfinite(d["SAP_FLUX"]) & \
                    np.isfinite(d["MOM_CENTR1"]) & np.isfinite(d["MOM_CENTR2"]) & np.isfinite(d["POS_CORR1"]) & \
                    np.isfinite(d["POS_CORR2"])
                out = dict(t=np.asarray(d["TIME"][ok], float) - 2167.0, f=np.asarray(d["SAP_FLUX"][ok], float),
                           cx=np.asarray(d["MOM_CENTR1"][ok], float), cy=np.asarray(d["MOM_CENTR2"][ok], float),
                           p1=np.asarray(d["POS_CORR1"][ok], float), p2=np.asarray(d["POS_CORR2"][ok], float), mission="K2")
    else:
        for prod in ("SPOC_120s", "TESS-SPOC_*s"):
            fs = sorted((LCDIR / str(tic)).glob(f"s{sec:03d}_{prod}.csv.gz"))
            if not fs:
                continue
            df = pd.read_csv(fs[0])
            if "mom_centr1" not in df:
                continue
            df = df[(df.quality.fillna(0) == 0) & np.isfinite(df.time) & np.isfinite(df.sap_flux) &
                    np.isfinite(df.mom_centr1) & np.isfinite(df.mom_centr2)].sort_values("time")
            if len(df) > 50:
                out = dict(t=df.time.to_numpy(float), f=df.sap_flux.to_numpy(float), cx=df.mom_centr1.to_numpy(float),
                           cy=df.mom_centr2.to_numpy(float), mission="TESS", prod=prod.split("_")[0])
                break
    _cent_cache[key] = out
    return out


def centroid_epoch(c, T, D):
    B = max(0.75, 2 * D)
    w = np.abs(c["t"] - T) < D / 2 + B
    x = c["t"][w] - T
    oot = np.abs(x) > D / 2 + 0.02
    inn = np.abs(x) < 0.35 * D
    if oot.sum() < 20 or inn.sum() < 2 or (x[oot] < 0).sum() < 5 or (x[oot] > 0).sum() < 5:
        return None
    cols = [np.ones(w.sum()), x, x ** 2]
    if c["mission"] == "K2":
        p1, p2 = c["p1"][w], c["p2"][w]
        cols += [p1, p2, p1 ** 2, p2 ** 2, p1 * p2]
    A = np.vstack(cols).T
    res = {}
    fmed = np.median(c["f"][w][oot])
    for name, y in (("f", c["f"][w] / fmed - 1), ("cx", c["cx"][w]), ("cy", c["cy"][w])):
        coef, *_ = np.linalg.lstsq(A[oot], y[oot], rcond=None)
        r = y - A @ coef
        ro = r[oot]
        sw = 1.4826 * np.median(np.abs(ro - np.median(ro)))
        cad = float(np.median(np.diff(c["t"][w])))
        beta = sa.red_noise_beta(ro, max(1, int(round(0.7 * D / cad))))
        b = 1.0 if not np.isfinite(beta) else max(1.0, beta)
        res[name] = (float(np.mean(r[inn])), b * sw / np.sqrt(inn.sum()))
    d, sd = -res["f"][0], res["f"][1]
    if d <= 0 or d / sd < 3:
        return None
    k = (1 - d) / d
    return dict(d=d, ox=-res["cx"][0] * k, oy=-res["cy"][0] * k, sx=res["cx"][1] * k, sy=res["cy"][1] * k)


def check_X(sig, ep):
    per = {"K2": [], "TESS": []}
    for e in ep:
        c = cent_data(sig["EPIC"], sig["TIC"], e["sec"])
        if c is None:
            continue
        T = e["res"][rank_key(e["sec"], e["res"])]["T"]
        r = centroid_epoch(c, T, sig["D"])
        if r is not None:
            per[c["mission"]].append(r)
    out = {}
    for mis, rs in per.items():
        if not rs:
            out[f"X_{mis}"] = "no-data"
            continue
        sc = PIX[mis]

        def wm(v, s):
            w = 1 / np.array(s) ** 2
            return float(np.sum(w * np.array(v)) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))

        ox, sx = wm([r["ox"] for r in rs], [r["sx"] for r in rs])
        oy, sy = wm([r["oy"] for r in rs], [r["sy"] for r in rs])
        off = np.hypot(ox, oy)
        s = float(np.sqrt(((ox * sx) ** 2 + (oy * sy) ** 2) / off ** 2)) if off > 0 else max(sx, sy)
        smax = max(sx, sy)
        if off >= 5:  # added 2026-11-24 after the first full run (README): larger than a K2 aperture
            v = "excursion"
        elif off / s >= 3 and off >= 0.5:
            v = "off-target"
        elif smax <= 0.5:
            v = "on-target"
        else:
            v = "uninformative"
        out[f"X_{mis}"] = v
        out[f"X_{mis}_detail"] = f"offset {off * sc:.1f}\" +- {s * sc:.1f}\" (dx {ox * sc:+.1f}, dy {oy * sc:+.1f}; n={len(rs)})"
    return out


# ---------------------------------------------------------------- C control
def check_C(sig, rep, per):
    if sig["kind"] == "posthoc":
        real, ctl = sig["stat"], sig["stat_control"]
    elif sig["kind"] == "BLS":
        real = sig["stat"]
        ctl = max((s["sde"] for b in per["bls_control"] for s in b["signals"]), default=0.0)
    else:
        real = sig["stat"]
        evs = [e for e in rep["events"] if any(abs(e["T_fit"] - T) < 0.05 for T in sig["T_events"])]
        smin = min((e["snr"] for e in evs), default=np.inf)
        ctl = sum(1 for e in rep["control_events"] if sp.eligible(e) and e["snr"] >= smin)
    return dict(C="control-level" if ctl >= 0.8 * real else "ok", C_detail=f"{sig['stat_name']} {real:.4g} vs own control {ctl:.4g}")


# ---------------------------------------------------------------- catalogues
def catalogues():
    cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv", dtype={"EPIC": str}).drop_duplicates("EPIC").set_index("EPIC")
    vz = pd.read_csv(REV / "results" / "literature" / "vizier_targeted.csv", dtype={"EPIC": str})
    nss = pd.read_csv(REV / "results" / "literature" / "gaia_nss.csv", dtype={"EPIC": str})
    f = HERE / "catalog_lookup.csv"
    look = pd.read_csv(f, dtype={"EPIC": str}) if f.exists() else pd.DataFrame(columns=["EPIC", "status"])
    look = look[look.status == "ok"].set_index("EPIC")
    eb = pd.read_csv(ROOT / "data" / "raw" / "tess_eb_prsa2022.tsv", sep="\t", comment="#", skiprows=[0, 1] , dtype=str)
    eb_tics = set(pd.to_numeric(eb.iloc[:, 0].str.strip(), errors="coerce").dropna().astype(int))
    allc = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv", dtype={"EPIC": str})
    cat_epochs = allc.groupby("EPIC").epoch_btjd.apply(list).to_dict()
    return cat, vz, nss, look, eb_tics, cat_epochs


def check_cat(sig, cats):
    cat, vz, nss, look, eb_tics, cat_epochs = cats
    e = sig["EPIC"]
    c = cat.loc[e]
    eb, flags, var_periods = [], [], []
    com = str(c.comment) if isinstance(c.comment, str) else ""
    if "secondary eclipse" in com.lower():
        eb.append("catalogue comment: secondary eclipse")
    rs = float(c.Rstar) if np.isfinite(c.Rstar) else np.nan
    rp = np.sqrt(max(sig["depth"], 0)) * rs / RJ_RSUN if np.isfinite(rs) else np.nan
    if np.isfinite(rp) and rp > 2:
        eb.append(f"Rp {rp:.1f} R_J > 2")
    if sig["TIC"] in eb_tics:
        eb.append("TESS EB catalogue")
    if e in look.index:
        L = look.loc[e]
        t = str(L.get("vsx_type")) if pd.notna(L.get("vsx_type")) else ""
        if t:
            if t.startswith("E"):
                eb.append(f"VSX {t}")
            else:
                flags.append(f"VSX {t}")
            if pd.notna(L.get("vsx_period")):
                var_periods.append(float(L.vsx_period))
        g = str(L.get("gaia_vclass")) if pd.notna(L.get("gaia_vclass")) else ""
        if g:
            (eb if g.startswith("ECL") else flags).append(f"Gaia class {g}")
        if pd.notna(L.get("RUWE")) and float(L.RUWE) > 1.4:
            flags.append(f"binary host: RUWE {float(L.RUWE):.2f}")
        if pd.notna(L.get("NSS")) and int(L.NSS) > 0:
            flags.append(f"binary host: Gaia non_single_star {int(L.NSS)}")
        if pd.notna(L.get("nearest_nb_sep")):
            flags.append(f"Gaia neighbour {L.nearest_nb_sep}\" dG {L.nearest_nb_dG:+.1f} ({int(L.n_nb30)} within 30\")")
    else:
        flags.append("catalog_lookup missing")
    for r in nss[nss.EPIC == e].itertuples():
        flags.append(f"binary host: Gaia NSS {r.solution} P={r.period_d:.1f} d")
    for r in vz[vz.EPIC == e].itertuples():
        if r.catalog.startswith("B/vsx") and isinstance(r.periods, str):
            m = re.search(r"Period=([\d.]+)", r.periods)
            if m:
                var_periods.append(float(m.group(1)))
    # any catalogued period (VizieR targeted list, Gaia NSS) within 1% of an allowed period, x1/2, x1 or x2
    known = []
    cper = []
    for r in vz[vz.EPIC == e].itertuples():
        if isinstance(r.periods, str):
            for name, val in re.findall(r"(\w+)=(-?[\d.]+(?:e-?\d+)?)", r.periods):
                if name.lower().startswith("freq"):
                    continue
                try:
                    v = float(val)
                except ValueError:
                    continue
                if v > 0.05:
                    cper.append((r.catalog, name, v))
    cper += [("Gaia NSS", r.solution, float(r.period_d)) for r in nss[nss.EPIC == e].itertuples()]
    for P in sig["P_list"] or []:
        for catn, name, v in cper:
            for k in (0.5, 1, 2):
                if abs(P / (k * v) - 1) < 0.01:
                    known.append(f"{catn} {name}={v:g} (P={P:.4f}, x{k})")
    sysv = False
    if sig["P"] is not None:
        for vp in set(var_periods):
            for k in (0.5, 1, 2):
                if abs(sig["P"] / (k * vp) - 1) < 0.01:
                    sysv = True
                    flags.append(f"P = {k} x variable period {vp}")
    ce = cat_epochs.get(e, [])
    hit = [T for T in sig["T_events"] if any(abs(T - x) <= max(0.1, sig["D"] / 2 + 0.1) for x in ce)]
    return dict(includes_catalogued_event=len(hit), Rp_RJ=rp, EB_cat="; ".join(eb), flags="; ".join(flags), period_is_variability=sysv, comment=com,
                known_period="; ".join(sorted(set(known))[:5]))


# ---------------------------------------------------------------- verdict
def verdict(r):
    why = []
    if r["R"] == "fail" or r["period_is_variability"]:
        why.append(f"R fail: {r['R_detail']}" if r["R"] == "fail" else "period matches host variability")
        return "systematic", why
    if r["C"] == "control-level":
        return "control-level", [r["C_detail"]]
    ebw = [x for x in [r["EB_cat"]] if x]
    if np.isfinite(r["oddeven_sigma"]) and r["oddeven_sigma"] >= 3:
        ebw.append(f"odd/even {r['oddeven_sigma']:.1f} sigma")
    if np.isfinite(r["sec_sigma"]) and r["sec_sigma"] >= 3:
        ebw.append(f"phase-0.5 dip {r['sec_depth_ppm']:.0f} ppm, {r['sec_sigma']:.1f} sigma")
    if ebw:
        return "EB", ebw
    xs = [m for m in ("K2", "TESS") if r.get(f"X_{m}") == "off-target"]
    if xs:
        return "off-target", [f"{m} centroid {r[f'X_{m}_detail']}" for m in xs]
    on = [m for m in ("K2", "TESS") if r.get(f"X_{m}") == "on-target"]
    binary = "binary host" in r["flags"]
    onebls = r["kind"] == "BLS" and np.isfinite(r["dominance"]) and r["dominance"] >= 0.6
    if r["R"] == "pass" and on and not binary and not onebls:
        return "on-target candidate", [f"R pass; centroid on target ({', '.join(on)})"]
    why = []
    if r["R"] != "pass":
        why.append(f"reductions {r['R']}")
    if not on:
        why.append("centroid not informative")
    if binary:
        why.append("binary host")
    if onebls:
        why.append(f"one-event BLS (dominance {r['dominance']:.2f})")
    return "unclear", why


ORDER = ["on-target candidate", "unclear", "off-target", "EB", "control-level", "systematic"]


def main():
    cats = catalogues()
    rows = []
    sigs = signals()
    cur = None
    for sig in sigs:
        e = sig["EPIC"]
        if cur != e:
            prods = load(e, sig["TIC"])
            rep = json.loads((DATA / "per_star" / f"{e}.json").read_text())
            per = json.loads((DATA / "periodic" / f"{e}.json").read_text())
            cur = e
        r = dict(EPIC=e, TIC=sig["TIC"], sig_id=sig["sig_id"], kind=sig["kind"],
                 P=sig["P"], n_allowed=sig["n_allowed"], P_shortest=min(sig["P_list"]) if sig["P_list"] else np.nan,
                 P_longest=max(sig["P_list"]) if sig["P_list"] else np.nan, D_h=round(sig["D"] * 24, 2),
                 depth_ppm=round(sig["depth"] * 1e6, 1), n_events=len(sig["T_events"]) if sig["kind"] != "BLS" else np.nan,
                 T_events=" ".join(f"{t:.4f}" for t in sig["T_events"]) if sig["kind"] != "BLS" else f"T0 {sig['T0']:.4f}")
        R = check_R(sig, prods)
        ep = R.pop("_epochs")
        r.update(R)
        r.update(check_C(sig, rep, per))
        r.update(check_E_lc(sig, prods, ep))
        r.update(check_X(sig, ep))
        r.update(check_cat(sig, cats))
        v, why = verdict(r)
        r["verdict"], r["reasons"] = v, " | ".join(why)
        print(f"{sig['sig_id']:14s} {sig['kind']:8s} {v:20s} {r['reasons'][:150]}", flush=True)
        rows.append(r)
    S = pd.DataFrame(rows)
    S.to_csv(HERE / "signals.csv", index=False, float_format="%.6g")
    star = []
    for e, g in S.groupby("EPIC", sort=True):
        g = g.assign(o=g.verdict.map(ORDER.index)).sort_values("o")
        b = g.iloc[0]
        star.append(dict(EPIC=e, TIC=b.TIC, verdict=b.verdict, best_signal=b.sig_id, best_kind=b.kind,
                         reasons=b.reasons, n_signals=len(g),
                         all_signals=" ; ".join(f"{x.sig_id}:{x.verdict}" for x in g.itertuples()),
                         posthoc=e in __import__("leads").POSTHOC, comment=b["comment"], flags=b["flags"], known_period=b["known_period"]))
    T = pd.DataFrame(star)
    T.to_csv(HERE / "triage.csv", index=False)
    print(T.verdict.value_counts().to_string())
    print("signals:", S.verdict.value_counts().to_dict())


if __name__ == "__main__":
    main()
