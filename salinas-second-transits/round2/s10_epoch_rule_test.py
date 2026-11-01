#!/usr/bin/env python3
"""Compare epoch exclusion rules on stars with known periods, injections and the real candidates.

st2core.epoch_verdicts has used half3s since this test. Rules compared, with
sigma = hypot(epoch noise, reference depth error):

  current  depth < max(depth0 - 3 sigma, 0.5 depth0)  (the pipeline's rule until 2026-11-01)
  both     depth < min(depth0 - 3 sigma, 0.5 depth0)
  half3s   depth < 0.5 depth0 - 3 sigma

Period truth: the 57 'multi' rows in candidates.csv carry BLS periods from one catalogue, rounded and
mostly unconfirmed (one is known to be twice the planet's period), so they are a weak test. The main
test uses TOIs dispositioned CP or KP on ExoFOP, chosen before looking at their light curves.

Modes (run from salinas-second-transits/):
  download    fetch and detrend every sector of a star set (--set toi or multi)
  known       real transits of that set: does the true period survive, how many false aliases do
  inject      box transits added to the cached light curves of the other stars, by epoch SNR
  candidates  the same rules on the current alias events of the real candidates
  summary     print and write results/round2/epoch_rule_test_summary.csv
"""
import argparse
import glob
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402

WORK = Path("work")
OUT = Path("results/round2")
RULES = ("current", "both", "half3s")
PIPELINE_RULE = "half3s"
SNRS = [4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50]
SNR_BINS = [0, 4, 6, 8, 10, 15, 20, 30, 1e9]
SAMPLE_PER_BIN = 30
TIMING_SIG_D = 15 / 1440  # scatter added to the two known times, about what the box search gives


def thresholds(d0, tot):
    return dict(current=np.maximum(d0 - 3 * tot, 0.5 * d0),
                both=np.minimum(d0 - 3 * tot, 0.5 * d0),
                half3s=0.5 * d0 - 3 * tot)


def verdicts(depth, frac, edge_d, sig, d0, d0_err):
    """Clean-epoch mask, sigma, and per rule excludes/supports masks."""
    err = sig / np.sqrt(np.clip(frac, 1e-3, None))
    tot = np.hypot(err, d0_err)
    clean = (frac >= core.CLEAN_FRAC) & (edge_d >= core.CLEAN_EDGE_D)
    thr = thresholds(d0, tot)
    out = {}
    for r in RULES:
        ex = clean & (depth < thr[r])
        out[r] = (ex, clean & ~ex & (np.abs(depth - d0) < 3 * tot))
    return clean, tot, out


def aliases(t1, t2, dur_d, d0, d0_err, t, f, sig, edges, p_floor=1.0, n_max=3000, pmin=np.nan):
    """Alias status under each rule, same logic as core.alias_table but one box_stats call."""
    dT = t2 - t1
    tmin, tmax = np.nanmin(t), np.nanmax(t)
    ns = np.arange(1, min(n_max, int(np.floor(dT / p_floor))) + 1)
    if not ns.size:
        return pd.DataFrame()
    Ps = dT / ns
    idx, ep = [], []
    for i, (n, P) in enumerate(zip(ns, Ps)):
        k = np.arange(int(np.floor((tmin - t1) / P)) - 1, int(np.ceil((tmax - t1) / P)) + 2)
        k = k[(k != 0) & (k != n)]
        idx.append(np.full(k.size, i))
        ep.append(t1 + k * P)
    idx, ep = np.concatenate(idx), np.concatenate(ep)
    depth, frac, _ = core.box_stats(t, f, ep, dur_d)
    _, _, v = verdicts(depth, frac, core.edge_distance(edges, ep), sig, d0, d0_err)
    df = pd.DataFrame(dict(n=ns, period_d=Ps))
    for r in RULES:
        ex = np.bincount(idx[v[r][0]], minlength=ns.size) > 0
        su = np.bincount(idx[v[r][1]], minlength=ns.size) > 0
        st = np.where(ex, "excluded", np.where(su, "supported", "allowed")).astype(object)
        st[(Ps < pmin) & (st == "allowed")] = "allowed_below_pmin"
        df[r] = st
    return df


def load_star(tic, work=WORK):
    """Cached detrended sectors of one star, without touching MAST."""
    T, F = [], []
    for p in sorted(glob.glob(str(work / "lc" / f"{tic}_s*.npz"))):
        d = np.load(p)
        T.append(d["t"])
        F.append(d["f"])
    if not T:
        return None, None
    t, f = np.concatenate(T), np.concatenate(F)
    o = np.argsort(t)
    return t[o], f[o]


def star_set(name):
    """Stars for the known-period test: tic, period_d, t0, duration_h, depth_ppm, depth_err."""
    cols = ["tic", "period_d", "t0", "duration_h", "depth_ppm", "depth_err"]
    c = pd.read_csv("candidates.csv")
    if name == "multi":
        return c.loc[c["type"] == "multi", cols].assign(period_err=np.nan)
    path = OUT / "epoch_rule_toi_sample.csv"
    if path.exists():
        return pd.read_csv(path)
    t = pd.read_csv(WORK / "exofop_toi.csv", comment="#", low_memory=False)
    t = t[t["TFOPWG Disposition"].isin(["CP", "KP"]) & ~t["TIC ID"].isin(c["tic"])]
    s = pd.DataFrame(dict(tic=t["TIC ID"], toi=t["TOI"], period_d=t["Period (days)"], period_err=t["Period (days) err"],
                          t0=t["Epoch (BJD)"] - 2457000, duration_h=t["Duration (hours)"], depth_ppm=t["Depth (ppm)"],
                          depth_err=t["Depth (ppm) err"], tmag=t["TESS Mag"], planet=t["Planet Name"]))
    s = s.dropna(subset=["period_d", "t0", "duration_h", "depth_ppm", "depth_err"]).drop_duplicates("tic")
    s = s[s["period_d"].between(3, 40) & s["duration_h"].between(1.5, 12) & (s["depth_err"] > 0)]
    # fixed draw per depth range, so shallow (low SNR) planets are not crowded out by hot Jupiters
    s["depth_bin"] = pd.cut(s["depth_ppm"], [500, 2000, 5000, 25000]).astype(str)
    s = s[s["depth_bin"] != "nan"]
    s = s.groupby("depth_bin", group_keys=False).sample(n=SAMPLE_PER_BIN, random_state=12).sort_values("tic")
    OUT.mkdir(parents=True, exist_ok=True)
    s.to_csv(path, index=False)
    return s


def refine_ephemeris(t, f, P0, t0, dur_d, half=None):
    """Best stacked box over a period grid around the starting period."""
    width = dur_d / 8
    base = np.ptp(t)
    half = max(0.002 * P0, 0.003) if half is None else half
    step = width / max(base / P0, 1)
    best = None
    r = f - 1
    nb = 8  # bins in one duration
    for P in np.arange(P0 - half, P0 + half + step, step):
        ph = (t - t0 + P / 2) % P - P / 2
        m = np.abs(ph) < 1.5 * dur_d
        b = np.floor((ph[m] + 1.5 * dur_d) / width).astype(int)
        n = np.bincount(b, minlength=int(3 * nb) + 1)
        s = np.bincount(b, r[m], minlength=int(3 * nb) + 1)
        cn, cs = np.concatenate([[0], np.cumsum(n)]), np.concatenate([[0], np.cumsum(s)])
        ni, si = cn[nb:] - cn[:-nb], cs[nb:] - cs[:-nb]
        with np.errstate(invalid="ignore", divide="ignore"):
            dep = -si / ni
            q = dep * np.sqrt(ni)
        q[(ni < 5) | ~np.isfinite(q)] = -np.inf
        j = int(np.argmax(q))
        if best is None or q[j] > best["q"]:
            best = dict(q=float(q[j]), P=float(P), depth=float(dep[j]),
                        t0=float(t0 + (j + nb / 2) * width - 1.5 * dur_d))
    best["at_grid_edge"] = abs(best["P"] - P0) > 0.95 * half
    best["P_offset_d"] = best["P"] - P0
    return best


def step_download(args):
    m = star_set(args.set)
    (WORK / "lc").mkdir(exist_ok=True)

    def one(r):
        try:
            t, f, secs, _ = pl.all_lc_for_tic(int(r.tic), WORK, r.duration_h / 24)
            return dict(tic=int(r.tic), n_sectors=len(secs), n_points=t.size, error="")
        except Exception as e:
            return dict(tic=int(r.tic), n_sectors=0, n_points=0, error=str(e)[:200])

    with ThreadPoolExecutor(args.workers) as ex:
        rows = []
        for i, row in enumerate(ex.map(one, m.itertuples())):
            rows.append(row)
            print(i + 1, row, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / f"epoch_rule_download_log_{args.set}.csv", index=False)


def step_known(args):
    rng = np.random.default_rng(10)
    stars, epochs = [], []
    for r in star_set(args.set).itertuples():
        row = dict(tic=int(r.tic), period_cat_d=r.period_d, depth_ppm=r.depth_ppm, depth_err_ppm=r.depth_err,
                   duration_h=r.duration_h, skipped="")
        t, f = load_star(r.tic)
        if t is None:
            stars.append({**row, "skipped": "no light curve"})
            continue
        dur_d, d0, d0_err = r.duration_h / 24, r.depth_ppm * 1e-6, r.depth_err * 1e-6
        half = None if args.set == "multi" else max(5 * np.nan_to_num(r.period_err), 2e-4 * r.period_d)
        eph = refine_ephemeris(t, f, r.period_d, r.t0, dur_d, half)
        row.update(period_d=eph["P"], period_offset_d=eph["P_offset_d"], stack_depth_ppm=eph["depth"] * 1e6,
                   baseline_d=float(np.ptp(t)))
        row["at_grid_edge"] = eph["at_grid_edge"]
        # the BLS periods have to prove themselves first; CP and KP periods are not judged by this data
        if args.set == "multi" and eph["at_grid_edge"]:
            stars.append({**row, "skipped": "period at edge of refinement grid"})
            continue
        if args.set == "multi" and eph["depth"] < 0.5 * d0:
            stars.append({**row, "skipped": "stacked depth under half the catalogue depth"})
            continue
        P, t0 = eph["P"], eph["t0"]
        sig, edges = core.duration_scale_noise(t, f, dur_d), core.gap_edges(t)
        k = np.arange(int(np.floor((t.min() - t0) / P)) - 1, int(np.ceil((t.max() - t0) / P)) + 2)
        ep = t0 + k * P
        depth, frac, _ = core.box_stats(t, f, ep, dur_d)
        clean, tot, v = verdicts(depth, frac, core.edge_distance(edges, ep), sig, d0, d0_err)
        if clean.sum() < 3 or np.ptp(ep[clean]) < 300:
            stars.append({**row, "skipped": "under 3 clean transits or under 300 d apart"})
            continue
        i1, i2 = np.where(clean)[0][[0, -1]]
        t1, t2 = ep[i1] + rng.normal(0, TIMING_SIG_D), ep[i2] + rng.normal(0, TIMING_SIG_D)
        n_true = int(round((ep[i2] - ep[i1]) / P))
        at = aliases(t1, t2, dur_d, d0, d0_err, t, f, sig, edges)
        for name, rows in (("transit", np.where(clean)[0]),):
            for i in rows:
                if i in (i1, i2):
                    continue
                epochs.append(dict(tic=row["tic"], kind=name, epoch=ep[i], depth=depth[i], d0=d0, tot=tot[i],
                                   **{f"{u}_excludes": bool(v[u][0][i]) for u in RULES}))
        # away from transit: quarter-ish phases, so a secondary eclipse at 0.5 stays out
        off = np.concatenate([ep + 0.3 * P, ep + 0.7 * P])
        dep_o, fr_o, _ = core.box_stats(t, f, off, dur_d)
        cl_o, tot_o, v_o = verdicts(dep_o, fr_o, core.edge_distance(edges, off), sig, d0, d0_err)
        for i in np.where(cl_o)[0]:
            epochs.append(dict(tic=row["tic"], kind="flat", epoch=off[i], depth=dep_o[i], d0=d0, tot=tot_o[i],
                               **{f"{u}_excludes": bool(v_o[u][0][i]) for u in RULES}))
        multiple = n_true % at["n"] == 0  # P_true x m: every epoch is a real transit, data can't object
        truth = at["n"] == n_true
        row.update(n_true=n_true, t1=t1, t2=t2, n_aliases=len(at), n_clean_transits=int(clean.sum()),
                   epoch_snr_median=float(np.median(d0 / tot[clean])), sig_ppm=sig * 1e6)
        for u in RULES:
            row[f"true_{u}"] = at.loc[truth, u].iloc[0] if truth.any() else "below p_floor"
            row[f"false_survivors_{u}"] = int(((at[u] != "excluded") & ~multiple).sum())
            row[f"false_supported_{u}"] = int(((at[u] == "supported") & ~multiple).sum())
            row[f"multiples_excluded_{u}"] = int(((at[u] == "excluded") & multiple & ~truth).sum())
        stars.append(row)
        print(row["tic"], {u: (row[f"true_{u}"], row[f"false_survivors_{u}"]) for u in RULES}, flush=True)
    pd.DataFrame(stars).to_csv(OUT / f"epoch_rule_known_stars_{args.set}.csv", index=False)
    pd.DataFrame(epochs).to_csv(OUT / f"epoch_rule_known_epochs_{args.set}.csv", index=False)


def step_inject(args):
    rng = np.random.default_rng(11)
    c = pd.read_csv("candidates.csv")
    rows = []
    for tic in sorted(c.loc[c["type"] != "multi", "tic"].unique()):
        t0, f0 = load_star(tic)
        if t0 is None or np.ptp(t0) < 300:
            continue
        cs = c[c["tic"] == tic]
        dur_d = float(cs["duration_h"].median()) / 24
        keep = np.ones(t0.size, bool)
        for tc in cs["t0"]:  # the star's own known events
            keep &= np.abs(t0 - tc) > dur_d
        t, f = t0[keep], f0[keep]
        sig = core.duration_scale_noise(t, f, dur_d)
        edges = core.gap_edges(t)
        for _ in range(args.draws):
            P = float(np.exp(rng.uniform(np.log(5), np.log(60))))
            ph = rng.uniform(0, P)
            k = np.arange(int(np.floor((t.min() - ph) / P)) - 1, int(np.ceil((t.max() - ph) / P)) + 2)
            ep = ph + k * P
            depth, frac, _ = core.box_stats(t, f, ep, dur_d)
            ok = (frac >= core.CLEAN_FRAC) & (core.edge_distance(edges, ep) >= core.CLEAN_EDGE_D)
            ok[~np.isfinite(depth)] = False
            if ok.sum() < 3 or np.ptp(ep[ok]) < 300:
                continue
            i1, i2 = np.where(ok)[0][[0, -1]]
            n_true = int(round((ep[i2] - ep[i1]) / P))
            t1, t2 = ep[i1] + rng.normal(0, TIMING_SIG_D), ep[i2] + rng.normal(0, TIMING_SIG_D)
            for snr in SNRS:
                d0 = snr * sig
                for scale in (1.0, 0.7):
                    # the two known events keep full depth; the other transits are scaled
                    g = f.copy()
                    for j, e in enumerate(ep):
                        m = np.abs(t - e) < dur_d / 2
                        g[m] -= d0 * (1.0 if j in (i1, i2) else scale)
                    at = aliases(t1, t2, dur_d, d0, 0.05 * d0, t, g, sig, edges)
                    multiple = n_true % at["n"] == 0
                    truth = at["n"] == n_true
                    row = dict(tic=int(tic), period_d=P, snr=snr, depth_scale=scale, n_true=n_true,
                               n_aliases=len(at), n_clean_transits=int(ok.sum()))
                    for u in RULES:
                        row[f"true_{u}"] = at.loc[truth, u].iloc[0] if truth.any() else "below p_floor"
                        row[f"false_survivors_{u}"] = int(((at[u] != "excluded") & ~multiple).sum())
                    rows.append(row)
        print(tic, len(rows), flush=True)
    pd.DataFrame(rows).to_csv(OUT / "epoch_rule_injections.csv", index=False)


def step_candidates(args):
    passing = pl.load_passing(WORK, args.event_times)
    hits = pl.alias_events(passing)
    cov = pd.read_csv(WORK / "coverage.csv").set_index("tic")
    hyp = pd.read_csv(WORK / "aliases" / "hypotheses_summary.csv")
    rows, hrows, mismatches = [], [], 0
    for h in hits.itertuples():
        t, f = load_star(h.tic)
        dur_d, d0, d0_err = h.dur_ref_h / 24, h.depth_ref_ppm * 1e-6, h.depth_err_ppm * 1e-6
        info = cov.loc[h.tic] if h.tic in cov.index else {}
        mstar = float(info.get("mstar", np.nan)) if len(info) else np.nan
        rstar = float(info.get("rstar", np.nan)) if len(info) else np.nan
        pmin = core.pmin_from_duration(dur_d, mstar, rstar, d0)
        sig, edges = core.duration_scale_noise(t, f, dur_d), core.gap_edges(t)
        at = aliases(h.t1, h.t2, dur_d, d0, d0_err, t, f, sig, edges, pmin=pmin)
        # the pipeline's rule here must reproduce its own table
        stored = pd.read_csv(WORK / "aliases" / f"{h.cand_id}.csv")
        diff = int((stored["status"].to_numpy() != at[PIPELINE_RULE].to_numpy()).sum()) if len(stored) == len(at) else -1
        mismatches += diff != 0
        row = dict(cand_id=h.cand_id, tic=int(h.tic), dT_d=h.t2 - h.t1, n_aliases=len(at),
                   epoch_snr=d0 / np.hypot(sig, d0_err), n_differ_from_pipeline=diff)
        for u in RULES:
            for s in ("allowed", "supported", "excluded"):
                row[f"{s}_{u}"] = int((at[u] == s).sum())
        rows.append(row)
        for q in hyp[hyp["cand_id"] == h.cand_id].itertuples():
            # the quoted period itself, which need not sit on this event's alias grid
            P = q.period_used_d
            k = np.arange(int(np.floor((t.min() - h.t1) / P)) - 1, int(np.ceil((t.max() - h.t1) / P)) + 2)
            ep = h.t1 + k * P
            ep = ep[(np.abs(ep - h.t1) >= dur_d / 2) & (np.abs(ep - h.t2) >= dur_d / 2)]
            depth, frac, _ = core.box_stats(t, f, ep, dur_d)
            _, _, v = verdicts(depth, frac, core.edge_distance(edges, ep), sig, d0, d0_err)
            st = {u: "excluded" if v[u][0].any() else "supported" if v[u][1].any() else "allowed" for u in RULES}
            hrows.append(dict(cand_id=h.cand_id, tic=int(h.tic), period_used_d=P, **st, pipeline_status=q.status_now,
                              n_excl_current=int(v["current"][0].sum()), n_excl_both=int(v["both"][0].sum()),
                              n_excl_half3s=int(v["half3s"][0].sum())))
    pd.DataFrame(rows).to_csv(OUT / "epoch_rule_candidates.csv", index=False)
    pd.DataFrame(hrows).to_csv(OUT / "epoch_rule_candidate_periods.csv", index=False)
    print(f"{len(rows)} alias events; {mismatches} differ from the pipeline's stored alias table")


def step_summary(args):
    lines = []
    for name in ("toi", "multi"):
        s = pd.read_csv(OUT / f"epoch_rule_known_stars_{name}.csv")
        used = s[s["skipped"].isna()]
        lines.append(f"KNOWN STARS, set {name}: {len(s)} total, {len(used)} used")
        lines.append(s["skipped"].value_counts().to_string())
        lines.append(f"median epoch SNR quartiles {used['epoch_snr_median'].quantile([.25, .5, .75]).round(1).tolist()}")
        for u in RULES:
            vc = used[f"true_{u}"].value_counts().to_dict()
            fs = used[f"false_survivors_{u}"]
            lines.append(f"{u:8s} true period {vc}  false survivors median {fs.median():.0f} mean {fs.mean():.1f} "
                         f"max {fs.max()}  of median {used['n_aliases'].median():.0f} aliases")
        e = pd.read_csv(OUT / f"epoch_rule_known_epochs_{name}.csv")
        e["snr"] = e["d0"] / e["tot"]
        e["snr_bin"] = pd.cut(e["snr"], SNR_BINS)
        e["ratio"] = e["depth"] / e["d0"]
        lines.append("epoch rates on real data (fraction of clean epochs that exclude)")
        g = e.groupby(["kind", "snr_bin"], observed=True)
        tab = g[[f"{u}_excludes" for u in RULES]].mean().round(4)
        tab["n"] = g.size()
        tab["n_stars"] = g["tic"].nunique()
        tab["ratio_p05"] = g["ratio"].quantile(0.05).round(2)
        tab["ratio_p50"] = g["ratio"].median().round(2)
        lines.append(tab.to_string())
        lines.append("median across stars, so stars with many epochs don't dominate")
        per = e.groupby(["tic", "kind"])[[f"{u}_excludes" for u in RULES]].mean().groupby("kind")
        lines.append(per.median().round(4).to_string() + "\n")
    lines.append("\nINJECTIONS")
    i = pd.read_csv(OUT / "epoch_rule_injections.csv")
    for scale, gi in i.groupby("depth_scale"):
        rows = []
        for snr, gs in gi.groupby("snr"):
            row = dict(snr=snr, n=len(gs))
            for u in RULES:
                row[f"true_excl_{u}"] = round((gs[f"true_{u}"] == "excluded").mean(), 3)
            for u in RULES:
                row[f"false_surv_{u}"] = round(gs[f"false_survivors_{u}"].median(), 1)
            rows.append(row)
        lines.append(f"other transits at {scale} x reference depth:")
        lines.append(pd.DataFrame(rows).to_string(index=False))
    c = pd.read_csv(OUT / "epoch_rule_candidates.csv")
    lines.append(f"\nCANDIDATES: {len(c)} alias events, {(c['n_differ_from_pipeline'] != 0).sum()} differ from the "
                 f"pipeline's stored table under {PIPELINE_RULE}")
    for u in RULES:
        ok = c[f"allowed_{u}"] + c[f"supported_{u}"]
        lines.append(f"{u:8s} allowed+supported median {ok.median():.0f} mean {ok.mean():.1f}; "
                     f"events with none left {(ok == 0).sum()}; with exactly 1 {(ok == 1).sum()}")
    p = pd.read_csv(OUT / "epoch_rule_candidate_periods.csv")
    lines.append(f"\nQUOTED PERIODS: {(p[PIPELINE_RULE] != p['pipeline_status']).sum()} of {len(p)} differ from the "
                 f"pipeline's status under {PIPELINE_RULE}")
    lines.append(p.to_string(index=False))
    text = "\n".join(lines)
    print(text)
    (OUT / "epoch_rule_test_summary.txt").write_text(text + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["download", "known", "inject", "candidates", "summary"])
    ap.add_argument("--set", choices=["toi", "multi"], default="toi", help="star set for download and known")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--draws", type=int, default=2)
    ap.add_argument("--event-times", default="results/round2/qlp_s94_timing.csv")
    args = ap.parse_args()
    globals()[f"step_{args.mode}"](args)


if __name__ == "__main__":
    main()
