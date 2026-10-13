"""
AAS82207 revision, editor Q6: audit of the 34 stars whose catalogued CTOI event was not
recovered by the pipeline's event check (pilot.py -> alias_pipeline.test_window at epoch_btjd).

For each star:
  - which sector(s) cover the catalogued epoch, which products MAST has there (fresh search,
    search_products.py), which product the pipeline loaded;
  - measure the event in EVERY product and flux column (SAP and PDCSAP etc.) and in our own
    TESScut aperture photometry: pipeline-equivalent test (quality==0, +-0.03 d), the same with
    quality-flagged cadences kept, and the strongest box dip within +-1 d (time offset);
  - a single-event box search over the whole epoch sector (strongest dip, number of dips of the
    catalogued depth) in the pipeline product and in TESScut;
  - the CTOI's ExoFOP record (cached target JSON in data/raw/exofop_targets);
  - a diagnostic plot of all products +-1 d around the epoch.

Usage:  python singles/revision/audit_34.py [TIC ...]
Output: singles/revision/audit_measurements.csv (one row per star x sector x product x column),
        singles/revision/audit_sector_scan.csv, singles/revision/plots/<TIC>_audit.png,
        singles/revision/audit_meta.json
The per-star reason categories are assigned in classify_audit.py from these outputs.
"""
import json
import sys
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REV = Path(__file__).resolve().parent
sys.path.insert(0, str(REV))
import hlsp_io  # noqa: E402
import revlib as rl  # noqa: E402

ROOT = rl.ROOT
PLOTS = REV / "plots"
PLOTS.mkdir(exist_ok=True)

# which column the pipeline used: 'flux' of the stored CSV = PDCSAP (SPOC, TESS-SPOC), SAP (QLP)
PIPE_COL = {"SPOC": "pdcsap_flux", "TESS-SPOC": "pdcsap_flux", "QLP": "sap_flux", "OWN": "flux"}
PLOT_COLS = {"T16": ["tfa1", "irm1"], "CDIPS": ["tfa1", "irm1"], "TASOC": ["flux_corr", "flux_raw"]}


def sector_ranges():
    o = pd.read_csv(ROOT / "data" / "raw" / "tess_orbit_times.csv", comment="#")
    for c in ["Start of Orbit", "End of Orbit"]:
        o[c] = (pd.to_datetime(o[c]) - pd.Timestamp("2000-01-01 12:00")).dt.total_seconds() / 86400 + 2451545 - 2457000
    s = o.groupby("Sector").agg(start=("Start of Orbit", "min"), end=("End of Orbit", "max"))
    s.index = s.index.astype(int)
    return s


def sector_of(T, rng, pad=0.0):
    m = rng[(rng.start - pad <= T) & (rng.end + pad >= T)]
    return [int(x) for x in m.index]


def load_products(tic, sectors, search):
    """[(sector, author, exptime, DataFrame, [flux cols], source)] for the given sectors."""
    loc = rl.local_products(tic)
    out, seen = [], set()
    # local CSVs (what the pipeline had)
    for (sec, author, exp), f in sorted(loc.items()):
        if sec not in sectors:
            continue
        df = pd.read_csv(f)
        if author == "OWN":
            cols = ["flux"]
        else:
            cols = [c for c in rl.FLUX_COLS.get(author, []) if c in df.columns] or ["flux"]
        out.append((sec, author, exp, df, cols, "local:" + f.name))
        seen.add((sec, author, exp))
    # anything else MAST has for those sectors
    for _, r in search.iterrows():
        try:
            sec = int(str(r.mission).split()[-1])
        except ValueError:
            continue
        exp = int(round(float(r.exptime)))
        if sec not in sectors or (sec, r.author, exp) in seen:
            continue
        if r.author in ("SPOC", "TESS-SPOC", "QLP", "TGLC", "GSFC-ELEANOR-LITE") and (sec, r.author, exp) in seen:
            continue
        try:
            p = hlsp_io.download(r.dataURI, tic)
            df, cols = hlsp_io.read_any(p)
        except Exception as e:  # noqa: BLE001
            print(f"   {r.author} S{sec}: read failed {e}")
            continue
        if "flux_corr" in cols and abs(np.nanmedian(df.flux_corr)) < 1e3:  # TASOC: ppm about 0
            df["flux_corr"] = 1 + df.flux_corr * 1e-6
        if r.author in PLOT_COLS:
            cols = [c for c in PLOT_COLS[r.author] if c in cols] or cols[:2]
        tag = f"{r.author}" + (f"_{Path(r.dataURI).name.split('_')[-1].replace('.fits', '')}" if r.author == "TASOC" else "")
        key = (sec, tag, exp)
        if key in seen:
            continue
        seen.add(key)
        out.append((sec, tag, exp, df, cols, "mast:" + Path(r.dataURI).name))
    return out


def audit_star(row, cl_row, rng):
    tic = int(row.TIC)
    epoch, D, delta = float(row.epoch_btjd), row.dur_h / 24, row.depth_ppm * 1e-6
    search = pd.read_csv(REV / "data" / "search" / f"{tic}_lc.csv")
    cut = pd.read_csv(REV / "data" / "search" / f"{tic}_tesscut.csv")
    cut_secs = sorted({int(m.split()[-1]) for m in cut.mission if int(m.split()[-1]) < 200})
    lc_secs = sorted({int(str(m).split()[-1]) for m in search.mission if int(str(m).split()[-1]) < 200})
    ep_secs = sector_of(epoch, rng, pad=1.0)
    obs_secs = [s for s in ep_secs if s in cut_secs or s in lc_secs]
    meta = dict(TIC=tic, CTOI=row.CTOI, user=row.User, epoch_btjd=epoch, dur_h=row.dur_h, depth_ppm=row.depth_ppm,
                Tmag=row.Tmag, epoch_sector=" ".join(map(str, ep_secs)), observed_epoch_sectors=obs_secs,
                tesscut_sectors=cut_secs, lc_sectors=lc_secs,
                pipeline_product={s: rl.pipeline_product(tic, s) for s in obs_secs},
                products_at_epoch_sector=sorted({f"{r.author}_{int(r.exptime)}" for _, r in search.iterrows()
                                                 if int(str(r.mission).split()[-1]) in obs_secs}))
    # ExoFOP record
    xf = ROOT / "data" / "raw" / "exofop_targets" / f"{tic}.json"
    if xf.exists():
        x = json.loads(xf.read_text())
        pp = [p for p in x["planet_parameters"][1:] if str(row.CTOI)[:len(str(tic)) + 3] in str(p.get("name"))]
        meta["exofop_ctoi"] = [dict(ctoi=c["ctoi"], notes=c.get("cnotes"), date=c.get("cdate"), disp=c.get("cdisp"))
                               for c in x["ctois"]]
        meta["exofop_params"] = [{k: p.get(k) for k in ("name", "epoch", "epoch_e", "per", "dep_p", "dur", "snr",
                                                        "pnotes", "puser", "pdate")} for p in pp]
        meta["exofop_files"] = [f.get("fname") or f.get("file_name") or str(f)[:120] for f in x.get("files", [])]
        meta["exofop_tois"] = [t.get("toi") for t in x.get("tois", [])]
    rows, scans, series = [], [], []
    for sec in obs_secs:
        prods = load_products(tic, [sec], search)
        try:
            own, info = rl.own_photometry(tic, sec, row.RA, row.Dec)
        except Exception as e:  # noqa: BLE001
            own, info = None, None
            print(f"   own photometry S{sec} failed: {e}")
        if own is not None:
            prods.append((sec, "TESScut-own", 0, own.rename(columns={"flux": "aper_flux"}), ["aper_flux"],
                          f"tesscut {Path(info['path']).name} n_ap={info['n_ap']}"))
        pipe_prod = meta["pipeline_product"].get(sec)
        for sec_, author, exp, df, cols, src in prods:
            q = df["quality"].to_numpy() if "quality" in df else np.zeros(len(df))
            qok = (q == 0)
            if author == "TESScut-own":
                qok = df["good"].to_numpy()
            pkey = f"SPOC_{exp}" if author == "SPOC" else author
            for c in cols:
                if c not in df:
                    continue
                m = rl.measure(df.time, df[c], epoch, D, delta, qok=qok, qlp_trim=(author == "QLP"))
                is_pipe = (pkey == pipe_prod and c == PIPE_COL.get(author.split("_")[0], "flux"))
                rows.append(dict(TIC=tic, sector=sec_, product=author, exptime=exp, column=c, source=src,
                                 pipeline_used=is_pipe, **m))
                series.append((f"S{sec_} {author}{'' if not exp else f' {exp}s'} {c}", df.time.to_numpy(),
                               df[c].to_numpy(dtype=float), qok, exp))
                if is_pipe or (author == "TESScut-own"):
                    s = rl.sector_scan(df.time, df[c], D, delta, q_ok=qok)
                    scans.append(dict(TIC=tic, sector=sec_, product=author, column=c, **s))
    if not obs_secs:
        # epoch not observed: search every sector the pipeline had for a dip of the catalogued
        # depth/duration (typo or wrong-sector hypothesis)
        lc = rl.ap.load_lightcurves(tic, products=rl.ap.PRODUCT_RANK + ["OWN"])
        for sec, g in lc.groupby("sector"):
            s = rl.sector_scan(g.time, g.flux, D, delta)
            scans.append(dict(TIC=tic, sector=int(sec), product=g["product"].iloc[0], column="flux",
                              note="epoch-not-observed", **s))
    plot_star(tic, row, series, epoch, D, delta)
    return meta, rows, scans


def plot_star(tic, row, series, epoch, D, delta):
    n = max(1, len(series))
    fig, ax = plt.subplots(figsize=(8, 1.2 + 0.75 * n))
    step = None
    offs = []
    for i, (lab, t, f, qok, exp) in enumerate(series):
        m = np.isfinite(t) & np.isfinite(f) & (np.abs(t - epoch) < 1.0)
        if m.sum() < 5:
            offs.append(None)
            continue
        ff = f[m] / np.nanmedian(f[m & qok]) if (m & qok).sum() > 5 else f[m] / np.nanmedian(f[m])
        # display only: divide by a quadratic fitted to good data outside the event window
        tt0 = t[m] - epoch
        fitm = qok[m] & (np.abs(tt0) > D / 2 + 0.05) & (np.abs(ff - 1) < max(0.3, 5 * delta))
        if fitm.sum() > 10:
            ff = ff / np.polyval(np.polyfit(tt0[fitm], ff[fitm], 2), tt0)
        good = qok[m] & np.isfinite(ff)
        sd = 1.4826 * np.nanmedian(np.abs(np.diff(ff[good]))) / np.sqrt(2) if good.sum() > 5 else delta
        step = max(step or 0, 2.5 * delta, 8 * sd)
        offs.append((t[m], ff, qok[m], lab, exp))
    step = step or 3 * delta
    k = 0
    for o in offs:
        if o is None:
            continue
        tt, ff, qq, lab, exp = o
        y0 = -k * step
        x = (tt - epoch) * 24
        ax.plot(x[qq], (ff[qq] - 1) + y0, ".", ms=1.5, color="0.55")
        ax.plot(x[~qq], (ff[~qq] - 1) + y0, "x", ms=2.5, color="tab:orange")
        if exp and exp < 600:
            nb = max(1, int(round(1800 / exp)))
            ok = qq
            xs, ys = x[ok], ff[ok]
            if len(xs) > nb:
                bins = np.arange(-24, 24.01, 0.5)
                idx = np.digitize(xs, bins)
                bx = [xs[idx == j].mean() for j in np.unique(idx) if (idx == j).sum() >= nb // 2]
                by = [ys[idx == j].mean() for j in np.unique(idx) if (idx == j).sum() >= nb // 2]
                ax.plot(bx, np.array(by) - 1 + y0, "o", ms=2.5, color="k")
        ax.plot([-D * 12, -D * 12, D * 12, D * 12], [y0, y0 - delta, y0 - delta, y0], "r-", lw=1)
        ax.text(-23.8, y0 + 0.25 * step, lab, fontsize=6.5, va="bottom")
        k += 1
    ax.axvspan(-D * 12, D * 12, color="r", alpha=0.06)
    ax.set_xlim(-24, 24)
    ax.set_ylim(-(k - 0.3) * step - delta, step * 0.9)
    ax.set_xlabel(f"hours from catalogued epoch BTJD {epoch:.4f}")
    ax.set_ylabel("relative flux (offset)")
    ax.set_title(f"TIC {tic} CTOI {row.CTOI:.2f} ({row.User}): depth {row.depth_ppm:.0f} ppm, D {row.dur_h:.2f} h"
                 + (" (display: quadratic baseline removed; orange x = quality-flagged)" if series else "  -- NO DATA AT EPOCH"), fontsize=8)
    if not series:
        ax.text(0, 0, "no TESS data in any product within +-1 d of the epoch", ha="center")
    fig.tight_layout()
    fig.savefig(PLOTS / f"{tic}_audit.png", dpi=110)
    plt.close(fig)


def main():
    cl = pd.read_csv(ROOT / "singles" / "results" / "final_classification_2026-10-11_rnaas.csv")
    bad = cl[cl.cls == "event-not-recovered"]
    tics = [int(x) for x in sys.argv[1:]] or bad.TIC.tolist()
    cs = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv")
    rng = sector_ranges()
    metas, allrows, allscans = {}, [], []
    mf = REV / "audit_meta.json"
    if mf.exists():
        metas = json.loads(mf.read_text())
    for tic in tics:
        crow = cl[cl.TIC == tic].iloc[0]
        row = cs[(cs.TIC == tic)].iloc[0]
        print(f"TIC {tic}", flush=True)
        try:
            meta, rows, scans = audit_star(row, crow, rng)
        except Exception:  # noqa: BLE001
            print(traceback.format_exc())
            continue
        metas[str(tic)] = meta
        allrows += rows
        allscans += scans
    mf.write_text(json.dumps(metas, indent=1, default=str))
    out = REV / "audit_measurements.csv"
    new = pd.DataFrame(allrows)
    if out.exists() and len(sys.argv) > 1:
        old = pd.read_csv(out)
        new = pd.concat([old[~old.TIC.isin(tics)], new], ignore_index=True)
    new.to_csv(out, index=False)
    out2 = REV / "audit_sector_scan.csv"
    new2 = pd.DataFrame(allscans)
    if out2.exists() and len(sys.argv) > 1:
        old = pd.read_csv(out2)
        new2 = pd.concat([old[~old.TIC.isin(tics)], new2], ignore_index=True)
    new2.to_csv(out2, index=False)


if __name__ == "__main__":
    main()
