"""
One vetting sheet per event (prereg_catalogue.yaml, vetting.sheet): a PNG and a JSON of the numbers.

The sheet carries no star name, no period and no alias table. Its file name is a hash of the event.
Panels, the same for every event:
  A  the reference and the event, raw (as loaded, quality-masked, 10-min bins) and detrended, one scale
  B  the event's times in every other product of that sector
  C  SAP beside the pipeline flux of the product used, where the file has both
  D  the event's sector: the whole light curve, flagged cadences (quality != 0, with the momentum-dump
     bit named), and the distance of the event from the sector's edges and gaps
  E  the event's times on up to five nearest stars of similar brightness (Tmag within 1) on the same
     camera and CCD, same sector and same product type
  F  a difference image if a TESScut cut-out of that sector is already on disk; else "no difference
     image"
  G  depth, duration and SNR of the event against the reference

`inject` (calibration events only): a function (sector, author, t, f) -> f applied to every light
curve of the target shown on the sheet, so that an injected transit appears in every product.

Inputs  : MAST light curves and the TIC (network), cut-outs under CUTOUT_ROOTS (never downloaded here)
Outputs : <out>/<sheet_id>.png, <out>/<sheet_id>.json
Usage   : imported by make_sheets.py
"""
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "catalogue"))
import catrun  # noqa: E402

pl, core = catrun.mc.pl, catrun.core
CUTOUT_ROOTS = [Path("/Users/mcbee/Documents/GitHub/coracle/duo-transits-cache/raw/tesscut"),
                Path("/Users/mcbee/Documents/GitHub/coracle/salinas-blind-periods/data/tesscut")]
N_NEIGH, NEIGH_TRIES, NEIGH_RADIUS_DEG, NEIGH_DMAG = 5, 20, 1.0, 1.0
MOMENTUM_DUMP_BIT = 32


def sheet_id(event_id, salt="catalogue-v0"):
    return hashlib.sha1(f"{salt}|{event_id}".encode()).hexdigest()[:12]


def _window(t, f, tc, half):
    m = np.abs(t - tc) <= half
    return t[m] - tc, f[m]


def _bin10(t, f):
    """10-min bins, as the pipeline bins faster cadences."""
    if t.size < 3 or np.nanmedian(np.diff(t)) >= 9 / 1440:
        return t, f
    idx = np.floor((t - t.min()) / (10 / 1440)).astype(int)
    n = np.bincount(idx)
    ok = n > 0
    return (np.bincount(idx, t) / np.maximum(n, 1))[ok], (np.bincount(idx, f) / np.maximum(n, 1))[ok]


def _local(x, y, dur):
    """Divide by the median outside 1.5 durations of the marked time (a product with its own trend)."""
    out = np.abs(x) > 1.5 * dur
    m = np.nanmedian(y[out]) if out.sum() >= 5 else np.nanmedian(y) if y.size else 1.0
    return y / m if np.isfinite(m) and m != 0 else y


def _raw_file(sr, row, cache, **kw):
    with pl.DOWNLOAD_LOCK:
        return sr[row["idx"]].download(download_dir=str(cache), **kw)


def neighbours(ra, dec, tmag, tic, sector, author, cam_ccd, cache):
    """Up to N_NEIGH (tic, sep_arcmin, tmag, t, f) of the nearest stars within NEIGH_DMAG of the
    target's Tmag that have a light curve of the same author in the sector, on the same camera/CCD."""
    import lightkurve as lk
    from astroquery.mast import Catalogs
    r = Catalogs.query_criteria(catalog="Tic", coordinates=f"{ra} {dec}", radius=NEIGH_RADIUS_DEG,
                                Tmag=[tmag - NEIGH_DMAG, tmag + NEIGH_DMAG]).to_pandas()
    r["dstArcSec"] = 3600 * np.hypot((r.ra.astype(float) - ra) * np.cos(np.radians(dec)), r.dec.astype(float) - dec)
    r = r[r.ID.astype(int) != int(tic)].sort_values("dstArcSec")
    out, tried = [], 0
    for x in r.itertuples():
        if len(out) >= N_NEIGH or tried >= NEIGH_TRIES:
            break
        tried += 1
        try:
            sr = lk.search_lightcurve(f"TIC {int(x.ID)}", mission="TESS", sector=int(sector), author=author)
            if len(sr) == 0:
                continue
            with pl.DOWNLOAD_LOCK:
                lc = sr[0].download(quality_bitmask="default", download_dir=str(cache))
            if cam_ccd and (lc.meta.get("CAMERA"), lc.meta.get("CCD")) != cam_ccd:
                continue
            if author == "QLP":
                lc = lc[lc.quality == 0]
            lc = lc.remove_nans().normalize()
            out.append((int(x.ID), float(x.dstArcSec) / 60, float(x.Tmag),
                        np.asarray(lc.time.value, float), np.asarray(lc.flux.value, float)))
        except Exception:
            continue
    return out, tried


def cached_cutout(ra, dec, sector):
    for root in CUTOUT_ROOTS:
        if root.exists():
            hits = sorted(root.glob(f"**/tess-s{int(sector):04d}-*_{ra:.6f}_{dec:.6f}_*astrocut.fits"))
            if hits:
                return hits[0]
    return None


def make_sheet(d, ev, ref, sectors, sr, prod_rows, info, out_dir, cache, inject=None, salt="catalogue-v0"):
    """d: definition; ev: the event (catrun.search_row); ref: dict(depth, depth_err, dur_d, snr, measured,
    sector) of the reference; sectors: the star's loaded sectors [(sector, author, t, f)] (already
    injected, if inject is given); sr, prod_rows: pl.lc_products(tic); info: pl.tic_info(tic)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sid = sheet_id(ev["event_id"], salt)
    dur, tc, s_ev = ref["dur_d"], ev["tc"], int(ev["sector"])
    half = max(1.0, 4 * dur)
    by = {s: (a, t, f) for s, a, t, f in sectors}
    notes = {}
    fig = plt.figure(figsize=(16, 15))
    gs = fig.add_gridspec(4, 4, hspace=0.42, wspace=0.28)
    depth_scale = max(ref["depth"], ev["depth"], 1e-5)
    ylim = (1 - 2.5 * depth_scale, 1 + 1.5 * depth_scale)

    def lcpanel(ax, t, f, t0, title, det_dur=None, local=False):
        if det_dur is not None:
            try:
                t, f = pl.detrend(t, f, det_dur)
            except Exception as e:
                ax.set_title(f"{title}: detrend failed")
                return
        x, y = _window(t, f, t0, half)
        if local:
            y = _local(x, y, dur)
        ax.plot(x, y, ".", ms=3, color="k")
        ax.axvspan(-dur / 2, dur / 2, color="C1", alpha=0.2)
        ax.axhline(1, color="0.6", lw=0.6)
        ax.axhline(1 - ref["depth"], color="C0", lw=0.6, ls="--")
        ax.set_xlim(-half, half)
        ax.set_ylim(*ylim)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("days from the marked time", fontsize=8)

    # A: reference and event, raw and detrended
    a_ref, t_ref, f_ref = by[ref["sector"]]
    a_ev, t_e, f_e = by[s_ev]
    lcpanel(fig.add_subplot(gs[0, 0]), t_ref, f_ref, d["t_cat"], f"A1 REFERENCE raw ({a_ref})")
    lcpanel(fig.add_subplot(gs[0, 1]), t_ref, f_ref, d["t_cat"], "A2 REFERENCE detrended", dur)
    lcpanel(fig.add_subplot(gs[0, 2]), t_e, f_e, tc, f"A3 EVENT raw ({a_ev})")
    lcpanel(fig.add_subplot(gs[0, 3]), t_e, f_e, tc, "A4 EVENT detrended", dur)

    # B: other products of the event's sector
    best = pl.best_per_sector(prod_rows).get(s_ev)
    others = [r for r in prod_rows if r["sector"] == s_ev and r is not best]
    seen, shown = set(), []
    for r in sorted(others, key=lambda r: (pl.AUTHOR_RANK[r["author"]], r["exptime"])):
        key = (r["author"], r["exptime"])
        if key in seen or len(shown) >= 2:
            continue
        seen.add(key)
        try:
            t, f = pl.load_sector_lc(sr, r, cache)
            if inject:
                f = inject(s_ev, r["author"], t, f)
            shown.append((f"{r['author']} {r['exptime']:.0f} s", t, f))
        except Exception as e:
            notes.setdefault("other_product_errors", []).append(str(e)[:100])
    for i in range(2):
        ax = fig.add_subplot(gs[1, i])
        if i < len(shown):
            lcpanel(ax, shown[i][1], shown[i][2], tc, f"B{i + 1} EVENT in another product: {shown[i][0]}", local=True)
        else:
            ax.text(0.5, 0.5, "no other product for this sector" if not shown else "no further product",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"B{i + 1}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
    notes["n_other_products"] = len(shown)

    # C and D need the file itself (all cadences, all columns)
    ax_c, ax_d = fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[2, 0:2])
    cam_ccd, flagged = None, []
    try:
        lc = _raw_file(sr, best, cache, quality_bitmask="none")
        cam_ccd = (lc.meta.get("CAMERA"), lc.meta.get("CCD"))
        t_all = np.asarray(lc.time.value, float)
        q = np.asarray(lc.quality, int)
        near = np.isfinite(t_all) & (np.abs(t_all - tc) <= 2 * dur)
        flagged = sorted(set(int(v) for v in q[near] if v != 0))
        notes.update(quality_values_within_2_durations=flagged,
                     n_flagged_within_2_durations=int((q[near] != 0).sum()),
                     momentum_dump_within_2_durations=bool(((q[near] & MOMENTUM_DUMP_BIT) != 0).any()))
        col = "sap_flux" if "sap_flux" in lc.columns else None
        if col:
            sap = np.asarray(lc[col].value, float)
            ok = np.isfinite(t_all) & np.isfinite(sap) & (q == 0)
            ts, fs = t_all[ok], sap[ok] / np.nanmedian(sap[ok])
            if inject:
                fs = inject(s_ev, best["author"], ts, fs)
            x, y = _window(*_bin10(ts, fs), tc, half)
            ax_c.plot(x, _local(x, y, dur), ".", ms=3, color="C3", label="SAP (10-min bins)")
            x2, y2 = _window(t_e, f_e, tc, half)
            ax_c.plot(x2, y2 + 2.0 * depth_scale, ".", ms=3, color="k", label="pipeline flux (offset up)")
            ax_c.axvspan(-dur / 2, dur / 2, color="C1", alpha=0.2)
            ax_c.set_xlim(-half, half)
            ax_c.set_ylim(1 - 2.5 * depth_scale, 1 + 3.5 * depth_scale)
            ax_c.legend(fontsize=7, loc="lower left")
            notes["sap_shown"] = True
        else:
            ax_c.text(0.5, 0.5, "file has no SAP column", ha="center", va="center", transform=ax_c.transAxes)
            notes["sap_shown"] = False
        dump_t = t_all[np.isfinite(t_all) & ((q & MOMENTUM_DUMP_BIT) != 0)]
        bad_t = t_all[np.isfinite(t_all) & (q != 0)]
        ax_d.plot(t_e, f_e, ".", ms=2, color="k")
        lo = np.nanpercentile(f_e, 0.5)
        ax_d.plot(bad_t, np.full(bad_t.size, lo), "|", color="C3", ms=8, label="flagged cadence")
        for x in dump_t:
            ax_d.axvline(x, color="C2", lw=0.5, alpha=0.5)
        ax_d.axvline(tc, color="C1", lw=1.5, label="event")
        if ref["sector"] == s_ev:
            ax_d.axvline(d["t_cat"], color="C0", lw=1.5, ls="--", label="reference")
        ax_d.legend(fontsize=7, loc="upper right")
    except Exception as e:
        notes["file_error"] = str(e)[:150]
        ax_c.text(0.5, 0.5, "file could not be read", ha="center", va="center", transform=ax_c.transAxes)
        ax_d.plot(t_e, f_e, ".", ms=2, color="k")
        ax_d.axvline(tc, color="C1", lw=1.5)
    ax_c.set_title("C SAP beside the pipeline flux, at the EVENT", fontsize=9)
    edges = core.gap_edges(t_e)
    notes.update(days_from_sector_start=float(tc - t_e.min()), days_to_sector_end=float(t_e.max() - tc),
                 days_to_nearest_gap_or_edge=float(core.edge_distance(edges, np.array([tc]))[0]))
    ax_d.set_title("D EVENT's sector, raw (green lines: momentum dumps; red ticks: flagged cadences)", fontsize=9)
    ax_d.set_xlabel("BTJD", fontsize=8)

    # E: neighbours
    ax_e = fig.add_subplot(gs[2, 2:4])
    try:
        nb, tried = neighbours(info["ra"], info["dec"], info["tmag"], d["tic"], s_ev, best["author"], cam_ccd, cache)
    except Exception as e:
        nb, tried = [], 0
        notes["neighbour_error"] = str(e)[:150]
    for i, (_, sep, tm, t, f) in enumerate(nb):
        x, y = _window(*_bin10(t, f), tc, half)
        if x.size:
            ax_e.plot(x, _local(x, y, dur) + i * 2.5 * depth_scale, ".", ms=2, color=f"C{i}",
                      label=f"{sep:.1f}' away, Tmag {tm:.1f}")
    ax_e.axvspan(-dur / 2, dur / 2, color="C1", alpha=0.2)
    ax_e.set_xlim(-half, half)
    ax_e.set_ylim(1 - 2.5 * depth_scale, 1 + (max(len(nb), 1) * 2.5 + 1) * depth_scale)
    if nb:
        ax_e.legend(fontsize=7, loc="upper right", ncol=2)
    ax_e.set_title(f"E EVENT's times on {len(nb)} neighbouring stars, same camera and CCD (offset)",
                   fontsize=9)
    notes.update(n_neighbours=len(nb), n_neighbours_tried=tried)

    # F: difference image, only from a cut-out already on disk
    ax_f = fig.add_subplot(gs[1, 3])
    hit = cached_cutout(info["ra"], info["dec"], s_ev) if not inject else None
    notes["difference_image"] = False
    if hit is not None:
        try:
            import lightkurve as lk
            sys.path.insert(0, str(ROOT / "salinas-second-transits"))
            import tesscut_phot as tp
            tpf = lk.read(hit)
            diff, _, n_in, n_out = tp.difference_image(tpf, tc, dur)
            if diff is not None:
                ax_f.imshow(diff, origin="lower", cmap="RdBu_r")
                c, r = tp.target_pixel(tpf, info["ra"], info["dec"])
                ax_f.plot(c, r, "k+", ms=14)
                notes["difference_image"] = True
        except Exception as e:
            notes["difference_image_error"] = str(e)[:150]
    if not notes["difference_image"]:
        ax_f.text(0.5, 0.5, "no difference image", ha="center", va="center", transform=ax_f.transAxes)
        ax_f.set_xticks([]); ax_f.set_yticks([])
    ax_f.set_title("F difference image (out minus in; + is the target)", fontsize=9)

    # G: the numbers
    num = dict(
        sheet_id=sid, event_sector=s_ev, event_product=str(ev["author"]), reference_sector=int(ref["sector"]),
        reference_measured=bool(ref["measured"]),
        reference=dict(depth_ppm=ref["depth"] * 1e6, depth_err_ppm=ref["depth_err"] * 1e6,
                       duration_h=ref["dur_d"] * 24, snr=ref["snr"]),
        event=dict(depth_ppm=ev["depth"] * 1e6, depth_err_ppm=ev["depth_err"] * 1e6, duration_h=ev["dur_d"] * 24,
                   snr=ev["snr"], depth_ratio=ev["depth_ratio"], duration_ratio=ev["dur_ratio"],
                   sector_box_snr=ev.get("sec_box_snr"), strongest_brightening_snr=ev["anti_snr"],
                   coverage=ev["frac"]),
        event_at_another_catalogued_transit_on_this_star=bool(ev.get("at_other_row", False)),
        target_tmag=info.get("tmag"), notes=notes)
    ax_g = fig.add_subplot(gs[3, :])
    ax_g.axis("off")
    r_, e_ = num["reference"], num["event"]
    txt = (f"G  sheet {sid}    target Tmag {info.get('tmag', float('nan')):.1f}\n"
           f"REFERENCE: depth {r_['depth_ppm']:.0f} +- {r_['depth_err_ppm']:.0f} ppm, duration {r_['duration_h']:.2f} h, "
           f"SNR {r_['snr']:.1f}" + ("" if ref["measured"] else "   (NOT MEASURED in the light curve: catalogue values)") + "\n"
           f"EVENT:     depth {e_['depth_ppm']:.0f} +- {e_['depth_err_ppm']:.0f} ppm, duration {e_['duration_h']:.2f} h, "
           f"SNR {e_['snr']:.1f}; depth ratio {e_['depth_ratio']:.2f}, duration ratio {e_['duration_ratio']:.2f}; "
           f"SNR against same-length boxes of its sector {e_['sector_box_snr']:.1f}; "
           f"strongest brightening in the sector SNR {e_['strongest_brightening_snr']:.1f}; coverage {e_['coverage']:.2f}\n"
           f"Flagged cadences within 2 durations: {notes.get('n_flagged_within_2_durations', 'unknown')} "
           f"(quality values {notes.get('quality_values_within_2_durations', 'unknown')}); "
           f"momentum dump within 2 durations: {notes.get('momentum_dump_within_2_durations', 'unknown')}\n"
           f"Event is {notes['days_from_sector_start']:.2f} d after its sector's first point, "
           f"{notes['days_to_sector_end']:.2f} d before its last, {notes['days_to_nearest_gap_or_edge']:.2f} d from the nearest gap or edge\n"
           f"Other products shown: {notes['n_other_products']}; SAP shown: {notes.get('sap_shown', False)}; "
           f"neighbours shown: {notes['n_neighbours']}; difference image: {notes['difference_image']}\n"
           f"Event coincides with another catalogued transit on this star: "
           f"{num['event_at_another_catalogued_transit_on_this_star']}\n"
           "Orange band: one reference duration at the marked time. Dashed blue line: the reference depth.")
    ax_g.text(0.0, 1.0, txt, va="top", ha="left", fontsize=10, family="monospace", transform=ax_g.transAxes)
    fig.savefig(out_dir / f"{sid}.png", dpi=80, bbox_inches="tight")
    plt.close(fig)
    (out_dir / f"{sid}.json").write_text(json.dumps(num, indent=1, default=float))
    return sid
