"""Difference imaging from TESScut cutouts.  Usage: s13_diffimg.py TIC

For each fully covered transit (ephemeris from out/ephem_<tic>.json):
  * per pixel, a linear-in-time baseline is fitted to out-of-transit frames
    (0.65-1.6 durations from mid-transit, both sides) and evaluated at the in-transit frames
    (|t-tc| < 0.3 duration). D = baseline - in-transit mean (positive = flux lost in transit);
    O = baseline evaluated at tc (the out-of-transit image).
  * Pixel errors for D come from the scatter of the out-of-transit residuals.
  * A scene model (Gaussian PSF at each Gaia DR3 star with G < 18.5, flux from G, positions from
    the cutout WCS) is fitted to O for a global shift (dx, dy), PSF width and scale.
  * D is then tested against "star k is the source" for every star within 200": PSF at star k
    (with the O-fitted shift/width) times a free amplitude plus a constant. chi2_k is reported.
  * D is also fitted with a free-position PSF: the offset of the transit source from the target's
    (shifted) WCS position, in arcsec, with a covariance-based error.
Outputs: out/diffimg_<tic>.txt, out/diffimg_<tic>.png
Option --exclude E[,E...]: epochs reported but left out of the combined offset, the chi2 sum and the
own-pixel test (e.g. a transit overlapped by a neighbouring eclipsing binary's eclipse).
"Key neighbour" = the nearest Gaia star that could mimic the dip (needed eclipse < 100%).
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import erf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TARGETS, OUT, DATA

TIC = int(sys.argv[1])
EXCL = {int(x) for x in sys.argv[sys.argv.index("--exclude") + 1].split(",")} if "--exclude" in sys.argv else set()
dur = TARGETS[TIC]["dur"]
eph = json.load(open(os.path.join(OUT, f"ephem_{TIC}.json")))
gaia = json.load(open(os.path.join(OUT, f"gaia_{TIC}.json")))
SID = int(gaia["source_id"])
NB_SID = min((n for n in gaia["neighbours"] if n["need"] < 1.0), key=lambda n: n["sep"])["source_id"]
P, T0 = eph["P"], eph["T0_E0"]
stars = [dict(source_id=SID, ra=TARGETS[TIC]["ra"], dec=TARGETS[TIC]["dec"], G=gaia["target"]["G"], sep=0.0)]
stars += [n for n in gaia["neighbours"] if np.isfinite(n["G"]) and n["G"] < 18.5]
PIX = 21.0  # arcsec / TESS pixel


def psf(nx, ny, x0, y0, s):
    """Pixel-integrated circular Gaussian, unit total flux."""
    xe = np.arange(nx + 1) - 0.5
    ye = np.arange(ny + 1) - 0.5
    fx = np.diff(0.5 * (1 + erf((xe - x0) / (np.sqrt(2) * s))))
    fy = np.diff(0.5 * (1 + erf((ye - y0) / (np.sqrt(2) * s))))
    return np.outer(fy, fx)


def load(fn):
    with fits.open(fn) as h:
        sec = int(h[0].header["SECTOR"])
        d = h[1].data
        t = np.array(d["TIME"], float)
        cube = np.array(d["FLUX"], float)
        q = np.array(d["QUALITY"])
        w = WCS(h[2].header)
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    t, cube = t[good], cube[good]
    med = np.nanmedian(cube, axis=0)
    thr = np.nanpercentile(med, 30)
    bkg = np.nanmedian(cube[:, med <= thr], axis=1)
    return sec, t, cube - bkg[:, None, None], w


def diff_image(t, cube, tc):
    dtt = t - tc
    inn = np.abs(dtt) < 0.3 * dur
    oot = (np.abs(dtt) > 0.65 * dur) & (np.abs(dtt) < 1.6 * dur)
    cad = np.median(np.diff(t))
    # full transits only: >= 80% of the expected in-transit frames, OOT on both sides
    if inn.sum() < 0.8 * 0.6 * dur / cad or (oot & (dtt < 0)).sum() < 5 or (oot & (dtt > 0)).sum() < 5:
        return None
    A = np.vstack([np.ones(oot.sum()), dtt[oot]]).T
    Y = cube[oot].reshape(oot.sum(), -1)
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    res = Y - A @ coef
    sig = res.std(axis=0)
    base_in = (np.vstack([np.ones(inn.sum()), dtt[inn]]).T @ coef)
    D = (base_in - cube[inn].reshape(inn.sum(), -1)).mean(axis=0)
    O = coef[0]
    eD = sig * np.sqrt(1 / inn.sum() + 1 / oot.sum())
    sh = cube.shape[1:]
    return D.reshape(sh), eD.reshape(sh), O.reshape(sh), inn.sum(), oot.sum()


def star_pix(w):
    xy = []
    for s in stars:
        x, y = w.all_world2pix(s["ra"], s["dec"], 0)
        xy.append((float(x), float(y)))
    return np.array(xy)


def fit_scene(O, xy, flux, box):
    ny, nx = O.shape

    def model(p):
        dx, dy, s, sc, b = p
        m = np.full(O.shape, b)
        for (x, y), f in zip(xy, flux):
            if -3 < x < nx + 3 and -3 < y < ny + 3:
                m += sc * f * psf(nx, ny, x + dx, y + dy, s)
        return m

    r = least_squares(lambda p: (O - model(p))[box] / np.sqrt(np.abs(O[box]) + 25), [0, 0, 0.9, O.sum(), 0],
                      bounds=([-1.5, -1.5, 0.4, 0, -np.inf], [1.5, 1.5, 3, np.inf, np.inf]))
    return r.x, model(r.x)


def run():
    files = sorted(glob.glob(os.path.join(DATA, "tesscut_diff", "**", f"*{TARGETS[TIC]['ra']:.6f}*.fits"), recursive=True) +
                   glob.glob(os.path.join(DATA, "tesscut", "**", f"*{TARGETS[TIC]['ra']:.6f}*.fits"), recursive=True))
    gflux = np.array([10 ** (-0.4 * (s["G"] - stars[0]["G"])) for s in stars])
    lines, panels, comb, own = [], [], [], {}
    nb_idx = [i for i, s in enumerate(stars) if s["source_id"] == NB_SID][0]
    NBL = f"G{stars[nb_idx]['G']:.1f}@{stars[nb_idx]['sep']:.1f}\""
    for fn in files:
        sec, t, cube, w = load(fn)
        if len(t) < 100:  # same rule as s03
            print(f"S{sec}: only {len(t)} good cadences - skipped")
            continue
        ny, nx = cube.shape[1:]
        xy = star_pix(w)
        # same off-silicon rule as s03: skip if the target's 3x3 peak is < 50 e/s above background
        ix, iy = int(round(xy[0, 0])), int(round(xy[0, 1]))
        med = np.nanmedian(cube, axis=0)
        peak = med[max(0, iy - 1):iy + 2, max(0, ix - 1):ix + 2].max() if (0 <= ix < nx and 0 <= iy < ny) else 0
        if peak < 50:
            print(f"S{sec}: target not on science pixels (peak {peak:.1f} e/s) - skipped")
            continue
        yy, xx = np.mgrid[:ny, :nx]
        box = np.hypot(xx - xy[0, 0], yy - xy[0, 1]) < 4.5
        n = np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1)
        for k in n:
            tc = T0 + k * P
            r = diff_image(t, cube, tc)
            if r is None:
                continue
            D, eD, O, nin, noot = r
            (dx, dy, s, sc, b), Om = fit_scene(O, xy, gflux, box)
            # hypothesis tests on D
            chis = []
            for i, (x, y) in enumerate(xy):
                if stars[i]["sep"] > 200 or not (0 <= x + dx < nx and 0 <= y + dy < ny):
                    continue
                pm = psf(nx, ny, x + dx, y + dy, s)
                A = np.vstack([pm[box], np.ones(box.sum())]).T / eD[box][:, None]
                c, *_ = np.linalg.lstsq(A, D[box] / eD[box], rcond=None)
                chis.append((np.sum((D[box] / eD[box] - A @ c) ** 2), i, c[0]))
            chis.sort()
            chi_t = [c for c in chis if c[1] == 0][0]
            chi_n = [c for c in chis if c[1] == nb_idx][0]
            # free-position fit of D

            def rr(p):
                return ((p[2] * psf(nx, ny, p[0], p[1], s) + p[3]) - D)[box] / eD[box]

            q = least_squares(rr, [xy[0, 0] + dx, xy[0, 1] + dy, chi_t[2], 0])
            J = q.jac
            cov = np.linalg.pinv(J.T @ J) * max(1.0, 2 * q.cost / (box.sum() - 4))
            ep = int(round((tc - T0) / P))
            # offsets on sky (arcsec) via WCS from the shifted target position
            ra_d, de_d = w.all_pix2world(q.x[0] - dx, q.x[1] - dy, 0)
            dra = (float(ra_d) - stars[0]["ra"]) * np.cos(np.radians(stars[0]["dec"])) * 3600
            dde = (float(de_d) - stars[0]["dec"]) * 3600
            pos_err = np.sqrt(np.diag(cov)[:2]) * PIX
            # separation from neighbour in pixel units
            off_t = np.hypot(q.x[0] - (xy[0, 0] + dx), q.x[1] - (xy[0, 1] + dy)) * PIX
            off_n = np.hypot(q.x[0] - (xy[nb_idx, 0] + dx), q.x[1] - (xy[nb_idx, 1] + dy)) * PIX
            depth_pix = D[box].sum() / O[box].sum()
            best = stars[chis[0][1]]
            lines.append(
                f"S{sec:3d} E{int(round((tc - T0) / P)):4d} tc={tc:.3f} n_in={nin} n_oot={noot} psf_sigma={s:.2f}px "
                f"shift=({dx:+.2f},{dy:+.2f})px | D source offset from target: dRA={dra:+.1f}\" dDec={dde:+.1f}\" "
                f"(|off|={off_t:.1f}\" ±{np.hypot(*pos_err):.1f}\"), from {NBL} neighbour: {off_n:.1f}\" | "
                f"chi2 target={chi_t[0]:.1f} neighbour {NBL}={chi_n[0]:.1f} (dof {box.sum()-2}); best-fitting star: "
                f"{'TARGET' if chis[0][1] == 0 else best['source_id']} (sep {best['sep']:.0f}\")" + ("  [EXCLUDED from combined]" if ep in EXCL else ""))
            if ep in EXCL:
                panels.append((sec, tc, O, D, eD, xy + [dx, dy], q.x[:2]))
                continue
            comb.append((sec, int(round((tc - T0) / P)), chi_t[0], chi_n[0], off_t, np.hypot(*pos_err),
                         dra, dde, np.max(pos_err), np.max(pos_err)))
            panels.append((sec, tc, O, D, eD, xy + [dx, dy], q.x[:2]))
            # own-pixel test for resolved neighbours bright enough to mimic the dip: if star k were the
            # source, its 3x3 pixels would lose (total lost flux) x (its PSF fraction in that box).
            # ratio = measured / expected; ~1 if k is the source, ~0 if not.
            lost = D[box].sum()
            for i, st in enumerate(stars):
                if i == 0 or st["sep"] < 30 or st.get("need", 1.0) >= 0.5:
                    continue
                x, y = xy[i, 0] + dx, xy[i, 1] + dy
                ix, iy = int(round(x)), int(round(y))
                if not (1 <= ix <= nx - 2 and 1 <= iy <= ny - 2):
                    continue
                sl = (slice(iy - 1, iy + 2), slice(ix - 1, ix + 2))
                frac = psf(nx, ny, x, y, s)[sl].sum()
                exp_ = lost * frac
                own.setdefault(i, []).append((D[sl].sum() / exp_, np.sqrt((eD[sl] ** 2).sum()) / exp_))
    # combined offset (weighted mean in pixel x/y)
    a = np.array([c[6:] for c in comb])
    wx, wy = 1 / a[:, 2] ** 2, 1 / a[:, 3] ** 2
    mx, my = np.sum(wx * a[:, 0]) / wx.sum(), np.sum(wy * a[:, 1]) / wy.sum()
    ex, ey = 1 / np.sqrt(wx.sum()), 1 / np.sqrt(wy.sum())
    nbsep = stars[nb_idx]["sep"]
    dchi = sum(c[3] - c[2] for c in comb)
    lines.append("")
    nb = stars[nb_idx]
    nra = (nb["ra"] - stars[0]["ra"]) * np.cos(np.radians(stars[0]["dec"])) * 3600
    nde = (nb["dec"] - stars[0]["dec"]) * 3600
    lines.append(f"Combined over {len(comb)} full transits{' (excluded epochs: ' + ','.join(map(str, sorted(EXCL))) + ')' if EXCL else ''}: transit-source offset from target "
                 f"dRA = {mx:+.1f} ± {ex:.1f}\", dDec = {my:+.1f} ± {ey:.1f}\" "
                 f"(|offset| = {np.hypot(mx, my):.1f}\"). The key neighbour {NBL} is at dRA = {nra:+.1f}\", "
                 f"dDec = {nde:+.1f}\" ({nbsep:.1f}\"), i.e. {np.hypot(mx-nra, my-nde)/np.hypot(ex, ey):.0f} sigma "
                 f"from the measured source position.")
    lines.append("Per-transit errors are least-squares position errors (residual-scaled); systematic PSF/WCS "
                 "errors of a few arcsec are not included.")
    lines.append(f"Sum over transits chi2(neighbour) - chi2(target) = {dchi:.1f}  (positive favours the target)")
    # scatter-based check of the combined offset
    lines.append(f"Unweighted mean offset dRA = {a[:, 0].mean():+.1f} ± {a[:, 0].std(ddof=1)/np.sqrt(len(a)):.1f}\", "
                 f"dDec = {a[:, 1].mean():+.1f} ± {a[:, 1].std(ddof=1)/np.sqrt(len(a)):.1f}\" (error from the scatter "
                 f"between transits)")
    lines.append("")
    lines.append("Own-pixel test for resolved Gaia stars (sep 30-200\") that could mimic the dip (needed eclipse < 50%):")
    lines.append("  ratio = flux lost in the star's own 3x3 pixels / flux it would lose if it were the source (1 = source)")
    for i, v in sorted(own.items(), key=lambda kv: stars[kv[0]]["sep"]):
        r = np.array(v)
        m, e, _ = (np.sum(r[:, 0] / r[:, 1] ** 2) / np.sum(1 / r[:, 1] ** 2), 1 / np.sqrt(np.sum(1 / r[:, 1] ** 2)), 0)
        st = stars[i]
        lines.append(f"  Gaia {st['source_id']} G={st['G']:.2f} sep={st['sep']:.0f}\" needed eclipse {st['need']*100:.1f}%: "
                     f"ratio = {m:+.3f} ± {e:.3f} over {len(r)} transits -> {(1-m)/e:.0f} sigma from being the source")
    txt = "\n".join(lines)
    open(os.path.join(OUT, f"diffimg_{TIC}.txt"), "w").write(txt + "\n")
    print(txt)
    # plot
    nr = len(panels)
    fig, axs = plt.subplots(nr, 3, figsize=(11, 3.4 * nr), squeeze=False)
    for row, (sec, tc, O, D, eD, xy, src) in zip(axs, panels):
        for ax, img, ttl in zip(row, [O, D, D / eD], ["out-of-transit", "difference (OOT - in)", "difference S/N"]):
            im = ax.imshow(img, origin="lower", cmap="viridis")
            plt.colorbar(im, ax=ax, fraction=0.046)
            NBL = f"G{stars[nb_idx]['G']:.1f} @{stars[nb_idx]['sep']:.1f}\""
            ax.plot(xy[0, 0], xy[0, 1], "r+", ms=14, mew=2, label="target")
            ax.plot(xy[nb_idx, 0], xy[nb_idx, 1], "wx", ms=10, mew=2, label=NBL)
            others = [i for i in range(1, len(stars)) if i != nb_idx and stars[i]["G"] < 16.5]
            ax.plot(xy[others, 0], xy[others, 1], "o", mfc="none", mec="orange", ms=6)
            if "difference" in ttl:
                ax.plot(src[0], src[1], "c*", ms=10)
            ax.set_xlim(-0.5, img.shape[1] - 0.5); ax.set_ylim(-0.5, img.shape[0] - 0.5)
            ax.set_title(f"S{sec} tc={tc:.2f}: {ttl}", fontsize=8)
        row[0].legend(fontsize=6, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"diffimg_{TIC}.png"), dpi=70); plt.close(fig)


if __name__ == "__main__":
    run()
