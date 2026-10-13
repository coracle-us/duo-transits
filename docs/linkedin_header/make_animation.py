"""
LinkedIn header animation for the AAS82207 work: TIC 159614322 (CTOI 159614322.01).
Left: the host star with the planet crossing it, planet drawn to scale (Rp/R* = 0.055, b = 0, from
singles/writeups/tic159614322/analysis.json). Right: the real TESS light curves of the three transits
(S14 TESS-SPOC, S21 TESS-SPOC SAP partial, S41 SPOC 2-min binned to 30 min), revealed as the planet moves,
over the fitted transit model. Data loading follows singles/writeups/tic159614322/make_figures.py.
Outputs (next to this script): header.mp4 (1920x1080), header.gif (960x540), header_still.png, header_mid.png
Run:  python docs/linkedin_header/make_animation.py
"""
import json
import sys
from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WU = ROOT / "singles" / "writeups" / "tic159614322"
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vc  # noqa: E402

TIC = 159614322
EPH = json.load(open(WU / "ephemeris.json"))
SH = json.load(open(WU / "analysis.json"))["best_b_P90"]
TR = {r["sector"]: r for r in EPH["transits"]}
P8 = EPH["ephemerides"]["n8"]["P"]
GAP21 = 1883.975                     # start of the S21 orbit-1 downlink gap (BTJD)
D = 8.773 / 24

BG, FG, MUTED = "#0a1220", "#e8eef5", "#8a9bb0"
COL = {14: "#ffb347", 41: "#4fc3f7", 21: "#f06292"}
LAB = {14: "Jul 2019 (Sector 14)", 41: "Jul 2021 (Sector 41)", 21: "Jan 2020 (Sector 21): partial, cut by a data gap"}
FPS, SWEEP_S, HOLD_S = 30, 9.0, 4.0
XLIM = (-14, 14)


def load():
    prods = vc.all_products(TIC)
    get = lambda s, pr: next(q for q in prods if q["sector"] == s and q["product"] == pr)  # noqa: E731
    out = {}
    for s, pr in ((14, "TESS-SPOC"), (41, "SPOC_120")):
        Tm = TR[s]["T"]
        t, f = vc.window(get(s, pr), Tm, D)[:2]
        x, y = (t - Tm) * 24, (f - 1) * 1e3
        if s == 41:                   # 2-min cadence: 30-min bins
            bins = np.arange(-14, 14.01, 0.5)
            idx = np.digitize(x, bins)
            u = [i for i in np.unique(idx) if (idx == i).sum() > 5]
            x, y = np.array([x[idx == i].mean() for i in u]), np.array([y[idx == i].mean() for i in u])
        out[s] = (x, y)
    Tm = TR[21]["T"]
    p = get(21, "TESS-SPOC")
    m = (p["time"] > Tm - 14 / 24) & (p["time"] < GAP21) & np.isfinite(p["sap_flux"])
    base = m & (p["time"] < 1883.55)
    out[21] = ((p["time"][m] - Tm) * 24, (p["sap_flux"][m] / np.median(p["sap_flux"][base]) - 1) * 1e3)
    tt = np.linspace(*XLIM, 800)
    model = (vc.model(TR[41]["T"] + tt / 24, TR[41]["T"], P8, SH["rp"], SH["b"], SH["aR"], 120) - 1) * 1e3
    return out, tt, model


def star_image(n=700, ext=1.75):
    g = np.linspace(-ext, ext, n)
    X, Y = np.meshgrid(g, g)
    r = np.hypot(X, Y)
    mu = np.sqrt(np.clip(1 - r ** 2, 0, 1))
    u1, u2 = 0.45, 0.2
    disk = np.where(r <= 1, 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2, 0)
    taper = np.clip((ext - r) / 0.45, 0, 1) ** 2          # fade the glow to zero before the image edge
    glow = np.where(r > 1, (0.35 * np.exp(-(r - 1) / 0.12) + 0.08 * np.exp(-(r - 1) / 0.5)) * taper, 0)
    return disk + glow, ext


def main():
    data, tt, model = load()
    img, ext = star_image()
    cmap = LinearSegmentedColormap.from_list("star", [BG, "#5a1e0a", "#c4561b", "#ffb15c", "#ffe9c4", "#ffffff"])

    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=BG)
    ax_s = fig.add_axes([0.03, 0.12, 0.38, 0.68])
    ax_s.imshow(img, extent=[-ext, ext, -ext, ext], cmap=cmap, vmin=0, vmax=1.02, origin="lower", interpolation="bilinear")
    ax_s.set_xlim(-ext, ext)
    ax_s.set_ylim(-ext, ext)
    ax_s.axis("off")
    planet = Circle((-3, SH["b"]), SH["rp"], color="#05080d", zorder=5)
    ax_s.add_patch(planet)
    ax_s.text(0, -1.62, "star and planet to scale", color=MUTED, fontsize=15, ha="center")

    ax = fig.add_axes([0.47, 0.17, 0.5, 0.6], facecolor=BG)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=15)
    ax.set_xlim(*XLIM)
    ax.set_ylim(-5.2, 3.8)
    ax.set_xlabel("hours from mid-transit", color=MUTED, fontsize=17)
    ax.set_ylabel("starlight (parts per thousand)", color=MUTED, fontsize=17)
    ax.axvspan((GAP21 - TR[21]["T"]) * 24, XLIM[1], ymin=0, ymax=1, color=COL[21], alpha=0.06, lw=0)
    gap_txt = ax.text((GAP21 - TR[21]["T"]) * 24 + 0.3, -4.9, "Jan 2020 data gap", color=COL[21], alpha=0.0, fontsize=14)
    ax.plot(tt, model, color="#ffffff", lw=2.2, alpha=0.18, zorder=1)
    mline, = ax.plot([], [], color="#ffffff", lw=2.2, alpha=0.75, zorder=4)
    sc = {s: ax.scatter([], [], s=34 if s != 21 else 22, color=COL[s], alpha=0.95, lw=0, zorder=3) for s in (14, 21, 41)}
    cursor = ax.axvline(XLIM[0], color="#ffffff", lw=1, alpha=0.25)
    for k, s in enumerate((14, 21, 41)):
        ax.text(0.01, 0.97 - 0.065 * k, "● " + LAB[s], transform=ax.transAxes, color=COL[s], fontsize=15, va="top")

    fig.text(0.03, 0.93, "Three transits of a planet candidate, seen by NASA's TESS", color=FG, fontsize=36, weight="bold")
    fig.text(0.03, 0.875, "TIC 159614322: a Jupiter-sized object dims its star by 0.3% for about 10 hours",
             color=MUTED, fontsize=21)
    end1 = fig.text(0.47, 0.085, "The partial 2020 transit narrows the orbit to 90 or 180 days", color=FG, fontsize=22,
                    weight="bold", alpha=0)
    end2 = fig.text(0.47, 0.045, "Next predicted transit: 28 Dec 2026, 06:40 UTC (±45 min)", color="#ffd27a",
                    fontsize=20, alpha=0)
    fig.text(0.03, 0.03, "Data: NASA TESS, Sectors 14, 21, 41  ·  Analysis: Jon McBee, Research Note submitted to the AAS",
             color=MUTED, fontsize=13)

    vx = 2 * np.pi * SH["aR"] / (P8 * 24)       # planet speed, stellar radii per hour
    n_sweep, n_hold = int(FPS * SWEEP_S), int(FPS * HOLD_S)
    writer = imageio.get_writer(HERE / "header.mp4", fps=FPS, codec="libx264", quality=9, pixelformat="yuv420p",
                                macro_block_size=8)
    gif_frames = []
    for i in range(n_sweep + n_hold):
        h = XLIM[0] + (XLIM[1] - XLIM[0]) * min(i / (n_sweep - 1), 1.0)
        planet.center = (vx * h, SH["b"])
        cursor.set_xdata([h, h])
        for s, (x, y) in data.items():
            m = x <= h
            sc[s].set_offsets(np.column_stack([x[m], y[m]]) if m.any() else np.empty((0, 2)))
        mm = tt <= h
        mline.set_data(tt[mm], model[mm])
        a = np.clip((i - n_sweep) / (FPS * 0.8), 0, 1)
        end1.set_alpha(a)
        end2.set_alpha(np.clip((i - n_sweep - FPS * 0.6) / (FPS * 0.8), 0, 1))
        gap_txt.set_alpha(0.8 if h > (GAP21 - TR[21]["T"]) * 24 else 0)
        cursor.set_alpha(0.25 if i < n_sweep else 0)
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        writer.append_data(frame)
        if i % 2 == 0:
            gif_frames.append(frame[::2, ::2])
        if i == int(n_sweep * 0.5):
            imageio.imwrite(HERE / "header_mid.png", frame)
        if i == n_sweep + n_hold - 1:
            imageio.imwrite(HERE / "header_still.png", frame)
    writer.close()
    imageio.mimsave(HERE / "header.gif", gif_frames, duration=2 / FPS, loop=0)
    print("wrote header.mp4, header.gif, header_still.png")


if __name__ == "__main__":
    main()
