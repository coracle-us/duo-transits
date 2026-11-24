"""
AAS82225 revision: what causes the 0.2405 d SPOC TCE on EPIC 211489484 (TIC 6951741)?
1. SPOC 2-min PDCSAP light curves of the run's sectors (S7, S44, S45, S46; downloaded 2026-11-23 into
   data/literature/lc_tic6951741/), our eclipses masked: Lomb-Scargle periodogram, and a fold at P and 2P per sector.
2. Neighbours within 60" in VizieR (CfA mirror): Gaia DR3, VSX, ZTF, ATLAS and Gaia DR3 variability, to find a
   variable star at the DV difference-image offset (mean 24.5 +- 3.5" from the target, to the south).
Output: results/literature/tce024_check.txt
Usage:  python -I k2singles/revision2/tce024_check.py
"""
import glob
import io
import os
import time
import urllib.parse
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.io.votable import parse
from astropy.timeseries import LombScargle

HERE = os.path.dirname(os.path.abspath(__file__))
LC = os.path.join(HERE, "data", "literature", "lc_tic6951741")
P_TCE, T_TCE = 0.24051265, 1491.6647469
P_EB, T_EB, PHI2 = 82.25843, 1508.757, 0.194  # our eclipse ephemeris (S7 primary), secondary phase
RA, DEC = 126.022711, 12.987866


def eclipse_mask(t, half=0.25):
    m = np.zeros(t.size, bool)
    for ph in (0.0, PHI2):
        d = ((t - T_EB - ph * P_EB + P_EB / 2) % P_EB) - P_EB / 2
        m |= np.abs(d) < half
    return m


def fold_depth(t, f, P, T0, half_h=0.25):
    """Mean flux within +-half_h of phase 0, minus the median elsewhere, in ppm, with its error."""
    d = ((t - T0 + P / 2) % P) - P / 2
    w = np.abs(d) < half_h / 24
    o = ~w
    sd = np.std(f[o]) / np.sqrt(w.sum())
    return 1e6 * (np.mean(f[w]) - np.median(f[o])), 1e6 * sd


print("## 1. SPOC 2-min PDCSAP, eclipses masked\n")
T, F = [], []
for f in sorted(glob.glob(os.path.join(LC, "*_lc.fits"))):
    d = fits.getdata(f, 1)
    sec = fits.getheader(f)["SECTOR"]
    g = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
    t, y = d["TIME"][g].astype(float), d["PDCSAP_FLUX"][g].astype(float)
    y = y / np.nanmedian(y)
    m = ~eclipse_mask(t)
    t, y = t[m], y[m]
    T.append(t); F.append(y)
    fr, pw = LombScargle(t, y).autopower(minimum_frequency=1 / 10, maximum_frequency=1 / 0.05, samples_per_peak=20)
    top = np.argsort(pw)[::-1]
    peaks = []
    for i in top:
        if all(abs(fr[i] - fr[j]) > 0.05 for j in peaks):
            peaks.append(i)
        if len(peaks) == 3:
            break
    a1, e1 = fold_depth(t, y, P_TCE, T_TCE)
    a2, e2 = fold_depth(t, y, 2 * P_TCE, T_TCE)
    a3, e3 = fold_depth(t, y, 2 * P_TCE, T_TCE + P_TCE)
    print(f"S{sec}: {g.sum()} good cadences, {m.sum()} after masking eclipses, rms {1e6 * np.std(y):.0f} ppm")
    print("   LS peaks (period d, power): " + ", ".join(f"{1 / fr[i]:.5f} ({pw[i]:.3f})" for i in peaks))
    print(f"   fold at P = {P_TCE} d: {a1:+.0f} +- {e1:.0f} ppm;  at 2P, phase 0: {a2:+.0f} +- {e2:.0f}, "
          f"phase 0.5: {a3:+.0f} +- {e3:.0f}")
t, y = np.concatenate(T), np.concatenate(F)
a1, e1 = fold_depth(t, y, P_TCE, T_TCE)
print(f"all four sectors: fold at P: {a1:+.0f} +- {e1:.0f} ppm")
fr, pw = LombScargle(t[t > 2400], y[t > 2400]).autopower(minimum_frequency=1 / 2, maximum_frequency=1 / 0.1,
                                                         samples_per_peak=50)
i = np.argmax(pw)
print(f"S44-46 combined LS peak: {1 / fr[i]:.6f} d (power {pw[i]:.3f}); 2x = {2 / fr[i]:.6f} d")

print("\n## 2. Neighbours within 60\" (VizieR, CfA mirror)\n")
cats = "I/355/gaiadr3 B/vsx J/ApJS/249/18 J/AJ/156/241 I/358/vclassre I/358/veb"
url = "https://vizier.cfa.harvard.edu/viz-bin/votable?" + urllib.parse.urlencode(
    {"-source": cats, "-c": f"{RA} {DEC:+}", "-c.rs": 60, "-out.max": 50, "-out.all": 1, "-out.add": "_r,_RA,_DE"})
cache = os.path.join(HERE, "data", "literature", "tce024_neighbours.vot")
raw = open(cache, "rb").read() if os.path.exists(cache) else None
for k in range(20):  # the mirror's database goes down, and some replies are cut: retry until complete
    if raw is not None:
        break
    try:
        x = urllib.request.urlopen(url, timeout=600).read()
    except Exception:  # noqa: BLE001
        x = b""
    tx = x.decode("utf-8", "replace")
    if "not currently reachable" not in tx and tx.rstrip().endswith("</VOTABLE>") and "I/355/gaiadr3" in tx:
        raw = x
        open(cache, "wb").write(raw)
    else:
        time.sleep(min(300, 30 * (k + 1)))
if raw is None:
    print("VizieR did not return a complete reply in 20 tries; not checked")
else:
    v = parse(io.BytesIO(raw), verify="ignore")
    for res in v.resources:
        for tb in res.tables:
            tab = tb.to_table(use_names_over_ids=True)
            print(f"{tb.name}: {len(tab)} rows")
            cra = next(c for c in tab.colnames if c.startswith("_RA"))
            cde = next(c for c in tab.colnames if c.startswith("_DE"))
            for r in tab:
                dra = (float(r[cra]) - RA) * 3600 * np.cos(np.radians(DEC))
                dde = (float(r[cde]) - DEC) * 3600
                keep = [c for c in tab.colnames if c in ("Source", "Gmag", "BPmag", "RPmag", "Name", "Type", "Period",
                                                         "Per", "per", "Amp", "max", "min", "ID", "Class", "Freq",
                                                         "VarFlag", "PS", "RUWE", "Plx")]
                print(f"   r = {float(r['_r']):5.1f}\"  dRA {dra:+6.1f}\"  dDec {dde:+6.1f}\"  "
                      + "  ".join(f"{c}={r[c]}" for c in keep))

print("\n## 3. Variability catalogues at the neighbour (Gaia DR3 651067327828913280), 3\"\n")
NB = (126.0233740, 12.9806440)  # filled from the Gaia row printed above, approximately; refined below
try:
    v = parse(io.BytesIO(raw), verify="ignore")
    g = next(tb for res in v.resources for tb in res.tables if tb.name == "I/355/gaiadr3").to_table(use_names_over_ids=True)
    row = g[g["Source"] == 651067327828913280][0]
    NB = (float(row["RA_ICRS"]), float(row["DE_ICRS"]))
except Exception as e:  # noqa: BLE001
    print("could not read the neighbour's position from the Gaia row:", e)
url = "https://vizier.cfa.harvard.edu/viz-bin/votable?" + urllib.parse.urlencode(
    {"-source": "I/355/gaiadr3 I/358 B/vsx J/ApJS/249/18 J/AJ/156/241 II/366", "-c": f"{NB[0]} {NB[1]:+}", "-c.rs": 3,
     "-out.max": 5, "-out.all": 1, "-out.add": "_r"})
cache2 = os.path.join(HERE, "data", "literature", "tce024_neighbour_vari.vot")
raw2 = open(cache2, "rb").read() if os.path.exists(cache2) else None
for k in range(20):
    if raw2 is not None:
        break
    try:
        x = urllib.request.urlopen(url, timeout=600).read()
    except Exception:  # noqa: BLE001
        x = b""
    tx = x.decode("utf-8", "replace")
    if "not currently reachable" not in tx and tx.rstrip().endswith("</VOTABLE>") and "I/355/gaiadr3" in tx:
        raw2 = x
        open(cache2, "wb").write(raw2)
    else:
        time.sleep(min(300, 30 * (k + 1)))
if raw2 is None:
    print("VizieR did not return a complete reply in 20 tries; not checked")
else:
    v = parse(io.BytesIO(raw2), verify="ignore")
    for res in v.resources:
        for tb in res.tables:
            tab = tb.to_table(use_names_over_ids=True)
            if not len(tab):
                continue
            print(f"{tb.name}: " + "; ".join(f"{c}={tab[c][0]}" for c in tab.colnames
                                             if str(tab[c][0]) not in ("--", "", "nan")))

print("\n## 4. Shape of the 0.2405 d signal in the target's light curve (all four sectors, eclipses masked)\n")
ph = ((t - T_TCE) / P_TCE + 0.5) % 1 - 0.5
edges = np.linspace(-0.5, 0.5, 25)
print("phase bin centre: mean flux - 1 (ppm) +- error")
for a, b in zip(edges[:-1], edges[1:]):
    w = (ph >= a) & (ph < b)
    print(f"   {0.5 * (a + b):+.3f}: {1e6 * (np.mean(y[w]) - 1):+7.0f} +- {1e6 * np.std(y[w]) / np.sqrt(w.sum()):.0f}")
ls = LombScargle(t, y)
f0 = 1 / P_TCE
fr = np.linspace(1 / 2, 1 / 0.05, 200000)
pw = ls.power(fr)
p0 = ls.power(np.array([f0, 2 * f0, 3 * f0, f0 / 2]))
print(f"LS power at f, 2f, 3f, f/2 (f = 1/{P_TCE} d): " + ", ".join(f"{x:.2e}" for x in p0)
      + f";  99.9th percentile over 0.05-2 d: {np.percentile(pw, 99.9):.2e}; max {pw.max():.2e} at {1 / fr[np.argmax(pw)]:.5f} d")

print("\n## 5. Periodogram to 50 c/d (the section 1 periodogram stopped at 20 c/d)\n")
fr = np.linspace(2, 50, 400000)
for lab, sel in [("S7", t < 1600), ("S44", (t > 2490) & (t < 2525)), ("S45", (t > 2525) & (t < 2551)),
                 ("S46", (t > 2551) & (t < 2600)), ("S44-46", t > 2490), ("all four", t > 0)]:
    ls = LombScargle(t[sel], y[sel])
    pw = ls.power(fr)
    i = np.argmax(pw)
    amp = np.sqrt(np.sum(ls.model_parameters(fr[i])[1:] ** 2))
    print(f"{lab:8s}: peak {fr[i]:.5f} c/d (P = {1 / fr[i]:.6f} d), power {pw[i]:.2e} "
          f"(median {np.median(pw):.1e}), semi-amplitude {1e6 * amp:.0f} ppm;  5 x f_TCE = {5 / P_TCE:.5f} c/d")
