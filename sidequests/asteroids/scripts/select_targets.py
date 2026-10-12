"""
Which members of a population cross a TESS sector, for how long, and how bright?

Pre-filter (fast, approximate): two-body (Keplerian) propagation of JPL SBDB osculating elements
to a 0.5-d time grid across the sector (sector dates from cache/TESS_orbit_times.csv, fetched
from https://tess.mit.edu/public/files/TESS_orbit_times.csv), geocentric RA/Dec (Earth from
astropy's built-in ephemeris; TESS-Earth parallax, <~2' at 4 au, ignored), then tess-point
(tess_stars2px, fixed sector) to test whether the position is on silicon.
Days on silicon and the predicted V (H-G, G=0.15) are reported; the exact track, r, Delta,
alpha and V come later from JPL Horizons (fetch_ephem in mtphot.py), and TESScut decides the
actual cutouts.

Usage:
  python select_targets.py count  POP SEC_LO SEC_HI       # members on silicon per sector
  python select_targets.py list   POP SECTOR [HMAX]       # per-object table for one sector
POP: trojan | hilda
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

HERE = Path(__file__).resolve().parents[1]
CACHE = HERE / "cache"
DATA = HERE / "data"
K = 0.01720209895  # Gaussian gravitational constant (au^1.5 / day)
EPS = np.deg2rad(23.43928)  # J2000 obliquity


def elements(pop):
    cls = {"trojan": "TJN", "hilda": "OMB"}[pop]
    f = CACHE / f"sbdb_elem_{cls}.json"
    if not f.exists():
        url = ("https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=pdes,name,H,epoch,e,a,i,om,w,ma"
               f"&sb-class={cls}&sb-kind=a&sb-ns=n")
        with urllib.request.urlopen(url, timeout=600) as r:
            f.write_bytes(r.read())
    j = json.loads(f.read_text())
    d = pd.DataFrame(j["data"], columns=j["fields"])
    for c in ["H", "epoch", "e", "a", "i", "om", "w", "ma"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if pop == "hilda":
        d = d[(d.a >= 3.70) & (d.a <= 4.20) & (d.e < 0.30) & (d.i < 20)]
    return d.dropna(subset=["H", "a"]).reset_index(drop=True)


def helio_ecl(el, jd):
    """Heliocentric ecliptic J2000 positions, shape (n_obj, n_t, 3)."""
    a, e = el.a.values[:, None], el.e.values[:, None]
    i, om, w = (np.deg2rad(el[c].values)[:, None] for c in ["i", "om", "w"])
    n = K / a ** 1.5
    M = np.deg2rad(el.ma.values)[:, None] + n * (jd[None, :] - el.epoch.values[:, None])
    E = M.copy()
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1 - e ** 2) * np.sin(E)
    cw, sw, co, so, ci, si = np.cos(w), np.sin(w), np.cos(om), np.sin(om), np.cos(i), np.sin(i)
    x = (co * cw - so * sw * ci) * xp + (-co * sw - so * cw * ci) * yp
    y = (so * cw + co * sw * ci) * xp + (-so * sw + co * cw * ci) * yp
    z = (sw * si) * xp + (cw * si) * yp
    return np.stack([x, y, z], axis=-1)


def earth_helio_ecl(jd):
    t = Time(jd, format="jd", scale="tdb")
    e = get_body_barycentric("earth", t).xyz.to("au").value
    s = get_body_barycentric("sun", t).xyz.to("au").value
    v = (e - s).T  # equatorial ICRS
    x, y, z = v[:, 0], np.cos(EPS) * v[:, 1] + np.sin(EPS) * v[:, 2], -np.sin(EPS) * v[:, 1] + np.cos(EPS) * v[:, 2]
    return np.stack([x, y, z], axis=-1)


def radec(el, jd):
    P = helio_ecl(el, jd)
    E = earth_helio_ecl(jd)[None, :, :]
    g = P - E
    # ecliptic -> equatorial
    x, y, z = g[..., 0], np.cos(EPS) * g[..., 1] - np.sin(EPS) * g[..., 2], np.sin(EPS) * g[..., 1] + np.cos(EPS) * g[..., 2]
    delta = np.sqrt(x * x + y * y + z * z)
    ra = np.rad2deg(np.arctan2(y, x)) % 360
    dec = np.rad2deg(np.arcsin(z / delta))
    r = np.linalg.norm(P, axis=-1)
    re = np.linalg.norm(E, axis=-1)
    cosa = (r ** 2 + delta ** 2 - re ** 2) / (2 * r * delta)
    alpha = np.rad2deg(np.arccos(np.clip(cosa, -1, 1)))
    return ra, dec, r, delta, alpha


def hg_mag(H, r, delta, alpha, G=0.15):
    a = np.deg2rad(alpha)
    p1 = np.exp(-3.33 * np.tan(a / 2) ** 0.63)
    p2 = np.exp(-1.87 * np.tan(a / 2) ** 1.22)
    return H + 5 * np.log10(r * delta) - 2.5 * np.log10((1 - G) * p1 + G * p2)


def sector_span(sector):
    o = pd.read_csv(CACHE / "TESS_orbit_times.csv", comment="#")
    o = o[o.Sector == sector]
    t0 = Time(pd.to_datetime(o["Start of Orbit"]).min().to_pydatetime(), scale="utc").jd
    t1 = Time(pd.to_datetime(o["End of Orbit"]).max().to_pydatetime(), scale="utc").jd
    return t0, t1


def on_silicon(ra, dec, sector):
    """tess-point test, after a cheap pre-cut: only positions within 17 deg of one of the four
    camera boresights (camera FOV is 24x24 deg, half-diagonal 17 deg) are passed to tess-point."""
    from tess_stars2px import TESS_Spacecraft_Pointing_Data, tess_stars2px_function_entry as t2p
    sc = TESS_Spacecraft_Pointing_Data(trySector=sector)
    cra, cdec = np.deg2rad(np.ravel(sc.camRa)), np.deg2rad(np.ravel(sc.camDec))
    r, d = np.deg2rad(ra.ravel()), np.deg2rad(dec.ravel())
    cosd = (np.sin(d)[:, None] * np.sin(cdec)[None] + np.cos(d)[:, None] * np.cos(cdec)[None]
            * np.cos(r[:, None] - cra[None]))
    near = np.where((cosd > np.cos(np.deg2rad(17.5))).any(axis=1))[0]
    hit = np.zeros(ra.size, bool)
    cam = np.zeros(ra.size, int)
    if len(near) == 0:
        return hit.reshape(ra.shape), cam.reshape(ra.shape)
    res = t2p(near, ra.ravel()[near], dec.ravel()[near], trySector=sector, scInfo=sc)
    ok = res[3] == sector
    hit[res[0][ok]] = True
    cam[res[0][ok]] = res[4][ok]
    return hit.reshape(ra.shape), cam.reshape(ra.shape)


def sector_table(pop, sector, hmax=99.0, step=0.5):
    el = elements(pop)
    el = el[el.H <= hmax].reset_index(drop=True)
    t0, t1 = sector_span(sector)
    jd = np.arange(t0, t1, step)
    ra, dec, r, delta, alpha = radec(el, jd)
    hit, cam = on_silicon(ra, dec, sector)
    days = hit.sum(1) * step
    V = hg_mag(el.H.values[:, None], r, delta, alpha)
    rows = []
    for k in np.where(days > 0)[0]:
        h = hit[k]
        rows.append(dict(pdes=el.pdes[k], name=el.name[k], H=el.H[k], a=el.a[k], days_on_silicon=days[k],
                         camera=int(np.bincount(cam[k][h]).argmax()),
                         V_med=round(float(np.median(V[k][h])), 2), r=round(float(np.median(r[k][h])), 3),
                         delta=round(float(np.median(delta[k][h])), 3),
                         alpha=round(float(np.median(alpha[k][h])), 2),
                         jd_first=jd[h][0], jd_last=jd[h][-1]))
    return pd.DataFrame(rows)


def main():
    mode, pop = sys.argv[1], sys.argv[2]
    if mode == "count":
        lo, hi = int(sys.argv[3]), int(sys.argv[4])
        out = []
        for s in range(lo, hi + 1):
            try:
                tab = sector_table(pop, s, hmax=13.5, step=1.0)
            except Exception as e:  # noqa: BLE001
                print(s, "ERR", e)
                continue
            n10 = int((tab.days_on_silicon >= 10).sum()) if len(tab) else 0
            nb = int(((tab.days_on_silicon >= 10) & (tab.V_med < 18.0)).sum()) if len(tab) else 0
            out.append(dict(sector=s, n_any=len(tab), n_ge10d=n10, n_ge10d_V_lt18=nb))
            print(out[-1], flush=True)
        DATA.mkdir(exist_ok=True)
        pd.DataFrame(out).to_csv(DATA / f"sector_counts_{pop}_H13.5.csv", index=False)
    else:
        s = int(sys.argv[3])
        hmax = float(sys.argv[4]) if len(sys.argv) > 4 else 99.0
        tab = sector_table(pop, s, hmax=hmax)
        DATA.mkdir(exist_ok=True)
        tab.to_csv(DATA / f"crossing_{pop}_s{s}.csv", index=False)
        print(tab.sort_values("V_med").to_string(index=False))


if __name__ == "__main__":
    main()
