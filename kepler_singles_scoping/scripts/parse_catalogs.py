"""
Parse Kepler single/double-transit catalogues from their arXiv LaTeX (and VizieR) into one event list.

Sources (raw files in kepler_singles_scoping/data/raw/):
  KM19  Kawahara & Masuda 2019, AJ 157, 218 (arXiv:1904.04980; VizieR J/AJ/157/218)
        tables 'clean' (table.tex = VizieR table3) and 'flagged' (table_rm.tex)
  W15   Wang et al. 2015, ApJ 815, 127 (Planet Hunters VIII; arXiv:1512.02559)
        orbital_params_1.tex (1 transit), orbital_params_2.tex (2 transits)
  FM16  Foreman-Mackey et al. 2016, AJ 152, 206 (arXiv:1607.08237), Table 'catalog'
  HZW19 Herman, Zhu & Wu 2019, AJ 157, 248 (arXiv:1901.01974), param_tab.tex
  U16   Uehara et al. 2016, ApJ 822, 2 (arXiv:1602.07848), Table 'list' (28 STEs)
  HC24  KOI 1271.02 single (arXiv:2403.03427), entered by hand from the text

Times: BKJD = BJD - 2454833; BTJD = BJD - 2457000 = BKJD - 2167.
Depth: (Rp/R*)^2 (no limb darkening; scoping only). Duration: catalogue value if given (FM16 'duration',
U16 fit 'T'), otherwise circular-orbit T14 from P, a/R* (or M*, R*) and b.

Output: data/events_by_source.csv (one row per source x event).
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
RAW = HERE / "data" / "raw"
RSUN_REARTH = 109.1
RJ_RSUN = 0.10045


def first_num(s):
    """Leading central value of a LaTeX cell: '$1.7^{+1.8}_{-0.5}$', '$(1.7^{+1.8}_{-0.5})\\times10^3$', '1203.8'."""
    if s is None:
        return np.nan
    s = s.strip()
    if "nodata" in s or s == "":
        return np.nan
    m = re.search(r"-?\d+\.?\d*(?:\s*\\times\s*10\^\{?-?\d+\}?)?", s)
    if not m:
        return np.nan
    v = float(re.match(r"-?\d+\.?\d*", m.group(0)).group(0))
    mul = re.search(r"\\times\s*10\^\{?(-?\d+)\}?", s)
    if mul:
        v *= 10 ** int(mul.group(1))
    return v


def errs(s):
    """(+err, -err) from ^{+a}_{-b}, scaled by \\times10^n if present."""
    up = re.search(r"\^\{\+?([\d.]+)\}", s or "")
    lo = re.search(r"_\{-?([\d.]+)\}", s or "")
    mul = re.search(r"\\times\s*10\^\{?(-?\d+)\}?", s or "")
    f = 10 ** int(mul.group(1)) if mul else 1
    return (float(up.group(1)) * f if up else np.nan, float(lo.group(1)) * f if lo else np.nan)


def rows_between(tex, start=r"\startdata", end=r"\enddata"):
    body = tex.split(start, 1)[1].split(end, 1)[0]
    body = re.sub(r"(?m)%.*$", "", body)
    out = []
    for line in body.split("\\\\"):
        line = line.strip()
        if "&" not in line:
            continue
        f = [x.strip() for x in line.split("&")]
        m = re.match(r"^(\d{6,9})", f[0])
        if m:
            out.append((int(m.group(1)), f[0], f[1:]))
    return out


def t14_hours(P, aRs, b, k):
    """Circular-orbit total duration (hours)."""
    if not np.all(np.isfinite([P, aRs, b, k])) or aRs <= 1:
        return np.nan
    inc = np.arccos(min(b / aRs, 1))
    x = np.sqrt(max((1 + k) ** 2 - b ** 2, 0)) / (aRs * np.sin(inc))
    return P / np.pi * np.arcsin(min(x, 1)) * 24


def aRs_from(P, M, R):
    G = 2.959122e-4  # AU^3 / (Msun d^2)
    a_au = (G * M * P ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return a_au * 215.032 / R


def parse_km19():
    out = []
    for fn, tab in [("table.tex", "clean"), ("table_rm.tex", "flagged")]:
        tex = (RAW / "arxiv_1904.04980" / fn).read_text(encoding="utf-8", errors="ignore")
        for kic, c0, f in rows_between(tex):
            koi = f[0]
            Kp, M, R, r = (first_num(x) for x in f[1:5])
            t0 = first_num(f[5])
            P = first_num(f[6])
            Pup, Plo = errs(f[6])
            b = first_num(f[7])
            k = r / RSUN_REARTH / R
            # KM19: 'P' with tiny errors = two observed transits (P = interval); flag in daggers not reliable
            n_tr = 2 if (np.isfinite(Pup) and Pup < 1.0) else 1
            aRs = aRs_from(P, M, R)
            out.append(dict(source="KM19", table=tab, KIC=kic, KOI=re.sub(r"\$.*?\$", "", koi).strip(),
                            Kp=Kp, Rstar=R, Mstar=M, Rp_re=r, t0_bkjd=t0, P_cat=P, P_up=Pup, P_lo=Plo,
                            b=b, depth_ppm=k ** 2 * 1e6, dur_h=t14_hours(P, aRs, b, k), dur_src="calc",
                            n_kepler_transits=n_tr, flags=re.sub(r"^\d+", "", c0)))
    return pd.DataFrame(out)


def parse_w15():
    out = []
    for fn, ntr in [("orbital_params_1.tex", 1), ("orbital_params_2.tex", 2)]:
        tex = (RAW / "arxiv_1512.02559" / fn).read_text(encoding="utf-8", errors="ignore")
        for kic, c0, f in rows_between(tex):
            if ntr == 1:
                P = first_num(f[1]); prange = f[2]
                aRs, inc, k, rp, ep = (first_num(x) for x in f[3:8])
                com = f[-1]
            else:
                P = first_num(f[1]); prange = ""
                aRs, inc, k, rp, ep = (first_num(x) for x in f[2:7])
                com = f[-1]
            b = aRs * np.cos(np.radians(inc))
            out.append(dict(source="W15", table=f"{ntr}tr", KIC=kic, KOI=f[0], Rp_re=rp, t0_bkjd=ep, P_cat=P,
                            P_range=prange, b=b, depth_ppm=k ** 2 * 1e6, dur_h=t14_hours(P, aRs, b, k),
                            dur_src="calc(W15 a/R*,i)", n_kepler_transits=ntr, flags=com))
    return pd.DataFrame(out)


def parse_fm16():
    tex = (RAW / "arxiv_1607.08237" / "ms.tex").read_text(encoding="utf-8", errors="ignore")
    seg = tex.split(r"\label{tab:catalog}", 1)[1]
    out = []
    for kic, c0, f in rows_between(seg):
        Teff, R, Kp, Pyr, t0, rad, dur, b, Teq, prp = (first_num(x) for x in f[:10])
        Pup, Plo = errs(f[3])
        two = "c" in c0.split("$")[1] if "$" in c0 else False
        k = rad * RJ_RSUN / R  # R_J -> R_sun
        out.append(dict(source="FM16", KIC=kic, KOI=f[10] if len(f) > 10 else "", Kp=Kp, Rstar=R,
                        Rp_re=rad * 11.209, t0_bkjd=t0, P_cat=Pyr * 365.25, P_up=Pup * 365.25, P_lo=Plo * 365.25,
                        b=b, depth_ppm=k ** 2 * 1e6, dur_h=dur, dur_src="FM16", n_kepler_transits=2 if two else 1,
                        flags=f"Pr_planet={prp}"))
    return pd.DataFrame(out)


def parse_hzw19():
    tex = (RAW / "arxiv_1901.01974" / "param_tab.tex").read_text(encoding="utf-8", errors="ignore")
    out = []
    for kic, c0, f in rows_between(tex):
        R, Kp, t0, P, rad, b, snr, fap = (first_num(x) for x in f[:8])
        Pup, Plo = errs(f[3])
        k = rad * 0.10045 / R
        aRs = aRs_from(P, R ** 1.0, R)  # M ~ R for dwarfs (scoping)
        ntr = 2 if (np.isfinite(Pup) and Pup < 1) else 1
        out.append(dict(source="HZW19", KIC=kic, KOI=f[8] if len(f) > 8 else "", Kp=Kp, Rstar=R, Rp_re=rad * 11.209,
                        t0_bkjd=t0, P_cat=P, P_up=Pup, P_lo=Plo, b=b, depth_ppm=k ** 2 * 1e6,
                        dur_h=t14_hours(P, aRs, b, k), dur_src="calc(M=R)", n_kepler_transits=ntr,
                        flags=f"FAP={fap}"))
    return pd.DataFrame(out)


def parse_u16():
    tex = (RAW / "arxiv_1602.07848" / "main.tex").read_text(encoding="utf-8", errors="ignore")
    seg = tex.split(r"\label{tab:list}", 1)[1].split(r"\end{tabular}", 1)[0]
    seg = re.sub(r"(?m)%.*$", "", seg)
    out, last = [], {}
    for line in seg.split("\\\\"):
        f = [x.strip() for x in line.split("&")]
        if len(f) < 9 or not re.match(r"^\d{6,9}$", f[0].split()[-1] if f[0] else ""):
            continue
        kic = int(f[0].split()[-1])
        Kp = first_num(f[3]); Kp = Kp if np.isfinite(Kp) else last.get(kic, np.nan)
        last[kic] = Kp
        out.append(dict(source="U16", KIC=kic, KOI=re.sub(r"\\tablenotemark\{.\}", "", f[1]), Kp=Kp,
                        t0_bkjd=first_num(f[4]), depth_ppm=first_num(f[5]), P_min_cat=first_num(f[8]),
                        dur_h=np.nan, dur_src="", n_kepler_transits=1, flags=f[2]))
    d = pd.DataFrame(out)
    # durations from the U16 MCMC table (T, hours) where available, matched on Tc
    fit = tex.split(r"\label{tab:list_fit}", 1)[1].split(r"\end{tabular}", 1)[0]
    fit = re.sub(r"(?m)%.*$", "", fit)
    for line in fit.split("\\\\"):
        f = [x.strip() for x in line.split("&")]
        if len(f) < 12:
            continue
        tc, T = first_num(f[4]), first_num(f[9])
        m = (d.t0_bkjd - tc).abs() < 0.05
        d.loc[m, "dur_h"] = T
        d.loc[m, "dur_src"] = "U16"
    return d


def hc24():
    # KOI 1271.02 (KIC 8631160), single transit at ~679 BKJD, Rp 5.32 Re, R* 1.48 Rsun (arXiv:2403.03427 text)
    k = 5.32 / RSUN_REARTH / 1.48
    return pd.DataFrame([dict(source="HC24", KIC=8631160, KOI="1271.02", Rstar=1.48, Rp_re=5.32, t0_bkjd=679.0,
                              depth_ppm=k ** 2 * 1e6, dur_h=np.nan, dur_src="", n_kepler_transits=1,
                              flags="t0 approximate (text: 'day 679 BKJD')")])


def main():
    parts = [parse_km19(), parse_w15(), parse_fm16(), parse_hzw19(), parse_u16(), hc24()]
    d = pd.concat(parts, ignore_index=True)
    d["t0_btjd"] = d.t0_bkjd - 2167.0
    d.to_csv(HERE / "data" / "events_by_source.csv", index=False)
    print(d.groupby("source").agg(rows=("KIC", "size"), stars=("KIC", "nunique"),
                                  duos=("n_kepler_transits", lambda x: (x == 2).sum())))
    print(d[["source", "KIC", "t0_bkjd", "P_cat", "depth_ppm", "dur_h", "n_kepler_transits"]].round(2).to_string())


if __name__ == "__main__":
    main()
