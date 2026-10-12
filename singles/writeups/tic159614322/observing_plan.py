"""
Observing plan for the predicted TIC 159614322 transit of 2026-12-28 from Providence, RI
(41.82 N, 71.41 W), with a comparison for 2027-03-28 (only the 89.965-d period predicts that one).

Inputs: predictions.csv (two-transit ephemeris, planning sigma = formal x sqrt(chi2) of the
three-transit fit; ephemeris.py). T14 = 9.7 h, depth 3.3 ppt (joint fit, analysis.json).

Noise model per exposure, for 8/11/14-inch SCTs (D = 20.3/27.9/35.6 cm, 34% central obstruction):
  photons   V = 11.8 star (Gaia G 11.63, BP-RP 0.87 -> V via Riello+2021); unfiltered CMOS,
            1.05e6 e/s/cm2 at V = 0 (bandwidth ~3000 A x total throughput ~0.35; uncertain by ~2x)
  sky       18.5 mag/arcsec2 (moonlit suburban sky), photometric aperture radius 8"
  scint.    Young (1967) with the Osborn+2015 empirical factor 1.5, sea level:
            0.09 D^(-2/3) X^1.75 / sqrt(2 t) * 1.5   (D in cm)
  comps     ensemble of 4 comparison stars of similar brightness: noise x sqrt(1 + 1/4)
  60-s exposures, 70-s cadence.
White-noise only. The real limit for a 9.7-h event is red noise: airmass-dependent
differential extinction, focus and guiding drifts, flat-field errors. This is discussed in the
plan; it is not modelled here.

Output: observing_plan.json, fig_providence.png, observing_plan_providence.md
Run:    python singles/writeups/tic159614322/observing_plan.py
"""
import json
from pathlib import Path

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time

HERE = Path(__file__).resolve().parent
TGT = SkyCoord(220.00977 * u.deg, 65.484041 * u.deg)
LOC = EarthLocation(lat=41.82 * u.deg, lon=-71.41 * u.deg, height=20 * u.m)
T14_H, DEPTH = 9.7, 3.3e-3
V = 11.8
SCOPES = {'8"': 20.3, '11"': 27.9, '14"': 35.6}
EXP, CAD = 60.0, 70.0
SKY_MAG, AP_R = 18.5, 8.0
NCOMP = 4


def sigma_exp(D, X):
    area = np.pi * (D / 2) ** 2 * (1 - 0.34 ** 2)
    star = 1.05e6 * 10 ** (-0.4 * V) * area * EXP
    sky = 1.05e6 * 10 ** (-0.4 * SKY_MAG) * area * np.pi * AP_R ** 2 * EXP
    phot = np.sqrt(star + sky) / star
    scint = 0.09 * D ** (-2 / 3) * X ** 1.75 / np.sqrt(2 * EXP) * 1.5
    return np.sqrt(phot ** 2 + scint ** 2) * np.sqrt(1 + 1 / NCOMP)


def night(mid_utc, sig_min):
    tm = Time(mid_utc, scale="utc")
    tt = tm + np.arange(-8, 8.01, 0.25) * u.hour
    fr = AltAz(obstime=tt, location=LOC)
    alt = TGT.transform_to(fr).alt.deg
    sun = get_sun(tt).transform_to(fr).alt.deg
    moon = get_body("moon", tt, LOC)
    malt = moon.transform_to(fr).alt.deg
    msep = moon.separation(TGT, origin_mismatch="ignore").deg
    elong = get_sun(tm).separation(get_body("moon", tm, LOC), origin_mismatch="ignore").rad
    illum = 0.5 * (1 - np.cos(elong))
    X = 1 / np.sin(np.radians(np.clip(alt, 5, 90)))
    rows = []
    for t_, a, x, s, ma, ms in zip(tt, alt, X, sun, malt, msep):
        dt_h = (t_ - tm).to(u.hour).value
        state = ("in transit" if abs(dt_h) < T14_H / 2 - sig_min / 60 else
                 "out of transit" if abs(dt_h) > T14_H / 2 + sig_min / 60 else
                 ("ingress window" if dt_h < 0 else "egress window"))
        row = dict(UTC=t_.iso[5:16], EST_UTCminus5=(t_ - 5 * u.hour).iso[11:16], alt=round(float(a), 1), airmass=round(float(x), 2),
                   sun_alt=round(float(s), 1), moon_alt=round(float(ma), 1), moon_sep=round(float(ms)), state=state, dt_h=dt_h)
        for name, D in SCOPES.items():
            row[f"sig10min_{name}_ppt"] = round(float(sigma_exp(D, x) / np.sqrt(600 / CAD)) * 1e3, 2)
        rows.append(row)
    df = pd.DataFrame(rows)
    dark = df.sun_alt < -12
    summ = dict(mid_UTC=mid_utc, sigma_plan_min=sig_min, moon_illum=round(float(illum), 2),
                ingress_UTC=(tm - T14_H / 2 * u.hour).iso[11:16], egress_UTC=(tm + T14_H / 2 * u.hour).iso[11:16],
                dark_UTC=[df[dark].UTC.iloc[0] if dark.any() else None, df[dark].UTC.iloc[-1] if dark.any() else None])

    def step_snr(D, t_c, side=1.5):
        """SNR of the 3.3-ppt step at time t_c (h from mid), comparing 'side' h before and after,
        dark time only, white noise only."""
        out = []
        for lo, hi in ((t_c - side, t_c), (t_c, t_c + side)):
            m = dark & (df.dt_h >= lo) & (df.dt_h < hi)
            if m.sum() < 2:
                return None
            n = m.sum() * 900 / CAD
            out.append(np.mean([sigma_exp(D, x) for x in df[m].airmass]) / np.sqrt(n))
        return round(float(DEPTH / np.hypot(*out)), 1)

    def transit_snr(D):
        """In-transit (dark, fully inside) vs all dark out-of-transit, white noise only."""
        inn = dark & (df.dt_h.abs() < T14_H / 2 - sig_min / 60)
        oot = dark & (df.dt_h.abs() > T14_H / 2 + sig_min / 60)
        if inn.sum() < 2 or oot.sum() < 2:
            return None
        si = np.sqrt(np.sum([sigma_exp(D, x) ** -2 * 900 / CAD for x in df[inn].airmass])) ** -1
        so = np.sqrt(np.sum([sigma_exp(D, x) ** -2 * 900 / CAD for x in df[oot].airmass])) ** -1
        return round(float(DEPTH / np.hypot(si, so)), 1), round(float(oot.sum() * 0.25), 2)

    for name, D in SCOPES.items():
        summ[name] = dict(sigma_60s_ppt_at_X1_2=round(float(sigma_exp(D, 1.2)) * 1e3, 2),
                          sigma_60s_ppt_at_X3=round(float(sigma_exp(D, 3.0)) * 1e3, 2),
                          ingress_step_snr=step_snr(D, -T14_H / 2), egress_step_snr=step_snr(D, T14_H / 2),
                          transit_snr_white_and_oot_hours=transit_snr(D))
    return df, summ


def main():
    pr = pd.read_csv(HERE / "predictions.csv")
    e1 = pr[(pr.alias.str.contains("n=8")) & (pr.mid_UTC.str.startswith("2026-12-28"))].iloc[0]
    e2 = pr[(pr.alias.str.contains("n=8")) & (pr.mid_UTC.str.startswith("2027-03-28"))].iloc[0]
    out = {}
    fig, axs = plt.subplots(2, 1, figsize=(9, 6.2))
    for ax, e in zip(axs, (e1, e2)):
        df, s = night(e.mid_UTC, int(e.sigma_plan_min))
        out[e.mid_UTC[:10]] = dict(summary=s, table=df.drop(columns="dt_h").to_dict("records"))
        x = df.dt_h
        ax.plot(x, df.alt, "-", color="#1d2b3a", label="target altitude")
        ax.plot(x, df.moon_alt, ":", color="#5d6b7c", label=f"Moon ({s['moon_illum']:.0%} lit)")
        ax.fill_between(x, 0, 90, where=df.sun_alt > -12, color="#e9edf2", step="mid", label="Sun above -12 deg")
        h = T14_H / 2
        sg = int(e.sigma_plan_min) / 60
        ax.axvspan(-h + sg, h - sg, color="#d9682b", alpha=0.15, label="certainly in transit")
        for c in (-h, h):
            ax.axvspan(c - sg, c + sg, color="#d9682b", alpha=0.35)
        ax.axvline(0, color="#d9682b", lw=0.8)
        ax.set_ylim(0, 90)
        ax.set_xlim(-8, 8)
        ax.set_ylabel("altitude (deg)", fontsize=8.5)
        ttl = "both periods" if e.mid_UTC.startswith("2026") else "89.965-d period only"
        ax.set_title(f"Providence, RI: transit mid {e.mid_UTC} UTC +/- {int(e.sigma_plan_min)} min ({ttl})", fontsize=9)
        tick = np.arange(-8, 8.01, 2)
        ax.set_xticks(tick)
        off, zone = (4, "EDT") if e.mid_UTC[5:7] in ("04", "05", "06", "07", "08", "09", "10") or e.mid_UTC.startswith("2027-03-28") else (5, "EST")
        ax.set_xticklabels([(Time(e.mid_UTC) + t * u.hour - off * u.hour).iso[11:16] for t in tick], fontsize=8)
        ax.set_xlabel(f"local time ({zone})", fontsize=8.5)
        ax.legend(fontsize=6.5, loc="upper left", ncol=3)
    fig.tight_layout()
    fig.savefig(HERE / "fig_providence.png", dpi=150)
    (HERE / "observing_plan.json").write_text(json.dumps(out, indent=1, default=str))
    for k, v in out.items():
        print(k, json.dumps(v["summary"], indent=1))
        t = pd.DataFrame(v["table"])
        print(t[t.index % 4 == 0][["UTC", "EST_UTCminus5", "alt", "airmass", "sun_alt", "moon_alt", "state", 'sig10min_8"_ppt', 'sig10min_14"_ppt']].to_string(index=False))


if __name__ == "__main__":
    main()
