"""
Saved inputs for numbers quoted in the K2 RNAAS draft (fixes review_fable3 items 2 and 10):
  1. Gaia DR3 host parameters (gaia_source + astrophysical_parameters) for the stars whose stellar
     parameters the text quotes -> k2singles/data/gaia_hosts.csv
  2. Predicted transit/eclipse times from the fitted ephemerides, as BJD_TDB and as geocentric UTC
     (the barycentric-to-geocentric light-travel time is removed with astropy)
     -> k2singles/rnaas/predictions.csv
Usage:  python k2singles/scripts/hosts_and_predictions.py
"""
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
from gaia_flags import query  # noqa: E402

EPH = ROOT / "k2singles" / "results" / "ephemerides"
HOSTS = {  # EPIC: Gaia DR3 source_id (from k2singles/data/lead_status.csv)
    211924561: 660640844292494336,
    211503363: 650964454772110976,
    246194159: 2633217515438236544,
    211892898: 660920219027943552,
    211390677: 602455410504038016,
}
# (label, EPIC, period divisor n, phase of the event, earliest BTJD)
PREDICT = [
    ("TOI-5560.01 transit, both periods", 211924561, 1, 0.0, 4325.0),
    ("TOI-5560.01 transit, 70.25 d only", 211924561, 2, 0.0, 4400.0),
    ("HD 220967 secondary eclipse", 246194159, 1, None, 4325.0),
    ("HD 220967 primary eclipse", 246194159, 1, 0.0, 4325.0),
]


def main():
    rows = []
    for epic, sid in HOSTS.items():
        g = query(f"""SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.parallax, g.ruwe, g.radial_velocity_error,
                      g.rv_nb_transits, g.rv_chisq_pvalue, g.rv_renormalised_gof, g.rv_amplitude_robust,
                      a.teff_gspphot, a.logg_gspphot, a.radius_gspphot, a.radius_flame, a.mass_flame
                      FROM gaiadr3.gaia_source g LEFT JOIN gaiadr3.astrophysical_parameters a
                      ON g.source_id = a.source_id WHERE g.source_id = {sid}""")
        r = g.iloc[0].to_dict()
        r["source_id"] = str(sid)
        r["EPIC"] = epic
        rows.append(r)
    hosts = pd.DataFrame(rows)
    hosts.to_csv(ROOT / "k2singles" / "data" / "gaia_hosts.csv", index=False)
    print(hosts[["EPIC", "phot_g_mean_mag", "radius_flame", "mass_flame", "logg_gspphot", "teff_gspphot"]].to_string())

    geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)
    out = []
    for label, epic, n, phase, after in PREDICT:
        j = json.loads((EPH / f"{epic}.json").read_text())
        P, sP, T0, sT0 = j["P"] / n, j["sigma_P"] / n, j["T0"], j["sigma_T0"]
        if phase is None:
            phase = j["joint"]["phase_sec"]
        k = np.ceil((after - T0) / P - phase)
        if n == 2:                                  # odd cycles of P/2 occur only for the short period
            k = k if k % 2 == 1 else k + 1
        t = T0 + (k + phase) * P
        s = float(np.hypot(sT0, (k + phase) * sP))
        h = hosts[hosts.EPIC == epic].iloc[0]
        c = SkyCoord(h.ra * u.deg, h.dec * u.deg)
        tb = Time(t + 2457000, format="jd", scale="tdb", location=geo)
        tgeo = tb - tb.light_travel_time(c, kind="barycentric")
        out.append(dict(label=label, EPIC=epic, BTJD=round(t, 5), BJD_TDB=round(t + 2457000, 5), sigma_d=round(s, 5),
                        sigma_min=round(s * 1440, 1), utc_geocentric=tgeo.utc.iso[:16],
                        ltt_min=round(float(tb.light_travel_time(c, kind='barycentric').to(u.min).value), 1)))
    pr = pd.DataFrame(out)
    pr.to_csv(ROOT / "k2singles" / "rnaas" / "predictions.csv", index=False)
    print(pr.to_string())


if __name__ == "__main__":
    main()
