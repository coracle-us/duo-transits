"""Gaia DR3 for a target (usage: s12_gaia.py TIC): the target's astrophysical parameters and all sources within
200" (9.5 TESS pixels). For each neighbour, the eclipse depth it would need to mimic the observed
dip if it were the source: depth_needed = depth_obs * (F_target + F_nbr) / F_nbr (T ~ G - small
colour term ignored; G-band flux ratio used)."""
import sys, os, json, warnings; warnings.filterwarnings("ignore")
import numpy as np
from astroquery.gaia import Gaia
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import TARGETS
from astroquery.mast import Catalogs
TIC = int(sys.argv[1])
SID = int(Catalogs.query_criteria(catalog="TIC", ID=TIC)[0]["GAIA"])
RA, DEC = TARGETS[TIC]["ra"], TARGETS[TIC]["dec"]
DEPTH = round(json.load(open(os.path.join(ROOT, "out", f"shape_{TIC}.json")))["depth"], 4)  # fitted (deblended) depth
_lines = []
_print = print
def print(*a):
    _lines.append(" ".join(str(x) for x in a)); _print(*a)
q = f"""SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.parallax_error,
 g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.non_single_star, g.phot_variable_flag,
 DISTANCE(POINT(g.ra, g.dec), POINT({RA}, {DEC}))*3600 AS sep
 FROM gaiadr3.gaia_source AS g
 WHERE 1=CONTAINS(POINT(g.ra, g.dec), CIRCLE({RA}, {DEC}, {200/3600}))
 ORDER BY sep"""
nb = Gaia.launch_job(q).get_results()
ap = Gaia.launch_job(f"""SELECT source_id, teff_gspphot, logg_gspphot, radius_gspphot, radius_flame,
 radius_flame_lower, radius_flame_upper, mass_flame, mass_flame_lower, mass_flame_upper,
 lum_flame, age_flame, evolstage_flame, distance_gspphot
 FROM gaiadr3.astrophysical_parameters WHERE source_id = {SID}""").get_results()
out = {}
print("Target astrophysical parameters (Gaia DR3):")
for c in (ap.colnames if len(ap) else []):
    v = ap[c][0]
    print(f"  {c:22s} {v}")
    out[c] = None if np.ma.is_masked(v) else float(v)
t = nb[nb["source_id"] == SID][0]
gt = float(t["phot_g_mean_mag"])
print(f"\nTIC {TIC}  Gaia DR3 {SID}  depth used {DEPTH*1e6:.0f} ppm\nTarget: G={gt:.3f} BP-RP={t['bp_rp']:.3f} plx={t['parallax']:.4f}±{t['parallax_error']:.4f} "
      f"RUWE={t['ruwe']:.3f} ipd_multi_peak={t['ipd_frac_multi_peak']} ipd_gof_harm={t['ipd_gof_harmonic_amplitude']:.4f} "
      f"non_single_star={t['non_single_star']} variable={t['phot_variable_flag']}")
out["target"] = dict(G=gt, ruwe=float(t["ruwe"]), plx=float(t["parallax"]), bp_rp=float(t["bp_rp"]),
                     ipd_frac_multi_peak=int(t["ipd_frac_multi_peak"]), non_single_star=int(t["non_single_star"]))
ft = 10 ** (-0.4 * gt)
print(f"\nNeighbours within 200\" ({len(nb)-1}); depth needed if the neighbour were the eclipsing source:")
rows = []
for r in nb:
    if r["source_id"] == SID:
        continue
    g = float(r["phot_g_mean_mag"]) if not np.ma.is_masked(r["phot_g_mean_mag"]) else np.nan
    fn = 10 ** (-0.4 * g)
    need = DEPTH * (ft + fn) / fn
    rows.append(dict(source_id=int(r["source_id"]), ra=float(r["ra"]), dec=float(r["dec"]), sep=float(r["sep"]),
                     G=g, dG=g - gt, need=need))
    if need < 1.0 or float(r["sep"]) < 40:
        print(f"  {r['source_id']} sep={r['sep']:6.1f}\" G={g:6.2f} dG={g-gt:+5.2f} needed depth={need*100:6.1f}%"
              f"{'  <-- could mimic' if need < 0.5 else ''}")
out["neighbours"] = rows
out["source_id"] = SID
json.dump(out, open(os.path.join(ROOT, "out", f"gaia_{TIC}.json"), "w"), indent=1)
open(os.path.join(ROOT, "out", f"gaia_{TIC}.txt"), "w").write("\n".join(_lines) + "\n")
