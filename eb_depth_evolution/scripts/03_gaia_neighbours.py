"""Gaia DR3 neighbours within 3 arcmin, plus NSS / wide-companion search within 60 arcsec with matching parallax/PM."""
import warnings; warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
q = """SELECT source_id, ra, dec, phot_g_mean_mag AS G, bp_rp, parallax, parallax_error, pmra, pmdec, ruwe,
 radial_velocity, non_single_star, phot_variable_flag, has_epoch_photometry,
 DISTANCE(POINT(97.98012,25.48210), POINT(ra,dec))*3600 AS sep
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(97.98012,25.48210,0.05)) ORDER BY sep"""
t = Gaia.launch_job(q).get_results()
t.write("data/gaia_neighbours.csv", overwrite=True)
print(t[t["G"]<18][["source_id","G","bp_rp","parallax","parallax_error","pmra","pmdec","ruwe","sep","non_single_star","has_epoch_photometry"]])
