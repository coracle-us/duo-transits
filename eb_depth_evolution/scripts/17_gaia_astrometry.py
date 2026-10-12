"""Gaia DR3 astrometric quality of the target vs. comparison stars of similar magnitude/colour within 2 deg; NSS tables;
wide co-moving companions within 10 arcmin; LAMOST MRS binary catalogue J/ApJS/258/26."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
sid = 3383995253529149056
t = Gaia.launch_job(f"""SELECT source_id, phot_g_mean_mag, bp_rp, parallax, parallax_over_error, ruwe, astrometric_excess_noise,
 astrometric_excess_noise_sig, astrometric_n_good_obs_al, astrometric_chi2_al, ipd_gof_harmonic_amplitude, ipd_frac_multi_peak,
 visibility_periods_used, rv_amplitude_robust, non_single_star, phot_variable_flag, duplicated_source
 FROM gaiadr3.gaia_source WHERE source_id={sid}""").get_results()
for c in t.colnames: print(f"  {c} = {t[c][0]}")
ref = Gaia.launch_job("""SELECT ruwe, astrometric_excess_noise AS aen, astrometric_excess_noise_sig AS aens FROM gaiadr3.gaia_source
 WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(97.98,25.48,2.0)) AND phot_g_mean_mag BETWEEN 10.5 AND 11.5 AND bp_rp BETWEEN 0.2 AND 0.7""").get_results()
ruwe = np.array(ref["ruwe"]); aen = np.array(ref["aen"]); aens = np.array(ref["aens"])
print(f"comparison stars (G 10.5-11.5, BP-RP 0.2-0.7, within 2 deg): N={len(ref)}")
print(f"  RUWE median {np.nanmedian(ruwe):.3f}, 16-84% {np.nanpercentile(ruwe,16):.3f}-{np.nanpercentile(ruwe,84):.3f}; target percentile {np.mean(ruwe < t['ruwe'][0])*100:.0f}")
print(f"  AEN median {np.nanmedian(aen):.3f} mas, 84% {np.nanpercentile(aen,84):.3f}; target AEN {t['astrometric_excess_noise'][0]:.3f} -> percentile {np.mean(aen < t['astrometric_excess_noise'][0])*100:.0f}")
print(f"  AEN_sig median {np.nanmedian(aens):.1f}; target {t['astrometric_excess_noise_sig'][0]:.1f} -> percentile {np.mean(aens < t['astrometric_excess_noise_sig'][0])*100:.0f}")
for tab in ["gaiadr3.nss_two_body_orbit", "gaiadr3.nss_acceleration_astro", "gaiadr3.nss_non_linear_spectro", "gaiadr3.vari_eclipsing_binary"]:
    r = Gaia.launch_job(f"SELECT COUNT(*) AS n FROM {tab} WHERE source_id={sid}").get_results()
    print(f"  {tab}: {r['n'][0]} rows")
w = Gaia.launch_job("""SELECT source_id, phot_g_mean_mag AS G, parallax, parallax_error, pmra, pmdec,
 DISTANCE(POINT(97.98012,25.48210), POINT(ra,dec))*3600 AS sep FROM gaiadr3.gaia_source
 WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(97.98012,25.48210,0.1667)) AND parallax BETWEEN 1.2 AND 1.85
 AND ABS(pmra+2.03) < 1.5 AND ABS(pmdec+3.99) < 1.5 ORDER BY sep""").get_results()
print("co-moving candidates within 10 arcmin (plx 1.2-1.85, |dPM| < 1.5 mas/yr):"); print(w)
v = Vizier(columns=["**"], row_limit=20)
c = SkyCoord(97.98012, 25.48210, unit="deg")
for cat in ["J/ApJS/258/26", "J/ApJS/256/14", "V/162"]:
    try:
        r = v.query_region(c, radius=5*u.arcsec, catalog=cat)
        print(cat, "tables:", list(r.keys()))
        for k in r.keys():
            if "258" in k: print(r[k])
    except Exception as e: print(cat, "ERR", e)
