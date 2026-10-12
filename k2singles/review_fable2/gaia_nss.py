"""Pull Gaia DR3 NSS two-body-orbit solutions and gaia_source RV stats for the review targets."""
from astroquery.gaia import Gaia
import pandas as pd, sys
ids = {660920219027943552: "EPIC 211892898", 602455410504038016: "EPIC 211390677", 604962365734255744: "EPIC 211411112",
       3358696792346046080: "EPIC 202072917", 3377970891880909440: "EPIC 202071902", 3383995253529149056: "EPIC 202126877"}
q = f"""SELECT source_id, nss_solution_type, period, period_error, t_periastron, t_periastron_error, eccentricity, eccentricity_error,
 arg_periastron, arg_periastron_error, semi_amplitude_primary, semi_amplitude_primary_error, center_of_mass_velocity, center_of_mass_velocity_error,
 inclination, inclination_error, goodness_of_fit, significance, obj_func, efficiency, flags, bit_index, corr_vec
 FROM gaiadr3.nss_two_body_orbit WHERE source_id IN ({','.join(map(str, ids))})"""
j = Gaia.launch_job(q); t = j.get_results().to_pandas()
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(t.T)
t.to_csv("k2singles/review_fable2/gaia_nss.csv", index=False)
q2 = f"""SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, ruwe, non_single_star, radial_velocity, radial_velocity_error, rv_nb_transits, rv_amplitude_robust,
 rv_chisq_pvalue, rv_renormalised_gof, rv_template_teff, rv_template_logg, parallax, parallax_error, astrometric_excess_noise, astrometric_excess_noise_sig, phot_variable_flag, teff_gspphot, logg_gspphot
 FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(map(str, ids))})"""
t2 = Gaia.launch_job(q2).get_results().to_pandas(); t2["name"] = t2.source_id.map(ids)
print(t2.T); t2.to_csv("k2singles/review_fable2/gaia_source.csv", index=False)
q3 = f"""SELECT * FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN ({','.join(map(str, ids))})"""
t3 = Gaia.launch_job(q3).get_results().to_pandas(); print("vari_eclipsing_binary rows:", len(t3)); print(t3.T if len(t3) else "")
q4 = f"""SELECT source_id, best_class_name, best_class_score FROM gaiadr3.vari_classifier_result WHERE source_id IN ({','.join(map(str, ids))})"""
t4 = Gaia.launch_job(q4).get_results().to_pandas(); print(t4)
