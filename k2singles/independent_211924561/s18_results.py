"""Collect key numbers into results.json."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(H, f)))
fr, al, bp, dc, ap, sl = L("fit_results.json"), L("alias_table.json"), L("b_profile.json"), L("density_check.json"), L("aperture_tests.json"), L("secondary_limit_P70.json")
j = fr["joint4_primary"]; ja = fr["joint4_alt"]
surv = [p for p in al["periods"] if p["status"] == "survives"]
res = dict(
    target=dict(epic=211924561, tic=238606535, ra=132.652546, dec=19.091365, gaia_dr3=660640844292494336, G=10.696,
                note="RA/Dec given in the task (131.0849, +20.3497) are wrong; header and Gaia position is 132.6525, +19.0914",
                flame_R=1.739, flame_M=1.424, teff_gspphot=6683, rho_gcc=dc["rho_star_gcc"]),
    events_BTJD=[dict(src=s, tc=t[0], err=0.5 * (t[1] + t[2])) for s, t in zip(j["events"], j["tc"])],
    events_alt_reductions_BTJD=[dict(src=s, tc=t[0], err=0.5 * (t[1] + t[2])) for s, t in zip(ja["events"], ja["tc"])],
    blind_prediction="3-event ephemeris (C5,S46,S72) predicted C16 epoch at 1163.3555; observed 1163.3554+-0.0013",
    shape_joint=dict(T14_h=j["T14_h"], b=j["b"], k=j["k"], model_depth_ppt=j["depth_ppt"], T23_h=j["T23_h"]),
    b_profile_dchi2={str(p["b"]): round(p["dchi2"], 2) for p in bp["b_profile"] if p["b"] <= 1.1},
    surviving_periods=[dict(P=p["P"], sP=p["sP"], T0=p["T0"], sT0=p["sT0"], N=p["N"], epochs=p["epochs"],
                            rho_circ_pct_2p5_16_50_84_97p5=bp["per_P"][f"{p['P']:.5f}"]["rho_circ_pct"] if f"{p['P']:.5f}" in bp["per_P"] else None,
                            e_min_pct=bp["per_P"][f"{p['P']:.5f}"]["e_min_pct"] if f"{p['P']:.5f}" in bp["per_P"] else None) for p in surv],
    density_tied_fits=dc["per_period"],
    excluded_periods_count=sum(p["status"] == "excluded" for p in al["periods"]),
    excluded_shortest_aliases=[dict(P=p["P"], n_excl_epochs=len(p["excluding_epochs"]), first=[(e["fam"], e["tc"]) for e in p["excluding_epochs"][:4]]) for p in al["periods"][2:5]],
    centroid_tests=dict(S46=ap["TESS_S46"]["centroid_offset_pix"], S72=ap["TESS_S72"]["centroid_offset_pix"], K2_full_aperture_depth=ap["K2_C5"]["all r<=2.3"]),
    secondary_P70_phase05=dict(mean_ppt=sl["mean_ppt"], err_ppt=sl["err_ppt"], upper3sig_ppt=sl["upper3sig_ppt"]),
)
json.dump(res, open(os.path.join(H, "results.json"), "w"), indent=1)
print(json.dumps(res["surviving_periods"], indent=1))
