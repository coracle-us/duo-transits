"""Collate results into results.json (verdict table + key numbers) and list upcoming predicted times."""
import json, os
import numpy as np
from astropy.time import Time
OUT = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(OUT, f)))
fk, ep, et, s21, ing, gap, hde, pix = (L("fit_known.json"), L("epochs.json"), L("epoch_tests.json"), L("vet_s21.json"),
                                       L("ingress_fap.json"), L("gap_epochs_ownphot.json"), L("highn_density_ephem.json"), L("pixel_level.json"))
cen = L("centroid_neighbours.json")
P = fk["primary"]
verdict = {
 1: ("allowed by coverage (no third epoch in data); inconsistent with the S21 dip if that dip is a transit", []),
 2: ("allowed by coverage (no third epoch in data); inconsistent with the S21 dip if that dip is a transit", []),
 3: ("excluded", [1943.735, 2663.455, 3383.174]),
 4: ("supported (S21 dip at predicted time; no other epoch in data)", [1883.758]),
 5: ("allowed by coverage (no third epoch in data); inconsistent with the S21 dip if that dip is a transit", []),
 6: ("excluded", [1943.735, 2663.455, 3383.174]),
 7: ("excluded", [1909.463, 2629.182, 3348.902]),
 8: ("supported, degenerate with n=4 (its extra epochs all fall in gaps)", [1883.758]),
 9: ("excluded", [1943.735, 2663.455, 3383.174]),
 10: ("excluded", [1919.744, 2639.464, 3359.184]),
 11: ("excluded", [1900.116, 2619.835]),
 12: ("excluded (despite the S21 epoch matching)", [1883.758, 1943.735, 2663.455, 3383.174]),
}
table = []
for n in range(1, 13):
    eps = sorted({round(r["tc"], 3) for r in ep["epochs"] if r["n"] == n and not r["known"]})
    table.append(dict(n=n, P=ep["DT"] / n, predicted_epochs_in_or_near_data=eps, verdict=verdict[n][0], decisive_epochs=verdict[n][1]))
upcoming = {}
for n in (4, 8):
    e = hde["ephemeris"][str(n)]
    ks = np.arange(int((3400 - e["T0"]) / e["P"]) + 1, int((4500 - e["T0"]) / e["P"]) + 1)
    upcoming[n] = [dict(btjd=round(e["T0"] + k * e["P"], 3), sig_d=round(float(np.hypot(e["sig_T0"], k * e["sig_P"])), 3),
                        utc=Time(e["T0"] + k * e["P"] + 2457000, format="jd", scale="tdb").utc.iso[:16]) for k in ks]
res = dict(
 target="TIC 159614322",
 known_events=dict(tc1=P["tc1"], tc2=P["tc2"], dT=P["dT"], depth_rp2=P["depth_rp2"], T14_hr=P["T14_hr"], b=P["b"],
                   note="joint fit S14 TESS-SPOC 30-min + S41 SPOC 2-min PDCSAP; values [median, -err, +err]; QLP-only and TESS-SPOC-only fits agree",
                   qlp_fit=dict(tc1=fk["qlp"]["tc1"], tc2=fk["qlp"]["tc2"]), tessspoc_fit=dict(tc1=fk["tessspoc"]["tc1"], tc2=fk["tessspoc"]["tc2"])),
 sectors_with_data=list(ep["spans"].keys()), future_sectors_tess_point=[120, 121, 122, 126, 127, 128],
 alias_table=table,
 high_n_multiples_of_4={k: v["verdict"] for k, v in hde["high_n"].items()},
 s21_event=dict(predicted_tc=s21["tc_pred"], predicted_sigma=s21["sig_pred"], measured_tc=hde["tc_S21"],
                depth_fits_at_predicted_time=s21["fits"], known_event_depths_own_phot=s21["known_depths_own"],
                ingress_side_timing=ing["ingress_timing"],
                coverage="ingress and transit floor covered; data end 1883.971 (= predicted 4th contact 1883.962); last ~1.5 h on orbit-end scattered-light ramp; QLP flags from 1883.85",
                empirical_false_alarms=ing["fap_search"],
                centroid_inout_mpix=cen, pixel_level=pix),
 gap_epochs_own_photometry=gap,
 pipeline_epoch_fits=et,
 duration_density_check=hde["density"],
 ephemerides=hde["ephemeris"],
 upcoming_transits=upcoming,
 conclusion="n=4 (P=179.93 d) or n=8 (P=89.96 d), indistinguishable with existing TESS data, IF the partial S21 event at BTJD 1883.77 is a transit of the same object; n=3,6,7,9,10,11,12 (and n=16..48 multiples of 4) excluded; n=1,2,5 survive only if the S21 dip is rejected.",
)
json.dump(res, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
for n in (4, 8):
    print(n, upcoming[n][:6])
