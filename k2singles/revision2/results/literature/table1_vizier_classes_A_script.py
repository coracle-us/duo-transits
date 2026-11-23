import pandas as pd
d = pd.read_csv("data/literature/table1_vizier_rows.csv")
pairs = d[["EPIC", "catalog"]].drop_duplicates()
C = {  # (EPIC, catalog): (class, period_d, kind, note)
 (202071902, "J/A+A/579/A19/table3"): ("LISTED", "", "", "Armstrong+2015 K2 variable catalogue, Type=EB, no period"),
 (202071902, "J/ApJS/264/12/param"): ("LISTED", "", "", "LAMOST activity survey"),
 (202071902, "J/ApJS/278/46/table1"): ("PERIOD", "22.89", "orbital", "SB2 orbit, e_Per 0.395"),
 (202071902, "J/MNRAS/452/3561/table2"): ("LISTED", "", "", "LaCourse+2015 single-eclipse list, no period"),
 (202071902, "J/MNRAS/527/521/tableb1"): ("LISTED", "", "", "SB2 candidate list"),
 (202071902, "J/MNRAS/527/521/tablec1"): ("ANALYSIS", "", "", "dynamical mass ratio"),
 (202071902, "J/MNRAS/527/521/tablec2"): ("PERIOD", "22.85881;1.207985", "orbital", "Pkep 22.85881(153) Keplerian fit; Psin 1.208 sinusoid fit"),
 (202071902, "J/other/RAA/21.292/table2"): ("LISTED", "", "", ""),
 (202071902, "J/other/RAA/21.292/table3"): ("PERIOD", "8.6973", "eclipsing", "Type EA, Plc from VSX (settled 2026-11-22)"),
 (202071902, "J/other/RAA/21.292/tablea2"): ("LISTED", "", "", ""),
 (202071902, "J/other/RAA/21.292/tablea4"): ("LISTED", "", "", ""),
 (202072917, "B/vsx/vsx"): ("PERIOD", "3.842372", "eclipsing", "VSX Type E"),
 (202072917, "J/A+A/691/A242/obaf-eb1"): ("PERIOD", "23.25528;34.74438", "eclipsing", "EB subset; SPOC y4 and QLP y4 periods"),
 (202072917, "J/A+A/691/A242/tessobaf"): ("PERIOD", "23.25528;34.74438", "eclipsing", "same row as obaf-eb1"),
 (202072917, "J/AJ/170/287/table13"): ("PERIOD", "34.74574;34.74165;1007.648;15.17308;14.19501", "TCE", "SPOC TCEs, ExoMiner FP"),
 (202072917, "J/AJ/170/287/table5"): ("PERIOD", "5.578824", "TCE", "labelled TCE"),
 (202072917, "J/ApJS/264/17/table1"): ("LISTED", "", "", "LAMOST spectra of TESS objects"),
 (202072917, "J/ApJS/264/17/table3"): ("LISTED", "", "", "LAMOST; VarRV=Y"),
 (202072917, "J/MNRAS/452/3561/table2"): ("LISTED", "", "", "LaCourse+2015, no period"),
 (202126877, "B/vsx/vsx"): ("PERIOD", "3.877092", "eclipsing", "VSX Type E"),
 (202126877, "J/A+A/691/A242/obaf-eb1"): ("PERIOD", "21.78076;21.77555", "eclipsing", "EB subset; SPOC y4, QLP y4"),
 (202126877, "J/A+A/691/A242/tessobaf"): ("PERIOD", "21.78076;21.77555", "eclipsing", "same row as obaf-eb1"),
 (202126877, "J/ApJS/259/63/table1"): ("LISTED", "", "", "LAMOST CP-star list"),
 (202126877, "J/MNRAS/452/3561/table1"): ("LISTED", "", "", "LaCourse+2015, no period"),
 (211411112, "J/AJ/161/190/table1"): ("ANALYSIS", "", "", "RVs"),
 (211411112, "J/AJ/161/190/table3"): ("PERIOD", "62.593", "orbital", "SB2 orbit, Geller+2021"),
 (211411112, "J/AJ/162/184/table1"): ("ANALYSIS", "", "", "APOGEE SB2, RV amplitudes, no period"),
 (211411112, "J/MNRAS/459/1060/table2"): ("ANALYSIS", "21.175186", "rotation", "Gonzalez 2016 periodogram period, class BM, not orbital"),
 (211489484, "J/AJ/170/287/table13"): ("PERIOD", "106.9361;0.240513", "TCE", "SPOC TCEs, ExoMiner FP"),
 (211489484, "J/A+A/672/A176/table1"): ("LISTED", "", "", "contact-binary sample row, no period"),
 (211489484, "J/ApJS/244/11/stars"): ("LISTED", "", "", "Kruse+2019 star list"),
 (211503363, "J/MNRAS/495/1252/table2"): ("ANALYSIS", "14.0529", "rotation", "Prot"),
 (211503363, "J/MNRAS/495/1252/table3"): ("LISTED", "", "", ""),
 (211503363, "J/MNRAS/495/1252/table5"): ("ANALYSIS", "12.470412", "rotation", "Prot, flare table"),
 (211892898, "I/357/tbosb1"): ("PERIOD", "72.53358", "orbital", "Gaia DR3 SB1"),
 (211892898, "J/AJ/154/224/table2"): ("ANALYSIS", "8.51", "rotation", "Prot; plx 7.8 is the cluster value, Gaia DR3 plx 0.88 (background giant)"),
 (211892898, "J/ApJ/842/83/table3"): ("ANALYSIS", "", "", "Notes=eclipse, Bin?=cand; no period"),
 (211892898, "J/ApJ/842/83/table2"): ("LISTED", "", "", "photometric binary flag"),
 (211924561, "J/A+A/691/A242/tessobaf"): ("LISTED", "0.119355", "", "whole OBAF catalogue, not the EB subset; scoreSPOCy4=0; not an EB detection"),
 (211924561, "J/AJ/170/287/table13"): ("PERIOD", "31.26866;58.29026", "TCE", "SPOC S42-46 TCEs (MES 7.4, 11.5), ExoMiner FP"),
 (211924561, "J/ApJS/244/43/table1"): ("LISTED", "100.0", "placeholder", "Villanova placeholder copied (settled 2026-11-22)"),
 (211924561, "J/A+A/672/A176/table1"): ("LISTED", "", "", "contact-binary sample row"),
 (212012030, "J/AJ/156/241/table4"): ("LISTED", "1.257543;6.794257;0.721049", "", "ATLAS variable, Class=dubious; periodogram peaks, not orbital"),
}
rows = []
for r in pairs.itertuples():
    c = C.get((r.EPIC, r.catalog), ("LISTED", "", "", "survey/membership/stellar-parameter table, no period"))
    rows.append(dict(EPIC=r.EPIC, catalog=r.catalog, **dict(zip(["class", "period_d", "period_kind", "note"], c))))
o = pd.DataFrame(rows)
missing = set(C) - set(zip(pairs.EPIC, pairs.catalog))
assert not missing, missing
o.to_csv("data/literature/table1_vizier_classes_A.csv", index=False)
print(len(o), o["class"].value_counts().to_dict())
