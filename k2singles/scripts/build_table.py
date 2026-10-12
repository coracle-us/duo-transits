"""
Table for the K2-singles RNAAS note. Periods come from the dedicated eclipse-time fits
(k2singles/results/ephemerides/<EPIC>.json, scripts/eclipse_times.py). Allowed alias numbers come
from the blind re-derivations (independent_table/, independent_table2/) and the second Fable review.
Each row's note and prior-work flag cite their source.
Outputs: k2singles/rnaas/table.csv and k2singles/rnaas/table_rows.tex (pasted into rnaas.tex).
Usage:  python k2singles/scripts/build_table.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EPH = ROOT / "k2singles" / "results" / "ephemerides"
OUT = ROOT / "k2singles" / "rnaas"

# EPIC, TIC, group, alias divisors n (P = P_fit / n), events, prior period (with source), note
ROWS = [
    (211924561, 238606535, "planet", [1, 2], "C5, C16; S46, S72",
     "SPOC TOI: $\\leq702.47$\\,d", "TOI-5560.01; V-shaped (see text)"),
    (211503363, 20444116, "planet", [1, 2], "C5; S45, S72",
     "none (single; Kov\\'acs 2020)", "$R_2\\approx2.1$--$2.4\\,R_{\\rm J}$; evolved host"),
    (211411112, 437034781, "validation", [1], "C5, C16; S46",
     "62.59484(3)\\,d$^a$", "Sanders 617 (M67)"),
    (246194159, 301328232, "eb", [1], "C12; S70, S92",
     "none", "HD 220967; secondary at phase 0.867"),
    (211489484, 6951741, "eb", [1], "C5, C18; S7, S46",
     "TCE 106.94\\,d (excluded)", "secondary at phase 0.194"),
    (211892898, 184895645, "eb", [1], "C5, C16, C18",
     "Gaia SB1 72.53(3)\\,d", "red giant; eclipse at Gaia primary conjunction"),
    (248811085, 349987241, "eb", [1], "C14; S46, S72", "none", ""),
    (246089124, 49603732, "eb", [1], "C12; S70, S92", "none", "K2 secondary at phase 0.494"),
    (246696483, 436612436, "eb", [1], "C13; S43, S70", "none", "$P/2$ excluded by one S98 epoch"),
    (211390677, 366568535, "eb", [1, 2], "C5, C18; S45", "none", "red-clump giant; $M_2$ stellar"),
    (247967714, 268396693, "eb", [1, 2, 3], "C13; S44, S71", "none", ""),
    (212012030, 14394716, "eb", [1, 2, 3, 4, 5], "C5, C18; S71", "none", "grazing; $2P$ possible"),
    (202126877, 86491361, "known", [1], "C0; S43--45, S71--72",
     "21.78\\,d$^b$", "eclipses deepen: 6\\% (2014) to 32\\% (2025)"),
    (202072917, 308255000, "known", [1], "C0; S6, S43, S45",
     "34.744\\,d$^{b,c}$", "VSX 3.842\\,d is wrong"),
    (202071902, 57051340, "known", [1], "C0; S44, S71",
     "22.859(15)\\,d$^d$", "VSX 8.697\\,d is wrong"),
]


def fmt(p, s):
    """Period with the uncertainty in the last two significant digits, e.g. 72.5089(15)."""
    if s <= 0:
        return f"{p:.5f}"
    nd = max(0, -int(f"{s:e}".split("e")[1]) + 1)
    return f"{p:.{nd}f}({round(s * 10 ** nd):d})"


def main():
    rows, tex = [], []
    for epic, tic, grp, ns, events, prior, note in ROWS:
        j = json.loads((EPH / f"{epic}.json").read_text())
        P, sP = j["P"], j["sigma_P"]
        if epic == 246194159:                     # joint primary + secondary fit (see triage notes)
            P, sP = 73.245178, 0.000007
        if epic == 202126877:                     # eclipse-timing variations: quote the mean only
            ptxt = f"{P:.3f}"
        elif len(ns) > 2:
            ptxt = f"{fmt(P, sP)}/$n$, $n={ns[0]}$--${ns[-1]}$"
        else:
            ptxt = " or ".join(fmt(P / n, sP / n) for n in ns)
        rows.append(dict(EPIC=epic, TIC=tic, group=grp, events=events, P_fit=P, sigma_P=sP,
                         n_allowed=" ".join(map(str, ns)), periods=ptxt, T0_btjd=j["T0"], sigma_T0=j["sigma_T0"],
                         prior=prior, note=note))
        tex.append(f"{epic} & {tic} & {events} & {ptxt} & {prior} & {note} \\\\")
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "table.csv", index=False)
    (OUT / "table_rows.tex").write_text("\n".join(tex) + "\n")
    print("\n".join(tex))


if __name__ == "__main__":
    main()
