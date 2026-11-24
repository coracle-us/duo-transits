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

# EPIC, TIC, group, alias divisors n (P = P_fit / n), events, prior period (with source), blind recovery, note.
# Groups are the literature verdicts of the AAS82225 revision (k2singles/revision2/PLAN.md, "Per-star verdicts").
# Blind: "yes" = the Table 1 period is among the pre-registered blind search's allowed periods or BLS
# (revision2/results/compare_table1.csv); "post hoc" = only the post-hoc linked-event test (link_calibrated.csv);
# "no" = neither.
ROWS = [
    (248811085, 349987241, "none", [1], "C14; S46, S72", "", "yes", ""),
    (246696483, 436612436, "none", [1, 2], "C13; S43, S70", "", "yes", "$P/2$ disfavoured"),
    (211390677, 366568535, "none", [1, 2], "C5, C18; S45", "", "yes", ""),
    (247967714, 268396693, "none", [1, 2, 3], "C13; S44, S71", "", "yes", ""),
    (212012030, 14394716, "none", [1, 2, 3, 4, 5], "C5, C18; S71", "", "yes", ""),
    (211503363, 20444116, "estimate", [1, 2], "C5; S45, S72", "$\\gtrsim58$$^g$", "no", "planet cand."),
    (211924561, 238606535, "differs", [1, 2], "C5, C16; S46, S72", "TCEs$^h$", "no", "TOI-5560.01"),
    (211489484, 6951741, "differs", [1], "C5, C18; S7, S46", "TCEs$^h$", "post hoc", "sec.\\ 0.194"),
    # S70 TESS-SPOC TCEs, each on one eclipse (revision2/dv_tce_s70.py; moved from "none" 2026-11-24)
    (246194159, 301328232, "differs", [1], "C12; S70, S92", "TCE$^h$", "no", "HD 220967; sec.\\ 0.867"),
    (246089124, 49603732, "differs", [1], "C12; S70, S92", "TCE$^h$", "yes", "sec.\\ 0.494"),
    (211411112, 437034781, "agrees", [1], "C5, C16; S46", "62.59484(3)$^a$", "yes", "Sanders 617"),
    (211892898, 184895645, "agrees", [1], "C5, C16, C18", "72.53(3)$^e$", "yes", "Gaia SB1"),
    (202126877, 86491361, "agrees", [1], "S43--45, S71--72", "21.78$^b$", "yes", "sec.\\ 0.33; deepening$^f$"),
    (202072917, 308255000, "agrees", [1], "C0; S6, S43, S45", "34.744$^b$, 34.7457$^c$", "yes", ""),
    (202071902, 57051340, "agrees", [1], "C0; S44, S71", "22.859(15)$^d$", "yes", ""),
]


def fmt(p, s):
    """Period with the uncertainty in the last two significant digits, e.g. 72.5089(15)."""
    if s <= 0:
        return f"{p:.5f}"
    nd = max(0, -int(f"{s:e}".split("e")[1]) + 1)
    return f"{p:.{nd}f}({round(s * 10 ** nd):d})"


GROUPS = {"none": "No prior orbital period", "estimate": "Prior estimate only", "differs": "Published periods differ",
          "agrees": "Published period agrees"}


def main():
    rows, tex = [], []
    last = None
    for epic, tic, grp, ns, events, prior, blind, note in ROWS:
        if grp != last:
            tex.append(f"\\cutinhead{{{GROUPS[grp]}}}")
            last = grp
        j = json.loads((EPH / f"{epic}.json").read_text())
        P, sP = j["P"], j["sigma_P"]
        if epic == 211411112:
            # K2 PDCSAP removes the C5 eclipse core, so the ephemeris fit used an ingress-only time. Re-timed on the
            # full eclipse in EVEREST: k2singles/revision/audit21/t0_211411112.py (stdout, 2026-10-13).
            P, sP = 62.594821, 0.000027
        if epic == 202126877:                     # eclipse-timing variations: quote the mean only
            ptxt = f"{P:.3f}"
        elif len(ns) > 2:
            ptxt = f"{fmt(P, sP)}/$n$, $n={ns[0]}$--${ns[-1]}$"
        else:
            ptxt = ", ".join(fmt(P / n, sP / n) for n in ns)
        rows.append(dict(EPIC=epic, TIC=tic, group=grp, events=events, P_fit=P, sigma_P=sP,
                         n_allowed=" ".join(map(str, ns)), periods=ptxt, T0_btjd=j["T0"], sigma_T0=j["sigma_T0"],
                         prior=prior, blind=blind, note=note))
        tex.append(f"{epic} & {events} & {ptxt} & {prior} & {blind} & {note} \\\\")   # TIC IDs are in table.csv
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "table.csv", index=False)
    (OUT / "table_rows.tex").write_text("\n".join(tex) + "\n")
    print("\n".join(tex))


if __name__ == "__main__":
    main()
