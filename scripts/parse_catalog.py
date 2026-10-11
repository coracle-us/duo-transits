"""
Parse Hawthorn et al. (2024) appendix Tables A1-A3 from the arXiv LaTeX source into
data/catalog/candidates.csv, then diff against the unverified web transcription in
data/transcribed/ (LaTeX wins; differences written to data/catalog/transcription_diff.csv).

Source files (data/raw/paper_src/, per appendix.tex):
    table_1.tex  -> Table A1  host properties
    table_4.tex  -> Table A2  transit events
    table_5.tex  -> Table A3  planet properties / aliases

Usage:  python scripts/parse_catalog.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "raw" / "paper_src"
CAT = ROOT / "data" / "catalog"
TRANS = ROOT / "data" / "transcribed"

NUM = r"[-+]?\d+(?:\.\d+)?"


def table_rows(path):
    """Return the data rows of a longtable as lists of raw cell strings."""
    text = path.read_text(encoding="utf-8")
    text = "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("%"))
    body = text.split(r"\hline", 1)[1]
    body = body.split(r"\bottomrule", 1)[0]
    rows = []
    for chunk in body.split(r"\\"):
        chunk = chunk.strip()
        if not chunk:
            continue
        cells = [c.strip() for c in chunk.split("&")]
        rows.append(cells)
    return rows


def tic_of(cell):
    m = re.fullmatch(r"\\href\{[^}]*id=(\d+)\}\{(\d+)\}", cell)
    assert m and m.group(1) == m.group(2), f"bad TIC cell: {cell!r}"
    return int(m.group(2))


def val_err(cell):
    """'$ 451.71^{+20.73}_{-15.90} $' / '$ 27.0 \\pm 9.0 $' / '1.5' / '--' / ''
    -> (value, err_plus, err_minus); err_minus returned as a positive number."""
    s = cell.replace("$", "").replace(" ", "")
    if s in ("", "--", "-"):
        return (np.nan, np.nan, np.nan)
    m = re.fullmatch(rf"({NUM})\^\{{\+({NUM})\}}_\{{-({NUM})\}}", s)
    if m:
        return tuple(float(x) for x in m.groups())
    m = re.fullmatch(rf"({NUM})\\pm({NUM})", s)
    if m:
        return (float(m.group(1)), float(m.group(2)), float(m.group(2)))
    m = re.fullmatch(NUM, s)
    if m:
        return (float(s), np.nan, np.nan)
    raise ValueError(f"unparsed value cell: {cell!r}")


def num(cell):
    s = cell.replace("$", "").strip()
    if s in ("", "--", "-"):
        return np.nan
    return float(s)


def parse_a1():
    out = []
    for c in table_rows(SRC / "table_1.tex"):
        assert len(c) == 10, c
        d, dp, dm = val_err(c[8])
        out.append(dict(
            ID=int(c[0]), TIC=tic_of(c[1]),
            TOI=None if c[2] in ("--", "") else c[2],
            RA=num(c[3]), Dec=num(c[4]), Tmag=num(c[5]), Rstar=num(c[6]), Teff=num(c[7]),
            Dist_pc=d, Dist_errp=dp, Dist_errm=dm, Nsec=int(c[9]),
        ))
    return pd.DataFrame(out)


def parse_a2():
    out = []
    for c in table_rows(SRC / "table_4.tex"):
        assert len(c) == 6, c
        t = []
        for cell in c[2:4]:
            m = re.fullmatch(rf"({NUM})\s*\[(\d+)\]", cell)
            assert m, cell
            t.append((float(m.group(1)), int(m.group(2))))
        deps = re.findall(rf"({NUM})\s*\[(\d+)\]", c[4])
        assert len(deps) == 2, c[4]
        assert [int(s) for _, s in deps] == [t[0][1], t[1][1]], f"depth sectors mismatch: {c}"
        out.append(dict(
            ID=int(c[0]), TIC=tic_of(c[1]),
            Tc1=t[0][0], S1=t[0][1], Tc2=t[1][0], S2=t[1][1],
            Depth1_ppt=float(deps[0][0]), Depth2_ppt=float(deps[1][0]), Dur_h=num(c[5]),
        ))
    return pd.DataFrame(out)


def parse_a3():
    out = []
    for c in table_rows(SRC / "table_5.tex"):
        assert len(c) == 10, c
        rp, rpp, rpm = val_err(c[2])
        pmin, pmax = (float(x) for x in c[3].split(","))
        pm, pmp, pmm = val_err(c[5])
        b, bp, bm = val_err(c[6])
        out.append(dict(
            ID=int(c[0]), TIC=tic_of(c[1]),
            Rp_Re=rp, Rp_errp=rpp, Rp_errm=rpm, Pmin=pmin, Pmax=pmax, N_alias=int(c[4]),
            P_marg=pm, P_marg_errp=pmp, P_marg_errm=pmm, b=b, b_errp=bp, b_errm=bm,
            Gaia_NSS=int(c[7]), Gaia_RV_amp_kms=num(c[8]), SNR_total=num(c[9]),
        ))
    return pd.DataFrame(out)


def parse_monotools_draft():
    """latex_tab.tex: a draft MonoTools-fit table that is NOT included in the published paper
    (commented out of appendix.tex) but carries Tc uncertainties, which Table A2 lacks.
    Returns TIC, t0, t0_err, t1, t1_err, tD_d (duration, days). Errors are the mean of +/-."""
    out = []
    for c in table_rows(SRC / "latex_tab.tex"):
        if len(c) != 13 or not c[0].strip().isdigit():
            continue
        t0, t0p, t0m = val_err(c[3])
        t1, t1p, t1m = val_err(c[4])
        td, _, _ = val_err(c[6])
        out.append(dict(TIC=int(c[0]), mt_t0=t0, mt_t0_err=(t0p + t0m) / 2,
                        mt_t1=t1, mt_t1_err=(t1p + t1m) / 2, mt_tD_d=td))
    return pd.DataFrame(out)


def diff_transcribed(cat):
    """Compare against the web transcription. Returns a DataFrame of disagreements."""
    diffs = []

    def cmp(tr, mapping, label):
        tr = tr.set_index("ID")
        c = cat.set_index("ID")
        for i in c.index:
            if i not in tr.index:
                diffs.append(dict(table=label, ID=i, TIC=c.at[i, "TIC"], column="(row)",
                                  latex="present", transcribed="missing"))
                continue
            for tcol, ccol in mapping.items():
                lv, tv = c.at[i, ccol], tr.at[i, tcol]
                if isinstance(lv, float) or isinstance(tv, float):
                    try:
                        same = (pd.isna(lv) and pd.isna(tv)) or abs(float(lv) - float(tv)) < 1e-9
                    except (TypeError, ValueError):
                        same = False
                else:
                    same = str(lv) == str(tv)
                if not same:
                    diffs.append(dict(table=label, ID=i, TIC=c.at[i, "TIC"], column=tcol,
                                      latex=lv, transcribed=tv))

    a1 = pd.read_csv(TRANS / "hosts_tableA1.csv", dtype=str)
    # normalise the transcription's formats to comparable values
    a1["TIC"] = a1["TIC"].astype(int)
    for col in ("RA", "Dec", "Tmag", "Rstar", "Teff"):
        a1[col] = pd.to_numeric(a1[col], errors="coerce")
    a1["Nsec"] = pd.to_numeric(a1["Nsec"], errors="coerce")
    a1["TOI"] = a1["TOI"].where(~a1["TOI"].isin(["_", "--", "-"]), None)
    a1["ID"] = a1["ID"].astype(int)

    def dist_parts(s):
        s = str(s)
        m = re.fullmatch(rf"\s*({NUM})\s*\+({NUM})\s*-({NUM})\s*", s)
        return [float(x) for x in m.groups()] if m else [np.nan] * 3 if s in ("nan", "", "_") else [s, s, s]

    parts = a1["Distance"].map(dist_parts)
    a1["Dist_pc"], a1["Dist_errp"], a1["Dist_errm"] = zip(*parts)
    cmp(a1, {"TIC": "TIC", "TOI": "TOI", "RA": "RA", "Dec": "Dec", "Tmag": "Tmag", "Rstar": "Rstar",
             "Teff": "Teff", "Dist_pc": "Dist_pc", "Dist_errp": "Dist_errp", "Dist_errm": "Dist_errm",
             "Nsec": "Nsec"}, "A1")

    a2 = pd.read_csv(TRANS / "transits_tableA2_partial.csv")
    cmp(a2, {c: c for c in ("TIC", "Tc1", "S1", "Tc2", "S2", "Depth1_ppt", "Depth2_ppt", "Dur_h")}, "A2")
    return pd.DataFrame(diffs)


def main():
    a1, a2, a3 = parse_a1(), parse_a2(), parse_a3()
    for name, t in (("A1", a1), ("A2", a2), ("A3", a3)):
        assert len(t) == 85, (name, len(t))
        assert list(t["ID"]) == list(range(1, 86)), name
        assert t["TIC"].is_unique, name
    assert (a1["TIC"].values == a2["TIC"].values).all() and (a1["TIC"].values == a3["TIC"].values).all()

    cat = a1.merge(a2, on=["ID", "TIC"]).merge(a3, on=["ID", "TIC"])
    cat["dT"] = cat["Tc2"] - cat["Tc1"]
    # Sanity: P_max should equal Tc2 - Tc1 (n = 1) and P_min = dT / n for integer n
    cat["n_max_from_Pmin"] = (cat["dT"] / cat["Pmin"]).round(3)

    # Tc uncertainties from the unpublished MonoTools draft table, used only where its
    # epochs agree with the published Table A2 (within 0.01 d); otherwise left NaN.
    mt = parse_monotools_draft()
    cat = cat.merge(mt, on="TIC", how="left")
    ok = ((cat.mt_t0 - cat.Tc1).abs() < 0.01) & ((cat.mt_t1 - cat.Tc2).abs() < 0.01)
    cat["Tc1_err"] = cat.mt_t0_err.where(ok)
    cat["Tc2_err"] = cat.mt_t1_err.where(ok)
    cat["Dur_h_monotools"] = (cat.mt_tD_d * 24).where(ok).round(3)
    print(f"Tc errors from MonoTools draft table: {ok.sum()} of {len(cat)} candidates matched")
    cat = cat.drop(columns=[c for c in cat.columns if c.startswith("mt_")])
    CAT.mkdir(parents=True, exist_ok=True)
    cat.to_csv(CAT / "candidates.csv", index=False)
    print(f"wrote {CAT / 'candidates.csv'}  ({len(cat)} rows, {len(cat.columns)} cols)")

    d = diff_transcribed(cat)
    d.to_csv(CAT / "transcription_diff.csv", index=False)
    print(f"{len(d)} disagreements with data/transcribed/ -> {CAT / 'transcription_diff.csv'}")
    if len(d):
        print(d.to_string(index=False))


if __name__ == "__main__":
    main()
