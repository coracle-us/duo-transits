"""
Decision functions for the Gaia cross-match. No network, no file access.

Thresholds are arguments; the values used in a run come from prereg_dr4.yaml.
"""
import math


def _missing(v):
    return v is None or (isinstance(v, float) and math.isnan(v)) or v != v


def rv_variable(p, gof, n):
    """Gaia DR3 RV variability cut (Katz et al. 2023): p <= 0.01, gof > 4, n >= 10.
    A missing value means 'not tested' and gives False. p = 0 is the most variable case."""
    if any(_missing(v) for v in (p, gof, n)):
        return False
    return bool(p <= 0.01 and gof > 4 and n >= 10)


def phase_match(dt, dt_err, p, p_err, nsig=3.0, max_window=0.5):
    """Is the gap dt between two transits a whole number of Gaia periods p?

    dphi   = 1-sigma uncertainty of dt/p in cycles
    window = fraction of all phases that would be accepted (2 * nsig * dphi); it is also the
             chance that an unrelated period 'matches'
    verdict: 'uninformative' if window > max_window or an error is missing,
             else 'match' if |frac| <= nsig * dphi, else 'no_match'.
    """
    out = dict(n=None, frac=None, dphi=None, window=None, verdict="uninformative")
    if any(_missing(v) for v in (dt, p, p_err)) or p <= 0:
        return out
    dt_err = 0.0 if _missing(dt_err) else dt_err
    cycles = dt / p
    n = round(cycles)
    frac = cycles - n
    dphi = math.hypot(dt_err / p, dt * p_err / p ** 2)
    window = 2 * nsig * dphi
    out.update(n=int(n), frac=frac, dphi=dphi, window=window)
    if window > max_window or n < 1:
        return out
    out["verdict"] = "match" if abs(frac) <= nsig * dphi else "no_match"
    return out


def period_match(p_hyp, p, p_err, nsig=3.0, rel_floor=1e-3, harmonics=(0.5, 1.0, 2.0)):
    """Does the Gaia period p agree with a hypothesised period p_hyp, directly or at a harmonic?

    harmonic h means p * h is compared with p_hyp (h = 0.5: Gaia found twice the period).
    Tolerance: max(nsig * h * p_err, rel_floor * p_hyp). The floor guards against Gaia period
    errors that are quoted too small.
    """
    out = dict(verdict="no_match", harmonic=None, ratio=None)
    if any(_missing(v) for v in (p_hyp, p)) or p <= 0 or p_hyp <= 0:
        return out
    out["ratio"] = p / p_hyp
    p_err = 0.0 if _missing(p_err) else p_err
    for h in sorted(harmonics, key=lambda x: abs(x - 1.0)):
        tol = max(nsig * h * p_err, rel_floor * p_hyp)
        if abs(p * h - p_hyp) <= tol:
            out.update(verdict="match", harmonic=float(h))
            return out
    return out


def binarity_flags(row, ruwe_max=1.4):
    """Flags from one gaia_source row (a mapping)."""
    ruwe = row.get("ruwe")
    nss = row.get("non_single_star")
    f = dict(
        ruwe_flag=bool(not _missing(ruwe) and ruwe > ruwe_max),
        nss_flag=bool(not _missing(nss) and nss > 0),
        rv_flag=rv_variable(row.get("rv_chisq_pvalue"), row.get("rv_renormalised_gof"),
                            row.get("rv_nb_transits")),
    )
    f["any_flag"] = f["ruwe_flag"] or f["nss_flag"] or f["rv_flag"]
    return f


def exact_ids(values):
    """A pandas Int64 column of 19-digit Gaia source IDs from strings or integers, with no pass
    through float64 (which cannot hold them exactly). Blank, non-numeric and float inputs give <NA>."""
    import pandas as pd

    def one(v):
        if isinstance(v, (bytes, bytearray)):
            v = v.decode()
        if isinstance(v, str):
            v = v.strip()
            return int(v) if v.isdigit() else pd.NA
        if isinstance(v, bool) or v is None or v is pd.NA:
            return pd.NA
        if hasattr(v, "dtype") and getattr(v.dtype, "kind", "") in "iu":
            return int(v)
        return int(v) if isinstance(v, int) else pd.NA

    return pd.array([one(v) for v in values], dtype="Int64")
