import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core import binarity_flags, period_match, phase_match, rv_variable  # noqa: E402


def test_rv_variable_cut():
    assert rv_variable(0.001, 5.0, 12) is True
    assert rv_variable(0.0, 5.0, 12) is True          # p = 0 is the most variable case
    assert rv_variable(0.02, 5.0, 12) is False
    assert rv_variable(0.001, 3.9, 12) is False
    assert rv_variable(0.001, 5.0, 9) is False
    assert rv_variable(None, 5.0, 12) is False
    assert rv_variable(float("nan"), 5.0, 12) is False


def test_phase_match_exact_alias():
    # dt = 20 cycles of 37.31497 d, Gaia period known to 0.001 d
    r = phase_match(20 * 37.31497, 0.002, 37.31497, 0.001)
    assert r["n"] == 20
    assert abs(r["frac"]) < 1e-9
    assert r["verdict"] == "match"
    assert r["window"] < 0.5


def test_phase_match_off_alias():
    # 20.4 cycles: clearly not an integer, and the period is precise enough to say so
    r = phase_match(20.4 * 37.31497, 0.002, 37.31497, 0.001)
    assert r["n"] == 20
    assert math.isclose(r["frac"], 0.4, abs_tol=1e-9)
    assert r["verdict"] == "no_match"


def test_phase_match_uninformative_when_period_is_loose():
    # dphi = dt * p_err / p^2 = 746.3 * 0.5 / 37.31497^2 = 0.268 cycles; 2*3*dphi > 0.5
    r = phase_match(20 * 37.31497, 0.002, 37.31497, 0.5)
    assert r["verdict"] == "uninformative"
    assert r["window"] > 0.5


def test_phase_match_window_is_fraction_of_phase_accepted():
    r = phase_match(1000.0, 0.0, 50.0, 0.01)      # dphi = 1000*0.01/2500 = 0.004
    assert math.isclose(r["dphi"], 0.004, rel_tol=1e-9)
    assert math.isclose(r["window"], 0.024, rel_tol=1e-9)   # 2 * 3 * dphi


def test_phase_match_missing_error_is_uninformative():
    assert phase_match(1000.0, 0.0, 50.0, float("nan"))["verdict"] == "uninformative"
    assert phase_match(1000.0, 0.0, 50.0, None)["verdict"] == "uninformative"


def test_period_match_direct_and_harmonics():
    assert period_match(37.31497, 37.3152, 0.0005)["verdict"] == "match"
    assert period_match(37.31497, 37.3152, 0.0005)["harmonic"] == 1.0
    r = period_match(37.31497, 74.630, 0.001)     # Gaia found twice the period
    assert r["verdict"] == "match" and r["harmonic"] == 0.5
    r = period_match(37.31497, 18.6575, 0.0005)   # Gaia found half the period
    assert r["verdict"] == "match" and r["harmonic"] == 2.0
    assert period_match(37.31497, 40.0, 0.01)["verdict"] == "no_match"


def test_period_match_relative_floor():
    # error quoted as 0, off by 0.05%: inside the 0.1% floor
    assert period_match(100.0, 100.05, 0.0)["verdict"] == "match"
    assert period_match(100.0, 100.2, 0.0)["verdict"] == "no_match"


def test_binarity_flags():
    f = binarity_flags(dict(ruwe=3.79, non_single_star=0, rv_chisq_pvalue=None,
                            rv_renormalised_gof=None, rv_nb_transits=None))
    assert f == dict(ruwe_flag=True, nss_flag=False, rv_flag=False, any_flag=True)
    f = binarity_flags(dict(ruwe=1.0, non_single_star=2, rv_chisq_pvalue=0.5,
                            rv_renormalised_gof=1.0, rv_nb_transits=20))
    assert f["nss_flag"] is True and f["ruwe_flag"] is False and f["any_flag"] is True
    f = binarity_flags(dict(ruwe=float("nan"), non_single_star=float("nan")))
    assert f["any_flag"] is False


def test_exact_ids_keeps_all_19_digits():
    import numpy as np
    import pandas as pd
    from core import exact_ids
    out = exact_ids(["6427464123776727123", " 2947647054662200193 ", "", None, float("nan"),
                     np.int64(1669046625055549313), 6.427464123776727e18])
    assert out[0] == 6427464123776727123
    assert out[1] == 2947647054662200193
    assert out[5] == 1669046625055549313
    assert pd.isna(out[2]) and pd.isna(out[3]) and pd.isna(out[4])
    assert pd.isna(out[6])            # a float ID has already lost digits: refuse it
    assert str(out.dtype) == "Int64"
