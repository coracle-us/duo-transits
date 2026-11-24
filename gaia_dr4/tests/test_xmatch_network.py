"""Checks on the DR3 run's outputs against values already recorded in this repo.
Run after `python gaia_dr4/xmatch.py --schema gaiadr3`:  pytest -m network gaia_dr4/tests"""
from pathlib import Path

import pandas as pd
import pytest

OUT = Path(__file__).resolve().parents[1] / "results_dr3"
pytestmark = [pytest.mark.network,
              pytest.mark.skipif(not (OUT / "report.csv").exists(), reason="no DR3 run yet")]


@pytest.fixture(scope="module")
def ids():
    return pd.read_csv(OUT / "ids.csv", dtype={"source_id": "Int64"}).set_index("tic")


@pytest.fixture(scope="module")
def rep():
    return pd.read_csv(OUT / "report.csv", dtype={"source_id": "Int64"})


def test_known_source_ids(ids):
    # singles/TODO_gaia_dr4.md
    assert int(ids.loc[147035057, "source_id"]) == 2947647054662200192
    assert int(ids.loc[159614322, "source_id"]) == 1669046625055549312


def test_known_ruwe(rep):
    r = rep[rep.tic == 159614322].iloc[0]
    assert abs(r.ruwe - 1.51) < 0.01
    r = rep[rep.tic == 230133616].iloc[0]      # CLAUDE.md, Salinas round 3
    assert abs(r.ruwe - 3.79) < 0.01 and bool(r.ruwe_flag)


def test_known_gaia_binaries_have_matching_orbits(rep):
    # singles/scripts/gaia_flags.py: both TESS periods match Gaia DR3 SB1 orbits
    for tic in (231969683, 230130556):
        r = rep[(rep.tic == tic) & (rep.catalogue == "own_ctoi")]
        assert len(r) and (r.n_nss >= 1).all()
        assert (r.period_verdict == "match").any()


def test_every_host_is_accounted_for(ids):
    hosts = pd.read_csv(OUT.parent / "hosts_frozen.csv")
    assert set(hosts.tic) == set(ids.index)
    assert (ids.id_method != "").all()
