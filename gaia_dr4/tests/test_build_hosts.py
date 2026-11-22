import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_hosts as bh  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def test_every_source_row_is_kept_or_dropped():
    raw = bh.load_all()
    kept, dropped = bh.split(raw)
    assert len(kept) + len(dropped) == len(raw)
    n = raw.catalogue.value_counts()
    assert n["hawthorn"] == 85
    assert n["ctoi"] == len(pd.read_csv(ROOT / "singles/data/ctoi_singles.csv")) == 368
    assert n["salinas"] == 173
    assert n["k2"] == 164
    assert n["toi_single"] == 156


def test_kept_rows_are_clean():
    kept, dropped = bh.split(bh.load_all())
    assert not kept.duplicated(["catalogue", "cand_id"]).any()
    assert str(kept.tic.dtype) == "int64" and (kept.tic > 0).all()
    assert (dropped.reason != "").all()


def test_duotransits_carry_both_times():
    kept, _ = bh.split(bh.load_all())
    h = kept[kept.catalogue == "hawthorn"]
    assert h.t1_btjd.notna().all() and h.t2_btjd.notna().all()
    assert (h.t2_btjd > h.t1_btjd).all()


def test_split_drops_missing_tic_and_duplicates():
    df = pd.DataFrame(dict(tic=[1.0, float("nan"), 2.0, 2.0], catalogue=["a"] * 4,
                           cand_id=["x", "y", "z", "z"], t1_btjd=0.0, t2_btjd=0.0, period_d=0.0))
    kept, dropped = bh.split(df)
    assert kept.tic.tolist() == [1, 2]
    assert dropped.reason.tolist() == ["no TIC", "duplicate (catalogue, cand_id)"]


def test_hosts_one_row_per_tic():
    kept, _ = bh.split(bh.load_all())
    h = bh.hosts(kept)
    assert h.tic.is_unique and h.n_cand.sum() == len(kept)
    row = h[h.tic == 147035057].iloc[0]
    assert "ctoi" in row.catalogues.split(";")
