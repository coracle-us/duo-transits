"""
Difference image (scripts/diff_image.analyse) at specific event times, e.g. a candidate repeat that the pipeline
run did not image (only the first 3 repeats get one).
Usage:  python toi_singles_pilot/scripts/s09_diff_at.py SET TOI T:SECTOR [T:SECTOR ...]
Output: printed result + toi_singles_pilot/plots/<TIC>_diffimg_s<SECTOR>.png; appended to
        toi_singles_pilot/results/<SET>/diffs_<TOI>.json
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402
import s01_run  # noqa: E402


def main():
    which, toi = sys.argv[1], sys.argv[2]
    row = next(r for r in s01_run.targets(which) if r["TOI"] == toi)
    D = row["dur_h"] / 24
    gaia = common.di.gaia_field(row["ra"], row["dec"])
    tgt = gaia.iloc[0]
    f = common.PILOT / "results" / which / f"diffs_{toi}.json"
    out = json.loads(f.read_text()) if f.exists() else []
    for arg in sys.argv[3:]:
        T, s = arg.split(":")
        r = common.di.analyse(int(row["TIC"]), 1e6, float(T), D, int(s), tgt.ra, tgt.dec, gaia,
                              out_dir=common.PILOT / "plots")
        r = r or {}
        r["T"] = float(T)
        out.append(r)
        print({k: r.get(k) for k in ("T", "sector", "n_transits", "offset_from_target_arcsec", "offset_err_px",
                                     "offset_sigma", "snr_peak", "crowded_reference", "nearest_gaia_to_source")})
    f.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
