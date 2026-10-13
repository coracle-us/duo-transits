"""
Download the K2 HLSP light curves (EVEREST, K2SFF, K2SC) for the campaign of each of the 21 catalogued K2 events
that are not 'supported' in k2singles/revision/catalogued_event_check.csv. Files go to
k2singles/revision/audit21/data/ (lightkurve mastDownload tree). Failures are logged and skipped.
Output: k2singles/revision/audit21/data/hlsp_manifest.csv (EPIC, campaign, author, path or error).
Usage:  python k2singles/revision/audit21/fetch_hlsp.py
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = HERE / "data"
AUTHORS = ["EVEREST", "K2SFF", "K2SC"]


def targets():
    d = pd.read_csv(ROOT / "k2singles" / "revision" / "catalogued_event_check.csv")
    c = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    c["event_idx"] = c.groupby("EPIC").cumcount()
    bad = d[d.verdict != "supported"][["EPIC", "event_idx", "verdict"]]
    return bad.merge(c, on=["EPIC", "event_idx"])


def main():
    warnings.filterwarnings("ignore")
    import lightkurve as lk
    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in targets().iterrows():
        epic, camp = int(r.EPIC), int(r.campaign)
        try:
            sr = lk.search_lightcurve(f"EPIC {epic}", mission="K2", campaign=camp)
        except Exception as e:  # noqa: BLE001
            rows.append(dict(EPIC=epic, campaign=camp, author="search", path="", error=repr(e)[:200]))
            continue
        for a in AUTHORS:
            s = sr[[str(x) == a for x in sr.table["author"]]]
            if len(s) == 0:
                rows.append(dict(EPIC=epic, campaign=camp, author=a, path="", error="not in MAST"))
                continue
            for k in range(len(s)):
                try:
                    lc = s[k].download(download_dir=str(DATA), quality_bitmask="none")
                    rows.append(dict(EPIC=epic, campaign=camp, author=a, path=str(Path(lc.meta["FILENAME"])
                                                                                   .relative_to(ROOT)), error=""))
                except Exception as e:  # noqa: BLE001
                    rows.append(dict(EPIC=epic, campaign=camp, author=a, path="", error=repr(e)[:200]))
        print(epic, camp, [(x["author"], bool(x["path"])) for x in rows if x["EPIC"] == epic], flush=True)
    pd.DataFrame(rows).to_csv(DATA / "hlsp_manifest.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())
