"""
Known multiple-system list for cross-matching: TIC numbers harvested from the arXiv LaTeX
source of TESS triple / quadruple / ETV papers (see notes/eb_literature.md).

Every TIC ID in the source is taken: 'TIC 123', 'TIC~123', and table rows starting with a
6-10 digit number followed by '&'. Over-inclusive on purpose (a paper may mention a TIC that
is not a triple, e.g. a comparison star); use the flag as "already discussed in the
multiples literature", not as "confirmed triple".

Downloads go to sidequests/eb_timing/cache/arxiv/ (gitignored).
Output: sidequests/eb_timing/data/known_multiples.csv (TIC, papers)
"""
import io
import re
import tarfile
import time
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parents[1]
CACHE = HERE / "cache" / "arxiv"
OUT = HERE / "data" / "known_multiples.csv"

PAPERS = {
    "2402.01486": "Mitnyan+2024 NCVZ ETV triples (S14-60)",
    "2604.04004": "Borkovits+2026 NCVZ ETV revisit (S14-86)",
    "2403.07694": "Marcadon&Prsa 2024 CVZ FFI timing",
    "2503.20900": "Esmer&Daylan 2025 circumbinary ETV",
    "2309.14200": "Kostov+2024 101 quadruple candidates",
    "2510.27470": "Kostov+2025 52 quadruple candidates",
    "2510.04565": "Borkovits+2025 ten triply eclipsing",
    "2403.12041": "Rappaport+2024 seven triply eclipsing",
    "2301.13810": "Rappaport+2023 nine triply eclipsing",
    "2309.01639": "Czavalinga+2023 four triply eclipsing (Gaia)",
    "2212.02903": "Czavalinga+2023 Gaia DR3 triple candidates",
    "2002.06902": "Borkovits+2020 TIC 167692429/220397947",
    "2111.10074": "Borkovits+2022 northern triply eclipsing",
    "2006.10449": "Borkovits+2020 TIC 209409435",
    "2009.03378": "Mitnyan+2020 TIC 278825952",
    "2208.05789": "Powell+2022 TIC 114936199",
    "2410.01711": "Kostov+2024 TIC 290061484",
    "2504.12239": "Powell+2025 two quadruples",
    "2605.20080": "Powell+2026 TIC 295741342",
    "2609.25332": "Rappaport+2026 TIC 366240660",
    "2608.22935": "Zasche+2026 six doubly eclipsing quadruples",
    "2403.16657": "Zasche+2024 six eccentric EBs with third body",
    "2202.06964": "Eisner+2022 TIC 470710327",
    "2006.08979": "Rowden+2020 TIC 278956474",
    "2303.14275": "Kostov+2023 TIC 219006972",
    "2105.12586": "Kostov+2021 TIC 454140642",
    "2209.01142": "Hajdu+2022 CoRoT-TESS LTTE",
}
PAT = [re.compile(r"TIC[\s~\-_:]*(?:\,)?\s*(\d{5,10})"), re.compile(r"^\s*(\d{6,10})\s*&", re.M)]


def source_text(aid):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{aid}.bin"
    if not f.exists():
        r = requests.get(f"https://arxiv.org/e-print/{aid}", timeout=120,
                         headers={"User-Agent": "duotransit-research/0.1"})
        r.raise_for_status()
        f.write_bytes(r.content)
        time.sleep(3)
    data = f.read_bytes()
    try:
        with tarfile.open(fileobj=io.BytesIO(data)) as tf:
            return "\n".join(tf.extractfile(m).read().decode("utf-8", "replace") for m in tf.getmembers()
                             if m.isfile() and m.name.lower().endswith((".tex", ".txt", ".dat", ".csv")))
    except tarfile.ReadError:
        import gzip
        try:
            return gzip.decompress(data).decode("utf-8", "replace")
        except OSError:
            return data.decode("utf-8", "replace")


def main():
    hits = {}
    for aid, label in PAPERS.items():
        try:
            txt = source_text(aid)
        except Exception as e:  # noqa: BLE001
            print(f"{aid}: FAILED {e}")
            continue
        ids = {int(m) for p in PAT for m in p.findall(txt)}
        print(f"{aid} {label}: {len(ids)} TIC IDs")
        for i in ids:
            hits.setdefault(i, []).append(aid)
    df = pd.DataFrame({"TIC": list(hits), "papers": [";".join(v) for v in hits.values()]}).sort_values("TIC")
    df.to_csv(OUT, index=False)
    print(f"{len(df)} unique TIC -> {OUT}")


if __name__ == "__main__":
    main()
