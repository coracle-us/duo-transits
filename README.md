# duotransits

Recovering orbital periods for TESS duotransit exoplanet candidates
(Hawthorn et al. 2024, arXiv:2310.17268) using newer TESS sectors.
Target output: a Research Note of the AAS.

## Quick start (Windows, PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/fetch_data.py
```

See `CLAUDE.md` for the full project briefing and plan.
