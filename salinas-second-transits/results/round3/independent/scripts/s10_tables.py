"""Emit markdown tables for RESULTS.md from out/ephem_<tic>.json (no hand transcription)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import OUT, TARGETS
for tic in TARGETS:
    e = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
    sc = e["scale"]
    print(f"### TIC {tic}\n")
    print(f"- P = {e['P']:.6f} ± {e['sP']:.6f} d" + (f" (× √χ²ᵣ: ± {e['sP']*sc:.6f})" if sc > 1 else ""))
    print(f"- T0 (E=0, paper event) = {e['T0_E0']:.5f} ± {e['sT0_E0']:.5f} BTJD" + (f" (scaled ± {e['sT0_E0']*sc:.5f})" if sc > 1 else ""))
    print(f"- T0 (E={e['Epiv']}, least-correlated) = {e['T0_piv']:.5f} ± {e['sT0_piv']:.5f} BTJD" + (f" (scaled ± {e['sT0_piv']*sc:.5f})" if sc > 1 else ""))
    print(f"- χ² = {e['chi2']:.2f} for {e['dof']} dof, reduced χ² = {e['rchi2']:.2f}; N = {e['n']}\n")
    print("| E | Sector | Source | Tc (BTJD) | ± (min) | O−C (min) | O−C/σ | coverage in/ingress/egress |")
    print("|---:|---:|---|---|---:|---:|---:|---|")
    for (sec, src, ep, tc, err, cov, ing, egr) in sorted(e["used"], key=lambda u: u[2]):
        oc = tc - (e["T0_E0"] + ep * e["P"])
        flag = " **>3σ**" if abs(oc / err) > 3 else ""
        print(f"| {ep} | {sec} | {src} | {tc:.5f} | {err*1440:.1f} | {oc*1440:+.1f} | {oc/err:+.2f}{flag} | {cov:.2f}/{ing:.2f}/{egr:.2f} |")
    print()
