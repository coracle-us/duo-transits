# Opening prompt for the fresh session

**Setup (Jon):**
1. Create an empty repo, e.g. `C:\git\robovetter-autoresearch` (and on GitHub if you want it public).
2. Copy `WARMSTART.md` into it as `CLAUDE.md`.
3. Open a new Claude Code session in that folder and paste the prompt below.

---

Read CLAUDE.md. It's the brief for this project: an autoresearch loop to beat the Kepler DR25 Robovetter at
separating planets from false alarms at long periods and low MES, using NASA's simulated injection, inversion and
scrambling sets.

This session builds and validates the harness. Do NOT start the autoresearch loop yet. Work through the
"First-session checklist" in section 7:

1. Set up the repo, venv and gitignored `data/`.
2. Download and parse nasa/kepler-robovetter. Show me the per-set TCE counts and the counts in the
   100–500 d, MES ≤ 10 regime.
3. Compile the reference Robovetter and reproduce one published output file exactly.
4. Build the frozen `prepare.py` (KIC-hashed train/val/locked split, scr3 held out, SHA-256 manifest) and
   `evaluate.py` (completeness at matched effectiveness, bootstrap over stars). Compute the Robovetter
   baseline on val.
5. Run adversarial validation (obs vs inv/scr) and tell me which features leak "simulated-ness".
6. Write a simple `train.py` v0 (gradient boosting on the Robovetter metrics) and score it against the
   baseline, with σ.
7. Draft `program.md` for the loop. Keep it short. The rules: edit only `train.py`, keep a change only at > 2σ,
   never touch `locked`, log everything.

Before anything is downloaded, read KSCI-19114 enough to know what each column means. Verify every citation
you rely on.

When done, stop and give me a blunt report:
- counts, baseline, v0 and σ;
- leakage risk;
- whether the regime has enough injected planets for gains to be detectable;
- a go/no-go for an overnight loop, and how you'd run it safely (restricted permissions; edit only `train.py`).

Commit as you go. Ask before any push.
