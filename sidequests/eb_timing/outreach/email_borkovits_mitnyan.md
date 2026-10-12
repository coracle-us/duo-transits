# Draft email: southern CVZ eclipse timing (Borkovits / Mitnyan)

**To:** borko@bajaobs.hu (Tamás Borkovits, corresponding author)
**Cc:** mtibor@titan.physx.u-szeged.hu (Tibor Mitnyan)
These addresses come from the arXiv source of Borkovits et al. 2026 (arXiv:2604.04004, Universe 12, 107),
checked 2026-10-12.

**Subject:** Southern TESS CVZ eclipse timing: is a study in progress? (plus four triple candidates)

**Attach (optional):** `figs/32000625_oc.png`

---

Dear Dr Borkovits and Dr Mitnyan,

I'm an independent amateur. I have been working through TESS data with my own scripts, and I've
read your two northern CVZ eclipse-timing papers (Mitnyan et al. 2024; Borkovits et al. 2026,
Universe 12, 107) with great interest. Before going further, I wanted to ask you directly whether
you are already working on the southern CVZ.

As far as I can tell, nobody has published a systematic ETV study of the southern CVZ EBs using the
full baseline now available: Sectors 1–98, up to about 43 sectors over 7.4 years. If your group is
doing it, I'll stop, and I'm glad to hand over the candidates below for you to use however you
like, with no credit needed. If not, I'd plan a short Research Note of the AAS listing southern
triple candidates. I'd welcome your advice, or a collaboration if that interests you.

**What I've done so far.** I fit each eclipse with an empirical per-system template and build
separate O−C curves for primary and secondary eclipses. I then fit an LTTE orbit with a quadratic
term.

As a check, I ran it on systems from your NCVZ list on the same S14–86 baseline:
- TIC 424461577: P2 = 1008 d vs your 1013.7 d, amplitude 466 s vs 462 s.
- TIC 229771234: 359.3 d vs 359.1 d, amplitude 176 s vs 181 s.
- TIC 219900027: 1635 d vs 1631 d, amplitude 139 s vs 134 s.
- TIC 167692429 and TIC 220397947 (Borkovits et al. 2020): outer periods recovered to <1%. My
  pure-LTTE model does not capture their dynamical terms.

**A pilot on 20 long-baseline southern EBs** found four third-body candidates. As far as I can find,
none has a published ETV solution:

- TIC 32000625 (TYC 9155-22-1; P1 = 0.7018 d):
  - P2 ≈ 750 d (about ±5 d realistically), A_LTTE ≈ 370 s, e2 ≈ 0.43;
  - about 3.6 outer cycles covered, with primary and secondary moving together;
  - f(m3) ≈ 0.13 Msun;
  - Gaia DR3 RUWE 11.5.
- TIC 279323859 (P1 = 1.0800 d): P2 ≈ 2,700 d, A ≈ 790 s. About one cycle covered.
- TIC 278861826 (P1 = 1.0676 d): P2 ≳ 2,600 d, A ≳ 750 s. About one cycle covered.
- TIC 294327312 = HD 56910 (P1 = 0.9185 d; Gaia DR3 SB2): P2 ≈ 2,200 d, A ≈ 690 s. About one
  cycle covered.

**One request.** Your 2026 paper says the NCVZ times of minima are available on request. If you're
willing to share them, I would use them to test my timings point by point, especially their error
bars, which I know are underestimated where spots are present.

I'm happy to send the eclipse times, O−C plots and scripts for any of these systems. Corrections
are very welcome.

A disclosure: the analysis was done with substantial help from an AI assistant (Anthropic's
Claude). I have checked the results against the saved data and figures.

Best regards,
Jon McBee
jon.mcbee@gmail.com
