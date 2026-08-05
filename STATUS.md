# STATUS — single-page state of the program

*Last updated: increment 139 (2026-08-05). 117 commits, 6 days.*

## Where things stand

**Goldbach for large even N** reduces (Huang–Li 2022 + this program's
chain) to one measurable object: the dilate-averaged Möbius bound E1.
This repository contains:

1. **A complete proof program for E1** (PROOF_SKETCH_E1.md +
   paper/e1_proof.tex + paper/e1_transcription.md): ten transcription
   rows, all adjudicated; gates checked with margin δ < 0.17; every
   identity machine-verified (error exactly 0); every cancellation
   measured half-normal (three N values, moduli to 2×10⁷, 1200+ pairs).

2. **One remaining seam**: K ∈ (x^{0.3}, x^{1/3}] for the second-layer
   type I piece — exponent width 1/30. Its object is measured
   half-normal at two scales; its resistance is fully profiled (below
   3/5 unconditional, below 5/8 Selberg, GRH-√y weaker than trivial in
   its thin-progression regime); its interval shadow is already a
   theorem (zero-density). Three closing routes are named in
   paper/e1_proof.tex.

3. **The measurement corpus** (~50 scripts): the ln 2 constant,
   Conjecture P / Buchstab profile, structure laws, χ² ladders, the
   final-axiom landscape through the √N barrier, and the full engine
   stamps. One-shot reproduction: `python code/verify_all.py`.

## What is and is not claimed

- No theorem toward Goldbach is claimed.
- Claimed: the reductions (following published theorems), the exact
  identities, the measurements, and the gate arithmetic — all public
  and reproducible.
- The honest distance to a certified E1: expand the stated proof
  scopes line by line (transcription against Lichtman §5–§10), and
  close the 1/30-wide seam by one of its three routes.

## How to continue

- **Expert entry point**: paper/e1_proof.tex (the target theorem, the
  certified propositions, the seam's profile).
- **Verification entry point**: code/verify_all.py, then any engine.
- **History**: 139 dated increments with 24+ corrections and 14
  documented "teeth" (failure coordinates) — the map of every route
  that died, so no one repeats them.

*The wall, at final resolution, is one thin-progression cancellation
that nature performs everywhere we can measure and no current
technique can certify. It has an address, a profile, three roads, and
a complete reproducible dossier. — the program*
