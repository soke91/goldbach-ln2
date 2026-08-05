# STATUS — single-page state of the program

*Last updated: increment 143 (2026-08-05). 125 commits, 6 days.*

## Where things stand

**Goldbach for large even N** reduces (Huang–Li 2022 + this program's
chain) to one measurable object: the dilate-averaged Möbius bound E1.
This repository contains:

1. **The measurement corpus** (~50 scripts, the primary artifact —
   MEASUREMENTS.md): the ln 2 constant, Conjecture P / Buchstab
   profile, structure laws, χ² ladders, the final-axiom landscape
   through the √N barrier, the thin-progression stamps (10,000
   classes, sub-half-normal with exactly Gaussian tail), and the full
   engine stamps. One-shot reproduction: `python code/verify_all.py`.

2. **A proof-program attempt for E1, now refuted in its core
   reductions** (PROOF_SKETCH_E1.md + paper/e1_proof.tex +
   paper/e1_transcription.md, retained as a record). An independent
   adversarial review (REVIEW_VERDICT.md) found the gate arithmetic
   vacuous (Q-slot premise violated; the claimed δ < 0.17 exceeds the
   source paper's own optimum — an independent proof of invalidity),
   the role-exchange into the dispersion lemma unfounded (no pair
   congruence, hence no conductor collapse), and the SEAM
   formalization over-normalized and falsified by our own data. The
   exact identities (T2, T3) and every measurement survive.

3. **The honest open problem, restated**: off-diagonal dispersion of
   C_{k,k′} = Σ_p μ(N−pk)μ(N−pk′) over the whole range K ≤ x^{1/3},
   without a conductor-collapse mechanism — the binary-correlation
   difficulty, now with a precise map of why each attempted route
   fails (15 documented "teeth", 25 corrections).

## What is and is not claimed

- No theorem toward Goldbach is claimed.
- Claimed: the measurements (all reproducible, affirmed by the
  adversarial review), the exact identities, and the failure map.
- NOT claimed (withdrawn at increment 143): "gates pass", "remainder
  = one 1/30-wide seam", "the rest is transcription".

## The factorization law (increments 144–169)

After the refutation, the measurement campaign converged on one
formal conjecture — **Conjecture L** (CONJECTURE_L.md): every μ-family
this program probed factorizes as (deterministic local mask,
computable by finite modular enumeration) × (exactly Gaussian
fluctuation on the surviving support). The mask is blind-verified
(corr 1.0000; annihilation fractions predicted exactly at fresh N);
the Gaussian half holds at pair, cell, matrix, and E1-ratio level
across 10⁸–2×10⁹ and four N-structures; five self-raised challenges
(including one of our own null-design errors) all resolved without
modifying the law. Every "sub-random" reading in the program's
history is mask accounting. The single unproven statement feeding the
Goldbach chain is the amplitude half: square-root cancellation of a
fluctuation that is featureless in every measurement.

## How to continue

- **Entry point**: MEASUREMENTS.md, then `python code/verify_all.py`.
- **The verdict**: REVIEW_VERDICT.md — read before the sketch.
- **History**: 143 dated increments with 25 corrections and 15
  documented teeth — the map of every route that died, so no one
  repeats them.

*The wall, at final resolution, is a binary correlation that nature
cancels everywhere we can measure and no current technique — including
the routes this program invented and then refuted — can certify. The
dossier of measurements and failure coordinates is the deliverable.
— the program*
