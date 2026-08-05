# Adversarial Review Verdict (increment 143) — the proof program is refuted in its core reductions

*An independent adversarial review (fresh context, refutation mandate,
line-checked against Lichtman arXiv:2309.08522 §3–§5) returned the
following verdict. We publish it in full spirit: the failures are part
of the record.*

## Refuted (fatal)

1. **The gate arithmetic (T6/G4).** The dictionary placed the EH-layer
   modulus x^{1/2+δ} into Lemma 5.1's Q-slot, violating its premise
   Q < N·x^{−ε} ≤ x^{1/3−ε}; with Q that large the W-sum's
   off-diagonal is empty and the "gates pass" computation is vacuous.
   Independent reductio: the claimed δ < 0.17 exceeds the level
   0.5977 (δ ≤ 0.0977) that the lemma's own author extracts from the
   same inequalities — one cannot beat the source paper by
   substitution into its own lemma. Three inconsistent Q-assignments
   coexist across our documents.
2. **The role exchange (T5).** Lichtman's W carries the pair
   congruence n₁ ≡ n₂ (mod qd), which powers the conductor-collapse
   change of variable; our (k, k′) pairs are free — no congruence, no
   collapse, conductors remain ~K²·q (the thin regime we ourselves
   called void). Sourcing R, S from k's own factors contradicts the
   lemma's coprimality (n₁, r₁s₁) = 1; the route into W-form
   (his Prop 5.2) demands Siegel–Walfisz-type coefficients that our
   ξ_k do not satisfy.
3. **The SEAM formalization.** Over-normalized by √P: it demands more
   than the chain needs and is falsified by our own measurement
   (half-normal constant 0.717, not (log)^{−A}). Even corrected, with
   T5/T6 down it is not "the entire remainder."

## Refuted (major)

4. The h-range collapse does not survive the q-layer (excess K^{2/3}).
5. The second Vaughan layer is not gate-neutral — its divisor
   variables multiply the conductor (as our own type-I analysis had
   shown; the type-II variant is worse, and re-enters a binary
   equation).
6. The tex's Theorem E1 (L²) does not match the pipeline (signed L¹);
   the tex T1 identity conflates the two expansions.

## Survives

- The exact identities T2 (reparametrization/integrality) and T3
  (residue decomposition) — independently re-derived and confirmed.
- The m-face μ–BV core (conditionally, pending the Huang–Li weight
  check) and the seam-width arithmetic as arithmetic.
- **The entire measurement corpus** (MEASUREMENTS.md) — the reviewer
  explicitly affirms the numerics' integrity. Measurements suggest the
  truth of the underlying cancellations; they never substituted for
  the reductions, and the reductions failed.

## Honest restatement of the open problem

The unresolved remainder is not a 1/30-wide seam. It is the
off-diagonal dispersion over the whole range K ≤ x^{1/3} without a
conductor-collapse mechanism — that is, the binary-correlation
difficulty this program faced at the outset, now with a much more
precise map of why each route fails.

*This is the program's 25th documented correction and its largest.
The proof-program documents (PROOF_SKETCH_E1.md, paper/e1_proof.tex,
paper/e1_transcription.md) are retained as a record with this verdict
prepended in banner; MEASUREMENTS.md is unaffected.*
