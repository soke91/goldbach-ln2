# Adversarial Review Verdict — the proof program is refuted in its core reductions

*An independent adversarial review (fresh context, refutation mandate,
line-checked against Lichtman arXiv:2309.08522 §3–§5) returned this
verdict at increment 143; items decided on a magnitude or norm basis
were re-audited at increments 199–200 against the corrected E1 target,
and are stated below in their settled form. We publish the whole thing
in full spirit: the failures are part of the record.*

## Refuted — structural, and fatal

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
   T5/T6 down it is not "the entire remainder". *Note for anyone
   reviving it*: being over-strong and false is not evidence against
   the correctly normalized statement, which has not been derived.

**These three are independent of any target and settle the program.**
Items 4–6 below were budget or norm calls; none of them changes that.

## Re-audited — none of these revives anything

4. **The h-range collapse and the q-layer ("excess K^{2/3}").** The
   excess was computed against a target that was itself misstated (see
   CLOSURE_REAUDIT.md); at K = N^{1/3} it is N^{2/9}, far inside the
   actual margin N/K ≥ N^{2/3}, so it is probably no obstruction at
   all. It is not re-derived here, because items 1 and 2 kill the
   program that contains it.
5. **The second Vaughan layer.** Its divisor variables do multiply the
   conductor, and that half was likewise a budget call now in doubt;
   the surviving half — that the type-II variant re-enters a binary
   equation — is structural and stands. Same conclusion: immaterial to
   the verdict.
6. **Theorem E1 (L²) versus the pipeline (signed L¹).** The mismatch
   is real as an observation and **not fatal as an objection**. What
   the chain consumes is the signed sum Σ_k b_k D(k); the L² statement
   enters only through Cauchy–Schwarz, so L² is *sufficient* — the
   implication runs the safe way — and strictly stronger than needed.
   Measured price of that step: a factor ≍ √K (58× and 45× at
   K = 10³ and 3·10³, `code/norm_audit.py`), which the margin absorbs.
   The corollary matters more than the verdict: **Cauchy–Schwarz
   discards exactly the sign structure of b_k**, whose power is
   already on record — for b_k ≡ 1 the *unrestricted* double sum is
   exactly μ(N−1). The program has been aiming at something strictly
   harder than its own target.

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
conductor-collapse mechanism — the binary-correlation difficulty this
program faced at the outset, now with a much more precise map of why
each route fails, and with one route (C-III) still open.

*The proof-program documents (PROOF_SKETCH_E1.md, paper/e1_proof.tex,
paper/e1_transcription.md) are retained as a record with this verdict
prepended in banner; MEASUREMENTS.md is unaffected. Supersessions and
the corrections that produced them are recorded in one place:
CLOSURE_REAUDIT.md.*
