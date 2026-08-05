# Technique Forge — engineering the missing invariance (opened at increment 171)

*The adjudication (AMPLITUDE_ADJUDICATION.md) reduced "no technique
exists" to a precise design specification. This document is the forge:
candidate designs, fast kill-tests against the measurement corpus, and
the escalation rule. Discipline: every candidate gets a numeric
kill-test BEFORE any proof prose; a candidate that fails its kill-test
is recorded and closed (the cascade lesson: one-step R² = 0.68 looked
alive, two-step R²₂ = 0.19 killed it).*

## The design target (from the adjudication, verbatim spec)

Construct a transform/kernel/invariance T acting on the dilate family
{D(k)} such that the pair constraint of T-averaged correlations
becomes LINEAR in some character family (additive or multiplicative),
while T costs at most log-power in the E1 budget. Equivalently: what
translation is to the h-average and dilation-invariance of the log
measure is to Tao's average, T must be to the k-average.

## Assets unique to this program

- **A1. The exact dilation ladder** (identity error 0, verified):
  the p-divisible sub-sum of D(k) maps exactly onto −D^{(p∤)}(pk).
  The only exact invariance the family is known to possess.
- **A2. The determinant observation**: the pair constraint
  mu′ − m′u = N(m − m′) is a 2×2 determinant equation — the natural
  home of GL₂ spectral theory (Kloosterman/DI). Known blockage: our
  parameters violate the premises of the existing dispersion lemmas
  (REVIEW_VERDICT items 1–2); the forge question is whether a NEW
  weight/congruence structure can be manufactured, not borrowed.
- **A3. Conjecture L**: the target field is measured
  mask × unit-Gaussian at every level — any candidate technique may
  assume-and-verify against this exact statistical profile.

## Candidate designs and kill-tests

| # | Design | Kill-test (fast, numeric) | Status |
|---|---|---|---|
| K1 | **Multiplicative Fejér kernel on the ladder orbit**: represent D(k) as an optimally weighted combination of its orbit {D(sk) : s squarefree, s \| 30030} via A1; if the representation is near-exact, orbit-overlap across k linearizes the k-average | Full-orbit least-squares R² of D(k) on {D(sk)}: R² → 1 = alive; R² saturates < 0.9 = dead (cascade already suggests death; this is the definitive orbit-level version) | RUNNING (`code/e1_forge_kt1.py`) |
| K2 | **Determinant/Kloosterman route with manufactured congruence**: impose an artificial pair congruence k ≡ k′ (mod d) by splitting the k-average into progressions — recovers the conductor-collapse mechanism at cost d; question is whether some d-average wins a log-power | Measure the conductor-collapsed correlation strength as a function of d: does Σ_{k≡k′ (d)} C_{k,k′} gain more than the d-splitting cost? | queued |
| K3 | **Wishart/operator route**: E1 needs only the band ratio; the C-matrix is measured Wishart-clean. Design: prove the ratio bound from a moment-method bound on the Gram spectrum of the row family {μ(N−pk)}_p — rows are deterministic, but row-inner-products are the SAME correlations (circularity check needed at 4th-moment level) | Compute trace-moment ratios tr((RRᵀ)^j)/Wishart-prediction for j = 2,3,4: exact match = the moment method consumes only lower-order correlations (alive); deviation or circular dependence = dead | queued |
| K4 | **N-average descent**: the N-averaged theorem is provable (adjudication 5-(2)); design a descent that trades the N-average for a k-average within a single N via the ladder's N-independence | Identify the exact point where the N-average is consumed in Lichtman Lem 6.1 and test numerically whether the k-orbit supplies the same decorrelation | queued |

## Escalation rule

A design survives its kill-test → write the reduction in full,
adversarially review it fresh-context (the increment-143 protocol),
and only then attach it to the chain documents. A design fails → one
paragraph in this table, closed, never re-litigated.

## Ledger

- (increment 171) K1 kill-test launched.
