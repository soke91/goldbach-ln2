# Conjecture L — the factorization law (the program's final conjecture)

*Crystallized at increments 144–158, after the adversarial refutation
of the proof-program sketch (REVIEW_VERDICT.md). This document states
the one conjecture the measurement corpus now supports, its exact
relation to the Goldbach chain, and what is and is not evidence.*

## Statement

For an even integer N and the families probed in this repository —

- prime-indexed dilate pairs:
  C_{k,k′} = Σ_{p∼P} μ(N−pk) μ(N−pk′),
- integer-indexed dilates and their pairs:
  D(k) = Σ_m μ(m) μ(N−mk), T = Σ_m μ(N−k₁m) μ(N−k₂m),
- Möbius sums over thin progressions (moduli L > √y) —

**Conjecture L.** Each field factorizes as

> field = **M** × **G**,

where **M** is the deterministic local mask — the support pattern
computed exactly by finite modular enumeration (forced μ = 0 at
densities determined by the v_q-data of (N, k, k′); total
annihilation when the forced valuation reaches 2) — and **G** is,
on the surviving support, a fluctuation that is **exactly Gaussian
at half-normal scale**: |Σ|/√(support) → half-normal with mean
0.798, variance ratio 1, kurtosis 3, no mean field, no class
structure, and Wishart-consistent spectral behavior of the pair
matrix.

## Measured support (all in `code/`, one-shot subset in
`verify_all.py` V5)

| Level | Stamp |
|---|---|
| Pair statistics | free class 0.97–1.02 / kurt 2.99–3.03 (two N) |
| Exact cells | every viable (v₂,v₃)-cell 0.99–1.06 / 2.8–3.1 |
| Mask, blind | corr(pred, obs) = 1.0000, max err 0.027, amplitude 1.5% |
| Matrix | λ_max at the Wishart null dead center (z = −0.19) |
| Integer field | no mean field (0.00017 ± 0.00024), LOO r1 = 0.800 |
| Seam band | naive 0.450 vs predicted 0.798×(1−0.434) = 0.452 |
| Scale 10⁹, blind | 6000 pairs: m2 1.025 (pred 1.00±0.04 ✓), r1 0.810 (0.798±0.02 ✓), kurt 2.88 (3.0±0.3 ✓) — the 1500-pair marginal r1 resolved on enlargement |
| Fresh N, v₂(N)=2, combined blind | annihilation 0.458 vs 0.458 predicted (exact); an apparent viable-class dip (0.917, −2.4σ at 1600 pairs) resolved to 1.003 / kurt 2.96 at 6539 pairs; per-cell support fractions (0.831 / 0.417 / 0) match the mask's hand-derived densities exactly |
| E1 itself at 10⁹ | dyadic-band ratio Σ\|D(k)\|²/Σ support = 0.966 / 0.950 / 0.922 over K ~ 10³ / 3·10³ / 10⁴ (predicted ~1.0, all within 1σ) — no growth; the chain's consumable sits inside the unit-Gaussian budget |
| E1 across N (robustness) | four N (v₂ = 1 and 2), three bands each: all unit-consistent after high-power settlement — an apparent +3.7σ band at one N regressed to 1.078 ± 0.072 (bootstrap, 600 k); its candidate mechanism (parity aliasing) directly refuted (corr(D_even, D_odd) = +0.06) |
| E1 definitive grid | 4 N (v₂ ∈ {1,2}, v₃ ∈ {0,4}) × 2 bands × 300 k with bootstrap SE: mean z = +0.02 (globally unbiased), 7/8 cells within 1.25σ; the one suspect cell settled by near-census (800 of the band's 1000 k): 0.896, which is z ≈ −1.5 against the law's own μ-randomness budget — within budget; grid closed with no violation (`code/e1_grid_final.py`, `code/e1_settle_900.py`) |

**Challenge ledger** (all raised by our own measurements, all resolved
before publication): a −2.4σ viable-class dip (sampling, resolved at
4× sample); a z = 9 spectral excess (our null-design category error —
iid-entry Wigner null applied to a Gram matrix; correct Wishart null
sits dead center); a +3.7σ band elevation (heavy-tailed estimator
with understated SE; bootstrap settlement z = +1.09). The law has
survived every challenge unmodified.

Every historically "sub-random" reading of this program
(0.29–0.74 across all families and bands) is explained by **M**
alone; no measurement anywhere detects structure in **G**.

## Relation to Goldbach

By Huang–Li (arXiv:2005.03811) + this repository's reduction chain,
binary Goldbach for large even N needs only a (log-power saving)
**amplitude bound on G** for the dilate family — i.e., a provable
form of "G does not conspire across k". Conjecture L is *stronger*
than needed (it asserts exact Gaussianity; the chain needs only
square-root cancellation on average). What is missing is not
knowledge of the structure — the mask is an algorithm, and G is
featureless in every measurement — but a *proof technique* for the
amplitude of a featureless object: the binary-correlation difficulty
in its purest measured form.

## What this is not

- Not a theorem, and not evidence-by-authority: the adversarial
  review (REVIEW_VERDICT.md) refuted this program's earlier attempt
  to route the amplitude bound through existing dispersion machinery;
  that verdict stands.
- Measurements at 10⁸–10⁹ cannot exclude a conspiracy that onsets
  beyond the measured range; three blind extrapolation hits and one
  marginal reading (r1 at 10⁹, +1.6σ) are recorded honestly, not
  averaged away.

## Coverage adjudication (increment 170)

A fresh-context adjudication against the source papers' lemma
hypotheses (AMPLITUDE_ADJUDICATION.md) found: all five candidate
routes (shift→dilate substitution, entropy decrement, technique
rerun, Dirichlet-polynomial mean values, partial slices) are blocked
at named-lemma level; the only provable E1-shaped statement is the
N-averaged version (exceptional-set territory, not consumable by
Huang–Li); and the common obstruction is that the bilinear pair
constraint of μ(m)μ(N−mk) is diagonalized by no additive or
multiplicative character family, while the k-average supplies no
linearizing invariance.

*If Conjecture L (or just its amplitude half) is known, provable, or
refutable by current technology, the authors of this repository would
be grateful for a pointer — see paper/contact_drafts.md.*
