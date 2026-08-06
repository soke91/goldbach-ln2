# The Construction — a spectral representation of μ-pairs (opened at increment 181)

*The forge (TECHNIQUE_FORGE.md) closed every borrowing route and
identified the exact missing object. This document is the
construction site. Discipline unchanged: every candidate
representation must state its basis, its Parseval, and its lever;
every measurable consequence gets a pre-registered experiment before
any proof prose.*

## What a representation must provide

For the pair field t_m(k) = μ(m)μ(N−mk), a *spectral representation*
is a triple (B, P, L):

- **B (basis)**: a family {φ_λ} in which the field's generating
  object decomposes, with the pair constraint acting diagonally (or
  triangularly) on λ;
- **P (Parseval)**: an exact or asymptotic energy identity connecting
  Σ|t|²-type sums to Σ|coefficients|²;
- **L (lever)**: the coefficients must be computable, or boundable by
  information EXTERNAL to μμ-correlations (else circular — forge
  rounds 1–2 died exactly there).

## Candidate representation classes

| Class | Basis | Known analog | Risk |
|---|---|---|---|
| C-I. Abelian / rational spectrum | e(mα), peaks at a/q with Ramanujan-type amplitudes | Hardy–Littlewood major arcs | The mask already generates rational structure; the question is whether the field carries any EXCESS — if not, the abelian spectrum is mask-exact and the representation must be non-abelian |
| C-II. Quadric spherical basis | harmonics on the variety mu′−m′u = c under its SL₂ action, with μ⊗μ as test function | DI/Kuznetsov on determinant equations | R2 showed the raw Kloosterman phase is invisible; the full spherical basis is finer — but the lever must avoid re-consuming μμ input |
| C-III. Eisenstein-inverse / Motohashi type | spectral expansion of shifted convolutions of coefficients of 1/ζ-related objects | Motohashi's formula (d(n) = Eisenstein coefficients) | d(n) comes from a POSITIVE-weight object; μ's generating object 1/ζ has no known automorphic realization — this class asks for exactly the new mathematics |
| C-IV. Manufactured modularity | symmetry-search on Φ_N(z) = Σ_m t_m e(mz): does the generating function satisfy any approximate transformation law? | none (search, not borrow) | pure exploration; a discovered approximate symmetry would be the seed |

## Experiments

| # | Question | Pre-registered decision | Status |
|---|---|---|---|
| C1 | Does the abelian spectrum of the field exceed the mask? FFT of t_m at fixed k; rational-peak energy (q ≤ 32) real vs mask-null (same support, random signs), 8 draws, multiple k | Excess z ≥ 4 at two k's → C-I alive (handle found); else C-I closes and the construction goes non-abelian (C-II/III/IV) | **C-I CLOSED** — 0/6 hits; rational-peak energy consistent with the mask-null at every k (E_rat/E_tot ≈ 0.001, z ∈ [−4.9, +1.3]). One deficit outlier (k=3499, z=−4.86) noted-not-pursued: 8-draw null SE is imprecise and a deficit is not a handle. The abelian spectrum is mask-exact; the representation must be non-abelian (`code/e1_constr_c1.py`) |
| C2 | **Modular-inverse domain visibility** (threshold test for C-II): re-index the row field by p ↦ p̄ (mod k) and FFT over residue classes: h_k(r) = Σ_{p̄≡r} μ(N−pk), ĝ_k(a) = FFT(h). Kloosterman territory begins exactly where inverse-domain frequencies carry structure. Real profile vs mask-null (random signs, same p-support): special-a excess and profile-shape divergence | Excess z ≥ 4 at ≥ 2 of 300 k (Bonferroni-aware: threshold set for family size), or systematic profile divergence → C-II alive at row level; else C-II's row-level shadow is absent | **FIRED then DOWNGRADED on verification** — 4/300 hits (max z +10.97) collapsed under the 64-draw accurate null (+4.14 / +2.41 / +1.19 / +0.30: three evaporate, one is the selection-tail of a 300-family max), permutation null concurs, and second-N replication is 0/4. Pre-registered downgrade applies: the inverse-domain spectrum is ALSO mask-exact within measurement power; C-II's row-level shadow is absent. The 8-draw-null inflation is now a named pitfall (`code/e1_constr_c2.py`, `e1_constr_c2b.py`) |

## Ledger

- (increment 181) Construction opened; C1 launched.
