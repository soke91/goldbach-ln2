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
| C1 | Does the abelian spectrum of the field exceed the mask? FFT of t_m at fixed k; rational-peak energy (q ≤ 32) real vs mask-null (same support, random signs), 8 draws, multiple k | Excess z ≥ 4 at two k's → C-I alive (handle found); else C-I closes and the construction goes non-abelian (C-II/III/IV) | RUNNING (`code/e1_constr_c1.py`) |

## Ledger

- (increment 181) Construction opened; C1 launched.
