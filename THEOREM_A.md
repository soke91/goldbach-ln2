# Theorem A (candidate) and the permanent closure of the demand side

*Increment 191. Derived by fresh-context derivation from the
sign-preserving form (W) of the demand audit. **Status: candidate —
adversarial verification commissioned; not adopted until it returns.**
The program's history: three optimistic derivations refuted
(REVIEW_VERDICT, C3_REVIEW). No claim is made here beyond what the
verification survives.*

## The split

The divisor-switching collapse of (W) is legitimate and exact
(machine-verified to 1e−11 as a finite rearrangement identity), and
it does flip the inner sum into Λ-sums over APs with moduli
m < N^{1−θ′} < N^{1/2} — Bombieri–Vinogradov territory, with the
t-dependent truncation absorbed cleanly into BV's own max_y. But the
complete divisor sum the switch leaves behind splits the two weights:

| weight | complete sum Σ_{k|h} μ(h/k) w_k | consequence |
|---|---|---|
| w_k = 1 | **[h = 1]** — vanishes | **Theorem A** below |
| w_k = log k | **Λ(h)** (since μ ∗ log = Λ) | returns Σ_u Λ(N−u)Λ(h) = the binary Goldbach sum itself |

## Theorem A (candidate)

> For fixed θ′ > 1/2 and every A > 0,
> sup_{t<N} | Σ_{k<N^{θ′}, (k,N)=1} μ(k) E_μ(t; k, N mod k) |
> ≪_A N (log N)^{−A}.

Ingredients, all classical: (i) Bombieri–Vinogradov at level
N^{1−θ′}(log N)^C < N^{1/2−δ}; (ii) Σ_{m≤M,(m,N)=1} μ(m)λ(m)/m ≪
(N/φ(N)) e^{−c√log M} (PNT strength, Goldston–Yıldırım Lemma 2.1
family); (iii) Lemma 1 of Huang–Li. **No new ingredient.**

**Load-bearing line, independently verified here**: the main-term
density is c(m) = Ψ(m)/m, resting on
Σ_{g|m} μ(g)/(φ(m/g)·g·φ(g)) = **1/m** for squarefree m — exact
rational arithmetic confirms this for all squarefree m < 400 with
zero mismatches (`code/thmA_density_check.py`). Density exponent
exactly 1 is what yields PNT-strength cancellation; any other
exponent would give only (log M)^{−c} via Selberg–Delange and
Theorem A would be FALSE. A second consistency check: c(1) equals
Huang–Li's own constant A(N) of their (7).

**Consequence if it survives**: Huang–Li's E₄ / Lemma 4 consumption
becomes unconditional — one of the two EH_μ demands is free. This
advances the Goldbach chain by zero (E₄ carries no Goldbach content),
but it collapses the demand structure to a single scalar.

## The permanent closure (this is the real result)

For w_k = log k the same switch returns the binary Goldbach sum, and
in the subtracted piece the coefficient of C(t) = Σ Λ(n)μ(N−n) is
Σ_k μ(k) log k/φ(k) → **−𝔖(N) ≍ 1** (not the O(log²N) this program's
audit wrote — **that was an error, corrected here**). Both
observations combine into an *unconditional identity*:

> **E₃(α) = Σ_{n<N} Λ(n)Λ(N−n) − 𝔖(N)(N − C(N)) + O_A(N(log N)^{−A})**

— which is Huang–Li's own (22). Therefore:

> **(W)_log ⟺ binary Goldbach for large even N.**

The weakest sufficient form the demand-side audit was hunting is
*equivalent to the conclusion*. There was never any weakening to
find. The root cause is μ ∗ log = Λ — the very identity Huang–Li
start from — so the log k weight is the carrier of the Goldbach
content and **no divisor switch can fail to hand it back**. This is
closure at the level of identities, not of estimates: no choice of
θ′, truncation, or smoothing can evade it.

**Map coordinate**: the C-III verdict said "beyond the spectral door
waits the same wall." Its demand-side twin: **beyond the
divisor-switching door waits the binary Goldbach sum itself** — and
here the door is even free; the wall simply stands immediately
behind it.

## Verification tasks (before adoption)

1. Adversarial review of Theorem A (commissioned, increment 191),
   with the density computation flagged as the load-bearing line.
2. Explicit statement/proof of ingredient (ii) with the N-dependent
   Euler factor made explicit.
3. The Δ-lemma: the k-dependent truncation correction
   (k < (N−n)/α, dropped in the printed (18) of Huang–Li) — same
   BV + PNT machinery, separate statement required.
4. Map updates: WALL_AUDIT's crack-candidate section replaced by
   this verdict; the C(t)-coefficient error corrected (done above).
