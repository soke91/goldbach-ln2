# Theorem A (verified) and the permanent closure of the demand side

*Increments 190–192. Derived from the sign-preserving form (W) of the
demand audit, then put through the program's adversarial-review
protocol (the one that refuted three earlier derivations —
REVIEW_VERDICT.md, C3_REVIEW.md). **Verdict: SURVIVES**, with four
mandatory corrections incorporated below. This is the first
derivation of this campaign to survive its own review.*

## The split that decides everything

The divisor-switching collapse of (W) is legitimate and exact (finite
rearrangement identity; machine-checked to 5e−10), and it does flip
the inner sums into Λ-sums over APs with moduli m < N^{1−θ′} <
N^{1/2} — Bombieri–Vinogradov territory, with the t-dependent
truncation absorbed into BV's own max_y. What splits the two weights
is the complete divisor sum left behind:

| weight | complete sum Σ_{k\|u, (k,N)=1} μ(k) | consequence |
|---|---|---|
| w_k = 1 | **1_{rad(u)\|N}** — support N^{o(1)} | **Theorem A** |
| w_k = log k | **−Λ(u's N-coprime part)** (μ ∗ log = Λ) | returns the binary Goldbach sum itself |

## Theorem A

> **Theorem A.** Fix θ′ > 1/2. Then, uniformly in t < N,
> Σ_{k<N^{θ′}, (k,N)=1} μ(k) · E_μ(t; k, N mod k) ≪_{A,θ′}
> N e^{−c√log N} ≪ N (log N)^{−A} for every A > 0,
> where E_μ(t;k,a) = Σ_{n≤t, n≡a (k)} Λ(n)μ(N−n) −
> (1/φ(k)) Σ_{n≤t} Λ(n)μ(N−n).

(The review found the true bound is the stronger e^{−c√log N}.)

**Why it works, in one line**: after the substitution the Möbius on
the LONG variable cancels itself (μ(k)² = 1) and only μ² ≥ 0 remains
there; the surviving Möbius sits on the SHORT variable m ≤ N^{1/2−δ}
— the classical "Möbius on the short variable + BV" configuration.

**Ingredients — all classical, no new input**:
1. Bombieri–Vinogradov at level Q = N^{1/2−δ}(log N)^{4A+8}.
2. Huang–Li's Lemma 1 (Goldston–Yıldırım): Σ_{d≤R,(d,N)=1} μ(d)/φ(d)
   ≪ e^{−c√log R}.
3. Σ_{m≤x,(m,N)=1} μ(m)λ(m)/m ≪ (N/φ(N)) e^{−c√log x}, with
   λ(m) = ∏_{p|m, p∤N}(1 − 1/(p(p−1)))^{−1}. The N/φ(N) comes from
   the (m,N)=1 restriction (not from λ) and is harmless
   (≪ log log N). Density exponent is exactly 1 — 1/ζ appears to the
   first power, so the classical zero-free region gives
   e^{−c√log x}; a non-integral exponent would give only
   (log x)^{−c} by Selberg–Delange and **Theorem A would be false**.
4. Main-term density c(m) = A(N)·λ(m)/m for (m,N) = 1, resting on
   Σ_{g|m} μ(g)/(φ(m/g)·g·φ(g)) = 1/m (squarefree m) — verified here
   by exact rational arithmetic, zero mismatches for all squarefree
   m < 400 (`code/thmA_density_check.py`); c(1) = A(N) reproduces
   Huang–Li's own constant of their (7).
5. Degeneracy lemma: p | (k,N) forces n = p^j, contributing
   ≪ N^{1/2}(log N)³.

## The four mandatory corrections (from the review, incorporated)

1. **[fatal] Never Möbius-expand (k,N)=1.** Expanding over e | N
   blows the BV call count to N^ε and the budget with it. The
   condition must be discarded by the degeneracy argument (ingredient
   5) instead.
2. **[fatal, silent] Assign main terms only to (q,N)=1 classes.**
   Classes with (q,N) > 1 are degenerate (true value O(log N), not
   T/φ(q)). Missing this shifts the density by exactly N/φ(N) — and
   carried into the w = log k branch it would have produced an
   apparent *refutation of EH_μ*. **This is the most plausible
   false-positive trap this program has met**; recorded as such.
3. Truncation parameters must be A-dependent: D = (log N)^{A+2},
   E = (log N)^{2A+4} (the earlier fixed choice covers only A ≤ 6).
4. Statement fix: the complete sum is **1_{rad(u)|N}**, not [h = 1]
   (harmless — support N^{o(1)} — but the earlier wording was wrong;
   in the log-weight branch the analogous correction is material).

**Numerical confirmation** (θ′ = 0.56, N = 2.5·10⁴ … 8·10⁵):
switching identity exact to 5e−10; residual/N decreasing
0.1445 → 0.0483 (scaling as N^{1−θ′/2}, not N); corrected main term
0.1190/0.0964/0.0770/0.0616/0.0491 vs observed residual
0.1140/0.0965/0.0785/0.0618/0.0483 — 1–4% agreement, i.e. the
residual *is* the main term and the main term dies by Σμλ/m → 0
(`code/thmA_*.py`).

## Corollary (the actual gain)

> **Huang–Li's E₄ / Lemma 4 consumption is unconditional.**
> |E₄| ≪ (log N)·sup_t |Σ_k μ(k) E_μ(t;k,·)| ≪ N(log N)^{−A+1}, the
> partial summation matching exactly because the weight log(N−n) is
> k-independent. Hence their S₄(α) = O(N(log N)^{−A}) needs no EH_μ,
> and **the entire EH_μ demand collapses to the single scalar E₃**
> (plus the Δ-correction below, also unconditional in the
> Corollary-1 regime).

**Δ-correction (a genuine defect of the published paper).** The
printed (18) of Huang–Li drops the n-dependent constraint k <
(N−n)/α present in the definition of S₂(α); the missing piece is
Δ = −Σ_{2≤m≤α} μ(m)(log m) Σ_{k<K,(k,m)=1} μ²(k) Λ(N−mk), whose
trivial bound ≪ N log N exceeds the target, so it needs its own
lemma — which closes by the same BV machinery (μ again on the short
variable), main term ≪ N e^{−c√log N}.

## The permanent closure (the more important half)

For w_k = log k the same switch returns
Σ_{n<N} Λ(n)Λ(N−n) + O(N^{1/2+ε}), the residual main term does **not**
vanish (Σ_{m≤x} μλ(log m)/m → −G̃(1), with A(N)G̃(1) = 𝔖(N) exactly),
and the subtracted piece carries C(t) with coefficient
Σ_k μ(k)log k/φ(k) → −𝔖(N) ≍ 1 (my audit's "O(log²N)" was an error).
All three verified numerically. They combine into an **unconditional
identity**:

> **E₃(α) = Σ_{n<N} Λ(n)Λ(N−n) − 𝔖(N)(N − Σ_{n<N} Λ(n)μ(N−n))
> + O_A(N(log N)^{−A})**

— which is Huang–Li's own (22). Therefore

> **(W)_log ⟺ binary Goldbach for large even N.**

The weakest sufficient form the demand-side audit was hunting *is the
conclusion*. There was never a weakening to find, and no choice of
θ′, truncation, or smoothing can evade it: the root cause is
μ ∗ log = Λ, the identity Huang–Li start from, which makes the log k
weight the carrier of the Goldbach content, so every divisor switch
must hand it back. **Closure at the level of identities, not of
estimates.**

## Honest accounting

- **Net progress toward Goldbach: zero.** Theorem A removes only the
  demand that carries no Goldbach content.
- **What was gained**: the demand structure collapses to one scalar;
  that scalar is proved *equivalent* to the conclusion; the
  demand-side search for a weaker sufficient condition is
  permanently closed; and one defect of the published paper (the
  dropped constraint in (18)) is identified with its repair.
- **Map coordinate** (twin of the C-III verdict "beyond the spectral
  door waits the same wall"): **beyond the divisor-switching door
  waits the binary Goldbach sum itself — and here the door is free;
  the wall simply stands immediately behind it.**

*Adopted at increment 192 after adversarial review. Verification
code: `code/thmA_audit.py`, `thmA_scale.py`, `thmA_logw.py`,
`thmA_mtlog.py`, `thmA_fix.py`, `thmA_E3.py`,
`thmA_density_check.py`.*
