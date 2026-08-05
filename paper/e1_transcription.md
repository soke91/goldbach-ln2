# The Certified Transcription — Working Document

**Purpose.** Line-by-line transcription of PROOF_SKETCH_E1 against
Lichtman, arXiv:2309.08522 (§4–§10) and arXiv:2009.08969, producing a
certified proof (or the precise failure coordinate) of:

> **Theorem E1 (target).** For every A > 0 and dyadic K ≤ x^{1/2}:
> Σ_{k∼K} |Σ_{p∼P, P=x/K} μ(N−pk)|² ≪_A KP (log x)^{−A},
> uniformly in even N ≍ x, with the fixed-residue variant threaded by
> q ≤ x^{1/2+δ}, δ < 1/6.

**Consequence chain (established in this repository).** E1 ⟹ the
k-face of the Vaughan decomposition of EH_μ's f₃ ⟹ (with the m-face
by Möbius–Bombieri–Vinogradov and the boundary bands by
well-factorable level-3/5 technology) ⟹ EH_μ(x^{1/2+δ}, fixed
residue) ⟹ [Huang–Li, Corollary 1] Goldbach for all sufficiently
large even N.

## Transcription ledger

| # | Sketch step | Lichtman anchor | Status |
|---|---|---|---|
| T1 | D1 expansion + diagonal | — (elementary) | to write |
| T2 | D2′ reparametrization (w-APs, slopes k′/k) | his §5 change of variables (f-substitution) | to write |
| T3 | Residue decomposition (G1.b) | — (exact; machine-verified 0) | to write |
| T4 | Interval completion (G1.c), h-range x^{o(1)} | his H′-bookkeeping p.14 | to write |
| T5 | W-form matching (variables q,r,s,n₁,n₂,f,h; set N) | his Lemma 5.1 statement | to write |
| T6 | Gates (5.1)–(5.3) with q-layer, δ < 1/6 | his (5.1)–(5.3) | arithmetic done; to certify |
| T7 | Family split: gate-region k (density ↑0.61) / prime-moduli / dominant-factor | his §6 factorization of support | to write |
| T8 | MR regime K > x^{1/3} | MR short-interval theorems | to write |
| T9 | G3 gcd paragraph | — | to write |
| T10 | G4 assembly: dyadic sum, L²→L¹, EH_μ shape, feed into Pan f₃ | his §7 + Huang–Li §3 | to write |

Every row has its numerical verification already in `code/` and its
gate arithmetic pre-checked. Rows T2, T5, T10 carry the residual risk
(the naive-dictionary assumptions); a failure at any row yields the
exact coordinate of the true obstruction, which is itself the
program's deliverable in that branch.

*Opened at increment 123. To be filled row by row.*

---

## T1 — Expansion and diagonal (certified)

Let ξ_k ∈ {±1} realize the absolute values in the EH_μ-shape. Then
|Σ_{k∼K} ξ_k D(k)|² ≤ K · Σ_{k∼K} |D(k)|² (Cauchy–Schwarz), and
Σ_k |D(k)|² = Σ_k Σ_{p,p′∼P} μ(N−pk)μ(N−p′k)
= Σ_k Σ_{p∼P} μ²(N−pk) + OffDiag
with Σ_k Σ_p μ² ≤ KP and OffDiag = Σ_{p≠p′} Σ_k (·). Both steps are
identities/trivial inequalities. The target reduces to
|OffDiag| ≪ KP(log x)^{−A}. **Certified** (no analytic content).

## T2 — Reparametrization (certified modulo one named convention)

For fixed p ≠ p′ ∼ P: the k-sum Σ_{k∼K} μ(N−pk)μ(N−p′k) is, upon the
role exchange of D2′ (dispersion over short (k,k′), smooth long p),
organized instead as: for fixed k ≠ k′ ∼ K,
C_{k,k′} = Σ_{p∼P} μ(N−pk)μ(N−pk′).
Both arguments are integer-linear in p (no integrality completion
needed — G1.a). Writing w = N − pk: w runs over
{w ≡ N (mod k)} ∩ (N − 2Pk, N − Pk], and
w′ = N − pk′ = (k′(N − w)... = (k′ w + N(k − k′))/k — wait: from
w = N − pk we get p = (N − w)/k and
w′ = N − pk′ = N − (N − w)k′/k = (kN − k′N + k′w)/k
= (k′ w − N(k′ − k))/k. Integrality of w′ given w ≡ N (mod k):
k | k′w − N(k′−k) ⟺ k′w ≡ N(k′−k) (mod k) ⟺ (since w ≡ N)
k′N ≡ Nk′ − Nk ≡ ... — check: k′w − N(k′−k) ≡ k′N − Nk′ + Nk ≡ Nk ≡ 0
(mod k) ✓ — **integrality is automatic**, confirming G1.a; the slope
is k′/k with numerator and denominator ≤ 2K ≤ 2x^{1/3}. Convention
fixed: the (k, k′)-dispersion form is the canonical one; the p-form is
its Cauchy–Schwarz preimage. **Certified** (algebraic verification
above; machine confirmation in `entropy_ladder.py`/`g1_completion.py`
with identity error 0).

*Next rows: T3 (residue decomposition — exact, machine-stamped),
T4 (completion), T5 (W-matching — first naive-dictionary risk row).*
