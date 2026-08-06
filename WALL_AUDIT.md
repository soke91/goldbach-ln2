# The Demand-Side Audit (increment 189) — Huang–Li consumption, exact form

*A fresh-context lemma-level audit of what Huang–Li's proof actually
consumes, against the strongest unconditional supply. Full report
adopted; highlights below.*

## The wall, confirmed from the demand side

The distribution object Huang–Li require is not μ in APs but the
**fixed-sum correlation Λ(n)μ(N−n)** in APs (their §3: substituting
n ≡ N (mod k) gives μ(N−n) = μ(k)μ(m) with the modulus entering
INSIDE the arithmetic functions — the dilate coordinate). Dispersion
machinery is inapplicable at the level of definitions (the summand is
not a q-independent convolution), and in the absolute-value form the
demand already contains an open binary correlation **at the single
modulus k = 3**. The unconditional level of this object is **0**, not
1/2. This independently reconfirms, from the demand side, the same
bilinear wall the supply-side adjudication mapped.

Rewiring through well-factorable technology (BFI 4/7, Maynard 3/5,
Lichtman 66/107): **impossible** — those theorems are for Λ alone,
and the consumption weight λ(k) = μ(k)·{1, log k} is provably not
well-factorable (prime-modulus support obstruction: any
well-factorable weight vanishes on primes in (√Q, Q]; μ does not).

## The demand is weaker than this program assumed — three coordinates

1. **Fixed log-power only**: the chain needs (log N)^{3+ε}
   (absolute form) or (log N)^{4+ε} (sign-preserving form) — not
   arbitrary A. (Corollary 1's A = 2 + ε/2 threaded through E₃/E₄
   with one extra log each.)
2. **Sign-preserving weighted form suffices**: Huang–Li discard the
   weights μ(k) log k by triangle inequality on the first line; the
   proof structure permits keeping them. The weakest sufficient
   statement is
   > **(W)**  sup_{t<N} | Σ_{k<N^{θ′}, (k,N)=1} μ(k) w_k
   > E_μ(t; k, N mod k) | ≪ N (log N)^{−(4+ε)},  w_k ∈ {1, log k} —
   **inter-modulus cancellation is allowed**, a demand never examined
   because the signs die on line one of the published proof.
3. **y-uniformity relaxable** to an integral (L¹-type) form.

## VERDICT ON THE CRACK (increment 191) — see THEOREM_A.md

The derivation returned: the switch is exact and does flip the inner
sums into BV territory, but the complete divisor sum splits the two
weights — w = 1 vanishes ([h=1]) giving **Theorem A (candidate,
under verification)**, while w = log k returns **Λ(h)** (μ ∗ log = Λ)
and hence the binary Goldbach sum itself. **(W)_log is EQUIVALENT to
Goldbach**; the weakest sufficient form was the conclusion all along.
Closure at the level of identities.

**Error corrected**: the paragraph below states C(t) enters with
coefficient O(log²N). That is wrong — the coefficient is
Σ_k μ(k)log k/φ(k) → −𝔖(N) ≍ 1. The original text is retained below
for the record.

## The crack candidate (status: SUPERSEDED — retained as record)

In (W), μ(k)·E_μ carries μ(k)² = 1 on squarefree k, and the k-sum
telescopes by divisor-switching: u = mk with the truncation k ≤ N^{θ′}
pushes the surviving support to divisors k > N^{θ′}, whose cofactors
m = u/k are **small** (< N^{1−θ′} < N^{1/2}) — the inner sums flip
into Λ-sums over APs with SMALL moduli (Bombieri–Vinogradov
territory). Known hazards, explicitly flagged: (i) the collapse may
merely reproduce Huang–Li's own construction in reverse (circular);
(ii) t-dependent boundary coupling of the (m, k) ranges; (iii) main
term matching between the switched Λ-AP main terms and the
subtracted (1/φ(k))·C(t) pieces, where C(t) = Σ Λ(n)μ(N−n) is itself
an open correlation (though it appears multiplied by only O(log²N)
and would be consumed at √N-scale if half-normal — an assumption
that must NOT be smuggled in).

Derivation commissioned to a fresh context with a refutation-first
mandate (increment 190). No margin, no claim, until it returns.

*The audit's one-line summary: the wall is real on both sides of the
ledger, but the demand's weakest form (W) — with its signs alive —
is a coordinate the twenty-seven-correction map has never visited.*
