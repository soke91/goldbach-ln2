# What this campaign is trusting, and what happens if it is wrong

*Opened at increment 220, prompted by the right question: Huang–Li's
result might be wrong, and this program should not lean on it. We
already found one defect in it — the dropped n-dependent truncation in
equation (18) — so the paper is demonstrably not error-free, and
"published" is not the same as "verified".*

**What we have and have not checked.** We audited the two consumption
sites (E₃, E₄), the step at (18), and Lemma 1 as we use it. We have
**not** independently verified Theorem 1 or Corollary 1 end to end.
That is the honest state.

## Exposure, result by result

| result | depends on Huang–Li? | if their Theorem 1 is wrong |
|---|---|---|
| **Theorem A** (w = 1 sum is ≪_A N(log N)^{−A}) | **No.** Uses their Lemma 1, which is a Goldston–Yıldırım estimate checkable on its own, plus Bombieri–Vinogradov | stands |
| **Theorem D / D′** (no weight extracts C(N) by the switch, even under EH) | **No.** A statement about the divisor-switch method | stands |
| **Proposition D″** (polynomial weights) | **No.** μ ∗ log^D = Λ_D is classical | stands |
| **Proposition E** (circle method has zero margin on C(N)) | **No.** Parseval and the trivial bound | stands |
| **Conjecture L** and every measurement | **No.** Numerical facts about μ-fields | stand |
| **C(N) = √(𝔖N)·G** law; second-moment identity | **No.** Derived and verified here | stand |
| **Theorem C** (E₃ = ΣΛΛ − 𝔖(N−C) + O_A) | **No** — it is our own derivation. That it reproduces their (22) is a *consistency check on them* | stands |
| **Corollary B** (E₄ consumption unconditional) | **Partly.** The statement that E₄ is where they consume EH_μ is a reading of their proof | the corollary would lose its target but Theorem A survives |
| **"The wall is C(N) = o(N)"** | **Yes.** This is their Theorem 1's equivalence clause | the wall would have to be relocated |
| **"Goldbach follows from EH_μ(θ′ > 1/2)"** | **Yes.** Their Corollary 1 | the whole strategic frame goes |
| **"③ needs a fixed log power ≈ 6"** | **Yes.** Read off their Corollary 1 proof (A = 2 + ε/2) | the exponent is unfounded |

So: **the theorems this program proved do not rest on Huang–Li. The
strategic frame does.**

## What is actually conditional in the equivalence

Our own Theorem C gives, unconditionally,
E₃(α) = Σ Λ(n)Λ(N−n) − 𝔖(N)(N − C(N)) + O_A(N(log N)^{−A}).
The equivalence "r̃(N) ∼ 𝔖(N)N ⟺ C(N) = o(N)" therefore needs
E₃ = o(N), which is a hypothesis, not a theorem. Empirically the
discrepancy R(N) = r̃ − 𝔖(N − C) measures as ∼ N^{0.599} = o(N)
(MEASUREMENTS §9), so the equivalence looks true — but that is
measurement, not proof, and it is exactly the kind of thing this
campaign has learned not to treat as settled.

**Consequence for how we speak.** "The wall is C(N) = o(N)" should be
stated as *conditional on the Huang–Li reduction*, not as a fact about
Goldbach. Documents that assert it flatly are overstating.

## The task this creates

**Independent re-derivation of Huang–Li Theorem 1.** Not a reading of
their proof — a derivation from Pan's decomposition forward, in our own
bookkeeping, with the (18) repair included from the start. Two reasons
it is worth the effort:

1. we have already found one error in the chain we are standing on, so
   the prior that there is exactly one is not defensible;
2. every quantitative claim we make about what the chain *needs* — the
   exponent for ③, the level θ′ > 1/2, the collapse to E₃ — is read
   off their proof, so if the proof moves, those move.

Until that is done, the honest labels are: **our theorems, verified
here; their frame, used but unverified.**

## Re-derivation, first pass (increment 221)

Walked the chain in our own bookkeeping and checked every step by
brute force at N = 2000, 4000, 6000 (`code/hl_rederive.py`).

| step | statement | verdict |
|---|---|---|
| 1 | Λ(u) = −Σ_{d\|u} μ(d) log d for u > 1 | **exact**, 0 mismatches |
| 2 | r(N) = S₁(α) + S₂(α) after splitting the divisor sum at α | **exact**, diff ≤ 2·10⁻¹¹ |
| 3 | the μ² insertion | **legitimate**, cost 164 / 18 / 248 against a budget √N log²N of 2584 / 4351 / 5862 |
| 4 | the switch on squarefree u, k = u/d | **exact**, diff ≤ 10⁻¹² |

**A false alarm of ours, recorded.** The first version of step 3
applied μ² to S₂ *after* the split and reported a non-squarefree part
larger than S₂ itself — which would have been a second defect. It is
not. Reading their §3 shows the μ² goes into the **product**
Λ(n)Λ(N−n) *before* the split, and there it is justified because
Λ(N−n) is supported on prime powers, so Λ(N−n)(1 − μ²(N−n)) survives
only on N−n = p^ℓ with ℓ ≥ 2, of which there are O(√N log N). Our
misreading, not their error — and exactly why one checks before
claiming.

## Re-derivation, second pass: the S₁ evaluation (increment 222)

Their (13) claims S₁(α) = 𝔖(N)N + O(N(log N)^{−A}) under EH at level
α. This is the half of the chain that consumes EH for Λ — the half
Theorems A and C never touch — and it had not been checked here.
Computed directly (`code/hl_S1_check.py`), with
S₁ = −Σ_{d≤α} μ(d)log d · Σ_{n≡N (d)} Λ(n)μ²(N−n):

| N | S₁/N | 𝔖(N) | ratio | \|ratio−1\| | (log N)^{−1} |
|---|---|---|---|---|---|
| 5·10⁴ | 1.5508 | 1.7604 | 0.8809 | 0.1191 | 0.0924 |
| 10⁵ | 1.6805 | 1.7604 | 0.9546 | 0.0454 | 0.0869 |
| 2·10⁵ | 1.7048 | 1.7604 | 0.9684 | 0.0316 | 0.0819 |
| 4·10⁵ | 1.7290 | 1.7604 | **0.9822** | **0.0178** | 0.0775 |

**Consistent with (13)**: the ratio climbs monotonically toward 1 and
the discrepancy falls inside the (log N)^{−1} band throughout. Stated
with its limit: at these N the (log N)^{−A} band is wide, so this
detects a gross failure of (13) and nothing finer — it cannot confirm
the exponent A.

**Status of the re-derivation.** Verified independently so far: the
skeleton (inversion, split, μ² insertion, switch) and the leading
behaviour of the S₁ evaluation. Still not re-derived: the error term
in S₁ at the claimed exponent, the treatment of E₃ and E₄ under EH_μ,
and the assembly with the Δ-repair. The frame stays labelled "used but
unverified" with those two pieces exempted.

## Not to be confused

This is not a claim that Huang–Li is wrong. It is a statement about
what we have checked. Their Theorem 1 may well be correct; the defect
we found is repairable and we showed how. The point is that this
program has been treating the frame as bedrock while treating its own
derivations to adversarial review — an asymmetry with no justification.
