# The amplitude statement vs. existing machinery — adjudication (increment 170)

*A fresh-context adjudication (same discipline as REVIEW_VERDICT: no
optimism bias, verbatim comparison against the source papers' lemma
hypotheses) of whether any existing machine proves a nontrivial slice
of the distilled amplitude statement*

> **E1 (weak form).** For even N, D(k) = Σ_{√N<m≤N/k} μ(m)μ(N−mk):
> **Σ_{k∼K}|D(k)|² ≪ (log N)^{−2A−2} Σ_{k∼K} M_k²**, K ≤ x^{1/3}
> dyadic — the entirety of what the Huang–Li chain consumes (fixed
> log-power saving is the currency; o(1) or log log savings are not
> consumable).

> ### ⚠️ Correction #30 (increment 198) — this statement was wrong
>
> The line above originally read `≪ (log N)^{−A} Σ_{k∼K} M_k`, with
> **M_k where M_k² belongs**. The difference is a factor M ≍ N/K —
> a *power* of N, not a constant.
>
> **Why the printed form was false.** Σ_k M_k is the square-root scale
> (|D(k)| ≈ √M_k ⟹ Σ|D|² ≈ Σ M_k). So the printed target demanded
> *better than square-root cancellation by a log power* — and this
> repository's own measurements refute it: at N = 10⁸,
> **Σ|D|²/Σ M_k = 0.305 / 0.310 / 0.319** over K = 10³/3·10³/6·10³,
> against a demand of (log N)^{−4} = 8.8·10⁻⁶. False by four and a
> half orders of magnitude, at every N (`code/e1_target_audit.py`).
>
> **The correct derivation.** The wall is the type-II term
> T_II = Σ_{k∼K} b_k D(k) with |b_k| ≪ log N, needed at
> |T_II| ≪ N(log N)^{−A}. Cauchy–Schwarz in k gives
> |T_II| ≤ ‖b‖₂(Σ_k|D(k)|²)^{1/2} with ‖b‖₂² ≍ K(log N)², hence
> Σ_k|D(k)|² ≪ N²/(K(log N)^{2A+2}), and Σ_{k∼K}M_k² ≍ K(N/K)² = N²/K.
>
> **What actually changes.** The target is a fixed log-power saving
> over the **trivial** bound Σ M_k², not over the square-root scale.
> Nature's margin is then (N/K)/(log N)^{2A+2} → ∞ (since K ≤ N^{1/3}),
> though it is invisible at N = 10⁸, where (log N)^{10} = 4.5·10¹²
> dwarfs N/K = 10⁵ — the slack is asymptotic, and no computation at
> accessible N can display it.
>
> **What does not change.** All five route verdicts below stand: they
> are blocked at named-lemma level, structurally, not on margin. In
> particular Route 2's third obstruction survives — a triple-log saving
> is still not a fixed power of log. The correction sharpens the
> statement of the open problem ("any fixed log-power saving over
> trivial in an L² average") without reopening any route.

## Route verdicts (all checked at named-lemma level)

| Route | Verdict | Blocking coordinate |
|---|---|---|
| 1. MRT / Lichtman-2020 shift-averaged Chowla, shift→dilate | **Blocked** | Lemma 2.1's orthogonality factorization needs the pair constraint to be linear (h = m−n); the dilate pair constraint m′u − mu′ = N(m′−m) is bilinear — the detection phase carries a product of two active variables and factorizes under no fixed frequency. The h-average is a translation (Fejér kernel); the k-average is a dilation (hyperbolic kernel, no diagonalizing character family). |
| 2. Tao 2016 entropy decrement, k-averaged port | **Blocked ×3** | (i) Lemma 2.5 (approximate affine invariance of log measure) has no k-analog — uniform measure on [K,2K] is disjoint from its dilate; (ii) the circle-method step (Lemma 3.6) requires the sampling prime to enter the phase linearly; in dilate coordinates it enters multiplicatively (e(αpm′k)), the restriction-theorem sparsification fails, and the surviving input is the k-averaged binary correlation itself — circular; (iii) even if ported, the method's saving is triple-log (his footnote 2), below the consumption spec. |
| 3. Lichtman 2020 proof technique rerun in dilate coordinates | **Blocked** | The congruence coupling (p₁ \| N−pk ⟺ k ≡ Np̄⁻¹) is a genuine isomorphism but powers only the typical-set restriction, not the decoupling engine; the decoupling blocks at the same bilinear coordinate as Route 1, and the archimedean-uniformity slot (VK) would have to be replaced by uniformity over all moduli k ∼ K — an EH_μ-grade input, i.e., the target (circular). The Kloosterman/DI escape was already refuted (REVIEW_VERDICT items 1–2). |
| 4. Dirichlet-polynomial fourth moment / Perron | **Blocked** | Every log-saving mean-value pillar assumes coefficient multiplicativity in its own variable: Ramaré's identity (his Lemma 4.7) explicitly needs a_{mq} = b_m c_q; VK needs an Euler product / zero-free region — μ(N−u) has neither. The generic Montgomery–Vaughan mean value gives no log saving, and unfolding u = mk by Perron costs T ≳ N^{1−o(1)}, swallowing even the maximal VK saving. Opening the fourth moment reproduces the binary correlation (circular). |
| 5. Partial slices | **Partial** | (a) The type-I/linear slice: already consumed (the m-face Möbius–BV core — survives, not new). (b) **The N-averaged theorem — provable**: averaging over N ∼ X linearizes the constraint into shift-correlations covered by Lichtman 2020 Lemma 6.1/Thm 6.2, giving E1 for all even N outside a density-O((log X)^{−A/2}) exceptional set. But Huang–Li consumes E1 at every large even N, so this reproduces classical exceptional-set Goldbach — no new contribution to the chain. (c) Mask-annihilated cells: vacuous. Fixed-N nontrivial slices: **none** (any sub-range of K, any weakening of A — even o(1) — still requires the pair-object's cancellation, which is Conjecture L's own object). |

## The common obstruction (one sentence)

μ(m)μ(N−mk) couples its two variables simultaneously through the
product mk and the difference N−mk, so the pair-correlation constraint
is bilinear and is diagonalized by no single character family —
additive (Fourier/circle method: MRT, Tao's Lemma 3.6, Lichtman's
Lemma 2.1) or multiplicative (Mellin/Dirichlet polynomials: MR mean
values, VK) — while each route's decisive lemma consumes exactly that
diagonalization as a hypothesis, and the k-average, unlike the
h-average (translation) or log-average (dilation invariance), supplies
no invariance that linearizes the constraint.

## Meta

The repository's self-diagnosis in CONJECTURE_L.md — "what is missing
is a proof technique for the amplitude of a featureless object; the
binary-correlation difficulty in its purest measured form" — is
confirmed by verbatim comparison: E1 is out of reach of current
technology, and the unreachability is certified at the level of the
explicit hypotheses of the decisive lemmas of all five routes.

*Sources checked (all local): Tao 2016 (Thm 1.3, footnote 2, Prop 2.4,
Lemma 3.6), Lichtman 2020 (Lemmas 2.1, 4.5, 4.7, 6.1, Thm 6.2),
Lichtman 2023 (Lemma 5.1), Huang–Li (Thm 1 / Cor 1),
REVIEW_VERDICT.md, CONJECTURE_L.md.*
