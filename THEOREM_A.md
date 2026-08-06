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
> N (log N)^{−A} for every A > 0,
> where E_μ(t;k,a) = Σ_{n≤t, n≡a (k)} Λ(n)μ(N−n) −
> (1/φ(k)) Σ_{n≤t} Λ(n)μ(N−n).

Every ingredient other than Bombieri–Vinogradov (main term, mean term,
complete sum) gives an exponential saving N e^{−c√log N}; BV itself
yields only N(log N)^{−A} for each fixed A, its internal
Siegel–Walfisz range being q ≤ (log N)^B with B fixed. So the power of
log above is imposed by BV alone, and it is exactly what Huang–Li need.

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

## Theorem C is not numerically testable below ~10¹³ (increment 209)

Worth stating next to the theorem, because the obvious check fails and
a reader could mistake that for evidence against it.

Theorem C's right-hand side R(N) := r(N) − 𝔖(N)(N − C(N)) is
**α-free**, so the theorem predicts E₃(α) is α-independent up to its
error term. Sweeping θ ∈ [0.30, 0.70] with K = N^θ at
N = 4·10⁵, 8·10⁵, 1.2·10⁶ (`code/thmC_alpha_scan.py`), E₃/N instead
drifts monotonically — at N = 8·10⁵ from −0.021 to −2.648 — and
|E₃−R|/N reaches 2.64 against a pre-registered tolerance of 0.05.
**The test returns PROBLEM, and the diagnosis is that it has no power
here, not that the identity is in doubt.**

The reason is exact. For squarefree k,
μ(k)·μ(mk) = μ(m)μ(k)² = μ(m), so

> **Σ_{k<K} μ(k)·[AP sum at k] = Σ_{k<K} Σ_{(m,k)=1} μ(m)Λ(N−mk)**

— verified to machine precision (differences 7·10⁻¹², 0, 1.5·10⁻¹¹).
The k-weight in E₃ is therefore **log k with no sign**: E₃ gets no
cancellation in k. That is the μ(k)² = 1 mechanism of Theorem A
appearing a third time, here working against us. The resulting scale
is |E₃| ≲ 2√(NK)·log K, which bounds the measurements in all three
cases (ratios 0.21, 0.97, 0.18), so

> E₃/N ≍ N^{−(1−θ)/2}·log K.

This tends to 0 — Theorem C's error term is fine asymptotically, since
N^{(1+θ)/2} ≪_A N(log N)^{−A} — but only as a small power of N. At
θ = 0.7 one needs **N ≳ 10¹³** for E₃/N < 0.05. So the identity
cannot be confirmed or refuted by computation in the accessible range,
and increment 191's checks were of the constituent constants
(B(K) → −𝔖(N), the density identity) rather than of the identity
itself — which is the right thing to have checked.

## Theorem D — the demand side is empty (increment 195)

Theorems A and C treat the two weights that occur in Huang–Li,
w = 1 and w = log k. The obvious next question is whether some weight
*between* them works: complete divisor sum still cheap (as for w = 1),
main-term coefficient still of size 1 (as for w = log k). **No such
weight exists**, and the obstruction is quantitatively the √N barrier.

Write b := μ ∗ w, so w_k = Σ_{d|k} b_d — every weight has this form.
Two exact facts oppose each other:

- **Extraction.** B_w = Σ_{d<K,(d,N)=1} b_d·(μ(d)/φ(d))·ρ_{dN}(K/d)
  with ρ ≪ e^{−c√log x} (Huang–Li's Lemma 1). So B_w is controlled by
  the part of b sitting **at** the truncation point K = N^{θ′}: mass
  at d ≤ K^{1−ε} is damped by e^{−c√(ε log K)}.
- **BV-accessibility.** Expanding w inside the residual gives moduli
  m·d with m < N^{1−θ′}, so BV needs b supported on
  **d ≤ N^{θ′−1/2−δ}**.

The two thresholds are separated by exactly N^{1/2}.

> **Theorem D.** If b = μ∗w is supported in [1, N^{θ′−1/2−δ}] then
> ‖b‖₁/|B_w| ≫ exp(c√((1/2+δ)log N)), so the switch identity
> B_w·C(N) = Σ_u Λ(N−u)μ²(u)b_u − 𝓡_w yields at best
> |C(N)| ≪ exp(c√(½log N))·N(log N)^{−A} — **no saving of any power of
> log.** No weight extracts C(N) by divisor switching plus BV.

Clean by-product (Lemma: for squarefree u, Σ_{k|u}μ(k)w_k = μ(u)b_u):
the complete part is always **Σ_{u<N} Λ(N−u)μ²(u)b_u**. The two known
cases fall out as the two ends — w=1 ⟹ b=δ₁ (‖b‖₁=1, B_w→0);
w=log ⟹ b=Λ (‖b‖₁≍N, complete part *is* the Goldbach sum, B_w≍1).

**Scope, stated honestly**: this is a no-go for one precisely specified
method — divisor switching with BV as the only input — over that
method's entire weight space. It is not an obstruction to other
methods. Its value is that it closes the design space Theorems A and C
opened, instead of leaving it to be re-explored.

## Proposition V — the wall's exact scale, in closed form (increment 287)

Increment 283 established that Conjecture L's Gaussian half holds for
`C(N)` under the **exact** second moment
`V(N) = Σ_{v<N} μ²(v)Λ(N−v)²` and fails at `z = 98` under the fitted
`𝔖N`-based stand-in. That left `V(N)` exact but opaque. It is not.

`Λ(w)²` is supported on prime powers with weight `(log p)²`, so with
`w = N − v`,

> `V(N) = Σ_{p^k<N} (log p)²μ²(N−p^k) = Σ_{p<N}(log p)²μ²(N−p) + O(√N log²N)`,

a `(log p)²`-weighted count of **squarefree shifted primes**. Its local
density at a prime `q`:

- **`q | N`**: `q²|(N−p)` forces `q|(N−p)`, hence `q|p`, hence `p = q`
  — one term. **The primes dividing `N` impose no condition at all.**
- **`q ∤ N`**: the bad class is `p ≡ N (mod q²)`, and since `q ∤ N`
  that is a **unit** class, so it costs a density `1/φ(q²) = 1/(q(q−1))`.

> **Proposition V.** With `W(N) := Σ_{w<N}Λ(w)²` — a prefix sum
> independent of `N`'s factorisation — and
> `𝔄(N) := ∏_{q∤N}(1 − 1/(q(q−1))) = A_{Artin}/∏_{q|N}(1 − 1/(q(q−1)))`,
>
> `V(N) = W(N)·𝔄(N)·(1 + o(1))`,  hence  `V(N) ~ 𝔄(N)·N·log N`.

**This is Mirsky's 1949 theorem on squarefree shifted primes plus
partial summation. The derivation is recalled, not claimed.** What is
this program's is the identification of `V(N)` as the wall's scale
(283) and what follows next.

**Measured** (`code/lab_V_asymptotic.py`, every even `N ≤ 1.6·10⁷`).
Dividing by `W(N)` removes the analytic factor *exactly*, so no
asymptotic for `Σ(log p)²` is assumed and the arithmetic stands alone:

| band | mean `R/𝔄` | sd |
|---|---|---|
| `10⁵–2·10⁵` | 1.000023 | 0.00154 |
| `1.28·10⁷–1.6·10⁷` | **1.000000** | **0.000145** |

Per cell, `R/𝔄 = 1.00000` to five decimals for every radical from
`2\|N` through `2·3·5·7·11·13\|N` (`R = 0.748387, 0.898064, 0.945330,
0.968386, 0.977271, 0.983570`, matching `𝔄` in each).

### Proposition W — what the wall's excess cancellation is (increment 289)

Increment 288 measured `ρ(N) = Var C(N)/V(N)` rising 0.760 → 0.837 and
could not settle its limit: every parameter walked and `log log N`
barely moves over any reachable range. **More computing cannot settle
it. Identifying `ρ − 1` can.** Expanding the square, with `u = N − p`
and `h = p′ − p`,

> **Proposition W.**  `ρ − 1 = (1/V)·Σ_{h≠0} c(h)·S(h)`, where
> `c(h) = Σ_{p′−p=h}(log p)(log p′)` is a weighted prime-pair count and
> `S(h) = ⟨μ(u)μ(u−h)⟩` is the **binary Chowla correlation**.

So the wall's excess over square-root cancellation *is* a prime-pair-
weighted Chowla correlation. Chowla's conjecture gives `S(h) = o(1)`;
its averaged form over `h` is a theorem (Matomäki–Radziwiłł–Tao). Under
that input `ρ → 1`, and then

> **the wall is exactly square-root.** Nature over-delivers by a power
> of `log` and **no more** — sharper than the "nature over-delivers"
> this program has been repeating, and it closes off the hope that
> `C(N)` beats `√V`.

**Measured** (`code/lab_offdiag_chowla.py`, `X = 4·10⁶`, both
correlations computed exactly by FFT). The decomposition is an
identity, so checking *it* would be a check that cannot come out false
(#71); the content is quantitative and each test was given a bar
before the run.

| test | result | bar | verdict |
|---|---|---|---|
| (A) is `S(h)` at the random-sign floor? | rms `M(h)` / floor = **1.051, 1.054, 1.066, 1.067, 1.068** across five decades of shift | within 2× | **PASS** |
| (B) does the floor account for `ρ−1`? | reconstructed **−0.0976** vs measured −0.17…−0.19 | within 3× | **PASS** (0.54×) |
| (C) sign | **negative** | exact | **PASS** |

The floor in (A) is not `√X` but `√(0.32264·(X−h))` — the sum only
sees `n` with **both** `n` and `n+h` squarefree, density
`∏_q(1−2/q²)`. Against that the excess is a stable **5–7%**, which is
what `h` carrying small prime factors contributes and the generic
density does not.

**And the wall leans on the provable end of Chowla.** Gross mass by
shift (the ranges cancel, net `−2.2·10¹³` against gross `4.4·10¹³`, a
factor 2.0, so shares of the *net* would mislead and are not shown):

| shift | `h < 10³` | `10³–10⁴` | `10⁴–10⁵` | `10⁵–10⁶` | `>10⁶` |
|---|---|---|---|---|---|
| gross share | **1.1%** | 3.0% | 23.1% | **48.9%** | 23.8% |

Small shifts — where Chowla is hardest and the averaged theorem
weakest — carry **1.1%**. The mass sits at `h ≈ 10⁵–10⁶`, exactly where
the averaged theorem is strongest.

### The correction this forces, and it is the sharpest of the session

The campaign has used `𝔖(N) = 2C₂∏_{q|N,q>2}(q−1)/(q−2)` as *the* local
factor everywhere. For the wall's scale that is **the wrong function**.
Rescaling each predictor to `R`'s mean so only its *shape* in `N` is
judged:

| local factor | residual sd | max deviation |
|---|---|---|
| **`𝔄(N)`** (Mirsky) | **0.000323** | 0.0069 |
| `𝔖(N)` (used here) | **0.245235** | 0.617 |

A factor of **760**. So:

> **The singular series of the problem is not the local factor of the
> noise.** `𝔖(N)` is correct for the Goldbach *count* — Hardy–Littlewood
> — and was carried across to the variance of `C(N)`, where the right
> object is `𝔄(N)`. Both are products over `q | N`, which is why the
> substitution survived so long.

That single mis-identification explains three earlier findings at once:
increment 283's "93.7% of the variance of `(𝔖N)/V` is cell-explained"
(`𝔖` is the wrong cell function); increment 280's exponent `α` that
would not converge (a log-power was absorbing a wrong *arithmetic*
factor); and why `κ ≈ 0.465` still looked serviceable (`𝔖` and `𝔄`
correlate across `N`, so the mean fits and the fluctuation does not).

### Proposition D‴ — the no-go is monotone in the input (increment 278)

A no-go resting on an unconditional estimate invites one objection:
that it records our ignorance rather than an obstruction, and would
dissolve if the estimate improved. Theorem D rests on exactly one such
estimate, Huang–Li's Lemma 1, `ρ_n(x) ≪ e^{−c√log x}`. The objection is
answerable, and the answer runs the other way.

The proof of Theorem D passes through the inequality already displayed
in the paper (`theorem_A.tex`, l. 833):

> `|B_w| ≤ ‖b‖₁ · max_{d ≤ D} |ρ_{dN}(K/d)|`,

so the loss factor obeys

> **Proposition D‴.** `‖b‖₁/|B_w| ≥ 1 / max_{d≤D} |ρ_{dN}(K/d)|`.
> The right-hand side is **monotone decreasing in any upper bound for
> ρ**. Hence *every* improvement of Lemma 1 strengthens Theorem D; no
> improvement can weaken it.

The bound `e^{−c√log x}` is what the *classical zero-free region* gives,
and it is available because `1/ζ` enters to the first power
(ingredient 3 above). Under RH the same contour argument gives
`ρ_n(y) ≪ y^{−1/2+ε}`. Since b supported in `[1, N^{θ′−1/2−δ}]` forces
`K/d ≥ N^{1/2+δ}`, substitution replaces Theorem D's conclusion

> `exp(c√((1/2+δ)log N))`  by  `N^{(1/2+δ)(1/2−ε)}` — **a power of N.**

At N = 10⁵⁰ that is 3·10¹² in place of 44. The route is not blocked by
the weakness of what we can prove about ρ; it is blocked harder the
more we know.

**What is proof and what is measurement.** Proposition D‴ is the
one-line monotonicity, and it is a proof. The RH substitution is a
substitution into a published argument, not a new theorem, and is
labelled as such. Neither claims to identify ρ's true decay.
`code/lab_rho_decay.py` measures that decay directly for
`n ∈ {2, 6, 30, 2310, 30030}` up to `x = 2·10⁷`, with the discrimination
rule fixed before the run; `results/lab_rho_decay.txt`:

| n | RMS `H_exp` | RMS `H_pow` | verdict | fitted β |
|---|---|---|---|---|
| 2 | 0.763 | 0.728 | INDETERMINATE | 0.509 |
| 6 | 0.941 | 0.945 | INDETERMINATE | 0.527 |
| 30 | 0.979 | 1.034 | INDETERMINATE | 0.773 |
| 2310 | 0.202 | **0.081** | **H_pow** | 0.746 |
| 30030 | 0.209 | **0.089** | **H_pow** | 0.695 |

The split has an identified cause and it is not regime ambiguity: at
n = 2 and n = 6 **ρ changes sign** across the range (n=2 runs
−,−,−,+,−,+,+,−,+), so `log|ρ|` is ragged and neither straight line
fits. Sign oscillation is the signature of a tail governed by ζ-zeros;
a smooth power law does not oscillate. The clean `H_pow` fits are at
the *deep* n, where ρ is sign-stable over the whole accessible range —
which is plausibly pre-asymptotic, since a deep coprimality condition
thins the summand and delays the oscillation.

So the measurement is **not** evidence that β = 1/2 is the truth, and
is not offered as such. Its only load-bearing content is that the
observed decay is *at least* as fast as the unconditional bound, which
by Proposition D‴ can only help. **The robustness conclusion does not
depend on the fit**: it is the monotonicity, which holds whatever
regime ρ is really in.

### Proposition D″ — smooth weights: μ ∗ log^D = Λ_D (increment 197)

Theorem D assumes b = μ∗w is supported low enough for BV. That excludes
the most natural family of all — **w_k = f(log k), f polynomial** —
whose transform is spread over every scale. Nobody had looked past
D = 1. The structure there is explicit: **b = μ∗log^D = Λ_D**, the
generalized von Mangoldt function, which vanishes on integers with more
than D prime factors. So the complete part splits by ω(u):

| D | r = 1 | r = 2 |
|---|---|---|
| 1 | Σ_p Λ(N−p)log p — **the Goldbach sum** (Theorem C) | — |
| 2 | Σ_p Λ(N−p)log²p (Goldbach-type) | 2Σ_{pq} Λ(N−pq)log p log q (**Chen-type**) |

> **Proposition D″.** (i) D = 0 ⟹ B_w ≪ e^{−c√log K} — no extraction.
> (ii) For a monomial x^D, D ≥ 1, **every term is nonnegative**
> (Λ ≥ 0 and Λ_D ≥ 0), so CP_D ≍ N(log N)^{D−1} with fixed sign — never
> o(N). (iii) Cancelling across monomials requires tuning c₁ against
> c₂, and the ratio to be matched involves the asymptotics of the r=1
> piece — which at D=1 *is* the binary Goldbach sum. **Circular.**

**The two pieces are the same order** — the closure rests on
nonnegativity, not on separation of scales. Measured
(`code/thmD2_polyweight.py`):

| N | 10⁶ | 4·10⁶ | 1.6·10⁷ |
|---|---|---|---|
| r=1 / N | 22.51 | 25.04 | 27.46 |
| r=2 / N | 17.36 | 19.79 | 22.25 |
| **r₂/r₁** | 0.771 | 0.790 | **0.810** |
| CP₂/(N log N) | 2.886 | 2.949 | 2.997 |

The ratio drifts toward **1**, not toward 0 or ∞. Calibration: the D=1
column reproduces the Goldbach sum at 1.7565/1.7633/1.7614 N against
𝔖(10⁶) = 1.7604 ✓.

**The canonical tuning removes almost nothing.** There is exactly one
degree-2 tuning with an analytic justification: f(x) = x² − 2γx kills
the pole of ζ(s)W(s) at s = 1, i.e. makes b **mean-zero** — the best a
weight can do without knowing the correlations. Measured, it moves CP
from 39.88N/44.84N/49.71N to 37.85N/42.80N/47.67N — **about 5%**.
Mean-zero is not smallness: killing the pole removes the *average* of
b, not its *correlation* with Λ(N−·), and the correlation is the whole
problem.

### Theorem D′ — the no-go survives Elliott–Halberstam (increment 196)

The obvious objection is that BV is not the only input: Huang–Li's own
Theorem 1 *assumes* EH for Λ. Raising the level does not help.

> **Theorem D′.** If Λ has level of distribution θ_E ∈ (0,1), the
> residual is accessible only for b supported on d ≤ N^{θ_E−(1−θ′)},
> while extraction still needs mass at d ≍ K = N^{θ′}. The gap is
> **N^{1−θ_E}**, so ‖b‖₁/|B_w| ≫ exp(c√((1−θ_E)log N)) — beyond every
> power of log for each fixed θ_E < 1.

**The demand side stays closed even granting the full
Elliott–Halberstam conjecture**: the switch route needs not a stronger
level but a different mechanism.

#### Proposition D⁗ — where the boundary actually is (increment 279)

The theorem above is correctly hedged: *for each **fixed** θ_E < 1*.
The prose that stood here from increment 196 was not. It read

> ~~"Closing would require **θ_E = 1 exactly** — equidistribution of Λ
> to moduli of size N itself, where each progression holds O(1) terms
> and the statement carries no information."~~

and that skips a whole regime, because `θ_E = θ_E(N) → 1` **at a rate**
is neither a fixed level nor the boundary. Solving the inequality
instead of gesturing at it (`η := 1 − θ_E`, asking the loss to stay
below `(log N)^A`):

> `exp(c√(η log N)) ≤ (log N)^A ⟺ η ≤ (A/c)²·(log log N)²/log N.`

> **Proposition D⁗.** The switch route saves a power of log **iff**
> `1 − θ_E ≪ (log log N)²/log N`. That requirement is **strictly
> stronger than EH** (which fixes ε) and **strictly weaker than the
> vacuous θ_E = 1**.

**And the regime is not vacuous.** Moduli reach `N^{1−η}`, so a
progression retains `N^η = exp(C(log log N)²)` terms — which exceeds
`(log N)^A` for *every* fixed A, since the ratio is
`exp((log log N)(C log log N − A)) → ∞`. Not `O(1)`. The old sentence
was wrong on both halves: the boundary is not at 1, and what lies just
below it is not empty.

**Why the distinction is not pedantic.** "Vacuous" says the target is
meaningless; "stronger than EH" says it is meaningful and out of reach.
Only the second is true, and only the second leaves a well-posed
question — *how close to 1 can the level be pushed?* — instead of
closing one that was never closed.

**Increment 235 is not reopened.** Lichtman's `3/5` is a fixed
constant, and `results/lab_level_threshold.txt` column (C) confirms no
fixed level ever suffices; reason 1 there (object mismatch: `Λ(n)μ(N−n)`
is a correlation, not a single function) is independent and decisive.
The conclusion stands; the reason given for it did not.

#### What the same table says about every no-go here

Asked when Theorem D′ *acquires content*, the answer is sobering.
Solving `c√(η log N) = A log log N` for the crossover (c = 1, A = 3):

| η = 1 − θ_E | no-go bites from |
|---|---|
| 0.40 (Lichtman's 3/5) | **N ≈ 10⁴⁸⁰** |
| 0.10 | N ≈ 10³⁰⁷¹ |
| 0.01 | N ≈ 10⁵³⁷⁴⁴ |

Below those sizes the switch route loses nothing that matters and
Theorem D′ is silent. Two independent solves agree on this (the
required η at N = 10⁵⁰⁰ is 0.388, just under 0.40, and the 0.40
crossover is 10⁴⁸⁰, just under 10⁵⁰⁰). **Every no-go in this program
should be read this way**: they constrain *methods*, not any
computation anyone will run. This is stated because the theorems are
otherwise easy to read as though they bit at reachable N. They do not.

### Proposition E — the circle method has zero margin on C(N) (196)

Since the switch is closed over its whole design space, the exact
position of the *other* classical mechanism is worth recording.
With C(N) = ∫₀¹ S_Λ(α)S_μ(−α)e(−Nα)dα, both standard estimates lie
**at or above** the trivial bound |C(N)| ≤ ψ(N) ~ N:

- **(i) Cauchy–Schwarz**: ‖S_Λ‖₂‖S_μ‖₂ ~ (6/π²)^{1/2}·N(log N)^{1/2} —
  above trivial by a factor ≍ (log N)^{1/2}, which **grows**.
- **(ii) Pointwise × L¹**: sup_α|S_μ|·‖S_Λ‖₁ ≥ ‖S_μ‖₂‖S_Λ‖₁ ≫ N^{1/2}·
  N^{1/2} = N, **by Parseval** — an identity-level constraint. No
  improvement in Möbius exponential-sum technology can push
  sup_α|S_μ| below (6/π²)^{1/2}N^{1/2}.

So Davenport's uniform S_μ(α) ≪_A N(log N)^{−A} is useless here: it
saves against N, while the pairing needs the scale N^{1/2}, which
Parseval forbids. **This is the parity obstruction in circle-method
language — the binary problem sits exactly at the trivial bound with
no margin**, which is why the method that settles the ternary problem
cannot be pushed to the binary one by sharpening the Möbius input.

Measured (`code/circle_margin.py`, exact FFT on a 4N grid):

| N | 2¹⁴ | 2¹⁶ | 2¹⁸ | 2²⁰ |
|---|---|---|---|---|
| ‖S_μ‖₂/√N | 0.7798 | 0.7797 | 0.7797 | 0.7797 |
| sup\|S_μ\|/√N | 3.058 | 2.742 | 2.853 | 2.801 |
| ‖S_Λ‖₁/√N | 1.946 | 2.084 | 2.219 | 2.346 |
| (i) bound / N | 2.297 | 2.473 | 2.639 | 2.795 |
| **margin** N/(sup·‖S_Λ‖₁) | 0.168 | 0.175 | 0.158 | 0.152 |

The first row reproduces √(6/π²) = 0.7797 exactly (computation check);
the margin — which route (ii) would need to exceed 1 — sits near 0.16
and is **decaying**. Route (i) diverges like (log N)^{1/2} as
predicted. The object itself is small: C(N)/N = −0.0105, 0.0001,
0.0059, 0.0032.

**Numerical confirmation** (`code/thmD_tradeoff.py`, N = 99,999,998,
θ′ = 0.56, K = 30199): for w_k = [d₀|k] the brute-force B_w over all
k < K matches the factorised formula to every displayed digit, and
|B_w|·φ(d₀) is exactly |ρ(K/d₀)| — 1.000000 at K/d₀ = 1, 0.7500 at 7,
0.0667 at 13, 0.0034 at 143, 0.000529 at 30199. Order 1 only at the
truncation point; exponentially damped throughout the range BV admits
(here N^{θ′−1/2} = 3).

## Write-up (increment 193)

This file is the summary. **The full proof is `paper/theorem_A.tex`**:
the seven steps of Theorem A (divisor-switching identity, complete-sum
lemma, degeneracy lemma, unfolding of μ² and the coprimality with
truncations D = (log N)^{A+2} and E = (log N)^{2A+4}, τ₃-weighted BV,
the density lemma c(m) = A(N)λ(m)/m by an Euler-product computation,
and a t-uniform Abel summation); the exact partial-summation identity
behind Corollary B, E₄(α) = ∫₁^{N−1} T₁(t) dt/(N−t); the Δ-correction;
and the equivalence of Theorem C.

**Report for the authors: `paper/defect_report_18.md`** — the explicit
form of the term Δ missing from (18) and its repair, stating that
Theorem 1 and Corollary 1 stand, plus the two incidental results (E₄
unconditional; (W)_log equivalent to Goldbach). It reads standalone.
**Sending it is the user's decision and has not been done.**

Two points settled while writing out the proof, absent from the
increment-192 summary:

- **Σ_k μ(k) log k/φ(k) → −𝔖(N) is a direct consequence of Huang–Li's
  own Lemma 1**: subtract Σ μ/φ = O(e^{−c√log R}) from their
  Σ μ(d)/φ(d) · log(R/d) = 𝔖(N) + O(e^{−c√log R}). No separate lemma is
  needed. Independent check: f(s) = Π_{p∤N}(1 − p^{−s}/(p−1)) =
  ζ(s+1)^{−1} h(s) with h(0) = 2C₂ Π_{p|N,p>2}(p−1)/(p−2) = 𝔖(N).
- **Where the density c(m) actually comes from**: at p | m the local
  factor is 1/(p−1) − 1/(p(p−1)) − 1/(p²(p−1)) + 1/(p²(p−1)) = **1/p**.
  That four-term cancellation *is* the "density exponent exactly 1", and
  it is the same object as the load-bearing identity.
