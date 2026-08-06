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

*(Scope note, increment 262: the wall's own scalar C(N) is **not**
one of these families. It needs a location and a scale mask as well —
see "Scope, and the extension to C(N)" below. The statement as written
applies to the families listed above, where it has been re-verified.)*

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

## Scope, and the extension to C(N) (increments 237–261)

The families named in the statement above are the dilate pairs, the
integer-indexed dilates `D(k)`, and Möbius sums over thin
progressions. **The wall's own scalar `C(N)` is not among them**, and
the distinction now matters, because the extension to it needs a
different statement.

**Inside the original scope, nothing changes, and it has been
re-verified independently.** `sweep_A` reports 0 flags at |z| ≥ 4 out
of ~22 statistics on `D(k)` — and it counts the *mean* z as well as
the spread, so "no mean field" is tested and passes there. Increment
250 asked the question directly: the forcing that gives `C(N)` a mean
(n prime, `q | N` ⟹ `q ∤ N−n`) has no analogue for `D(k)`, because both
its factors are μ. Grouping 4000 values of k by how many small primes
divide k and not N gives no mean shift (all |z| < 1.2 against a
permuted control reaching 1.66). What `D(k)` carries instead is a pure
**support** mask, exactly as stated: `q² | gcd(k,N)` forces
`q² | N−mk` for every m, hence `D(k) = 0` identically — 1212 predicted
zeros, all observed.

**Outside it, the extension to `C(N)` is a different statement, and the
form below is the current one.** Its history — five superseded versions
between increments 236 and 289 — is in `CLOSURE_REAUDIT.md` (#36, #45,
#67, #68, #74, #75, #83, #84, #86, #87) and nowhere else.

> **The wall's law.**  `C(N) = m(N) + √(V(N))·G(N)`
>
> - `V(N) = Σ_{v<N} μ²(v)Λ(N−v)²`, the **exact** second moment
> - `m(N)`, the **location mask**, a deterministic term indexed by
>   which small primes divide `N`
> - `G(N)`, Gaussian — in the bulk **and in the tail**

Four things are now known about it, all measured on every even
`N ≤ 1.6·10⁷` with the mask removed by finite modular enumeration.

**1. The scale is closed-form** (Proposition V, increment 287).
`Λ(w)²` lives on prime powers with weight `(log p)²`, so `V(N)` is a
`(log p)²`-weighted count of **squarefree shifted primes**. Its local
density is `1` at `q | N` (there `q²|(N−p)` forces `p = q`) and
`1 − 1/(q(q−1))` at `q ∤ N` (a unit class mod `q²`). Hence, with
`W(N) = Σ_{w<N}Λ(w)²` — a prefix sum independent of `N`'s
factorisation —

> `V(N) = W(N)·𝔄(N)·(1+o(1))`,
> `𝔄(N) = ∏_{q∤N}(1 − 1/(q(q−1)))`, so `V(N) ~ 𝔄(N)·N·log N`.

Mirsky 1949 plus partial summation; recalled, not claimed. Verified to
`1.000000 ± 0.000145`, and per radical cell to five decimals.

**The local factor is `𝔄`, not `𝔖`.** Rescaled so only shape in `N` is
judged, the residual sd is **0.000323** for `𝔄` against **0.245235**
for `𝔖` — a factor of 760. `𝔖` is Hardy–Littlewood's, correct for the
Goldbach *count*; the *noise* has a different local factor, and the two
are both products over `q | N`, which is why the substitution survived
so long.

**2. The bulk is Gaussian, under that scale and no other**
(increment 283). Excess kurtosis **−0.0005 (z = −0.3)**, `E|X|/sd`
short of `√(2/π)` by 0.00018 (z = −0.8), on 6.3·10⁶ values, needing
cell **means** alone. Under an `𝔖N`-based scale the same data give
**+0.1704 at z = 98** — a normaliser that manufactures a heavy tail.

**3. The variance ratio, and what its excess is.** Write
`ρ(N) = Var C(N)/V(N)`, which is exactly `1` if the `μ(v)` on the
surviving support were independent signs. Measured, mask removed, it
**rises** 0.760 → 0.837 (increment 288); the raw ratio falls
1.006 → 0.858 and the two converge as the mask decays like `N^{−1/2}`.
So the wall beats a coin, by a margin that is **shrinking**. Whether
`ρ → 1` is not settled — the model comparison returns INDETERMINATE and
every parameter still walks.

What `ρ − 1` *is* was settled (Proposition W, increment 289). Expanding
`C(N)²` with `u = N−p`, `h = p′−p`:

> `ρ − 1 = (1/V)·Σ_{h≠0} c(h)·S(h)`, with `c(h) = Σ_{p′−p=h}(log p)(log p′)`
> a weighted prime-pair count and `S(h) = ⟨μ(u)μ(u−h)⟩` the **binary
> Chowla correlation**.

The wall's excess over square-root cancellation *is* a prime-pair-
weighted Chowla correlation. Chowla-type input therefore forces
`ρ → 1`, and then **the wall is exactly square-root** — nature
over-delivers by a power of `log` and no more. Measured: `S(h)` sits at
**1.051–1.068×** the random-sign floor `√(0.32264(X−h))` across five
decades of shift; the reconstruction gives `ρ−1 = −0.0976` against a
measured `−0.18`; the sign is negative. And the mass sits at
`h ≈ 10⁵–10⁶` — shifts below `10³`, where Chowla is hardest and the
averaged theorem weakest, carry **1.1%**.

**4. The tail is Gaussian too** (increment 290) — which is the half
that matters, since `C(N) = o(N)` constrains every `N` and a bulk can
match to five decimals while the tail is heavy. `max|Z|` tracks the
Gumbel law with mean deviation **+0.54 ± 0.45** from `E[max]` over
eight bands; aggregate tail counts give ratios **0.999** (`t=3`),
**0.997** (`t=4`), 0.878 (`t=5`); and the extremes are attained at
**generic** `N` (`2·8317`, `2·138917`, …), not at deep radicals, so the
mask removal is not leaking into the tail.

**The margin, at the extreme rather than at a typical `N`.**
`max|C| ≈ a_n√(𝔄(N)N log N)` gives `√N/(a_n√(𝔄 log N))` over the
requirement: `10^{4.4}` at `N = 10¹²`, **`10^{22.8}` at `10⁵⁰`**. The
requirement is not remotely tight; the whole difficulty is in *proving*
`o(N)`.

**5. The fluctuation is Gaussian in distribution but not phase-random
in `log N`** (increments 294–295). `C(N)` de-masked and `√V`-standardised
carries the **zeta ordinates**: a permutation test against 200
phase-randomised surrogates gives `z = +13.9`, with `M(N)/√N` — a known
zero sum — as positive control at `z = +15.8`. Regressing on
`cos(γ log N), sin(γ log N)` for the first ten ordinates:

> **`R² = 3.896·10⁻³`** (0.39% of the variance) against a surrogate mean
> of `2.49·10⁻⁶` and a 200-surrogate maximum of `5.09·10⁻⁶` — **1566×
> chance**, and every ordinate individually at `z ≥ 23.4` (`γ₁` at 209.5).

⚠️ **This does not overturn point 2.** Gaussian *in distribution* and
phase-random *in `log N`* are different statements, and a small
oscillatory component leaves kurtosis untouched — which is why 283's
`−0.0005` and this coexist. Conjecture L's "no class structure" is a
claim about the distribution; the spectrum is a second axis, and on it
the wall is **not** structureless.

⚠️ **0.39% is a floor, not the share.** The zeros do not stop at the
tenth; ten ordinates capture part of the zero-driven component and the
true fraction is larger by an amount not estimated here. Quoting it as
*the* share would be reading a truncated sum as a total. What is
settled is that the component is **real and subdominant**: the law needs
a spectral term **named**, not a new leading behaviour.

**And the wall is not the Mertens function** (increment 294), though
`Λ`'s mean of 1 suggests `C(N) ≈ Σ_v μ(v) = M(N)` and the exponents
nearly agree (`0.5272` against `0.5137`). `corr(C, M) ≈ 0.06`, flat
across every band and unchanged by removing the mask.

**Scope.** All five are measured at `N ≤ 1.6·10⁷`. Nothing here
constrains the sizes at which this program's own no-go theorems begin
to bite — `N ≈ 10⁴⁸⁰` and beyond.

### The original half, re-tested at that precision (increment 283 → 284)

The claim above was corrected for the **extension** to `C(N)`. The
**original** claim — for the μ-families `D(k) = Σ_m μ(m)μ(N−mk)` — was
stamped at "kurtosis 2.99–3.03" on ~12,000 pairs, i.e. `±0.045`. Retested
on **285,050** pairs (401 even `N`, 1000 values of `k`),
`code/lab_conjL_original_audit.py`:

| statistic | measured | Gaussian | z |
|---|---|---|---|
| excess kurtosis | **−0.0034** | 0 | **−0.4** |
| skewness | −0.0075 | 0 | −1.6 |
| `E\|Z\|/sd` | **0.79760** | 0.79788 | **−0.2** |
| variance ratio `E[D²]/support` | 0.99781 | 1 | — |

**It holds.** Exactly Gaussian at half-normal scale, now to `±0.009` in
excess kurtosis rather than `±0.045`. And **no class structure**, as
claimed: splitting by `gcd(k,N)` gives `|z| ≤ 2.5` across four classes
and variance ratios 0.978–1.020, a spread consistent with sampling.

**Why this half was safe and the other was not.** The same trap was
available here and the original work stepped over it. Repeating the
test with a **band-mean** support in place of the exact per-`k` count:

| normaliser | excess kurtosis | z | `E\|X\|/sd` z |
|---|---|---|---|
| exact per-`k` support (what §7 used) | −0.0034 | −0.4 | −0.2 |
| band-mean support — the stand-in | **+0.4645** | **+50.6** | **−13.1** |

Support varies by **40.4%** across `k` (range [416, 4619]), so a
constant stand-in manufactures a mixture exactly as `𝔖N` did for the
wall. **§7 counted the zero terms directly and was right; the `C(N)`
extension used a fitted `𝔖N`-based scale and was wrong at z = 98.** One
program, both halves, and the difference is whether the exact quantity
was available *and used*.

with `m(N)` a **location** mask — a deterministic mean indexed by
which small primes divide N, reaching 9 standard deviations below zero
for primorial N — and `√𝔖(N)` a **scale** mask. Both are finite
modular functions of N, so the *character* of the conjecture survives;
what fails is the clause **"no mean field"**, which is a property of
the original families and not of `C(N)`. With both masks applied, G
passes: every tail inside 3 SE and the extreme at z = +0.61
(LOCATION_MASK.md).

So the honest restatement is: **field = (support × scale × location
masks) × G**, with which masks are nontrivial depending on the family.
For the dilate families only the support mask is, and the original
statement stands as written. For `C(N)` all three are.

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

## What the chain consumes

The Goldbach chain does not need Conjecture L. It needs a bound on the
**signed** sum T_II = Σ_{k∼K} b_k D(k) ≪ N(log N)^{−A}, where b_k is
the arithmetic weight the Vaughan decomposition produces. Applying
Cauchy–Schwarz in k — a choice, not a requirement — gives the
sufficient L² form

> **Σ_{k∼K} |D(k)|² ≪ (log N)^{−2A−2} Σ_{k∼K} M_k²**,

a fixed log-power saving over the **trivial** bound Σ M_k² ≍ N²/K.
Square-root cancellation, which is what the measurements above record,
clears this with margin (N/K)(log N)^{−2A−2} → ∞. So the amplitude
half is not a strengthening of square-root cancellation; it is a
log-power saving over triviality, and nature supplies far more than is
asked. The whole difficulty is that no technique certifies any of it.

## Coverage: what is and is not reached

A fresh-context adjudication against the source papers' lemma
hypotheses (AMPLITUDE_ADJUDICATION.md) finds all five candidate routes
— shift→dilate substitution, entropy decrement, technique rerun,
Dirichlet-polynomial mean values, partial slices — blocked at
named-lemma level. The only provable E1-shaped statement is the
N-averaged version, which lands in exceptional-set territory and is
not consumable by Huang–Li. The common obstruction is that the
bilinear pair constraint of μ(m)μ(N−mk) is diagonalized by no additive
or multiplicative character family, while the k-average supplies no
linearizing invariance.

Beyond those five, the Technique Forge closed nine designs of its own
and the Construction closed three representation classes. **One route
is open**: C-III (Motohashi-type spectral realization). It needs a
legitimate transform, a classification covering the type-II region,
and quantitative averaged Chowla at fixed log-power strength; the
geometry and the current state of each requirement are in
CLOSURE_REAUDIT.md. So the accurate statement is *no coupling surface
has been found among the sources examined*, not that none exists in
any direction with a mathematical name.

*If Conjecture L (or just its amplitude half) is known, provable, or
refutable by current technology, the authors of this repository would
be grateful for a pointer — see paper/contact_drafts.md.*

*Supersessions and the corrections that produced them are recorded in
one place: CLOSURE_REAUDIT.md.*
