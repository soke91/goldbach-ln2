# The wall's location mask

*Increments 237–250. §1 is proved, §2 is derivation with the gap
stated, §3 is measurement, §4 is what failed, §5 the limits. Every
verdict here was reached against a pre-registered criterion; the ones
that were reached and then withdrawn are in CLOSURE_REAUDIT.md as
corrections #37–#43.*

The wall is `C(N) = Σ_{n<N} Λ(n)μ(N−n) = Σ_{v<N} μ(v)Λ(N−v)`. This
campaign carried the law `C(N) = √(𝔖(N)·N)·G(N)` with `G ∼ N(0,1)`,
i.e. a pure fluctuation with zero mean. **Both halves of that were
wrong**, and the corrections are of different kinds.

## 1. What is proved

**Proposition M.1 (the forcing).** *If n is prime and q | N with
n > q, then q ∤ (N−n). Hence*

> `C(N) = Σ_{v<N, (v, rad N)=1} μ(v)Λ(N−v) + O(ω(N) log N)`.

*Proof.* `q | N` and `q | N−n` give `q | n`, impossible for a prime
`n > q`. The exceptional terms are `n = q` for the `ω(N)` primes
dividing N, each contributing `O(log N)`. ∎

So the Möbius variable is confined to the integers coprime to rad(N).
That is not a neutral restriction: sieving out the small primes leaves
the rough numbers, among which the primes — where `μ = −1` with **no
sign variation at all** — are over-represented.

**Proposition M.2 (the identity).** *With `T_j(N)` the Λ-weighted count
of shifted primes at `ω = j`, i.e.
`T_j = Σ_{p<N, ω(N−p)=j} log p · μ²(N−p)`,*

> `C(N) = (Σ_j T_j) · R_A(N)`,  `R_A(N) := Σ_j (−1)^j T_j / Σ_j T_j`,

*where `Σ_j T_j = Σ_p log p·μ²(N−p) ∼ A(N)·N` is the trivial bound.*

*Proof.* `μ = (−1)^ω` on the squarefree integers; group the sum by
`ω(N−p)`. ∎

**So `R_A(N)` *is* the location mask**, up to the trivial factor. It is
not a fitted object.

**Proposition M.3 (the supply side has no analogue).** *For
`D(k) = Σ_{√N<m≤N/k} μ(m)μ(N−mk)`: if `q² | gcd(k,N)` then
`q² | N−mk` for every m, hence `D(k) = 0` identically.*

The wall's mask comes from **Λ**, a primality constraint, which forces
coprimality on the complementary variable. Both factors of `D(k)` are
μ, which imposes squarefreeness, so the analogue is a mask on the
**support**, not on the location. Verified: 1212 predicted zeros, all
observed (increment 250).

## 2. Derivation, and where it stops

Conditioning on `N−v` prime, `v` must avoid the class `v ≡ N (mod q)`,
leaving `q−1` classes of which exactly one is `v ≡ 0` when `q ∤ N`.
Hence

> `P(q | v) = 1/(q−1)` for `q ∤ N`,  `P(q | v) = 0` for `q | N`.

Note the sign of the effect: **a small prime factor of v *helps* N−v be
prime** when `q ∤ N`, because `v ≡ 0 ⟹ N−v ≡ N ≢ 0`. Treating the
indicators as independent,

> `R_A(N) ≈ Π_{q≤z, q∤N} (q−3)/(q−1)`.

**This is refuted quantitatively (correction #40).** Its sharpest
prediction — the factor at `q=3` is exactly 0, so `3∤N` kills the mask
— fails: the mask does not vanish but **changes sign** (mean `R_A` is
`−6.10·10⁻⁴` for `3|N` against `+2.45·10⁻⁴` for `3∤N`). Per-prime
factors measure 12.98, 4.70, 2.91, 2.62, 2.28, 2.16, 2.02 against a
predicted `(q−1)/(q−3)` = 2.00, 1.50, 1.25, 1.20, 1.143, 1.125, 1.10 —
decreasing as predicted but toward **2, not 1**, which no per-prime
effect can do.

> **The identified culprit is the independence of the divisibility
> indicators conditional on N−v being prime.** The local densities
> follow from counting classes and are almost certainly right; it is
> the independence that fails, and it fails by a factor rather than by
> a correction.

An intermediate form does better. With `M_P(x) = Σ_{v≤x,(v,P)=1} μ(v)`
and one fitted constant:

| model | free parameters | weighted R² on cell means |
|---|---|---|
| additive over prime indicators | 9 | +0.219 |
| Euler-multiplicative | 9 | **−714** |
| `M` alone | 1 | +0.420 |
| `κ(N)·M`, `κ = Π_{q\|N} q/(q−1)` **derived** | 1 | +0.569 |
| `𝔖(N)·M` **derived** | 1 | **+0.632** |
| `κ(N)^{4.75}·M` fitted power | 2 | +0.811 |

`κ(N)` is derived, not fitted: the numerator of the Λ-average is
`ψ(N) ∼ N` and the denominator is `N·Π_{q|N}(1−1/q)`. That the
shifted-prime factor `𝔖(N)` beats the coprime-density factor `κ(N)`
confirms which set matters. Neither is steep enough.

## 3. What is measured

**The mask is a finite modular function of N.** Enumerating one cell
per divisibility pattern of N over `q ∈ {3,…,23}` — which is what
Conjecture L means by "computable by finite modular enumeration" —
and removing it:

| | before | after |
|---|---|---|
| excess kurtosis | +0.4683 | **+0.0143** |
| tails t = 1..5 (z) | −22, −4, +38, +110, +387 | **−2.3, +2.0, +2.6, +1.7, +0.8** |
| max\|Z\| | 9.0034 (z = +16.78) | **5.1515 (z = +0.61)** |

**It survives to N = 10⁸** (targeted computation, 150 random even N as
the control, which itself measures N(0,1)):

| group | n | mean Z | control sd |
|---|---|---|---|
| shallow `N = 2q` | 42 | +0.104 | +0.10 |
| deep `k·510510` | 60 | −4.852 | **−5.62** |
| deeper `k·9699690` | 10 | −6.014 | **−6.96** |

Deeper is stronger, as the mechanism requires. In the scale-free form
`C(N)/√N` the deep family reads −28, −21.5, −19.4 at N ∼ 10⁶, ≤3·10⁷,
≤10⁸ — real, large, and shrinking; a power fit gives `(log N)^{−1.5}`
on three points, which does not determine a law.

**The sign balance is near-deterministic at depth.** Counting how many
of `{3,…,23}` divide N:

| # | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| P(C>0) | 0.502 | 0.520 | 0.494 | 0.364 | 0.193 | **0.109** |
| z | +3.1 | +36.3 | −7.1 | −85.6 | −72.4 | −27.4 |

**The enumerated mask is the whole of the deterministic part.** After
removing it the sign balance returns to within z ≈ 5 of ½, and that
residual is **skewness**, not a missing mask: raising the enumeration from
`q ≤ 13` to `q ≤ 37` barely moves it, while the Edgeworth term
`P(X>0) − ½ = −γ₁/(6√(2π))`, which has no fitted constant, predicts
+0.00258 against +0.00317 measured for `3|N` and −0.00103 against
−0.00233 for `3∤N`.

**And the variance law was wrong too** (correction #36): measured over
every even `N ≤ 4·10⁶`, `sd(G)` drifts 2.069 → 2.306 and `sd²/log N` is
flat to 0.5%, so

> `Var C(N) ≈ 0.465·𝔖(N)·N·log N` — **the recorded law was missing a
> factor `√(log N)`**. The normaliser needing no fit at all is the exact
> second moment `Σ_v μ²(v)Λ(N−v)²`, against which the measured variance
> runs 1.006 → 0.873.

That fit is superseded. It ran at increments 236–238, two increments
before this mask was found, and was never redone; `sd(G)` there removes
one band-wide mean, so `m(N)` sat inside the measured variance
throughout, supplying **14.1% of it at `N≈10⁵`** and 1.15% at
`1.6·10⁷`. Removing it moves the fitted exponent by 0.29 — and neither
exponent is a measurement, since both walk with the window (raw
0.83→1.02, de-masked 1.60→1.30). The wall's scale is the **exact**
second moment `V(N) = Σ_v μ²(v)Λ(N−v)²`, which needs no fit at all;
see CONJECTURE_L.md. History: CLOSURE_REAUDIT.md #36, #67, #68, #83.

### The mask's own scaling law (increment 280)

The quantity removed is `Var(m)/(𝔖N)`, and it has a clean law of its own
that nobody had measured:

> `Var_mask(Z) ∝ N^g`, **`g = −0.489 ± 0.005`**, walking `−0.443 → −0.489`
> across the window — i.e. `g → −1/2`, so **`m(N) ≍ √𝔖(N)·N^{1/4}`**.

The exact value is still drifting, but **the sign is not** (~100σ). So the
location mask is **lower order than the fluctuation** and does **not**
threaten `C(N) = o(N)`. Together with `E₃` cancelling 82% of it in the
Goldbach count, the mask is harmless twice over — a real feature of the
wall, not an obstruction to the conjecture.

It also explains the contamination mechanically: a term whose share of the
variance falls like `N^{−1/2}` must bias a fitted log-exponent downward by
an amount that itself shrinks with the range — which is exactly why the raw
and de-masked estimates converge toward each other (0.70 apart on three
bands, 0.29 apart on eight).

## 4. What failed, and it is instructive

Seven corrections came out of these twelve increments, and five are the
same species — **a statistic or a criterion that quietly measures
something other than what was asked**.

- **#37** The tail was read as heavy; it was a **location** effect. The
  outliers are primorials, every one with `C < 0`.
- **#40** The derived formula's exact zero at `q=3` came back as a sign
  flip.
- **#41** `sweep_B`'s B4 read the sign balance as "no signal" at
  z = −1.55 on ~1500 values. At 1.95·10⁶ values the two halves read
  z = −94.3 and +61.6. **Pooled over all N it reads −4.12**, because
  the halves have opposite signs and cancel: *a test that pools across
  a sign-flipping mask is nearly blind to a very large effect.*
- **#42** The Edgeworth term was first written with a **plus**, which
  inverted every verdict.
- **#43** `P(D>0) = 0.346` looked like a z = −12.7 refutation of
  Conjecture L's Gaussian half. It was **an atom at zero** counted as
  not-positive; conditioning on `D ≠ 0` gives 0.5044, z = +0.46.

Two more were method, not result: an exponent scan reporting a
**boundary hit** as an optimum, and a model comparison whose
`mean < −0.02` filter **selected negative fluctuations** and
manufactured a noise floor of R² = 0.79 where the true floor is
−0.0005.

## 5. Limits

- **No closed form.** The mask is an enumerated table with a derivation
  that reproduces its shape and not its size. Applying the table is
  enough for every use here; deriving it is open.
- **The deep cells are the least measured.** The cells with the largest
  effect are the rarest, and at N ≤ 4·10⁶ the deepest hold single
  digits. Increment 241's targeted computation exists precisely because
  a sweep cannot fix that.
- **The decay rate is undetermined.** `C/√N` for the deep family falls
  −28 → −19.4 across three scales; `(log N)^{−1.5}` fits three points
  and the consistency check between the `C/√N` and `Z` exponents misses
  by 0.35.
- **Nothing here threatens the wall.** `C(N)/N` for the deep family
  runs −0.0415, −0.0092, −0.0055, −0.0034, −0.0027, −0.0023 — plainly
  to zero. The deterministic term is of `√N` scale, hence `o(N)`.

## 6. What it means for Conjecture L

Conjecture L asserts that every μ-family factorises as a deterministic
local mask times an exactly Gaussian fluctuation. The mask half was
only ever tested **on the scale**.

> **The conjecture's shape is right and it was applied to the wrong
> moment.** The wall needs a mask on the **location** as well, and with
> both applied the Gaussian half survives: every tail inside 3 SE and
> the extreme at z = +0.61.

On the supply side no location mask exists (M.3) and the deterministic
structure is a support mask, which the conjecture already covers. So
Conjecture L stands on both sides, restated: **(support and location
and scale masks) × (Gaussian)**.
