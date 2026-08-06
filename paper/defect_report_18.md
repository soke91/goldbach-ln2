<!-- PRE-SEND CHECK, increment 277 (2026-08-06) -->
<!--
Asked before sending: might the authors already know?  Checked.

  arXiv:2005.03811 has TWO versions, v1 (2020-05-08) and v2
  (2022-08-28).  v2 is the latest, and its revision note reads
  "title changed, typos corrected" -- not a mathematical correction.
  The paper is published as a Springer book chapter, CANT 2020,
  doi 10.1007/978-3-030-67996-5_17.  No erratum was found, and no
  citing work noting the issue turned up in search.

  So the step at (18) stands uncorrected in the latest public version,
  about four years on, and this report has content.

  TWO CAVEATS, both real.  Absence of an erratum in a web search is
  weak evidence: a book chapter's errata are poorly indexed, and the
  authors may well know privately.  And the value here is the repair,
  not the defect -- the missing term is exhibited and shown harmless
  under hypotheses the paper already assumes, so nothing in their
  Theorem 1 is threatened.

  Sending remains the author's decision and nothing has been sent.
-->

# A missing term in equation (18) of Huang–Li, arXiv:2005.03811v2

*Prepared as a short, self-contained note that can be sent to the
authors as-is. It reports one defect and its repair; it does not
dispute Theorem 1 or Corollary 1, which stand.*

---

## Summary

In the passage from the definition of `S₂(α)` to the displayed equation
(18), an `n`-dependent constraint on the summation variable `k` is
replaced by an `n`-free one. The two ranges differ, and the difference
`Δ` is not negligible: its trivial bound is `≫ N(log N)²`, whereas the
target is `O(N(log N)^{-A})`.

`Δ` does close, by the same machinery already present in the paper and
under hypotheses already assumed. So the conclusions of Theorem 1 and
Corollary 1 are unaffected; what is missing is one lemma.

## The discrepancy

You define

```
S₂(α) = Σ_{n<N} Λ(n) · Λ̃_α(N−n),      Λ̃_α(u) = Σ_{d|u, d>α} μ(d) log(1/d),
```

and substitute `k = u/d = (N−n)/d`. Under that substitution the
constraint `d > α` becomes

```
k < (N − n)/α                                            (†)
```

which depends on `n`. The printed (18) reads

```
S₂(α) = Σ_{k < (N−1)/α} μ(k) Σ_{n<N, n≡N (mod k)} Λ(n) μ(N−n) log(k/(N−n)),
```

i.e. (†) has been replaced by the `n`-free bound `k < (N−1)/α`, which is
the supremum of the right side of (†) over `n ≥ 1`.

## The missing term

Write `N − n = m·k`, so that (†) reads `m > α`. The right-hand side of
(18) therefore contains, in addition to `S₂(α)`, exactly the terms with
`m ≤ α`. On those terms, `N − n` squarefree gives

```
μ(k)·μ(N−n) = μ(m)·μ²(k)·1_{(k,m)=1},      log(k/(N−n)) = −log m,
```

so that

```
S₂(α) = [RHS of (18)] + Δ,

Δ = Σ_{2 ≤ m ≤ α} μ(m) (log m) · Σ_{k < (N−1)/α, (k,m)=1} μ²(k) · Λ(N − mk).
```

(The term `m = 1` drops out because `log 1 = 0`.)

Trivially `Δ ≪ N (log N)²`, which exceeds the error term you are
tracking, so `Δ` cannot simply be absorbed.

## Why it is harmless

`Δ` has a favourable shape: the Möbius factor sits on the **short**
variable `m ≤ α`, while the long variable `k` carries only `μ² ≥ 0`.
Expanding `μ²(k) = Σ_{d²|k} μ(d)` and `1_{(k,m)=1} = Σ_{e|(k,m)} μ(e)`,
truncating at `D = (log N)^{A+2}` and `E = (log N)^{2A+4}` (the tails are
bounded trivially), reduces `Δ` to sums of `Λ` over arithmetic
progressions to moduli `≪ α (log N)^{4A+8}`, together with a main term

```
A(N) · Σ_{m ≤ α, (m,N)=1} μ(m) λ(m) (log m) T_m / m,
    λ(m) = Π_{p|m} (1 − 1/(p(p−1)))^{-1},
```

with `A(N)` your constant of (7). The main term is
`O(N e^{−c√log N})` by partial summation against
`Σ_{m≤x, (m,N)=1} μ(m)λ(m)/m ≪ (N/φ(N)) e^{−c√log x}` (the density
exponent is exactly 1, so `1/ζ` occurs to the first power and the
classical zero-free region applies).

Two important details in that argument:

1. **Main terms must be assigned only to classes with `(q,N) = 1`.**
   A class with `(q,N) = g > 1` is degenerate — `g | n` forces `n` to be
   a power of a prime dividing `N`, so its true contribution is
   `O(log² N)`, not `T_m/φ(q)`. Assigning a main term there shifts the
   density by the factor `N/φ(N)`.
2. **The condition `(k,N) = 1` should be discarded by the degeneracy
   argument, not expanded by Möbius over `e | N`** — the expansion
   multiplies the number of Bombieri–Vinogradov calls by `2^{ω(N)}`.

Consequently:

* **In the Corollary 1 regime** (`EH_μ(N^{θ′})`, `θ′ > 1/2`), one has
  `α = N^{1−θ′} < N^{1/2}`, and Bombieri–Vinogradov closes `Δ`
  **unconditionally**.
* **In the Theorem 1 regime**, `α ≍ N^{θ}` and `Δ` closes under
  `EH(N^{θ}(log N)^{2A+8})`, which Theorem 1 already assumes.

So the repair costs one lemma and no new hypothesis.

## Related unconditional observation

The same divisor-switching computation, applied to your `E₄(α)`, shows
that `E₄` does **not** require `EH_μ` at all. Precisely, for any fixed
`θ′ > 1/2` and `K = N^{θ′}`,

```
sup_{t<N} | Σ_{k<K, (k,N)=1} μ(k) · E_μ(t; k) | ≪_A N (log N)^{−A}   for every A > 0,
```

unconditionally, where `E_μ(t;k)` is the fixed-class discrepancy of
`Λ(n)μ(N−n)` in the class `n ≡ N (mod k)`. Since the weight `log(N−n)`
inside `E₄` does not depend on `k`, partial summation gives
`E₄(α) = ∫₁^{N−1} (that sum at t) dt/(N−t) ≪_A N (log N)^{1−A}`. Hence
your Lemma 4 is not needed for `E₄`, and the `EH_μ` demand of §3
collapses to `E₃(α)` alone.

The reason the signed sum is tractable while the absolute-value form is
not: after the switch, the Möbius on the long variable squares itself
away (`μ(u)μ(k) = μ(m)μ²(k)` for `u = mk` squarefree) and the surviving
Möbius sits on the short variable `m < N^{1−θ′} < N^{1/2}` — the
classical "Möbius on the short variable + BV" configuration.

For completeness, the corresponding statement for `E₃` (weight
`μ(k) log k`) is *not* a simplification: since `μ ∗ log = Λ`, the same
switch returns `Σ_{n<N} Λ(n)Λ(N−n)` itself, and one recovers
unconditionally

```
E₃(α) = Σ_{n<N} Λ(n)Λ(N−n) − 𝔖(N)(N − Σ_{n<N} Λ(n)μ(N−n)) + O_A(N(log N)^{−A}),
```

which is your own (22). So the bound `E₃(α) ≪_A N(log N)^{−A}` is
*equivalent* to binary Goldbach for large even `N`.

## Verification

The identities and constants above were checked numerically
(`θ′ = 0.56`, `N` from `2.5·10⁴` to `8·10⁵`, exact enumeration):
the switching identity holds to machine precision; the residual agrees
with the predicted main term to 1–4%; the `w = log k` branch stays of
size `≍ N` as the equivalence predicts; and the density identity
`Σ_{g|m} μ(g)/(φ(m/g)·g·φ(g)) = 1/m` was confirmed in exact rational
arithmetic for all squarefree `m < 400`.

Full write-up with proofs: `paper/theorem_A.tex` in the accompanying
repository; verification scripts: `code/thmA_*.py`.
