# v2 — the continuation

`v1` is frozen. Everything that follows lands here.

```
paper/      the v2 document
code/       new computations
results/    their outputs
```

## What v2 inherits

The four questions `v1` states and does not answer, in the order of how
much each would change the picture.

### 1. Does `rho -> 1`?

`rho = Var C / V` is `1` under random signs, and Proposition 16 makes
`rho - 1` a prime-pair-weighted Chowla correlation, so Chowla's
conjecture forces `rho -> 1` and the wall would be *exactly*
square-root. Measured, `rho ≈ 0.810` at `N ≈ 1e8`.

**The obstruction is not precision.** Lemma 18 (the coin control) kills
the natural estimator: replacing `mu` by random signs on the same
support leaves `V` and every other input byte-identical, and the
de-trended fit gives the same output on both. Any v2 attempt must first
exhibit an estimator that *distinguishes* `mu` from a coin, and show
that it does, before quoting a rate. With one realisation of `mu` this
may not be possible at all, and saying so would itself be a result.

### 2. The mask's decay law

Two things are open and they are different.

- **Which form.** `N^-a` against `(log N)^-b` is not separated over the
  factor `160` in `N` reachable so far. More range, or a derivation.
- **Why the exponent depends on depth.** The exponent rises
  monotonically as fewer small primes divide `N`, at 5 to 30 standard
  errors, and Proposition 22 rules out the mechanism that explains the
  mask's *size*: the singular-series excess is scale-invariant by
  construction and predicts an exponent of zero at every depth. So the
  decay is carried by something not yet identified. The frame handed
  over by that refutation is to fit `dm/se` rather than `dm`, with the
  error bar's own law taken from Proposition 21.

### 3. Is any part of the wall's spectrum `mu`'s?

The spectral measure is atomic on the rationals `j/q` with
Hardy-Littlewood weights `mu^2(q)/phi^2(q)` — which is `Lambda`'s
structure. `mu`'s own contribution is measured only through the
major-arc deficit (`8.40×` at `q=3`, `15.16×` at `q=5`). Whether
anything beyond that is `mu`'s is open.

### 4. Literature novelty of Theorems 5-6

The no-go over the whole weight space, and its survival under `EH`.
This needs a specialist or a literature search, not a computation, and
until it is done the theorems are stated as ours without a claim of
priority.

## The rules carry over

Pre-registration, adversarial review in fresh context, nulls before
thresholds, weights before comparisons, and a count is never an error
bar. See the root `README.md`.
