# Recall: the second pass scored against the first

This is the measurement `v1_verify2` exists for. `v1_verify` ran six
rounds and returned fifteen standing findings. This tree re-verified
`v1/paper/wall_v1.tex` from the statements without opening
`v1_verify/paper/` or `v1_verify/code/`, wrote its findings to
`FINDINGS.md`, and only then opened the first pass.

Both sealed hashes matched before and after, so the target did not move:

```
5d359db3…  v1_verify/paper/ADVERSARIAL_FINDINGS.md
454e92e0…  v1_verify/paper/wall_v1_corrected.tex
```

## Scoring rules, fixed before opening the first pass

- **HIT** — the second pass reports the same defect in the same
  statement, reached independently.
- **PARTIAL** — the same statement is flagged as defective, but the
  diagnosis differs or is weaker.
- **MISS** — the second pass did not report it.
- A miss is scored **in scope** if the second pass declared that area
  covered, and **out of scope** if `FINDINGS.md` listed it as uncovered
  before the comparison. Out-of-scope misses say nothing about recall;
  they say the second pass was smaller.

## Scorecard

| # | first pass's finding | second pass | note |
|---|---|---|---|
| 1 | `lem:MP` false as stated | **HIT** (M1) | same simplex-vs-box diagnosis, same repair, both verified numerically |
| 2 | "Chowla ⟹ `rho→1`" invalid | **PARTIAL** (A5) | both call the step invalid; first pass quantifies the amplification `Gamma ~ N/(A log N)` and shows *no* bound on `\|S(h)\|` suffices. Second pass only identified the shifted-prime restriction. Weaker. |
| 3 | reconstruction `-0.0976` divides by `W` not `V` | **MISS** (in scope) | second pass localised the error *to* the reconstruction ("a defect in that particular reconstruction") but never found the cause |
| 4 | "factor 0.54" withdrawn at `#106`; a coin reproduces `rho` | **MISS** (in scope) | the coin control on `rho` was never run, though the field was built. The withdrawal history needs `v1_log`, which this pass never opened |
| 5 | `prop:W` is not an identity | **HIT** (A5) | same |
| 6 | `prop:coh`'s derivation self-cancels | **MISS** (in scope) | second pass computed the correct three-term `Var` and used it throughout, and noted `prop:coh`'s constant is 2.5–6x off, but never noticed that the substitution makes `S - 2S + S = 0` |
| 7 | zeta ordinates: white null, per-ordinate claim dies | **HIT** (M13) | independent, and sharper on one point: the second pass *identified* the paper's null by matching a permutation surrogate's max `5.62e-6` to the quoted `5.09e-6`. Missed the `#110` withdrawal and the `\%` guard failure |
| 10 | "ten kill-tested" is nine | **HIT** (A3) | see correction below |
| 11 | `t=5` tail ratio has no error bar | **HIT** (M3) | same |
| 12 | `PROVENANCE.md` numbers off by one | **MISS** (out of scope) | second pass never opened `PROVENANCE.md` |
| 13 | `sec:R4` block ratios undefined, two readings | **HIT** (M8) | same two readings, same conclusion that the verdict survives on the reading the paper's own baseline picks |
| 14 | `sec:c3` `A_j`, `D_j`, `z` undefined | **HIT** (M9) | same. First pass is sharper: it shows the `z`-rounding convention moves the third decimal by `0.017` |
| 15 | K1's orbit is 13 of 63 live; open, not dead | **MISS** (out of scope) | K1 declared uncovered |
| 16 | R2's criteria are uncalibrated | **MISS** (out of scope) | R2 declared uncovered |
| 17 | R1's null is six draws | **MISS** (out of scope) | R1's null design declared uncovered |

## The recall number

| | count |
|---|---|
| first pass's standing findings | 15 |
| HIT | 7 |
| PARTIAL | 1 |
| MISS, in scope | 3 (#3, #4, #6) |
| MISS, out of scope | 4 (#12, #15, #16, #17) |

- **Raw recall, all 15: `7.5 / 15 = 50%`.**
- **Recall over the region the second pass covered: `7.5 / 11 = 68%`.**

Either way the answer to the question this tree was built to ask is:
**one pass of this kind catches roughly half to two thirds of what is
there.** That rate applies to everything the first pass declared
confirmed, and the section below shows it applying.

## What the second pass found that the first did not

Fifteen findings, four of them high severity. Three of these **contradict
things the first pass explicitly confirmed** — which is the part that
matters, because it means the first pass's "no defect found" list is not
safe either.

### Contradicting the first pass's confirmations

**M2 — `conj:wall` item 1.** The first pass lists item 1 under "Confirmed,
no defect found": *"With cells at depth 8 the excess kurtosis after
cell-mean removal is within 1.7 sigma of zero in every octave."*

But `sec:floor` defines cells as *"the number of `3,5,7,11,13` dividing
`N`"*, and item 1 reports a **pooled** figure on `6.3e6` values, not a
per-octave one. Measured under the paper's stated index and pooled as the
paper pools, on the field where `prop:V` reproduces:

| cell index | excess kurtosis | z |
|---|---|---|
| `{3,5,7,11,13}` — **as stated** | +0.0146 | **+8.4** |
| `{3,…,19}` — wider | +0.0076 | +4.4 |
| valuations, residue classes | +0.014 to +0.029 | +8.3 to +16.6 |
| per-cell **standardisation** at `{3,…,19}` | +0.0004 | +0.2 |

The first pass verified a **wider cell set, per octave**; the paper states
a narrower one, pooled. Under what the paper says, item 1 does not hold.
Pooling across octaves is also exactly the scale mixture the first pass
warns about in its own self-correction note — and the paper's figure is
the pooled one.

**M3 — `conj:wall` item 3.** The first pass's finding 11 addresses only
the `t=5` bar and states *"The `t=3` and `t=4` figures are fine (bars
±0.011 and ±0.07)"*. Under the paper's stated index and operation:

| t | expected | measured ratio | paper |
|---|---|---|---|
| 3 | 21,463 | 1.029 | 0.999 |
| 4 | 503.6 | **1.319** | 0.997 |
| 5 | 4.56 | **5.924** (27 events) | 0.878 |

`1.319` is roughly `4.5` of the first pass's own `t=4` bar. And item 3's
sub-claim that *"the extremes are attained at generic `N`, not at deep
radicals"* is false by a wide margin: of the 50 largest `|G|`, **47 sit at
depth ≥ 3**, which is 3.1% of the field, against 1.5 expected.

**Phase 2 resolved the disagreement, and it is not a contradiction.**
`CORRECTIONS_CHECK.md` sweeps the treatment on the identical field: the
first pass's corrected counts `21441/502/4` are reproduced by treatments
that normalise **within each octave** (`{3,5,7,11,13}` × octave,
standardised, gives `21425/506/5`), while the paper's stated treatment —
cells `{3,5,7,11,13}`, means only, pooled — gives `22086/664/27`. Both
passes are arithmetically right about different operations. What neither
the paper nor the corrected paper says is which one is being reported,
and item 3's own stated reason for caring (*"`C(N)=o(N)` constrains every
`N`"*) is a claim about the aggregate field, where the tail is heavy.

**M14 — `sec:floor`'s decay table.** The first pass's finding 6 closes
with *"Nothing downstream breaks: `lem:cellmom`'s bars are computed
exactly, so the mask significances are unaffected."*

They are affected. Using the exact floor — confirmed here against 60
independent-sign draws — the mask amplitude at **depths 0 and 2 never
reaches `|z| = 3` in any of nine octaves**, and depth 0 never exceeds
`1.0`. Yet the paper's table assigns those two cells significances of
`52.0` and `70.9`. Refitting the exponents with the errors the exact floor
supports:

- no step in the "monotone rise" exceeds `2.5` standard errors, and under
  the 32-cell reading one step **reverses**;
- `chi^2/dof` for a common exponent falls from `251` to `3.3` (32 cells)
  and to `1.9` under the six-depth reading the section itself states —
  where a common exponent is **not rejected at all**;
- the paper drops depth 1 as "not measurable" and keeps depths 0 and 2,
  but on the exact floor depth 1 is the *most* detectable of the three.

### New, not touching the first pass's confirmations

| | statement | finding |
|---|---|---|
| A1 | `conj:wall` item 2 | `se(gap)=0.0056` is smaller than either arm's `0.0293`/`0.0245`; the gap is `5.9 sigma`, not 40 |
| A2 | `sec:floor` | "steps at 5 to 30 s.e." — the paper's own table gives 4.3 to 19.8 |
| A4 | `sec:closures` | `1.29` s.e. gives `9.85%`, not `8.8%`; and **"C4" is defined nowhere** |
| A6 | `thm:A` | the one-line mechanism omits the divisor-sum completion, so `m < N^{1-theta'}` is asserted for a range where it is false — printed in both `wall_v1.tex` and `theorem_A.tex`'s overview. The first pass scanned `thm:A` line by line and reports no defect |
| A7 | `prop:E` | "by Parseval" does not cover the `>> N` step, which needs `\|\|S_Lambda\|\|_1 >> sqrt N` |
| A8 | `rem:rho` | the `1e8` conversion is a `3.83%` spread where the same remark claims `0.75%` at a smaller `N`, with the spread shrinking |
| A9 | `prop:E` | four values quoted at seven abscissae; the sequence is described as "decaying" and rises twice |
| A10 | `thm:D` | the loss exponent changes from `(1/2+delta)` to `1/2` between the statement and its consequence |
| M4 | `sec:coin` | the major-arc factors `8.40`/`15.16` swing over `2.06…9.66` and `1.78…7.26` across `N = 2^14…2^24`; `15.16` is not attained anywhere |
| M5 | `sec:coin` | the `1.051–1.068` autocorrelation excess is a **mis-specified null**: the squarefree-pair density depends on `h`, and against the correct floor the ratio is `1.006` (h odd) and `1.010` (`h≡2 mod 4`) |
| M7 | `sec:margin` | `N^{0.454}` against the paper's own formula `N^{0.336}` and its own measured trend `N^{0.30}` |
| M10 | `lem:placebo` | as worded ("a fixed permutation of the label set") the lemma is vacuous; the operation with power is a permutation of the assignment |
| M11 | `prop:coh` | confirmed, and stronger than stated: the exact floor is 6x to 160x wider than a count bar, not "about ten times" |
| M6 | `sec:coin` | the shift-mass table is specific to an unstated `X = 4e6`; at `1.6e7` the shares are entirely different |
| — | `conj:wall` item 1 | the stated sample size `6.3e6` matches neither `8.0e6` (every even `N ≤ 1.6e7`) nor `7.95e6` (the field where `prop:V` reproduces) |

## Corrections to the second pass, from the comparison

Recorded because the comparison found them and a pass that only reports
the other side's errors is not doing the job.

**Retracted: the second pass's original M6 was wrong.** It reported that
`sec:coin`'s five shift-mass percentages "do not reproduce", computing
gross mass as `sum_h |M(h)P(h)|` term by term. The paper's own totals
(net `-2.2e13` against gross `4.4e13`, a ratio of 2) show "gross" means
the sum of **absolute bucket nets**. Under that reading, at `X = 4e6`, all
five reproduce exactly — `1.1 / 3.0 / 23.1 / 48.9 / 23.8` — with net
`-2.2186e13` and gross/net `1.99`. The first pass had this right and the
second pass picked the wrong reading of an undefined statistic, which is
the failure mode both trees have a rule against. M6 has been rewritten in
`FINDINGS.md` as a reproduction plus the surviving `X`-dependence point.

**Narrowed: A3, the closure count.** The second pass wrote that the
abstract's "eighteen" is wrong because the body totals 17. The first pass
narrows this correctly: the three tables carry `5 + 9 + 3` rows **plus
C-III open = 18 rows**, so "eighteen closures" is defensible if the open
class is counted. What is not defensible is "**ten** kill-tested technique
designs" against a table of nine that is headed `(9)`. Only that half
stands.

**Weaker than the first pass on finding 2.** The second pass identified
that `Chowla ⟹ rho→1` needs an unstated hypothesis, but attributed it to
the shifted-prime restriction on `S(h)`. The first pass's diagnosis is
stronger and independent of that: the amplification factor
`Gamma = (sum_h c(h))/V ~ N/(A log N)` means `S(h) = o(1)` buys nothing,
the requirement is `S(h) = o(log N / N)`, and since `c(h) >= 0` no bound
on `|S(h)|` suffices at all — only signed cancellation across `h`, which
neither Chowla nor [MRT15] supplies. That is the better finding.

## What this says about the programme

The first pass's own headline was "fifteen standing findings, five of
them serious". The second pass adds fifteen more, four of them serious,
three of which contradict findings the first pass positively confirmed.

The two passes agree on the recall question's answer: **a single
adversarial pass of this design catches on the order of half of what is
there**, and its "confirmed, no defect found" list carries the same miss
rate as everything else. Neither pass is a certificate. The specific
lesson from the three contradictions is narrower and more useful: all
three turn on a statistic the paper leaves **undefined** — what a "cell"
is, and what field a figure was computed on — and in each case both
passes silently chose a reading, and chose differently.
