# Phase 2: the first pass's own corrections, checked

`v1_verify/paper/wall_v1_corrected.tex` rewrites nine passages of
`wall_v1.tex` and puts new numbers in them. Until now the first pass was
the only witness to its own repairs — the position `v1` was in before
`v1_verify` existed. This is the independent recomputation.

A correction that is itself wrong is worse than the error it replaced,
because it now carries a verification stamp. So the standard here is the
same one applied to `v1`: recompute from the statement, and say which
figures are reproduced and which are not.

## Summary

| correction | verdict |
|---|---|
| Lemma 13 — the truncated-convolution form | **CONFIRMED**, exact to 12 digits |
| Proposition 15 — `Gamma ~ N/(A log N)` | **CONFIRMED**, every quoted digit |
| §`sec:coin` — `-0.124` and `-0.052` | **CONFIRMED**, all four cells to 5 decimals |
| §`sec:coin` — the `W`-vs-`V` diagnosis | **CONFIRMED**; `W/V = 1.27080` against `1/A(N) = 1.27020` |
| Proposition 20 — `Var/(Q_cc/n_c^2)` between 2 and 60 | **CONFIRMED** independently at a different `N` |
| `conj:wall` item 3 — tail counts `21441/502/4` | **REPRODUCIBLE, BUT NOT UNDER THE STATED OPERATION** |
| §`sec:c3` — `0.933 … 0.977` | **DEFINITIONALLY DIFFERENT** from this tree's reading |
| `conj:wall` item 4 — the coin and rotation nulls | **CONFIRMED** in substance by this tree's own Phase 1 (M13) |
| §`sec:R4` — `z = -1.03 … -0.86` | not checked (sign-randomisation null not run here) |
| K1's row, R2's row, R1's row | not checked — outside this tree's coverage |

Eight of eleven checked; six confirmed outright, one confirmed in
substance, one qualified, three not checked.

## The confirmations

### Lemma 13's repair

The corrected form is: with `Chat(N) = sum_{n+v=N, n<=X, v<=X}
Lambda(n)mu(v)` and `N` running over its full support to `2X`,

```
sum_{N<=2X} Chat(N)^2 = sum_{|h|<X} M(h) P(h).
```

This tree derived the same repair independently in Phase 1 and verified
it at `X = 200, 400, 800, 1600` and again at `X = 1.6e7`, where the ratio
of the two sides is `1.000000000000`. Both passes also agree on the
diagnosis (simplex versus box) and on the size of the failure of the
published form.

### Proposition 15's amplification factor

`Gamma = (sum_{h!=0} c(h))/V = (psi(X)^2 - W(X))/V(X)`:

| `N` | first pass | here | `Gamma logN/N`, first pass | here |
|---|---|---|---|---|
| 10,000 | 1.549e3 | **1.5489e3** | 1.4266 | **1.4266** |
| 160,000 | 1.852e4 | **1.8517e4** | 1.3868 | **1.3868** |
| 4,000,000 | 3.580e5 | **3.5798e5** | 1.3605 | **1.3605** |

Every quoted digit. The consequence the correction draws — that
`S(h) = o(1)` buys nothing and the requirement is `S(h) = o(log N / N)`,
a gap of a factor `N/log N` — follows from these numbers directly.

This is also the finding this tree got **weaker** in Phase 1: it
identified that the `Chowla ⟹ rho→1` step needs an unstated hypothesis,
but attributed it to the shifted-prime restriction on `S(h)` rather than
to the amplification. The first pass's diagnosis is the better one, and
it is right.

### §`sec:coin`'s reconstruction, and the `W`-versus-`V` diagnosis

| `S(h)` reading | denominator | first pass | here |
|---|---|---|---|
| `M(h)/X` | `W` (v1's) | -0.09768 | **-0.09768** |
| `M(h)/X` | `V` (correct) | -0.12413 | **-0.12413** |
| `M(h)/(X-h)` | `W` | -0.04080 | **-0.04080** |
| `M(h)/(X-h)` | `V` | -0.05185 | **-0.05185** |

and the supporting quantities at `X = 4e6`:

| | first pass | here |
|---|---|---|
| `W(X) = sum_w Lambda(w)^2` | 5.67847e7 | **5.678467e7** |
| `V(X) = sum_v mu^2(v)Lambda(X-v)^2` | 4.46842e7 | **4.468418e7** |
| `W/V` | 1.27080 | **1.27080** |
| `1/A(N)` predicted | 1.27021 | **1.27020** |

So `v1`'s `-0.0976` does divide by `W` where the definition calls for
`V`, and the error is exactly the local factor `A(N)` — the object
`prop:V` of the same paper is about. The whole block is right.

*(A note on `A(N)` for this `X`: `4e6 = 2^8 · 5^6`, so both the `q=2` and
the `q=5` factors come out of the Artin constant. Using the even-`N`
value alone gives `1.337` and misses. This tree made that slip in a first
run of the check and it is recorded here because it is the same species
of error — a local factor evaluated at the wrong modulus — that the
correction is about.)*

### Proposition 20's self-cancellation

The correction's claim is that `prop:coh` computes `Q_cc/n_c^2`, which is
2x to 61x larger than the error bar `Var` it is about, and that the ratio
is a structural constant rather than noise. The first pass measured it at
`X = 2^22`. This tree's exact cell floor, built independently for
`sec:floor`, gives it at the top band of `X = 1.6e7`:

| depth | `Q_cc/n_c^2` here | `Var` here | ratio | first pass at `2^22` |
|---|---|---|---|---|
| 0 | 1.236e-1 | 1.398e-2 | **0.1131** | 0.113 |
| 1 | 1.212e-1 | 1.976e-3 | **0.0163** | 0.0163 |
| 3 | 1.842e-1 | 5.314e-2 | **0.2885** | 0.289 |
| 5 | 3.071e-1 | 1.694e-1 | **0.5514** | 0.557 |

Four for four, at an `N` four times larger. The correction's claim that
the ratio is structural is confirmed by the fact that it reproduces
across that gap.

### `conj:wall` item 4

The correction replaces the white null with coin and rotation nulls.
This tree reached the same place independently in Phase 1 (M13),
without seeing it: 200 frequency-redrawn surrogates give a maximum of
`3.64e-3` against `mu`'s `4.19e-3`, and 12 independent coin fields give
mean `3.04e-3` with `mu` at `z = +0.59` inside them, 3 of 12 above it.

One thing this tree adds. The first pass shows the paper's null is a
value permutation. This tree **identified it numerically**: 200
permutation surrogates give a maximum of `5.62e-6`, against the paper's
quoted `5.09e-6`. So the paper's number is not merely of the wrong kind;
it is that specific null, and the identification is checkable.

## The one qualified verdict

### `conj:wall` item 3's tail counts

The correction puts in `21441/21463` at `t=3`, `502/503.6` at `t=4`,
`4/4.6` at `t=5`.

**The expected values are right and they pin the field.** `21463`,
`503.6` and `4.6` are `2*Phi(-t)*n` at `n = 7,950,000`, which is exactly
`1e5 < N <= 1.6e7` — the field on which `prop:V`'s three figures
reproduce. This tree computed the same three expectations to the same
digits. So the two passes are measuring the same `N`.

**The observed counts are not reproducible under the operation the paper
states.** Sweeping the treatment on that field:

| cell index | operation | t=3 | t=4 | t=5 |
|---|---|---|---|---|
| **`{3,5,7,11,13}`** — as `sec:floor` states | **means only** | **22086** | **664** | **27** |
| `{3,5,7,11,13}` | standardised | 21601 | 529 | 6 |
| `{3,…,19}` | means only | 21837 | 576 | 6 |
| `{3,…,23}` ("depth 8") | means only | 21769 | 550 | 5 |
| `{3,…,23}` | standardised | 21468 | 518 | 3 |
| `{3,5,7,11,13}` × octave | means only | 21678 | 552 | 7 |
| **`{3,5,7,11,13}` × octave** | **standardised** | **21425** | **506** | **5** |
| `{3,…,23}` × octave | means only | 21485 | 505 | 3 |
| — | **first pass** | **21441** | **502** | **4** |

The first pass's counts sit among the rows that normalise **within each
octave** — `21425/506/5` and `21485/505/3` bracket `21441/502/4`. They
are nowhere near the row the paper's own words describe, which gives
`22086/664/27`, i.e. tail ratios `1.029 / 1.319 / 5.92`.

**Verdict.** The correction's numbers are sound *for the treatment it
used*, and per-octave normalisation is defensible — the first pass's own
self-correction note warns that pooling octaves is a scale mixture, and
it is right about that. But two things follow:

1. The corrected paper still reads as though "removing cell means alone",
   with cells indexed by `{3,5,7,11,13}`, produced these counts. It did
   not. The correction inherits the undefined-statistic problem it was
   written to fix — the same disease as the first pass's own findings 13
   and 14.
2. Normalising within octaves and within cells answers a **different
   question** from the one item 3 poses. Item 3's stated reason for
   caring is that "`C(N)=o(N)` constrains every `N`", which is a claim
   about the aggregate field. On the aggregate field, under the stated
   index, the tail is heavy: 27 exceedances at `t=5` against 4.6
   expected, and 47 of the 50 largest `|G|` at depth ≥ 3 against 1.5
   expected.

So this tree's Phase 1 finding M3 and the first pass's correction are
both arithmetically right. They are answers to different questions, and
neither the paper nor the corrected paper says which one it is asking.

## The definitional difference

### §`sec:c3`'s `0.933 … 0.977`

The first pass reads `A_j(a) = #{(m_1..m_j) : prod m_i = a, m_i <= z,
mu(m_i) != 0}` — with a squarefreeness condition on each factor. This
tree read it without that condition. The two give different tables:

| `M` | | `J=3` | `J=4` | `J=6` | `J=8` |
|---|---|---|---|---|---|
| 1e5 | first pass | 0.933 | 0.942 | 0.958 | 0.961 |
| 1e5 | here | 0.962 | 0.967 | 0.966 | 0.966 |
| 1e6 | first pass | 0.952 | 0.961 | 0.970 | 0.977 |
| 1e6 | here | 0.976 | 0.981 | 0.983 | 0.982 |

Both readings support the closure — "about 95% of the identity's weight
lies outside the region the classification covers, increasing with `M`" —
and both passes independently report that the paper's third decimal is
not supportable. The first pass's reading is the better one: Heath–Brown's
identity does carry `mu(m_i)`, so the absolute weight should count only
squarefree factors. This tree's reading is the looser of the two and its
numbers should be treated as the weaker check.

Neither pass disputes the other's conclusion here. What both establish is
that the paper's table cannot be reproduced to three decimals by a reader,
because `A_j`, `D_j` and `z` are undefined.

## Not checked

- §`sec:R4`'s `z = -1.03 … -0.86` — the sign-randomisation null was not
  run here. This tree ran a permutation null on the lag-1 statistic
  instead, and confirmed the paper's error bar there.
- K1's `13 of 63` and `0.904`, R2's `G_1 = 1.513 ± 0.187`, R1's six-draw
  versus forty-draw spread — all in the supply-side kill-test rows this
  tree declared uncovered in Phase 1. They remain single-witness.

**That is the residual risk.** Three of the first pass's corrections,
including the one that reopens K1 from "dead" to "open", still have
exactly one witness.

## Reproduction

```
python v1_verify2/code/verify/audit_first_pass_corrections.py
python v1_verify2/code/wall/audit_lem_mp.py
python v1_verify2/code/wall/audit_cell_floor.py
python v1_verify2/code/wall/audit_zero_coin_dist.py
```

---

## Is the corrected paper clean?

No. `code/verify/lint_corrected_vs_findings.py` matches each of this
tree's Phase 1 findings against the TeX source of
`wall_v1_corrected.tex`, in the LaTeX form the paper writes it in.

**Sixteen of seventeen survive unrepaired, five of them high severity.**

| id | sev | in the corrected paper |
|---|---|---|
| A1 | high | `$0.2238\pm0.0056$ --- forty standard errors` — verbatim |
| M2 | high | `Excess kurtosis $-0.0005$ ($z=-0.3$) … on $6.3\cdot10^6$ values, removing cell \emph{means} alone` — verbatim |
| M3 | high | the `extremes are attained at generic $N$` claim — verbatim (the tail *counts* were replaced; this sub-claim was not) |
| M14 | high | the whole decay table, `0.6289`/`0.0121`/`52.0` included |
| M14b | high | `A common exponent is rejected at $\chi^2/\mathrm{dof} = 251$` |
| A2 | med | `with the steps at $5$ to $30$ standard errors` |
| M7 | med | `the margin at $N=10^8$ is a factor $N^{0.454}$` |
| M4 | med | `8.40` at `q=3` and `15.16` at `q=5` |
| M5 | med | `1.051$--$1.068 times … $\sqrt{0.32264\,(X-h)}$` |
| A4 | low | `the C4 threshold of $0.5\times$ … $8.8\%$` |
| A6 | low | `squares itself away, and the surviving Möbius sits on the short variable` |
| A7 | low | `$\gg N$ by Parseval` |
| A8 | low | `0.841` |
| A9 | low | `$N=2^{14},\dots,2^{20}$ --- below $1$ and decaying` |
| M10 | low | `by a fixed permutation of the label set` |
| M6 | low | the five shift shares, still at an unstated `X` |
| **A3** | low | **repaired** — "ten kill-tested" is gone |

The one repair is the one finding both passes made.

This is the honest reading of `RECALL.md`'s recall figure, stated as a
property of the artifact rather than of the process: **the corrected
paper is `v1` with one pass's findings removed.** It is better than
`v1` — the fifteen repairs checked above are real and six of them are
confirmed exactly — but it is not a verified paper, and its verification
stamp covers only what one pass happened to look at.

Anything built on it should treat `conj:wall` items 1, 2 and 3 and
`sec:floor`'s decay table as open.
