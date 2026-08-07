# v1_verify2 — Phase 1, blind recall

Re-verification of `v1/paper/wall_v1.tex` written from the **statements**,
with `v1_verify/paper/ADVERSARIAL_FINDINGS.md` and `v1_verify/code/`
unopened throughout. `v1/` was not modified (`git status v1/`: clean).

## Sealed target

Recorded before starting and re-verified after, per `v1_verify2/README.md`:

```
5d359db38b00938e8b8de5e1e2537f4e80cdbb825165843a9621339e4db1ff0c
    v1_verify/paper/ADVERSARIAL_FINDINGS.md
454e92e034814525d7e5f3c85013c031888821e5b8aed356686ddadaa034964a
    v1_verify/paper/wall_v1_corrected.tex
```

Both match the README at both times. Scoring is against a fixed target.

## Coverage

This pass covers the wall, the demand side, the supply-side consumable
E1, `sec:R4`, `sec:c3`, and the Huang--Li source check. The residual
uncovered list is at the end and is part of the result: a recall figure
computed against this pass must be computed only over what this pass
actually looked at. The two substantial gaps are Table `tab:L`'s
dilate-pair rows and the five route adjudications.

## The findings in one place

**Where the paper is wrong in a way that matters**

| | statement | what is wrong |
|---|---|---|
| M1 | `lem:MP` | not an identity; RHS/LHS = 1.566. Repair verified |
| M2 | `conj:wall` item 1 | kurtosis `z = +8.4`, not `-0.3`, under the stated operation |
| M3 | `conj:wall` item 3 | `t=5` tail is `5.9x` Gaussian; extremes are at deep radicals, not generic `N` |
| M13 | `conj:wall` item 4 | `mu` sits `z=+0.59` inside the coin distribution; the null was computed on a shuffled field and is `648x` too tight |
| M14 | `sec:floor` decay table | depths 0, 1, 2 have no detectable mask; no step exceeds `2.5` s.e. |
| A1 | `conj:wall` item 2 | gap is `5.9 sigma`, not 40 |

**Where the paper is wrong in a way that does not change a conclusion**

M4 (major-arc factors), M5 (the autocorrelation null), M6 (shift-mass
shares), M7 (`N^{0.454}`), M8 (block-ratio weight), M9 (`M^{0.05}`),
M10 (`lem:placebo`'s wording), A2 (`5 to 30 s.e.`), A3 (closure count),
A4 (`8.8%`, and "C4"), A5 (`prop:W`'s Chowla step), A6 (`thm:A`'s
one-liner), A7 ("by Parseval"), A8 (`rem:rho` at `1e8`), A9, A10.

**Where the paper reproduces** — see the table below, and M11, M12, and
the Huang--Li check. The defect report against the published source
holds.

## What reproduced

Stated first, because a pass that only reports defects gives no evidence
that its own pipeline is right. Each of these was recomputed from
scratch and matches:

| Statement | Paper | Here |
|---|---|---|
| `prop:V` residual sd, candidate `A` | `0.000323` | `0.000323` |
| `prop:V` residual sd, candidate `S` | `0.245235` | `0.245235` |
| `prop:V` ratio of the two | `760` | `759.3` |
| `prop:V` mean of `(V/W)/A`, top octave | `1.000000` | `1.000000` |
| `prop:E` margin at `2^14, 2^16, 2^18, 2^20` | `0.168, 0.175, 0.158, 0.152` | `0.1679, 0.1749, 0.1576, 0.1521` |
| `sec:coin` `mu`-autocorrelation vs its floor | `1.051–1.068` | `1.051–1.069` (rms, five decades) |
| `sec:floor` `chi^2/dof` for a common exponent | `251` | `251.5` |
| `sec:floor` the `a_d/s.e.` column | six values | all six |
| `prop:coh` remark: `11.8`, `1.21`, `~10x`, `N^0.46` | — | all four |
| `sec:closures` R1 row (`z=-0.80`, `7.5%` at 3 s.e.) | — | both |
| supply side: the `8.8e-6` demand identifies `A=1` | — | `8.685e-6` |
| `prop:Dpp` piece ratio at `1e6, 4e6, 1.6e7` | `0.771, 0.790, 0.810` | `0.7712, 0.7903, 0.8103` |
| `prop:Dpp` `CP_2/(N log N)` | `2.886, 2.949, 2.997` | `2.8864, 2.9493, 2.9966` |
| `prop:Dpp` mean-zero tuning shift | "about five percent" | `-5.09%, -4.54%, -4.09%` |
| `sec:R4` exact identity `= mu(N-1)` | at 2 `N` | confirmed at 6 `N` |
| `sec:R4` lag-1 autocorrelation | `+0.0104, +0.0127`, s.e. `0.0112` | `+0.0055`, perm. null sd `0.0119` |
| `lem:cellmom` itself | stated | confirmed vs 60 sign draws |
| `conj:wall` item 2, `max_c\|z_c\|` top band | `8.4` | `9.1`, Bonferroni cleared in all 8 octaves |
| `prop:coh` `b` at the three shallowest depths | `0.0379, 0.0378, 0.0379` | `0.0395, 0.0397, 0.0394` |
| `sec:floor` Rarity: deepest-cell share / pooled | `0.94` / `0.018` | `0.913` / `0.011` |
| `sec:c3` weight outside `a <= M^{0.05}` | `~0.95` | `0.93`–`0.98` |
| `sec:c3` per-`j` share, `J=8`, `j in {6,7,8}` | `0.847` | `0.845` at `M=1e5` |
| Huang--Li eq. (18) defect report | stated | **confirmed against the source** |
| Huang--Li eq. (22), and Corollary 1 | quoted | both confirmed verbatim |

The **defect report is the paper's headline claim about the literature,
and it holds.** In `huangli_text.txt` the third display of the `S_2`
derivation carries the `n`-dependent constraint `k < (N-n)/alpha`, and
equation (18) two lines later carries `k < (N-1)/alpha`. The omitted
terms are exactly those with `d = (N-n)/k <= alpha`, which is the range
the paper's `Delta` sums over, with the sign coming from
`log(k/(N-n)) = -log m`. Huang--Li's equation (22) is verbatim the
identity of `thm:C`, and HL state the equivalence themselves; the paper
attributes it to them.

The `prop:V` block matching to **six significant figures on all three
numbers** is the load-bearing calibration: it establishes that `C(N)`,
`V(N)`, `W(N)`, `A(N)` and the field are computed here exactly as the
paper computes them. The divergences below are therefore not pipeline
differences.

One correction comes with it: those three figures reproduce on the field
`10^5 < N <= 1.6e7` (7,950,000 values) and **nowhere else**. The paper
says "over every even `N <= 1.6e7`". On the full field the `A` figure is
`0.000582` and the ratio is `421`, not `760`. The `S` figure is stable
across every cutoff, so the sentence's two numbers cannot both be placed
on the field it names.

---

# Findings

`A`-series: settled from the paper's own text and arithmetic.
`M`-series: settled by independent measurement.

## M1. `lem:MP` is not an identity

**Severity: high.** The paper calls it "exact and unconditional" and
cites it in the abstract and the Summary.

Expanding `C(N)^2` over `(n,v)` and `(n',v')` with `n+v = n'+v' = N` and
summing over `N <= X` leaves three free indices subject to `n + v <= X`.
The `Lambda`-pair `(n, n-h)` and the `mu`-pair `(v, v+h)` stay tied by
that **simplex** constraint. The proof's step "collecting the
`Lambda`-pairs and the `mu`-pairs separately gives the two factors"
replaces the simplex by the **box** `{n<X} x {v<X}`.

Brute force, `code/wall/audit_lem_mp.py`:

| X | LHS | RHS as stated | RHS/LHS |
|---|---|---|---|
| 200 | 1.0955e5 | 2.2635e5 | 2.0661 |
| 800 | 2.4826e6 | 3.8187e6 | 1.5382 |
| 3200 | 4.1930e7 | 6.5644e7 | 1.5656 |

Both readings of "below X" give the same numbers; the ratio settles near
**1.566**. The simplex-coupled form reproduces the left side to twelve
digits, so the diagnosis is the factorization step and not a coding error.

**Repair, verified to twelve digits at X = 200 … 1600 and again at
X = 1.6e7** (`audit_shift_mass.py`, ratio `1.000000000000`). With the
*truncated* convolution `Chat(N) = sum_{n+v=N, n<=X, v<=X} Lambda(n)mu(v)`,

```
sum_{N<=2X} Chat(N)^2 = sum_{|h|<X} M(h) P(h)
```

— the sum on the left running over **all** `N`, to `2X`. What is not
available is keeping `sum_{N<=X} C(N)^2` there.

The interpretation the paper draws from `lem:MP` (aggregate wall size is
a statement about binary Chowla and Hardy–Littlewood) survives the
repair intact. The word "exact" and the stated form do not.

## M2. `conj:wall` item 1 does not reproduce under the operation it states

**Severity: high.**

Quoted: excess kurtosis `-0.0005` (`z=-0.3`), `E|G|/sd(G)` short of
`sqrt(2/pi)` by `0.00018` (`z=-0.8`), "on `6.3e6` values, removing cell
**means alone**".

Measured on the field where `prop:V` reproduces (`10^5 < N <= 1.6e7`,
7,950,000 values), removing cell means alone:

| cell index | cells | excess kurtosis | z |
|---|---|---|---|
| divisibility by `{3,5,7,11,13}` | 32 | **+0.0146** | **+8.4** |
| valuations `v_q<=2` of `{3,5,7,11,13}` | 233 | +0.0145 | +8.3 |
| valuations `v_q<=3` of `{2,…,13}` | 1717 | +0.0146 | +8.4 |
| residue class `N mod 2310` | 1155 | +0.0288 | +16.6 |
| residue class `N mod 30030` | 15015 | +0.0160 | +9.2 |
| divisibility by `{3,…,19}` | 128 | +0.0076 | +4.4 |

No index tested returns `-0.0005`. Valuation-indexed cells change
nothing; residue classes make it worse. On the paper's literal field
(every even `N`, 8.0e6 values) it is `+0.0252`, `z = +14.6`.

The paper's figure **is** reproduced by a different operation: per-cell
**standardisation** — removing each cell's spread as well as its mean —
which gives `+0.0004`, `z = +0.2` at `{3,…,19}`. That is not "removing
cell means alone", and it matters beyond wording: item 1's headline is
"the bulk is Gaussian, **under this scale and no other**", and per-cell
standardisation is a second scale, applied after `sqrt(V)`.

Two further points on the same item:

- The stated sample size is wrong. Every even `N <= 1.6e7` is
  **8,000,000** values; the field where `prop:V` reproduces is
  **7,950,000**. Neither is `6.3e6`. (`6.3e6` is exactly the count for
  `N > 3.4e6`, which no statement in the paper names.)
- `z = 98` for the `S N`-scaled kurtosis is arithmetically inconsistent
  with `n = 6.3e6`, which gives `87.3`. It is consistent with
  `n = 8.0e6`. Measured here: `+0.1822` at `z = +105` on 8.0e6 — so the
  `S N` comparison itself is confirmed in substance, on a field that is
  not the one the same sentence names for the other z-score.

## M3. `conj:wall` item 3: the tail is not Gaussian, and the extremes are at deep radicals

**Severity: high.** This is the item the paper itself calls "the half
that matters, since `C(N)=o(N)` constrains every `N` and a bulk can
match to five decimals while the tail is heavy."

Quoted: "aggregate tail counts against the Gaussian expectation give
ratios `0.999` at `t=3`, `0.997` at `t=4`, `0.878` at `t=5`", and "the
extremes are attained at generic `N`, not at deep radicals, so the mask
removal is not leaking into the tail."

Measured, same field, cell means removed, `{3,5,7,11,13}`:

| t | expected | measured ratio | paper |
|---|---|---|---|
| 3 | 21,463 | **1.029** | 0.999 |
| 4 | 503.6 | **1.319** | 0.997 |
| 5 | 4.56 | **5.924** (27 events) | 0.878 |

27 events against 4.6 expected is not a Poisson fluctuation. And the
extremes claim fails badly. Of the 50 largest `|G|`:

| field | depths of the 50 largest `\|G\|` (d=0..5) | expected if generic |
|---|---|---|
| every even `N` | `[0, 2, 1, 26, 20, 1]` | `[19.2, 21.1, 8.2, 1.4, 0.1, 0.0]` |
| `10^5..1.4e7` | `[9, 6, 5, 14, 16, 0]` | `[19.2, 21.1, 8.2, 1.4, 0.1, 0.0]` |
| top octave | `[20, 13, 8, 4, 3, 2]` | `[19.2, 21.1, 8.2, 1.4, 0.1, 0.0]` |

Depths `>= 3` are 3.1% of the field. On the paper's literal field they
supply 47 of the 50 largest `|G|` against 1.5 expected. Even restricted
to the top octave — the most favourable window — they supply 9 against
1.5.

**Mechanism, and why removing means cannot fix it.** The mask has a
*variance* component, not only a location component. Top octave, after
its own mean is removed from each depth:

| depth | count | mean `Z` | sd `Z` | sd / pooled |
|---|---|---|---|---|
| 0 | 1,534,466 | +0.0059 | 0.9139 | 0.988 |
| 1 | 1,687,911 | +0.0314 | 0.9129 | 0.987 |
| 2 | 654,281 | -0.0072 | 0.9246 | 1.000 |
| 3 | 114,017 | -0.3572 | 1.0243 | 1.108 |
| 4 | 9,059 | -1.4121 | 1.2652 | **1.368** |
| 5 | 266 | -3.7294 | 1.1486 | 1.242 |

A scale mixture with unequal component sds has positive excess kurtosis
and a heavy tail by construction, and subtracting means does not touch
it. The paper's own "Rarity" paragraph states the principle exactly —
"the mask is rare, not small … pooling across cells is what makes it
look negligible" — for the location mask, and then reports the tail as
clean.

**The quoted `t=5` figure carries no information either way.** On a
field of this size the expected count at `t=5` is `4.6`, so the Poisson
spread of any `t=5` ratio is about `+-47%`. The paper's `0.878` is
`4` events against `4.6` — indistinguishable from `1.000`, and quoted to
three significant figures with no spread. That is the paper's own
rule 3 ("a bar quoted as an effect size is not a bar").

**Bonus, and it is a repair rather than a defect.** The mask does not
stop at 13. Widening the index from `{3,5,7,11,13}` to `{3,…,19}` cuts
the `t=5` count from 27 to 6 and the kurtosis `z` from `+8.4` to `+4.4`.
Whatever indexes the mask reaches at least 17 and 19.

## M4. `sec:coin`'s major-arc factors are single-`N` draws of a quantity that swings by 5x

**Severity: medium.**

Quoted: "the Möbius exponential sums are markedly smaller than a coin's
--- by a factor `8.40` at `q=3` and `15.16` at `q=5`."

The ratio `(coin scale)/|S_mu(j/q)|`, rms over `j` coprime to `q`,
scanned over `N = 2^14 … 2^24` (`code/wall/audit_majorarc_scan.py`):

| N | `q=3` | `q=5` |
|---|---|---|
| 2^16 | 9.22 | 5.17 |
| 2^18 | 3.75 | 1.99 |
| 2^20 | 2.40 | 2.76 |
| 2^21 | 9.66 | 5.26 |
| 2^22 | 5.18 | 7.26 |
| 2^24 | 4.50 | 2.73 |

At `q=3` the ratio ranges over `2.06 … 9.66` across eleven values of
`N`; at `q=5`, `1.78 … 7.26`. `8.40` is attainable at `q=3` at some `N`;
`15.16` at `q=5` is not attained anywhere in the range. The two figures
are quoted to three significant figures with no `N` and no spread, and
the coin convention is not stated either — rms and mean-modulus differ
by `sqrt(pi/2) = 1.2533`.

**The direction survives.** The ratio exceeds 1 at every `N` and every
small `q` measured, so "markedly smaller than a coin's" holds and the
major-arc mechanism behind `rho < 1` is not in question. Only the two
numbers are.

## M5. `sec:coin`'s `mu`-autocorrelation excess is an artifact of an `h`-independent floor

**Severity: medium.** This one both reproduces and dissolves.

Quoted: "The Möbius autocorrelation sits at `1.051`--`1.068` times the
random-sign floor `sqrt(0.32264 (X-h))` --- not `sqrt X`, since the sum
sees only `n` with both `n` and `n+h` squarefree --- stably across five
decades of shift."

Reproduced exactly. rms of `|M(h)|/sqrt(0.32264(X-h))` by decade:
`1.0514, 1.0584, 1.0656, 1.0689, 1.0690` over the five decades from
`10^2` to `8·10^6`.

But the floor is wrong, and in a way the paper's own parenthetical shows
it was half-seen. The density of `v` with `v` and `v+h` both squarefree
is `prod_{p^2 | h}(1-1/p^2) * prod_{p^2 not| h}(1-2/p^2)` — it depends
on `h`. `0.32264` is the value at generic `h` only. Splitting:

| | count | rms `\|M\|`/floor |
|---|---|---|
| `h` odd | 100,000 | **1.0055** |
| `h = 2 mod 4` | 50,000 | **1.0099** |
| `4 \| h` | 49,999 | **1.2265** |

`sqrt(1.5) = 1.2247`: for `4 | h` the true density is `1.5x` the generic
one, because the local factor at `p=2` is `1-1/4` rather than `1-2/4`.
Mixing the three classes reproduces the pooled `1.066` to three decimals.

So against the **correct** floor the Möbius autocorrelation sits at
`1.006` (h odd) and `1.010` (`h = 2 mod 4`) — essentially exactly at the
random-sign scale. The reported `5–7%` excess over random signs is the
mis-specified null, not a property of `mu`. The paper's rule ("null
before threshold") applies to its own floor.

## M6. `sec:coin`'s shift-mass table reproduces, at an `X` it does not state

**Severity: low.** This finding replaces an earlier, wrong one of mine;
the retraction is recorded at the end of this file.

Quoted: "`h<10^3` carries `1.1%`, `10^3-10^4` carries `3.0%`,
`10^4-10^5` carries `23.1%`, `10^5-10^6` carries `48.9%`, and above
`10^6` carries `23.8%` (the ranges cancel, net `-2.2e13` against gross
`4.4e13`)."

"Gross mass" is not defined, and the two readings differ by three orders
of magnitude. The paper's own totals settle it: with net `-2.2e13` and
gross `4.4e13` the ratio is `2`, so "gross" is the sum of the **absolute
bucket nets**, not the term-by-term `sum_h |M(h)P(h)|` (whose ratio to
the net is in the hundreds). Under that reading, at `X = 4e6`:

| bucket | measured | paper | per-`h` reading |
|---|---|---|---|
| `h < 1e3` | **1.1%** | 1.1 | 0.1% |
| `1e3 - 1e4` | **3.0%** | 3.0 | 0.6% |
| `1e4 - 1e5` | **23.1%** | 23.1 | 5.5% |
| `1e5 - 1e6` | **48.9%** | 48.9 | 45.1% |
| above `1e6` | **23.8%** | 23.8 | 48.8% |

with net `-2.2186e13` against the quoted `-2.2e13`, and gross/net
`1.99` against the quoted `2.0`. **Every figure reproduces.**

**What remains a finding.** The table is specific to `X = 4e6`, which
the paper never states, and the top bucket is unbounded above so its
share moves with `X`. At `X = 1.6e7` the same computation gives

| bucket | `h<1e3` | `1e3-1e4` | `1e4-1e5` | `1e5-1e6` | `>1e6` |
|---|---|---|---|---|---|
| share | 0.1% | 8.5% | 4.1% | **84.5%** | 2.7% |

— a different table entirely, and the cancellation ratio falls from
`1.99` to `1.06`. So the five percentages are a property of one
unstated `X`, and the reader cannot tell which.

**The conclusion is unaffected and holds at both `X`.** "Small shifts,
where Chowla is hardest and the averaged theorem weakest, carry almost
nothing" is true at `X = 4e6` (`1.1%` below `10^3`) and more so at
`X = 1.6e7` (`0.1%`).

One methodological note stands regardless: the paper apportions **gross**
mass in a sum that cancels, while `rho-1` is built from the **net**, and
the two disagree about which bucket dominates.

## M7. `sec:margin`: `N^{0.454}` is contradicted by measurement and by the paper's own formula

**Severity: medium.** Logged from the text as A7 before measuring; now
measured.

Quoted: "measured, `max|C|/N` falls from `0.056` to `0.0082` over the
range computed, and the margin at `N=1e8` is a factor `N^{0.454}`."

Measured `max|C|/N` by octave: `0.1135` at `3.1e4–6.3e4` falling to
`0.01007` at `8e6–1.6e7`. Over `N <= 1.4e7` the max of `|C|` is
`110,137`, i.e. `0.00787` of `1.4e7` — so the endpoint `0.0082` is
recognisably the paper's, and the fitted decay is `max|C|/N ~ N^{-0.433}`.

But the **margin** at the top of the measured range is `1/0.01007 = 99`,
which is `N^{0.277}`, and extrapolating the fitted decay to `N=1e8`
gives `249 = N^{0.30}`. The paper's own Gumbel formula gives
`490 = N^{0.336}`. `N^{0.454}` is `4285` — an order of magnitude above
all three. The same formula reproduces the paper at `N=1e12`
(`10^{4.50}` vs quoted `10^{4.4}`) and `N=1e50` (`10^{22.86}` vs quoted
`10^{22.8}`), so the formula and the reading of it are right and
`N^{0.454}` is the outlier.

## M14. `sec:floor`'s decay table reports exponents for cells with no detectable mask

**Severity: high.** The table carries the section's central claim ("the
exponent rises monotonically as the cell gets shallower"), the
`chi^2/dof = 251` rejection of a common exponent, and the paper's own
statement that the mask's decay is unexplained.

Using the exact floor of `lem:cellmom` — the floor the paper says it
used, and which this pass confirmed against 60 independent-sign draws —
here is `delta_c` and its `z` at every octave, cells = the six depths:

| octave top | d=0 | d=1 | d=2 | d=3 | d=4 | d=5 |
|---|---|---|---|---|---|---|
| 62,500 | +0.150 / **1.0** | +0.172 / 3.1 | -0.386 / -2.2 | -1.983 / -6.9 | -4.367 / -10.5 | — |
| 250,000 | +0.076 / **0.6** | +0.114 / 2.2 | -0.180 / -1.1 | -1.361 / -5.0 | -3.803 / -10.1 | -6.584 / -11.0 |
| 1,000,000 | +0.030 / **0.2** | +0.085 / 1.7 | -0.095 / -0.6 | -0.875 / -3.4 | -2.784 / -7.9 | -6.261 / -13.0 |
| 4,000,000 | +0.013 / **0.1** | +0.052 / 1.1 | -0.034 / -0.2 | -0.580 / -2.4 | -1.995 / -6.0 | -4.817 / -11.0 |
| 16,000,000 | +0.005 / **0.0** | +0.031 / 0.7 | -0.008 / -0.1 | -0.358 / -1.6 | -1.413 / -4.5 | -3.730 / -9.1 |

Octaves (of nine) in which `|z| >= 3`:

| depth | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| octaves with `\|z\| >= 3` | **0** | **1** | **0** | 5 | 9 | 7 of 7 |

**At depths 0 and 2 the mask is not detected at any scale**, and at
depth 1 only in the smallest octave. Depth 0's amplitude never exceeds
`1.0` standard errors of the paper's own exact floor.

Refitting `|delta_d| ~ N^{-a_d}` by the procedure `sec:floor` states —
weighted by the exact errors, exponent s.e. from the fit covariance —
under both cell readings:

| depth | `a_d` (6 cells) | s.e. | `a/se` | `a_d` (32 pooled) | s.e. | `a/se` | paper `a_d` | paper s.e. | paper `a/se` |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.1344 | 0.0265 | 5.1 | 0.1344 | 0.0265 | 5.1 | 0.1434 | 0.0155 | 9.3 |
| 4 | 0.1940 | 0.0264 | 7.3 | 0.2098 | 0.0151 | 13.9 | 0.2152 | 0.0065 | 33.1 |
| 3 | 0.2995 | 0.0570 | 5.3 | 0.2742 | 0.0238 | 11.5 | 0.2713 | 0.0040 | 67.8 |
| 2 | 0.5249 | 0.3544 | 1.5 | 0.2502 | 0.0590 | **4.2** | 0.3686 | 0.0052 | 70.9 |
| 1 | 0.2898 | 0.1250 | 2.3 | 0.3375 | 0.1623 | **2.1** | 0.0437 | 0.0556 | 0.8 |
| 0 | 0.5322 | 0.7688 | 0.7 | 0.5322 | 0.7688 | **0.7** | 0.6289 | 0.0121 | 52.0 |

Three things follow.

**(a) Where the mask exists, the paper's exponents are right.** At
depths 5, 4, 3 the central values agree closely — `0.134` vs `0.143`,
`0.210` vs `0.215`, `0.274` vs `0.271`. The disagreement is not in the
exponents but in their errors, which are `2x` to `100x` larger here.

**(b) Depths 0, 1, 2 are not measurable** (pre-registered bar:
`|a/se| >= 5` under some reading). The paper assigns depth 2 a
significance of `70.9` and depth 0 a significance of `52.0` — its second
and third largest — for cells whose amplitude never reaches `2.2` and
`1.0` standard errors respectively at any octave. An exponent fitted to
an amplitude consistent with zero is fitted to noise, whatever its
nominal error.

**(c) The monotonicity claim does not survive.** Recomputing the steps
the paper describes, along its own ordering with depth 1 excluded:

| step | 6 cells | 32 pooled |
|---|---|---|
| 5 -> 4 | +1.59 | +2.47 |
| 4 -> 3 | +1.68 | +2.28 |
| 3 -> 2 | +0.63 | **-0.38** |
| 2 -> 0 | +0.01 | +0.37 |

No step reaches `2.5` standard errors under either reading, and under
the 32-cell reading the `3 -> 2` step goes the **wrong way**. The
paper's "steps at 5 to 30 standard errors" (already wrong on its own
table, finding A2, where the true range is 4.3 to 19.8) becomes
`-0.4` to `+2.5` when the exponents are given the errors the exact floor
supports.

The `chi^2/dof = 251` rejection of a common exponent also goes:

| reading | weighted common exponent | `chi^2` (5 dof) | `chi^2/dof` | verdict |
|---|---|---|---|---|
| 32 pooled | 0.2128 | 16.61 | **3.3** | rejected, `p ~ 0.005` |
| 6 depth cells | — | 9.56 | **1.9** | **not rejected**, `p ~ 0.09` |
| paper | 0.3018 | 1257 | 251 | rejected |

So under the cell definition `sec:floor` itself states — the six depths —
a single common decay exponent for the mask **cannot be rejected at
all**.

**What survives.** The mask itself is real and large — `conj:wall`
item 2's Bonferroni claim reproduces (`max_c|z_c| = 9.1` in the top
band, cleared in all eight octaves), driven by depths 4 and 5. That the
mask **decays** is established at depths 3, 4, 5. `prop:scaleinv`'s
conclusion — that the size mechanism cannot explain the decay — is
untouched, since it concerns exponents of `0.13`--`0.27` that are
measured. What does not survive is the depth-ordering of the exponents
and the claim that the ordering is significant.

Also worth stating plainly: the paper drops depth 1 as "not measurable"
(`0.0437 +- 0.0556`) and keeps depths 0 and 2. On the exact floor,
**depth 1 is the more detectable of the three** — it is the only one of
`{0,1,2}` that reaches `|z| >= 3` in any octave.

## M13. `conj:wall` item 4 fails the paper's own coin control, and its null is the wrong null

**Severity: high.** Item 4 is the only **positive** claim in
`conj:wall`, and it is the one item for which the paper reports no coin
control — although `lem:coin`, stated two sections later in the same
paper, says in terms that an estimator reproduced under `mu -> eps` "is
not measuring `mu`", and the Methodology section lists the same rule.

Quoted: "Regressing `G` on `cos(gamma log N), sin(gamma log N)` for the
first ten ordinates gives `R^2 = 3.90e-3` against a 200-surrogate
maximum of `5.09e-6`, and every ordinate individually at `z >= 23`."

Measured on the field where `prop:V` reproduces
(`code/wall/audit_zero_spectrum.py`):

- **`R^2 = 4.19e-3`. The effect reproduces.** The paper's `3.90e-3` is
  right.
- **The surrogate maximum does not**, and the paper's null is
  identifiable. Two surrogates, 200 draws each, on the same field:

| surrogate | max | mean | sd |
|---|---|---|---|
| frequencies redrawn, **field intact** | **3.64e-3** | 1.87e-3 | 7.3e-4 |
| **field shuffled**, zeros kept | **5.62e-6** | 2.56e-6 | 7.7e-7 |
| iid expectation `2k/n` | 2.52e-6 | | |
| **paper's quoted maximum** | **5.09e-6** | | |

  The paper's `5.09e-6` sits inside the *permutation* surrogate's own
  spread and is `648x` below the frequency surrogate's. So the null was
  computed on a shuffled field. That is the wrong null here: `G(N)` and
  `G(N+2)` share almost every summand, and shuffling destroys exactly
  the correlation that lets a smooth regressor fit the field. Against a
  null that respects it, `R^2(mu)` sits at **1.15 times the surrogate
  maximum**, not 766 times.
- **`z >= 23` for every ordinate does not reproduce.** Individually:

| gamma | `R^2` | z | coin `R^2` | coin z |
|---|---|---|---|---|
| 14.134725 | 2.41e-3 | **+6.7** | 2.65e-3 | **+7.4** |
| 21.022040 | 8.26e-4 | +1.8 | 1.32e-4 | -0.3 |
| 25.010858 | 1.83e-4 | -0.1 | 2.45e-5 | -0.6 |
| 30.424876 | 2.30e-4 | +0.0 | 2.19e-4 | -0.0 |
| 32.935062 | 1.18e-4 | -0.3 | 2.79e-4 | +0.2 |
| 37.586178 | 1.80e-4 | -0.1 | 9.30e-5 | -0.4 |
| 40.918719 | 9.34e-5 | -0.4 | 4.94e-4 | +0.8 |
| 43.327073 | 3.23e-4 | +0.3 | 7.06e-4 | +1.5 |
| 48.005151 | 6.17e-5 | -0.5 | 2.20e-4 | -0.0 |
| 49.773832 | 4.14e-5 | -0.6 | 8.92e-5 | -0.4 |

  Nine of the ten sit within `2` standard errors of a
  frequency-redrawn null. Only the first ordinate carries anything — and
  **the coin carries more of it than `mu` does** (`+7.4` against `+6.7`).

**The coin control.** Replacing `mu` by random signs on the squarefree
support — `V(N)`, `Lambda` and the cells all untouched, exactly
`lem:coin`'s construction — over **12 independent sign fields**:

```
6.48e-3  6.75e-3  4.75e-3  1.56e-3  2.07e-3  1.65e-3
3.79e-3  2.86e-3  1.69e-3  1.67e-3  1.85e-3  1.36e-3

coin distribution: mean 3.04e-3, sd 1.96e-3, max 6.75e-3
R^2(mu)          : 4.19e-3
```

`R^2(mu)` sits at **`z = +0.59`** inside the coin distribution, and
**3 of 12 coin fields exceed it** (one-sided `p = 0.31`). By
`lem:coin`'s own words — "any estimator whose output is reproduced when
`mu` is replaced by `eps` is not measuring `mu`" — this statistic is not
measuring `mu`.

There is a mechanism for this and it is not subtle: `C(N) = sum_v mu(v)
Lambda(N-v)` contains `Lambda`, and the explicit formula puts
`e(gamma log N / 2pi)` oscillations into smooth averages of `Lambda`. A
regression of `C/sqrt(V)` on `cos(gamma log N)` can read the primes'
zeros without reading anything about `mu`. That is consistent with the
only surviving signal sitting at `gamma_1`, the largest term of the
explicit formula, and with the coin matching it.

**What survives.** The paper's cautious sentence — "*Gaussian in
distribution* and *phase-random in log N* are different statements" —
is fine, and item 4 feeds no other claim in the paper, so nothing else
moves. What does not survive is the evidence. `0.39%` is not a floor
established at `z >= 23` over ten ordinates against a `5.09e-6` null; it
is a `1.15x`-over-surrogate excess carried by a single ordinate, sitting
`0.6` standard deviations inside the distribution of what random signs
produce on the same field.

The claim item 4 makes — that `G` is *not* phase-random in `log N` — may
still be true. This pass shows only that the reported measurement does
not establish it.

## M8. `sec:R4`'s block ratio changes verdict with an unstated weight

**Severity: low — the conclusion survives.**

Quoted: "the block ratios are flat (`B=8`: `0.958` and `1.023` at two
`N`, against `B=1` baselines `0.980` and `0.979`)".

The block ratio is a ratio of two summaries and the paper does not say
its weight, which is the paper's own `H9`. The two natural weightings
disagree, measured at `N=1e8` over the whole non-empty band:

| B | ratio of sums | mean of ratios |
|---|---|---|
| 1 | **1.784** | **1.024** |
| 2 | 0.951 | 1.017 |
| 4 | 0.937 | 1.016 |
| 8 | 0.997 | 1.024 |
| 16 | 0.726 | 1.084 |
| 64 | 0.695 | 1.078 |
| 512 | 0.231 | 1.246 |

Under "ratio of sums" the sequence is not flat at all; under "mean of
ratios" it is flat and its `B=1` baseline is `1.02`, near the paper's
`0.980`. So the paper's own baseline identifies which weighting it used,
and **under that weighting the R4 conclusion reproduces**: flat from
`B=1` to `B=8`, degrading by `B=512` as the paper says.

Two smaller notes. The band at `N=1e8` is `k < sqrt N = 10,000`, i.e.
9,999 values of which 7,200 have non-empty support — not the "8000
values of `k`" the paper calls "the entire band"; `8000` corresponds to
`N = 6.4e7`, and the paper says "two `N`" without naming them.

And the lag-1 diagnostic is confirmed, including its error bar. Measured
`+0.0055` with a 400-draw permutation null of sd `0.01185`, against the
paper's `+0.0104/+0.0127` and quoted s.e. `0.0112`. My pre-registered
suspicion that `0.0112 = 1/sqrt(8000)` is a count-based bar and
therefore wrong was **refuted**: here the count-based bar is right,
because `D(k)` across `k` really are nearly independent. `prop:coh`'s
warning applies to cell means of `Z`, not to this statistic.

## M9. `sec:c3`: the `M^{0.05}` threshold is a near-trivial cut

**Severity: low.** The obstruction the paragraph describes is real and
is independently confirmed; the "95%" figure understates it.

Quoted: "the fraction with `a > M^{0.05}` is `0.939, 0.949, 0.960,
0.947` at `J = 3, 4, 6, 8`, rising to `0.961` and `0.969` at `M = 10^6`.
So about `95%` of the identity's weight lies outside the region the
classification covers."

Measured under the stated reading of `A_j` and `D_j`
(`code/supply/c3_hb_weight.py`): `0.937, 0.941, 0.928, 0.929` at
`M=1e4`; `0.962, 0.967, 0.966, 0.966` at `M=1e5`; `0.976, 0.981, 0.983,
0.982` at `M=1e6`. The paper's four figures sit between my `M=1e4` and
`M=1e5` rows and no `M` is stated for them, but the claim reproduces.

The concentration claim reproduces sharply: the top-`j` share at `J=3`
is `0.833` at `M=1e5` against the quoted `0.824`, and the
`j in {6,7,8}` share at `J=8` is `0.845` against the quoted `0.847`.

**What the fraction actually measures.** `M^{0.05}` is `1.585` at
`M=1e4` and `1.995` at `M=1e6`. So over the entire reachable range the
"covered region" `a <= M^{0.05}` is the single point `a = 1`, and the
quoted fraction is `1` minus the share of that one term. Its closeness
to 1 is a statement about the threshold, not about the identity. The
substantive obstruction — that the weight concentrates in the high-`j`
terms, where `a` may run to `x` — is carried entirely by the per-`j`
shares, which do establish it.

## M10. `lem:placebo` describes an operation that would make it vacuous

**Severity: low-medium** — the lemma is used, and used correctly, in
`sec:floor`; only its statement is wrong.

Quoted: "Permuting the cell labels across `N` --- assigning to each `N`
the label of some other `N`, **by a fixed permutation of the label
set** --- preserves every cell size and leaves the field `Z(N)`
byte-identical. Any cell-indexed statistic that survives the permutation
is a property of the cell sizes and not of the correspondence between
cells and arithmetic."

The sentence names two different operations:

- *a permutation of the label set* — renaming cell 1 to cell 2, and so
  on. This leaves the **partition** unchanged and only renames its
  parts. Every label-symmetric statistic (`max_c |z_c|`, the `chi^2`,
  the multiset of `delta_c`) is then invariant **by construction**, so
  under this reading the test has no power at all and the lemma proves
  that nothing is arithmetic.
- *assigning to each `N` the label of some other `N`* — permuting the
  assignment, `label'(N) = label(pi(N))`. This preserves the multiset of
  labels, hence every cell size, while destroying the cell-to-arithmetic
  correspondence. This is the operation with power, and it is the one
  `sec:floor`'s "Cause" paragraph actually applies.

Only the second makes the lemma non-vacuous, and the phrase "by a fixed
permutation of the label set" contradicts it.

Related, and worth stating since the same paragraph relies on it: the
"independent-sign value `(k-1)/n`" is the special case of `k` cells of
**equal** size, where `1/n_c - 1/n = (k-1)/n`. The depth cells here are
not remotely equal — in the top band they run from `1,534,466` down to
`266` — so `(k-1)/n` is not the independent-sign value for this
partition.

## M11. `prop:coh` is confirmed, and its consequence is larger than the paper says

**Severity: none — this strengthens the paper.**

`prop:coh`'s law is confirmed. Fitting `se ~ N^{-b}` to the exact floor
across eight octaves gives `b = 0.0395, 0.0397, 0.0394` at depths 2, 1,
0 against the paper's `0.0379, 0.0378, 0.0379`, and against the
prediction `1/(2<log N>) = 0.0358`. Depth 5 gives `0.0873` against the
paper's `0.1105`. So the error bar falls like `(log N)^{-1/2}` and not
like `n_c^{-1/2}`, exactly as stated.

The remark's "an interval built from a count is about ten times too
narrow at the top of this range" is a statement about the *shrink rate*
over the range and is correct as such. The *absolute* discrepancy is far
larger. In the top band, comparing the exact floor of `lem:cellmom`
against `sd(Z)/sqrt(n_c)`:

| depth | `n_c` | count-based s.e. | exact s.e. | ratio |
|---|---|---|---|---|
| 0 | 1,534,466 | 0.00074 | 0.1182 | **160x** |
| 1 | 1,687,911 | 0.00070 | 0.0445 | 63x |
| 2 | 654,281 | 0.00114 | 0.1414 | 124x |
| 3 | 114,017 | 0.00303 | 0.2305 | 76x |
| 4 | 9,059 | 0.01329 | 0.3169 | 24x |
| 5 | 266 | 0.07045 | 0.4115 | 5.8x |

which is what `prop:coh` predicts: the exact floor is `~1/log N`
independent of `n_c`, while a count bar falls like `1/n_c`, so the ratio
grows with the cell. A count-based interval is 6x to 160x too narrow
here depending on the cell, not "about ten times".

One caveat on the constant: measured `Q_cc/n_c^2` in the top band is
`0.124` to `0.307` by depth, against `prop:coh`'s estimate
`(6/pi^2)/(A log N) = 0.049`. The proposition claims only the form
`~1/log N`, which holds; the constant in its heuristic is low by a
factor of 2.5 to 6.

## M12. `prop:scaleinv` is sound as stated

**Severity: none.** Recorded because this pass checked it.

`D_c = E_{same,c}[S_2] - E_{all}[S_2]` depends only on the distribution
of the shift `h` in the residue classes mod `3,5,7,11,13` that the cell
fixes. Those classes are fixed by the cell and the shift distribution
within them does not depend on the scale of `N`, so `D_c` is
scale-invariant and predicts exponent zero. That is an argument about
the definition and it goes through; the measured predicted exponents
(`-0.0052 ... -0.0003`) are the numerical confirmation of a statement
that is already sound. Its conclusion — that the mechanism explaining
the mask's size cannot explain its decay — follows.

## A1. `conj:wall` item 2: the gap's standard error is incompatible with its arms

**Severity: high.** The paper singles this out as "the one figure here
that stakes a sign as well as a size", and `sec:floor` is built on it.

Quoted: deep arm `0.9476 +- 0.0293`, shallow arm `0.7238 +- 0.0245`,
gap `+0.2238 +- 0.0056`, "forty standard errors".

The gap value is right. For independent arms
`se(gap) = sqrt(0.0293^2 + 0.0245^2) = 0.0382` — the gap is **5.9
sigma, not 40**. `0.0056` is *smaller than either input*, which requires
the two arms to be correlated at `rho = 0.994`; the paper describes a
`300`-versus-`300` comparison and states no pairing. The other rescue
also fails: if the `+-` are SDs with `n=300` each, `se(gap) = 0.0022`.

One of the three numbers is wrong. The gap itself is not in doubt at
5.9 sigma; the quoted significance is off by a factor of about seven.

*(Not re-measured: the `N ~ 1e8` arm was outside this pass's compute.)*

## A2. `sec:floor`: "steps at 5 to 30 standard errors" — actual range 4.3 to 19.8

Recomputed from the paper's own table, along the paper's own ordering
(depth 1 excluded, as the paper excludes it):

| step | delta | s.e. | z |
|---|---|---|---|
| 5 → 4 | 0.0718 | 0.0168 | **4.27** |
| 4 → 3 | 0.0561 | 0.0076 | 7.35 |
| 3 → 2 | 0.0973 | 0.0066 | 14.83 |
| 2 → 0 | 0.2603 | 0.0132 | **19.76** |

Both ends are wrong. Including depth 1 gives `4.3, 7.4, 14.8, 5.8, 10.3`
— still nothing near 30.

This finding is only that the table is inconsistent with itself.
**M14 goes further**: refitting the exponents against the exact floor of
`lem:cellmom` makes no step exceed `2.5` standard errors and reverses
one of them, so the monotonicity claim does not survive either. The
sequence also has depth 1 removed from its middle, and by M14 the
dropped depth is the more detectable of the shallow three.

## A3. The abstract says "ten kill-tested technique designs"; the body has nine

Abstract: "eighteen pre-registered closures --- five route adjudications
…, **ten** kill-tested technique designs, and three representation-class
experiments". Body: `(5)` + `(9)`, tabulated as K1 K2 K3 K4 R1 R2 R3 R4
R5, + `(3)`.

**"Ten" is wrong**: the table has nine rows and its own heading says
`(9)`.

This finding was originally written as "the total is 17, so eighteen is
wrong". That half is **withdrawn**: the three tables carry `5 + 9 + 3`
rows *plus* C-III open, which is 18 rows, so "eighteen closures" is
defensible if the open class is counted. Only the kill-test count
stands.

## A4. `sec:closures`: the "8.8%" figure, and an undefined referent

"the C4 threshold of `0.5x` sits only `1.29` standard errors below its
null, so noise satisfies it `8.8%` of the time." The one-sided normal
tail at `z=-1.29` is **9.85%**; reproducing `8.8%` needs `z = 1.354`.
Separately, **"C4" is defined nowhere in the paper** — the classes are
C-I…C-IV and the kill-tests are K1–K4, R1–R5.

## A5. `prop:W`: the step from Chowla to `rho -> 1` needs a hypothesis not stated

**Severity: medium** — revised downward by measurement; see below.

`prop:W` writes `rho - 1 = (1/V) sum_{h != 0} c(h) S(h)` with
`c(h) = sum_{p'-p=h}(log p)(log p')` and `S(h) = <mu(u)mu(u-h)>`.
Collecting the off-diagonal of `C(N)^2` by `h = p'-p` forces `S(h)` to
be the average of `mu(u)mu(u-h)` over `{N-p : p, p+h both prime}` — a
shifted-prime-pair set, not `[1,X]`. Chowla, and the averaged theorem of
[MRT15], control the average over **all** `u`. Transferring either to
that sparse set is a statement of the same difficulty as `mu` over
shifted primes, not an application of the cited results.

**What measurement says about it.** I had flagged the paper's own
reconstruction gap ("`-0.0976` against a measured `-0.18`, a factor
`0.54`") as evidence for this. It is not. The exact aggregate identity
gives `rho - 1 = -0.1597` at `X = 1.6e7`, against a directly measured
band value of `-0.14` here and the paper's `-0.18`. So the aggregate
reconstruction **does** close, to within the band-versus-aggregate
difference, and the paper's factor-`0.54` shortfall is a defect in that
particular reconstruction rather than a symptom of the substitution.

The logical gap in the `Chowla => rho -> 1` inference stands on its own,
and is unaffected either way.

## A6. `thm:A`: the one-line mechanism is wrong as written, in two places

**Severity: low-medium — exposition only; the proof is correct.**

`wall_v1.tex`: "switching and writing `u = mk` gives, for squarefree `u`,
`mu(u)mu(k) = mu(m)mu^2(k)` --- … the surviving Möbius sits on the short
variable `m < N^{1-theta'}`."

As written this is false. The `k`-sum runs over `k < K`, so `m = u/k`
runs up to `~N`; for small `k` the surviving Möbius sits on the **long**
variable, the configuration the paper's own intro table says has "no
known machine". The bound holds only after the divisor sum is
*completed*: `sum_{k|u, k<K} mu(k) = 1_{rad(u)|N} - sum_{k|u, k>=K} mu(k)`,
and it is the complementary sum, where `k >= K`, that forces
`m = u/k < N^{1-theta'}`.

`v1/paper/theorem_A.tex` does the completion correctly in its Steps 1–3
(`lem:complete`, then "Write `u=mk` with `k>=K`"), so **the theorem is
unaffected**. But the note's own overview paragraph (its lines 215–227)
repeats the same defective one-liner. The sentence appears twice and
would not reproduce the theorem in either place.

This matters more than a typo because "Möbius on the short variable
versus the long variable" is the paper's central organizing contrast —
the intro table, `sec:R4`'s mirror paragraph, and constraint R4 all
turn on it.

## A7. `prop:E`: "by Parseval" does not cover the step it is attached to

"any bound `sup|S_mu| * ||S_Lambda||_1` is at least
`||S_mu||_2 ||S_Lambda||_1 >> N` by Parseval." Parseval gives
`sup|S_mu| >= ||S_mu||_2 = sqrt(6N/pi^2)`, which is the first inequality
and is correct. Concluding `>> N` additionally needs
`||S_Lambda||_1 >> sqrt(N)`, which is not Parseval — Parseval bounds the
`L^2` norm, and `||S_Lambda||_1 <= ||S_Lambda||_2` runs the wrong way. A
lower bound of the right order for the `L^1` norm of a prime exponential
sum is a nontrivial theorem and is uncited. The conclusion is supported
by the measured margins (which reproduce, above); the attribution is
what fails.

## A8. `rem:rho`: the `1e8` conversion contradicts the remark's own trend

The remark says the three `rho` conventions "agree to `0.75%` at
`N ~ 1.4e7`" and "differ by `10.3%` at `10^5`" — a spread that shrinks
with `N`. It then converts `0.810` at `N ~ 1e8` into `0.841`, a spread
of **3.83%**, five times larger at an `N` seven times bigger.

Measured here, the three conventions spread by `0.86%` at
`10^5..1.4e7`, `0.92%` at `7e6..1.4e7`, and `2.97%` at `1e6..2e6` — so
the remark's `0.75%` at `1.4e7` is confirmed in size and the shrinking
trend is real. That makes the `3.83%` at `1e8` the anomaly.

*(Incidental, not a defect: the stated ordering "as written" holds on
`7e6..1.4e7` here but not on `10^5..1.4e7`, where the first two
conventions swap.)*

## A9. `prop:E`: the margins do not match their stated abscissa

"`0.168, 0.175, 0.158, 0.152` at `N = 2^{14}, ..., 2^{20}`" names seven
exponents and gives four values. Measurement settles it: the four are
`2^14, 2^16, 2^18, 2^20`. The full seven are `0.1679, 0.1711, 0.1749,
0.1629, 0.1576, 0.1513, 0.1521` — so "below 1 and **decaying**" is also
wrong twice, at `2^14 -> 2^16` and at `2^19 -> 2^20`.

## A10. `thm:D`: the exponent changes between statement and consequence

The theorem gives `||b||_1/|B_w| >> exp(c sqrt((1/2+delta) log N))` and
concludes `|C(N)| << exp(c sqrt((1/2) log N)) N (log N)^{-A}`. The loss
factor in the conclusion is *smaller* than the lower bound just proved
for it. The intended reading is presumably `delta -> 0`, but as printed
the second display does not follow from the first.

---

## What this pass did not cover

Part of the result. A recall figure scored against this pass must be
computed only over the region above.

- `conj:wall` item 3's Gumbel arm was run and diverges (`+0.97 +- 0.26`
  over nine octaves here against `+0.54 +- 0.45` over eight), but the
  octave set is not stated so this is not recorded as a finding.
- **Table `tab:L` and `conj:L`.** The rows about the prime-indexed
  dilate-pair field `C_{k,k'} = sum_{p~P} mu(N-pk)mu(N-pk')` — pair
  statistics, the `(v_2,v_3)` exact cells, the blind mask stamp
  (`corr = 1.0000`, amplitude error `1.5%`), the Wishart pair matrix
  (`z = -0.19`) — were not measured; this pass never built that field.
  The E1 rows are consistent with what was measured here: the `1e9`
  band ratios `0.966/0.950/0.922` are the same statistic as the
  `sum|D|^2/sum supp` computed at `N=1e8`, which came out `0.96`–`1.17`.
- **The five route adjudications of `sec:closures`**, and the kill-test
  rows K1–K4, R1–R3, R5. These are readings of source papers against
  their verbatim lemma hypotheses, not measurements. Only the Huang--Li
  claims were checked against the source (and all hold).
- `v1/paper/e1_proof.tex` entirely; `theorem_A.tex` was read around the
  mechanism, the divisor-switch steps, and the `log k` branch (A6), but
  its numerics were not re-run.
- The `N ~ 1e8` arms of `conj:wall` item 2 and `rem:rho` (finding A1 is
  arithmetic on the quoted figures, not a re-measurement).
- `sec:margin`'s extrapolations to `N = 10^12` and `10^50` were checked
  against the paper's formula but are not measurable.

## Scripts

```
code/wall/audit_lem_mp.py               M1, and the repair
code/wall/lab_field_build.py            builds C, V, W, A, S to 1.6e7
code/wall/audit_propV_and_wall.py       prop:V, sec:margin, items 1 and 3
code/wall/audit_propV_readings.py       pins the field to N > 1e5
code/wall/audit_wall_gaussianity.py     M2, M3, and the rho cross-check
code/wall/audit_cell_definition.py      M2: six cell indices, two operations
code/wall/audit_shift_mass.py           M5, M6, repaired lem:MP at 1.6e7
code/wall/audit_majorarc_scan.py        M4
code/wall/audit_cell_floor.py           lem:cellmom vs Monte Carlo, prop:coh
code/wall/audit_mask_decay.py           M14
code/wall/audit_zero_spectrum.py        M13
code/wall/audit_zero_coin_dist.py       M13, the coin distribution
code/supply/e1_dilate_field.py          E1, M8, the R4 identity
code/supply/c3_hb_weight.py             M9
code/demand/audit_propE_majorarcs.py    prop:E margins, major arcs
code/demand/audit_propDpp_cp2.py        prop:Dpp
code/verify/audit_quoted_arithmetic.py  A1-A4, A8-A10
```

Each states its pre-registered decision rule and the prediction it was
written to test in its docstring. **Six predictions were refuted**, and
they are listed because that is the point of writing them down:

1. the `lem:MP` ratio was predicted near 2 and is 1.566;
2. `prop:V`'s figures were predicted to reproduce on the stated field
   and reproduce only on `N > 10^5`;
3. the `conj:wall` cell question was predicted to be resolved by
   valuation-indexed cells and is not;
4. `sec:coin`'s autocorrelation excess was predicted to be real and is a
   mis-specified null;
5. `max_c|z_c|` was predicted to come out much larger than the quoted
   `8.4` under the depth-cell reading, and came out `9.1`;
6. `sec:R4`'s quoted s.e. `0.0112` was predicted to be a wrong
   count-based bar, and a 400-draw permutation null gives `0.0119` —
   **the paper's bar is right there.**

Two predictions that were confirmed and matter: `conj:wall` item 4 was
predicted to fail the coin control (M13), and the shallow-depth decay
exponents were predicted to be unmeasurable (M14).

---

## Retractions

Recorded in one place, per the repository's convention that an
invalidated claim is rewritten in the body and its history kept once.

**M6, original form — withdrawn.** It reported that `sec:coin`'s five
shift-mass percentages do not reproduce, having computed "gross mass" as
`sum_h |M(h)P(h)|` term by term. The paper's own totals settle the
definition the other way: net `-2.2e13` against gross `4.4e13` is a ratio
of 2, and only the sum of **absolute bucket nets** gives that. Under that
reading all five percentages reproduce exactly at `X = 4e6`, with net
`-2.2186e13` and gross/net `1.99`. M6 has been rewritten as a
reproduction plus the surviving point, which is that the table is
specific to an `X` the paper never states.

The failure mode was this tree's own rule 4 — a statistic with no
definition and two defensible readings — applied to the reader rather
than the author.

**A3, second half — withdrawn.** "The body totals 17, so the abstract's
eighteen is wrong" does not stand: counting C-III open as a row gives 18.
The kill-test count ("ten" against a table of nine headed `(9)`) stands.

Scoring against the first pass is in `RECALL.md`.
