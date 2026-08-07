# Adversarial re-verification of `v1/paper/wall_v1.tex`

Five rounds, covering the wall sections, the demand side, the supply
side, the closure tables, and the paper's own bookkeeping. Every check
was written from the **statement** in `v1/paper/wall_v1.tex` or
`v1/paper/theorem_A.tex`, not from `v1/code/`, as `v1_verify/README.md`
requires. Where a number disagreed with `v1`, the disagreement was
chased until it resolved; `v1`'s code was read only after an
independent implementation existed and had produced a different number.

**Fifteen standing findings, five of them serious.** Four of my own
pre-registered predictions were refuted and two of my findings were
withdrawn after reading `v1_log/`; all are recorded as such. The demand
side survived every check intact.

**Round 6 cross-checked every finding against `v1_log/`, the program's
own internal record.** That changed the picture: two findings became
much more serious (the paper asserts two claims its own correction
registry **withdrew**), one was confirmed as introduced by the paper
alone, and two of mine were **withdrawn** because `v1_log` already
handled them correctly. See "Round 6" below.

| # | Statement | Verdict | Severity |
|---|---|---|---|
| 1 | Lemma 13 (`lem:MP`), the second-moment identity | **false as stated**; `v1_log` states it correctly — the paper introduced the truncation | high |
| 2 | "Under that input `rho -> 1`" (§`sec:coin`) | **invalid**, short by a factor `N/log N`; not noticed anywhere in `v1_log` | high |
| 3 | The reconstruction `-0.0976` | **wrong number**; divides by `W`, not `V` | high |
| 4 | "against a measured `-0.18`, a factor 0.54" | **withdrawn at `#106`** and re-asserted; a coin reproduces the target | **high** |
| 5 | "Proposition W is an identity" | **overstated**; `v1_log` itself calls it "the uniform-`u` approximation" | medium |
| 6 | Proposition 20 (`prop:coh`) derivation | **self-cancelling**; conclusion survives | medium |
| 7 | `conj:wall` item 4, the zeta ordinates | **withdrawn at `#110`** and re-asserted verbatim, including the withdrawn `0.39%`; the guard that exists to catch this misses it on a LaTeX `\%` | **high** |
| 10 | "ten kill-tested technique designs" | the paper's own table has **nine** | low |
| 11 | `conj:wall` item 3, the `t=5` tail ratio `0.878` | `v1_log` records the counts `4 vs 4.6`; the bar is ±0.43 | low |
| 12 | statement numbers in `v1/PROVENANCE.md` | **off by one** from `conj:L` onward; this tree and `v2/README.md` inherited it | medium |
| 13 | §`sec:R4`'s block ratios | the statistic is **not defined**, and two readings of it give opposite conclusions; neither carries power against its own null | medium |
| 14 | §`sec:c3`(2)'s Heath–Brown table | `A_j`, `D_j`, `z` **not defined**; the definition is recoverable and the conclusion holds, but the third decimal is a rounding convention | low |
| 15 | K1's "full 63-divisor orbit" | **13 of 63 columns are live**; the orbit is not computable on the object K1 addresses, and on the untruncated field the same statistic reads 0.904 — marginal, not dead | **high** |
| 17 | R1's quoted precision | the null spread is from **six draws** (a six-draw sd runs ~half a forty-draw one), and "excludes a 7.5% enhancement" is a fraction of a quantity that is 91% free parameters — against the genuine component it is 84% | low |
| 16 | R2's two criteria | the regression bar is "2× a control whose mean is zero" and does not decide; the coherent-gain bar sits **+5.38σ** of its own null while the measurement sits **+2.78σ**, so "determinant phases blind" is not what was measured | medium |
| ~~8~~ | ~~§3.1 "invisible at accessible N"~~ | **withdrawn by me** — `v1_log` quantifies it correctly | — |
| ~~9~~ | ~~K1's "half the field's energy"~~ | **withdrawn by me** — it rests on a separate residual-energy figure | — |

## Confirmed, no defect found

* **Theorem 1 (`thm:A`)**, tested directly rather than through the
  proof's decomposition: `sup_t |T_1(t)|` computed exactly for every
  `t`, at six `N` from 2e5 to 6.4e6. It falls like `N^{-0.346}`, and
  its size tracks the `e^{-c sqrt(log N)}` that the proof's non-BV
  terms predict (0.0242 at `N = 6.4e6` against `e^{-sqrt(log N)}
  = 0.019`). `v1`'s own numerics never compute `T_1`; they check the
  residual against its main term. This closes that gap.
* **`thm:C`, `thm:D`, `thm:Dprime`, `cor:B`** of `theorem_A.tex`
  (Theorems 3, 5, 6 and Corollary 2 in `wall_v1.tex`'s numbering):
  scanned line by
  line — the divisor switch, the degeneracy lemma, the density
  identity of `theorem_A.tex`'s Lemma 9 (local factor exactly `1/p`),
  the Abel
  summation of Corollary 2, the Step-4 truncation tails, the level
  computation `q <= N^{1-theta'}(log N)^{4A+8}`, and the
  `tau_3`-weighted Bombieri–Vinogradov. No defect.
* **Proposition 7 (`prop:E`)**: table reproduced (margins 0.1680,
  0.1749, 0.1580, 0.1521 against the quoted 0.168, 0.175, 0.158,
  0.152) **and shown to be grid-converged** — every entry moves less
  than 1% between a 4N and a 64N grid. My pre-registered prediction
  that `||S_Lambda||_1` was under-sampled on a 4N grid was **wrong**:
  it moves by 0.03%.
* **Proposition 8 (`prop:Dpp`)**: every quoted digit reproduced —
  22.51/17.36/0.771, 25.04/19.79/0.790, 27.46/22.25/0.810, the
  `CP_2/(N log N)` column, the calibration row, and the tuned column.
  The internal relation `CP_tuned = CP_2 - 2 gamma G` holds, so the
  tuning row is arithmetic on the other two rather than an independent
  measurement.
* **Proposition 11 (`prop:V`)**: confirmed as a by-product of finding 3
  — measured `W(X)/V(X) = 1.27080` against the predicted
  `1/A(N) = 1.27021`.
* **Conjecture 14 (`conj:wall`) item 1**: confirmed. With cells at
  depth 8 the excess kurtosis after cell-mean removal is within
  1.7 sigma of zero in every octave. See "Refuted objections" below.
* **Conjecture 10 (`conj:L`)** band ratios, kurtosis and half-normal
  stamps: corroborated at `N = 1e7` and `4e7` within their own error
  bars.

---

## 1. Lemma 13 (`lem:MP`) is false as stated

**Statement.** "For any `X`, `sum_{N<=X} C(N)^2 = sum_{|h|<X} M(h)P(h)`
... the inner sums being over the ranges that keep both arguments below
`X`." The paper calls it "exact and unconditional".

**Method.** `code/wall/audit_lemMP_identity.py`. Both sides by direct
enumeration — no FFT, no Parseval — at sizes where the enumeration is
the definition.

| `X` | LHS (paper) | RHS (paper) | LHS/RHS |
|---|---|---|---|
| 500 | 8.97506e5 | 1.51572e6 | 0.59213 |
| 2000 | 1.63769e7 | 2.61571e7 | 0.62610 |
| 8000 | 2.46004e8 | 3.98711e8 | 0.61700 |
| 16000 | 9.41243e8 | 1.57288e9 | 0.59842 |

The ratio sits near 0.6 and does not tend to 1.

**Why.** Summing over `N <= X` leaves three free variables `(n, n', v)`
subject to `n + v <= X`, a **simplex**. The product `M(h)P(h)` sums them
over a **box**, and so also counts every pair with `X < n+v <= 2X`. At
`X = 4000` that excluded block carries **36.8%** of the right-hand side.

**The repair, verified exact to machine precision.** With
`C_X(N) = sum_{n+v=N, n<=X, v<=X} Lambda(n) mu(v)` the truncated
convolution and `N` running over its full support to `2X`,
`sum_N C_X(N)^2 = sum_{|h|<X} M(h) P(h)`.

**`v1`'s own record.** `attack_wall_identity.py` carries the comment:
*"A first draft summed `C[:X+1]` against a complete h-sum and duly
missed by 46% — the identity was right and the check truncated one
side."* That 46% miss **is** the failure of the lemma as published; the
truncated side is the paper's left-hand side.

## 2. "Chowla's conjecture gives `S(h) = o(1)` ... under that input `rho -> 1`" is invalid

**Method.** `code/wall/audit_propW_chowla_gap.py`. `S(h)` is an average,
so `|S(h)| <= 1`. Feeding `|S(h)| <= eps` into the displayed formula
gives nothing better than `|rho-1| <= eps * Gamma(N)` with
`Gamma = (sum_{h!=0} c(h))/V = (theta(N)^2 - sum_p log^2 p)/V`.

| `N` | `Gamma` | `Gamma * log N / N` |
|---|---|---|
| 10 000 | 1.549e3 | 1.4266 |
| 160 000 | 1.852e4 | 1.3868 |
| 4 000 000 | 3.580e5 | 1.3605 |

`Gamma` grew 231x while `N` grew 400x. The pre-registered prediction
`Gamma ~ N/(A(N) log N)` is confirmed to three digits (predicted 1.359
at `N = 4e6` against 1.3605 measured, `1/A(N) = 1.2702` times the
`log N/(log N -1)` correction).

**Consequence.** The strength required is `S(h) = o(log N / N)`, not
`S(h) = o(1)` — a gap of a factor `N/log N`.

**Worse: no bound on `|S(h)|` suffices.** The absolute budget
`sum_h c(h)|S(h)|/V`, with `S` measured, is 13.5, 17.5, 23.3, 30.7 at
`X = 2e4 .. 1.6e5` — above 1 by an order of magnitude and **growing**.
Since `c(h) >= 0`, `rho -> 1` can only come from **signed cancellation
across `h`**. Chowla supplies smallness; `[MRT15]` bounds `sum_h |S(h)|`
— absolute values. Neither supplies cancellation.

This is the sentence the abstract compresses into "nature over-delivers
by a power of log and no more".

## 3. The quoted reconstruction `-0.0976` divides by the wrong quantity

**Method.** `code/wall/audit_propW_reconstruction.py`, at `X = 4e6`. The
numerator reproduces `v1` to 0.15% and the shift-mass table reproduces
`v1`'s **exactly** (1.1% / 3.0% / 23.1% / 48.9% / 23.8%, net
`-2.2186e13` against `v1`'s `-2.2152e13`, cancellation 2.0x). The
denominator does not. `lab_offdiag_chowla.py`:

```python
V = float((lam[: X + 1] ** 2 * (mu[: X + 1] != 0)).sum())
```

This evaluates `mu^2` **at the prime power itself**, not at the shift.
Since `Lambda(w) != 0` forces `w = p^k` and `mu^2(p^k) != 0` only for
`k = 1`, it is `sum_{p<=X}(log p)^2 = W(X)`, with no squarefree
condition anywhere.

| quantity | value at `X = 4e6` |
|---|---|
| `V(X) = sum_v mu^2(v) Lambda(X-v)^2` | 4.46842e7 |
| `W(X) = sum_w Lambda(w)^2` | 5.67847e7 |
| `v1`'s number, labelled `V(X)` | 5.6771e7 |

The error is exactly the factor `A(N)` — the local factor whose
identification is Proposition 11 of the same paper, whose §3.1
paragraph is headed *"The local factor is `A`, not `S`"*.

**Corrected figure.** `rho - 1 = -0.12413`, not `-0.0976`; the ratio to
the measured `-0.18` is 0.69, not 0.54.

Swept across `v1/code/`: the reproduction stamps `verify_all.py` and
`verify_deep.py` are **clean** (there `muv` already holds `mu(N-p)`), the
wall machinery builds `V` as a proper convolution and is clean, and
`demand/vaughan_pieces.py:95` has the same shape but is also clean. The
defect is local to `lab_offdiag_chowla.py`, the one place where `lam`
and `mu` are multiplied **at the same index**.

## 4. "a factor 0.54" is not stable under the paper's own definition of `S(h)`

The paper defines `S(h) = <mu(u) mu(u-h)>`, an **average**, whose
denominator is the number of terms `X - h`; the code divides by `X`.
Since 23.8% of the gross mass sits at `h > 1e6` with the **opposite
sign** to the bulk, and `X/(X-h)` reaches 4 there, the choice is not
cosmetic.

| | denominator `W` (`v1`) | denominator `V` (correct) |
|---|---|---|
| `S(h) = M(h)/X` (`v1`) | **-0.09768** | -0.12413 |
| `S(h) = M(h)/(X-h)` (a mean) | -0.04080 | **-0.05185** |

Ratio to the measured `-0.18`: 0.54, 0.69, 0.23, 0.29. The ratio
between the two normalisations is itself unstable in `X`: 0.147 at
`X = 1e6`, 0.418 at `X = 4e6`; at `X = 1e6` the fully corrected value is
`-0.0148`, a factor 0.08. The pre-registered "within a factor 3" bar
absorbs all of it, which is why the instability never surfaced.

**What survives.** Under the mean normalisation the shift shares become
0.8% / 2.2% / 17.0% / 38.4% / 41.6% — the top bin nearly doubles and the
cancellation factor goes from 2.0x to 6.6x — but small shifts still
carry almost nothing, so "the wall leans on the range where that
theorem is strongest" stands.

## 5. Proposition 15 (`prop:W`) is not an identity

§`sec:coin` says *"Proposition W is an identity and is untouched"*.
What is an identity is `C(N)^2 = V(N) + OffDiag(N)`. The step to
`sum_{h!=0} c(h)S(h)` is not: it replaces the values
`mu(N-p)mu(N-p-h)`, sampled on the sparse set of `p` with `p` and `p+h`
both prime, by the average of `mu(u)mu(u-h)` over **all** `u`. That is
an unproven decoupling.

`lab_offdiag_chowla.py`'s own pre-registration says so — *"a loose bar,
set loose deliberately because the derivation drops the prime-power
tail and **treats the u-range as uniform**"* — but the paper does not,
and the residual factor-of-two mismatch of finding 4 is what an
unproven decoupling produces. It matters because §4's coin control is
said not to touch Proposition W *because* it is an identity.

## 6. Proposition 20 (`prop:coh`) does not establish its conclusion — which is nonetheless true

Lemma 19 gives the error bar as **three** terms,
`Q_cc/n_c^2 - 2Q_ca/(n_c n) + Q_aa/n^2`. The substitution
`u_c(v) ~ n_c/sqrt(V)` that `prop:coh` makes is not specific to `c`;
applied consistently it gives all three the same value `S`, hence
`S - 2S + S = 0`. The approximation cancels itself exactly at the order
retained.

**Method.** `code/wall/audit_propcoh_cancellation.py`, all three terms
exact, cross-checked against `sum_v mu^2(v)(u_c/n_c - u_a/n)^2`
(agreement 1.5e-13). At `X = 2^22`:

| `d` | `n_c` | `Q_cc/n_c^2` | `-2Q_ca/(n_c n)` | `Q_aa/n^2` | Var | Var/`(Q_cc/n_c^2)` |
|---|---|---|---|---|---|---|
| 0 | 402252 | 1.355e-1 | -2.484e-1 | 1.283e-1 | 1.532e-2 | 0.113 |
| 1 | 442476 | 1.329e-1 | -2.590e-1 | 1.283e-1 | 2.166e-3 | **0.0163** |
| 3 | 29890 | 2.019e-1 | -2.719e-1 | 1.283e-1 | 5.826e-2 | 0.289 |
| 5 | 70 | 3.409e-1 | -2.792e-1 | 1.283e-1 | 1.900e-1 | 0.557 |

`prop:coh` computes a quantity **2x to 61x larger** than the error bar
it is about, and the ratio is a structural constant, not noise: it
reproduces to four significant figures across a factor 4 in `N`.

**The conclusion survives, independently confirmed.** `Var * log N`
across a factor 4 in `N` gives ratios 0.990 at every depth, against
0.250 for an `n_c^{-1/2}` law. Nothing downstream breaks: `lem:cellmom`'s
bars are computed exactly, so the mask significances are unaffected.

## 7. The zeta-ordinate regression: the effect is real, the quoted strength is not

**Statement.** `conj:wall` item 4: `R^2 = 3.90e-3` against a
"200-surrogate maximum of `5.09e-6`", and "every ordinate individually
at `z >= 23`".

**The null.** `lab_wall_spectral_share.py` draws surrogates as
`rng.permutation(n)` — a **value permutation**, which is white. Its own
docstring calls them "phase-randomised surrogates"; those are different
objects. The regressors sit in the lowest few dozen Fourier bins, which
is where a non-white series has its power.

`code/wall/audit_zeta_regression_null.py`, `X = 4e6`:

| null | mean `R^2` | max of 200 | measured / max |
|---|---|---|---|
| **P** value permutation (`v1`'s) | 1.029e-5 | 2.013e-5 | **230.6** |
| **F** phase randomisation | 1.846e-3 | 3.289e-3 | **1.41** |
| **C** circular rotation | 1.821e-3 | 3.320e-3 | **1.40** |

Null P reproduces `v1` (`20/n`); the correct null is **180x higher in
the mean**. And the per-ordinate claim does not survive at all:

| `gamma` | `z` (null P, `v1`) | `z` (null F) | `z` (null C) |
|---|---|---|---|
| 14.1347 | +99.4 | +3.52 | +3.48 |
| 30.4249 | +25.4 | +0.47 | +0.37 |
| 40.9187 | +14.0 | **-0.05** | -0.29 |
| 49.7738 | +10.6 | **-0.48** | -0.47 |

"Every ordinate individually at `z >= 23`" becomes `-0.48` to `+3.52`,
two of ten negative.

A local-background test —
`code/wall/audit_zeta_local_background.py`, using the comparison `v1`
itself adopted for the same question in `lab_E_zeta_spectrum.py`
(controls drawn from `gamma ± 4`, excluding draws within 1.0 of any
ordinate) — finds **6 of 10 ordinates above their own local 99th
percentile**, against a chance expectation of 0.1 of 10, and a joint
`R^2` above the maximum of 200 draws of ten frequencies from the same
neighbourhoods.

**Round 5 concluded from this that "the effect is real". That
conclusion was premature and is corrected in Round 6 below**: the
local test cannot attribute the effect, `v1_log` correction `#110`
withdrew the claim on a coin control, and re-running that control here
confirms `#110`'s aggregate reading while the local reading disagrees.
Whether the lines belong to `mu` or to `Lambda` is **open**. What is
settled either way is that the evidence as quoted — a ratio of 770
against a white null, and `z >= 23` per ordinate — is worthless.

This is the same species of error the paper's §Methodology records
catching once already: *"an i.i.d.-entry Wigner null applied to a Gram
matrix"*.

## 8. WITHDRAWN by me in round 6 — §3.1's "invisible at accessible N"

**This finding is withdrawn.** `v1_log/docs/AMPLITUDE_ADJUDICATION.md`
already states the point, and quantifies it. Kept below for the record.

The margin by which square-root cancellation clears E1 is
`(N/K)(log N)^{-2A-2}`. Computed at the parameters of the paper's own
supporting measurements:

| | `A = 1` | `A = 2` |
|---|---|---|
| `N = 1e8`, `K = 1e3` | 0.87 | 2.6e-3 |
| `N = 1e9`, `K = 1e3` | 5.4 | 1.6e-2 |
| `N = 1e9`, `K = 1e4` | 0.54 | 1.6e-3 |

The margin is not merely invisible; at Table 1's own top row
(`K ~ 1e4` at `1e9`) and at the `1e8` row it is **below 1**, and for
`A >= 2` it is below 1 by three orders of magnitude everywhere
accessible. The asymptotic statement is correct; the word "invisible"
is not the right word for a quantity that is less than one.

## 9. WITHDRAWN by me in round 6 — K1's `R^2 = 0.466`

**This finding is withdrawn.** "Half the field's energy" rests on a
separately measured residual energy of `0.499`, not on `1 - R^2`. Kept
below for the record.

`e1_forge_kt1.py` regresses `D(k)` on 63 orbit columns over **400**
values of `k`. Pure noise gives `R^2 = 63/400 = 0.157`; the adjusted
`R^2` is `1 - (1-0.466)(399/336) = 0.366`. The paper reads the raw
figure as "**Half** the field's energy is invisible to its entire
multiplicative orbit"; the orbit in fact reaches about a third, and
about two thirds is invisible. `v1`'s own `audit_killtest_nulls.py`
records K1 as carrying **no null at all**. The DEAD verdict is
unaffected — the pre-registered threshold was 0.9.

## 10. The kill-test count does not add up

The abstract says *"eighteen pre-registered closures --- five route
adjudications ..., **ten** kill-tested technique designs, and three
representation-class experiments"*, and the Summary repeats "the
eighteen pre-registered closures". §7's own headings and tables give
`(5)`, `(9)` — nine rows, K1–K4 and R1–R5 — and `(3, plus one open)`.
5 + 9 + 3 = **17**, and the kill-tested count is **nine**, not ten.

`v1_log/code/audit_quoted_numbers.py` passes (24 checks, 0 mismatches)
but checks only figures in `v1_log/` documents against `results/`;
`wall_v1.tex`'s own figures are not machine-checked, and its own LIMIT
note states the residual risk exactly: *"it does NOT check that it is
the right figure for the claim it supports"* — which is findings 2, 3,
4 and 7.

## 11. The `t = 5` tail ratio carries an error bar the paper does not quote

`conj:wall` item 3 gives "ratios 0.999 at `t=3`, 0.997 at `t=4`, 0.878
at `t=5`". At `t = 5` the Gaussian expectation over an octave is 1–2
counts, so the Poisson bar is ±0.77 to ±1.5 (`code/wall/
audit_gaussian_errorbars.py`, part 3). Three decimals on a ratio built
from a handful of exceedances is not a measurement. The `t=3` and `t=4`
figures are fine (bars ±0.011 and ±0.07).

---

## Refuted objections, recorded

The pre-registration discipline cuts both ways. Four of my own
predictions were wrong and are recorded here rather than dropped.

* **`prop:E`'s grid.** I predicted `||S_Lambda||_1` on a 4N grid was
  under-sampled by several percent. It moves 0.03% out to a 64N grid.
  Proposition 7's table is converged.
* **`conj:wall` item 1's attribution.** I suspected the paper
  misattributed `+0.1704 at z = 98` to an `S N`-based scale when
  `v1`'s file labels that row "centring only".
  `lab_gaussian_half_audit.py` line 173 reads
  `Z = C[Ns]/np.sqrt(S[Ns]*Ns)` — the whole script **is** on the `S N`
  scale, so the paper is right and I was wrong. Under the `V` scale the
  same protocol gives excess kurtosis within 1.7 sigma of zero in every
  octave (`code/wall/audit_conjwall_scalemask.py`).
* **Count-based error bars on the Gaussianity stamps.** I expected the
  count bar to be too narrow, by the paper's own `prop:coh` rule. A
  moving-block bootstrap inside the top octave gives
  `SE(kurtosis)/count SE = 0.94–1.14` at every block length from 1 to
  1e4 — the count bar is the right width — and 0.31–0.38 for
  `E|G|/sd`, i.e. three times too **wide**. Objection void.
* **Theorem 1's size.** I predicted `sup_t|T_1|/N` of order 1e-3. It is
  2.4e-2 at `N = 6.4e6`. I had misread `theorem_A.tex`'s 1–4% as
  bounding `T_1`; that figure bounds `R - MT`, while
  `T_1 = P - R - CB ~ -R`. The theorem is unaffected.

One further self-correction: the first version of
`audit_gaussian_errorbars.py` pooled every octave from 1e5 to 1.6e7 and
reported an excess kurtosis of +0.86. That is a scale mixture, not a
measurement — `v1`'s `lab_wall_tails.py` warns about exactly this — and
a second bug subtracted the cell mean of `C` instead of of
`Z = C/sqrt(V)`. Both are fixed; the script records the failure mode.

---

## Round 6: every finding against `v1_log/`

`v1_log/` is the program's internal record — `CLOSURE_REAUDIT.md` (183
numbered corrections), `STATUS.md`, `OPEN_QUESTIONS.md`,
`TECHNIQUE_FORGE.md`, `AMPLITUDE_ADJUDICATION.md`. Each finding above
was checked against it, to see whether the program had already
caught it, already answered it, or — the case that matters — already
**withdrawn** the claim the paper makes.

### 7 and 4 become the serious ones: two withdrawn claims, re-asserted

**Correction `#110`** (`CLOSURE_REAUDIT.md`), verbatim:

> "#94 and #96 read the zeta ordinates in `C(N)` as a property of **the
> wall**, against a **permutation** null | **the lines are `Λ`'s.**
> Replacing `μ` by a random `±1` on the same support, through the same
> `Λ` and the identical pipeline: real `R² = 3.896·10⁻³` against a coin
> mean of `2.994·10⁻³`, coin max `5.515·10⁻³`, and **6 of 20 coin draws
> at or above the real value**. The ratio is **1.30×**, not the `1566×`
> reported. … ⚠️ **Withdrawn**: #94's *"the wall's fluctuation is
> Gaussian in distribution but not phase-random in `log N`"* as a
> property of the wall — it remains literally true and is **empty** …;
> and #96's `0.39%` share, which is not attributable to `μ`."

`wall_v1.tex`, `conj:wall` item 4, states the withdrawn sentence and
the withdrawn number verbatim: *"**`G` is Gaussian in distribution but
not phase-random in `log N`**"* … *"The `0.39\%` is a floor, not the
share"*.

**Correction `#106`**, verbatim:

> "Substituting `ε(v) = ±1` at random on `{μ ≠ 0}` leaves `V`
> **identical** … so a coin must give `ρ = 1`. It gives 0.761, 0.780,
> 0.787, 0.792, 0.813, 0.828, 0.851, 0.860 — reproducing the real curve
> with `z` between −0.5 and +0.4 in every band. **The centred estimator
> cannot tell `μ` from a coin.** ⚠️ Withdrawn: … the *quantitative* half
> of **#86**, which compared a reconstruction against *"the measured
> −0.18"*."

`wall_v1.tex` §`sec:coin` states: *"Reconstructing `ρ−1` from
Proposition W gives `−0.0976` against a measured `−0.18`, a factor
`0.54`."* That is the withdrawn comparison, against the withdrawn
target.

**Both controls re-run here from the statement of Lemma 17**
(`code/wall/audit_coin_control_v1claims.py`), without reading `v1`'s
coin-control code:

| | real | coin mean | coin max | `z` of real |
|---|---|---|---|---|
| `ρ`, top band | 0.8323 | 0.8606 | — | **−0.30** |
| `ρ`, max over 6 bands | — | — | — | **&#124;z&#124; ≤ 0.94** |
| zeta `R²` | 4.640e-3 | 3.273e-3 | 6.655e-3 | 4 of 20 coins ≥ real |

`#106` is confirmed exactly: a coin reproduces `ρ` in every band. `#110`'s
aggregate comparison is confirmed exactly (ratio 1.42× against its
1.30×, 4 of 20 coins against its 6 of 20).

**Why the mechanical guard did not catch it.**
`v1_log/code/audit_withdrawn_forms.py` exists for precisely this
failure (correction `#127`: *"`paper/negative_map.tex` carried two of
them as live measured facts … nine increments of withdrawals never
reached it"*). It scans 38 `.md`/`.tex` files, **including
`v1/paper/wall_v1.tex`**, and has the `0.39%` form registered under
`#96/#110`. It reports one hit in the whole corpus — in *my* round-5
notes. It misses the paper because:

* the registered regex is `(?<![\d.])0\.39\s*%`, and the paper writes
  the LaTeX-escaped `$0.39\%$`. The backslash breaks the match. Zero
  hits in `wall_v1.tex`, one literal occurrence: `0.39\%$ is a flo`.
* the `#94/#110` form matches only `1566` and `z = +13.9`; the paper
  quotes the same withdrawn claim with *different* numbers
  (`5.09·10⁻⁶`, `z ≥ 23`), so nothing fires.
* **no form is registered at all** for `#106`'s withdrawal of the
  `−0.18` comparison. It survives in four live documents:
  `STATUS.md:206`, `CONJECTURE_L.md:151`, `negative_map.tex:572`, and
  `wall_v1.tex`.

The auditor's own closing note anticipates this: *"the form list is
hand-made, so a withdrawal nobody adds here stays unchecked."*

### The zeta effect: what I got wrong, and what is actually open

Round 5 concluded "the effect is real" from a local-background test
(6 of 10 ordinates above their own local 99th percentile). That test
cannot **attribute** the effect, and I reported it as settled. It is
not. Re-running the local test on eight independent coins
(`code/wall/audit_zeta_coin_local.py`):

| field | aggregate `R²` | ordinates above local p99 |
|---|---|---|
| real `μ` | 4.640e-3 | **6** of 10 |
| coins (8 draws) | mean 3.14e-3, max 5.32e-3 | 0, 3, 3, 0, 0, 1, 0, 2 — mean **1.12** |

So the two statistics genuinely disagree, and the reason is
identifiable: the aggregate `R²` at these frequencies is dominated by
the broadband low-frequency power that both fields share, which is why
`#110`'s statistic cannot separate them; the local test subtracts that
background. Real 6 is above every one of eight coins, but eight draws
give at best `p ≈ 0.11`.

**Verdict: open.** `#110`'s withdrawal rests on a statistic that cannot
see the lines; my local test suggests but does not establish that they
are `μ`'s. Neither `v1` nor this tree has settled whether the zeta
component belongs to `μ` or to `Λ`. What is settled is that the
paper's stated evidence — a ratio of 770 against a white null — is
worthless, and that the claim was withdrawn and re-asserted.

### 1 is the paper's own, introduced in transcription

`STATUS.md:175` states the identity **with no upper limit**:
`$\sum_N C(N)^2=\sum_h M(h)P(h)$`. Correction `#155` diagnoses the 46%
miss correctly — *"`μ*Λ` has support to `2X` and Parseval makes the
identity hold over the **full** range"*. The paper added *"For any
`X`"* and `\sum_{N\le X}`, which is exactly what makes it false.

### 2, 3, 5, 6 are not in `v1_log` either

* **2.** `#86` asserts *"Chowla-type input gives `ρ → 1`"*. A search of
  `v1_log` for any quantification of the required strength returns
  nothing. The amplification factor `Γ ~ N/(A log N)` is new here.
* **3.** No occurrence of the wrong denominator anywhere in `v1_log`.
* **5.** `v1_log` is itself split: `#86` calls the step *"the uniform-`u`
  approximation"*, while `#106` says *"Proposition W's identity is
  untouched — it is algebra"*. The paper takes the second reading. The
  first is right: the algebra is `C² = V + OffDiag`; the decoupling is
  not algebra.
* **6.** `#182` gives the identical derivation and `LOCATION_MASK.md:173`
  carries the three-term formula on the same page. The cancellation
  between them is not noticed anywhere.

### 11 sharpened

`#87` records the counts the paper drops: *"`t=5` **4 vs 4.6**
(0.878)"*. Four exceedances. The Poisson bar on the ratio is
`√4/4.6 = ±0.43`.

### 10 narrowed

`STATUS.md:4` reads *"18 recorded closures (13 standing)"*, and the
paper's three tables have 5 + 9 + 3 rows plus C-III open = 18 rows, so
"eighteen" is defensible. **"Ten kill-tested technique designs" is
not**: the table has nine, and its own heading says `(9)`.

### 8 and 9 withdrawn

* **8.** `AMPLITUDE_ADJUDICATION.md:33` already states it, and better
  than I did: *"The margin is asymptotic and invisible at accessible N
  — at `N = 10⁸`, `(log N)^10 = 4.5·10¹²` dwarfs `N/K = 10⁵` — so no
  computation can display it, and **none is offered as evidence for
  it**."* "Invisible" is the log's own word and the log quantifies it.
  Not a defect.
* **9.** `TECHNIQUE_FORGE.md` shows "half the field's energy is
  invisible" rests on a **separately measured** residual energy of
  `0.499` of the unit-Gaussian budget, not on `1 - R²`. My
  adjusted-`R²` point is immaterial: the pre-registered threshold was
  `0.9`, which `0.466` clears either way. Not a defect.

---

## 12. The statement numbers in `v1/PROVENANCE.md` are off by one

`wall_v1.tex` declares one shared counter --- every theorem,
conjecture, corollary, proposition, lemma **and remark** advances it ---
and the paper itself cites only by `
ef`, so it is safe. The
supporting documents write the numbers by hand.
`code/verify/lint_numbering.py` resolves the counter from source and
compares:

| label | `v1/PROVENANCE.md` says | it is |
|---|---|---|
| `conj:L` | Conjecture 11 | **Conjecture 10** |
| `prop:V` | Proposition 12 | **Proposition 11** |
| `prop:W` | Proposition 16 | **Proposition 15** |
| `lem:coin` | Lemma 18 | **Lemma 17** |
| `lem:placebo` | Lemma 19 | **Lemma 18** |
| `lem:cellmom` | Lemma 20 | **Lemma 19** |
| `prop:coh` | Proposition 21 | **Proposition 20** |

Seven of the thirty-three resolvable citations in the repository, all
in the one file, all in the same direction, all from `conj:L` onward.
`v2/README.md` inherited three of them and this tree inherited four
before the lint was written; both are now fixed and `v1` is frozen.
The lint is the fix that survives: it exits nonzero on any citation
that disagrees with the source.

## 13. §`sec:R4`'s block ratios are not defined, and both readings are inside the noise

The paper measures "block sums $S_B(j) = \sum_{k\in    ext{block}}D(k)$
against the $B$-independence that Conjecture 10 predicts" and reports
ratios `0.958` and `1.023` at `B=8` against `B=1` baselines `0.980` and
`0.979`. It does not say what the ratio is. Two natural readings both
equal `1` under independence:

```
r_supp(B)^2 = sum_j S_B(j)^2 / sum_j sum_{k in block j} supp(k)
rho(B)^2    = sum_j S_B(j)^2 / sum_k D(k)^2
```

and they disagree. `code/supply/audit_R4_switch.py` at
`N = 1.6e7`, over the full band `k < sqrt N`:

| `B` | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| `r_supp` | 1.75 | 1.18 | 0.91 | 0.78 | 0.68 |
| `rho` | 1.00 | 0.67 | 0.52 | 0.44 | 0.39 |

`r_supp` falls by a factor `2.6`, which is exactly the signature R4 was
built to detect and would read as a **surviving residue of the switch's
cancellation**. It is not one, and the reason is the paper's own rule
about nulls. Since `supp(k) ~ N/k`, the five smallest `k` carry `38%`
of `sum_k supp(k)`, so the effective number of independent blocks is a
handful and the sampling spread is enormous. Against the null that
holds every `|D(k)|` fixed and randomises only the signs across `k` ---
the one thing R4 is about --- `rho(B)` sits at

```
z = -1.03, -1.15, -1.19, -1.27, -1.19, -0.86   (B = 2 .. 64)
```

with null standard deviations of `0.28` to `0.46`. Nothing reaches two
standard errors at any block size, at either `N`.

So R4's **DEAD verdict stands**, but not on the block ratios: it stands
on the lag-1 autocorrelation, whose null spread is `n^{-1/2}` and which
reads `-0.0446` and `+0.0173` against `0.0264` and `0.0186` here, and
`+0.0104`/`+0.0127` against `0.0112` in `v1`. `v1`'s own
`reaudit_r4_errorbars.py` reached the same place from the other side
--- *"R4's DEAD is a 5-sigma statement at B = 8 and nearly nothing at
the block sizes where its own signature would appear"*.

Confirmed with no defect, in the same run: the full-range switch
identity `sum_k sum_m mu(m)mu(N-mk) = mu(N-1)`, exact at six `N`; and
the mirror identity putting `mu` on the long variable, exact at three.

---

## 14. §`sec:c3`(2)'s Heath–Brown weight table is not defined, and is reproducible only to two decimals

The paper measures "the absolute Heath–Brown weight
`W(a) = sum_j binom(J,j) A_j(a) D_j(M/a)`" and reports the fraction
with `a > M^0.05` as `0.939, 0.949, 0.960, 0.947` at `J = 3,4,6,8`.
`A_j`, `D_j` and `z` are never defined.

**The definition is recoverable.** Heath–Brown's identity of level `J`
with cut `z = M^{1/J}` splits each `n <= M`, in its `j`-th term, into a
Möbius side `a = m_1...m_j` with `m_i <= z` and a complementary side
`b = n_1...n_j`; taking absolute values, the only reading under which
`W` is a count of anything is

```
A_j(a) = #{(m_1..m_j) : prod m_i = a, m_i <= z, mu(m_i) != 0}
D_j(y) = sum_{b <= y} tau_j(b)
```

`code/supply/audit_c3_hb_weight.py` builds exactly that:

| `M` | `J=3` | `J=4` | `J=6` | `J=8` |
|---|---|---|---|---|
| `1e5`, here | 0.933 | 0.942 | 0.958 | 0.961 |
| `1e5`, `v1` | 0.939 | 0.949 | 0.960 | 0.947 |
| `1e6`, here | 0.952 | 0.961 | 0.970 | 0.977 |
| `1e6`, `v1` | — | — | 0.961 | 0.969 |

and the per-`j` concentration reproduces too: `0.811` at `j=3` when
`J=3` against the quoted `0.824`, and `0.840` at `j in {6,7,8}` when
`J=8` against `0.847`. So the reconstruction is the right object,
within about one percent.

**But the third decimal is a convention.** The rounding of `z` alone
moves the `J=8`, `M=1e5` entry from `0.9612` (`ceil`) to `0.9440`
(`floor`) — a swing of `0.017`, wider than the gap to `v1`'s `0.947`,
which sits between them. At `J = 3, 4` the convention moves nothing,
because `z` is large enough that one unit is immaterial.

**The closure is unaffected.** "About 95% of the identity's weight lies
outside the region the classification covers, increasing with `M`" is
confirmed at every `(J, M, z)` convention tested. What is not
supportable is the table at three decimals, and what was missing is the
definition that would let a reader tell.

---

## 15. K1's "full 63-divisor orbit" is thirteen divisors, and its verdict does not survive the full one

**Statement.** §7.2: *"K1 | multiplicative Fejér kernel on the exact
dilation ladder | **dead**: full 63-divisor orbit `R^2 = 0.466`. Half
the field's energy is invisible to its entire multiplicative orbit."*
Pre-registered: `R^2 >= 0.95` alive, `R^2 <= 0.9` dead, `(0.9, 0.95)`
marginal — *"repeat at a second `N` before deciding"*.

**The arithmetic.** E1's dilate field carries the type-II cut, so
`D(v) = sum_{sqrt N < m <= N/v} mu(m) mu(N - mv)` is **empty unless
`v < sqrt N`**. The orbit of a base point `k` is `{D(sk) : s | 30030}`,
so the column at `s` is identically zero unless `s k < sqrt N`. The
largest `s` is 30030, so a base point contributes its whole orbit only
when `N > (30030 k)^2` — which is `3.6e9` already at `k = 2`.

**Measured, at `v1`'s own configuration** (`N = 199,999,998`,
`k in [500,900]`, `code/supply/audit_K1_orbit_reach.py`):

| | |
|---|---|
| columns with any non-empty row | **13 of 63** |
| columns non-empty for every row | **10 of 63** |
| live entries in the design matrix | **17.3%** |
| `R^2` on all 63 columns | **0.4664** (`v1` quotes 0.466) |
| `R^2` on the 13 live columns | 0.4664 |
| pure-noise baseline, as quoted `63/400` | 0.1575 |
| pure-noise baseline, actual `13/400` | 0.0325 |

The fifty dead columns are `30, 33, 35, 39, 42, 55, ... 15015, 30030`
— every divisor of 30030 above 28.

**On the full orbit the verdict changes.**
`code/supply/audit_K1_orbit.py` runs the same design on the
untruncated field, where all 63 columns are live (`N = 1e8`, 400 base
points coprime to 30030):

| | `R^2` |
|---|---|
| real `mu` | **0.9040** (adjusted 0.8860) |
| coin null (Lemma 17), 3 draws | 0.44 |
| permutation null, 200 draws | 0.165 ± 0.051 |

`0.904` sits inside `(0.9, 0.95)` — the band the pre-registration
reserved for *"repeat at a second `N` before deciding"*. And the orbit
reaches far past its coin: `0.90` against `0.44`, so the ladder's reach
is a property of `mu`, which is what K1 was asking.

**Verdict: K1 is open, not dead.** The measurement that closed it
evaluated a fifth of the orbit it names, and the truncation runs in
exactly the direction that suppresses `R^2`. The two objects — type-II
field and untruncated field — are different, and the paper does not say
which K1 used; on the one where the design is computable it is
marginal. §`sec:R4` names the mechanism in another context: *"the
type-II cut makes every divisor sum incomplete."*

**Swept for the same failure elsewhere.** K1 is the only design whose
predictors are the field at **multiples** of the base point. K2, K3,
K4, R2 and C-II use prime-indexed pairs at `p ~ N/2K` with no
`sqrt N` cut; R1, C-I and C-IV use the type-II field at a single
fixed `k` and state their term counts. None of them truncates.

---

## 16. R2's thresholds are effect sizes, and one of them decides nothing

**Statement.** §7.2: *"R2 | determinant / Kloosterman phase | **dead**:
regression `R^2 = -0.0001/+0.0004` against controls `±0.0002` over
48,000 coprime pairs. The verdict rests on the measurement
(`-0.38` s.e.), not on the criterion."* And §7's closing summary:
*"determinant phases blind"*.

**Method.** `code/supply/audit_R2_criterion.py`, rebuilt from the
design's own description at `v1`'s parameters (`N = 199,999,998`,
`K in [2000,4000]`, 40 base `k`, 600 coprime partners each), with the
random-phase control run at 50 and 400 draws instead of 8.

**(a) The regression arm's criterion decides nothing.** The
pre-registration reads *"ALIVE iff mean `G >= 2` or regression capture
`>= 2x` random control"*. The implementation reads

```python
alive = (Gs.mean() >= 2) or (Gs_m.mean() >= 2) or         (mr > 0 and R2d >= 2*mr)
```

The guard `mr > 0` is not in the pre-registration. Measured here,
`mr = +0.000022 ± 0.000302` — positive, but **fourteen times smaller
than its own standard deviation**. So `2*mr = +0.000044`, and
"twice the control" is not a bar: the run the paper reports as
`R^2 = +0.0004` clears it. Whether that arm fired depends on the sign
and size of `mr`, which the paper does not report. **The published
numbers do not determine the published verdict.**

They do determine it on the measurement, which is what the paper says
it rests on: `R2_det = -0.000118` against the control, `z = -0.46`
(the paper quotes `-0.38` s.e.). **DEAD stands.**

**(b) The coherent-gain arm is a 5.4-sigma test, and the measurement
is not zero.**

| | value | `z` vs its own null |
|---|---|---|
| random-phase null, 400 draws | 0.994 ± 0.187 | — |
| `G_1`, phase `e(-N kbar'/k)` | **1.513** | **+2.78** |
| `G_2`, mirror phase | 1.184 | +1.02 |
| the pre-registered bar | 2 | **+5.38** |

The bar was set as an effect size against the *theoretical* null value
1, and the null's own spread was never computed. So the arm could only
have fired on an effect five and a half standard errors out, and the
measured `G_1` sits at `+2.78` — under the bar, but not blind. §7's
summary sentence "determinant phases blind" is not what the
measurement says.

**(c) The pair count.** The design is `40 × 600 = 24,000`, of which
`11,572` survive the design's own `nz > 100` filter and enter the
regression. The paper says "over 48,000 coprime pairs" — four times
the effective `n`.

**What stands.** R2's DEAD verdict, on the measurement. What does not
is the description: the criterion is uncalibrated in both arms, the
coherent gain is `+2.8σ` rather than absent, and the pair count is
overstated fourfold. The paper's own Remark on what a null verdict
costs restates exactly this failure mode for R4 and C4 — R2 belongs in
that list and is not in it.

---

## Reproduction

```
python v1_verify/code/wall/audit_lemMP_identity.py
python v1_verify/code/wall/audit_propW_chowla_gap.py
python v1_verify/code/wall/audit_propW_reconstruction.py
python v1_verify/code/wall/audit_propcoh_cancellation.py
python v1_verify/code/wall/audit_zeta_regression_null.py
python v1_verify/code/wall/audit_zeta_local_background.py
python v1_verify/code/wall/audit_conjwall_scalemask.py
python v1_verify/code/wall/audit_gaussian_errorbars.py
python v1_verify/code/wall/audit_coin_control_v1claims.py
python v1_verify/code/wall/audit_zeta_coin_local.py
python v1_verify/code/supply/audit_R4_switch.py
python v1_verify/code/supply/audit_c3_hb_weight.py
python v1_verify/code/supply/audit_K1_orbit.py
python v1_verify/code/supply/audit_K1_orbit_reach.py
python v1_verify/code/supply/audit_K3_wishart.py
python v1_verify/code/supply/audit_R2_criterion.py
python v1_verify/code/verify/verify_all.py
python v1_verify/code/verify/lint_numbering.py
python v1_verify/code/verify/lint_corrected_paper.py
python v1_verify/code/demand/audit_thmA_direct_sup.py
python v1_verify/code/demand/audit_propE_grid.py
python v1_verify/code/demand/audit_propDpp_table.py
python v1_verify/code/supply/audit_E1_norm_and_gauss.py
```

Outputs are in `v1_verify/results/`, one file per script, same
subdivision. Every script carries its pre-registration — including the
prediction it was written to test — in its docstring.
