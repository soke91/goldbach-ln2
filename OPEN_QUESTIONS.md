# Open questions — what the corrections left behind

*Increment 308.* `CLOSURE_REAUDIT.md` records what this program
withdrew and what replaced it. It does not record what each withdrawal
**left open**, and those holes were scattered across 126 rows as ⚠️
marks. This file collects them, in two registers that point in
opposite directions:

- **Register A** — a claim was withdrawn and *nothing replaced it*. The
  question it answered is now unanswered.
- **Register B** — a **closure** (a route declared dead) rested on a
  quantity that a later correction changed. A closure whose premise
  moved is not a closure until it is re-derived.
- **Register C** — a claim that still **stands** was *derived from* one
  that was withdrawn. It does not get withdrawn automatically; it just
  stops having support.

An earlier version of this file said a withdrawn *positive* claim
"costs nothing but a claim". That is wrong, and Register C is what it
missed. A positive claim used to derive something else takes the
derivation with it, silently, and nothing in this repo tracked that
(#135). `code/audit_withdrawn_forms.py` catches the *textual* case — a
document still printing a withdrawn number; `code/audit_withdrawn_dependents.py`
narrows the *logical* case to a reading list, and the reading is done
below rather than automated, because automating exactly this judgement
is what produced #134.

`code/audit_withdrawn_forms.py` enforces the other half of the rule —
that no live document still **asserts** a withdrawn form.

---

## Register A — holes with nothing in them

<!-- rho-ok: every rho in this table is the cancellation ratio
     Var C / V, not Huang-Li's Lemma 1 function; A6 mentions
     Theorem D nearby but does not use its rho. See NOTATION.md -->

| # | The question | What is known | What is not | What would close it |
|---|---|---|---|---|
| A1 | Does `ρ → 1` — is the wall **exactly** square-root? | `ρ = Var C / V` measured `0.857–0.878` at `N ≤ 1.6·10⁷` and **`0.810 ± 0.018` at `N ≈ 10⁸`** (#121). Proposition W makes `ρ − 1` a prime-pair-weighted Chowla correlation, exactly | Whether it tends to 1. Every measurement that claimed to confront it was withdrawn at #300: the de-trended estimator cannot tell `μ` from a coin | ⚠️ **Possibly nothing.** With one realisation of `μ` the quantity may not be estimable at all (#301): the unbiased estimator's spread `0.125` exceeds the deficit `0.094` |
| ~~A2~~ | The **size** of the cell floor | ✅ **CLOSED (inc. 333) to within ~6%.** It was a **weighting fault in the prediction**, found by hazard 9. The closed form averages same-cell pairs with weight `n_c/n`; increment 312 pooled `min(4n_c, 40000)` per cell, under-weighting the large shallow cells whose excess is smallest. The weighting alone inflates by `1.19 → 1.69` against an unexplained overshoot of `1.13 → 1.82`, the first five bands agreeing to `3%`, and the cap reaching 100% of the mass from band 6 — exactly where the curves part (#178) | — | A residual near `6%` at the top three bands is unaccounted for || A3 | The **law** of the mask's decay | ✅ **The exponent's VARIATION is now measured, not asserted (inc. 334).** Supplying the standard errors A3 never had: `a = 0.1434 ± 0.0155` (depth 5), `0.2152 ± 0.0065`, `0.2713 ± 0.0040`, `0.3686 ± 0.0052`, `0.0437 ± 0.0556` (depth 1, **not measurable**), `0.6289 ± 0.0121` (depth 0). A common exponent is rejected at `χ²/dof = 251`. ⭐ And it rises **monotonically as the cell gets shallower**, steps of 5–30 s.e. — the mask decays faster where fewer small primes divide `N` (#180) | ⚠️ **Still open: which law.** `N^{−a}` against `(log N)^{−b}` remains indistinguishable over a factor 160 in `N`, so the numbers above are exponents *of an assumed form* | A wider range in `N`, or a derivation || ~~A4~~ | The mask's **share of the variance** at large `N` | ✅ **ANSWERED (inc. 318).** There are two questions. Pooled `B` answers *how much of a **random** `N`'s fluctuation is the mask* — `0.176 → 0.018` across the octaves, a factor 9.6, and correctly small since a random `N` is in no deep cell. Per-cell `share_c = d_c²/(d_c²+ρ)` answers *how much for an `N` **in that cell*** — for `3·5·7·11·13 \| N` it runs `0.981 → 0.937`, a factor **1.05** | — | **The mask is not small, it is rare.** ⚠️ The debiasing carries a `+0.020` residual (#152), which moves the deep-cell shares by <0.2% and the count of cells above 50% by more || A5 | Is any part of the wall's spectrum **`μ`'s**? | `C(N)` carries the zeta ordinates, but so does a coin on the same support through the same `Λ`: ratio `1.30×` with 6 of 20 coin draws at or above the real value (#110). The lines are `Λ`'s | Whether `μ` contributes any spectral component at all. `#94`'s `z = +13.9` and `#96`'s share were measured against a permutation null that destroys everything | A statistic that isolates the `μ`–`Λ` correlation. Proposition W names it; A1 says it may be unmeasurable ❌ **And the obvious mechanism is refuted (inc. 335).** A2's closed form predicts exponents of `0.000` ± `0.0007` against measured `0.14–0.63`, because `D_c` is **scale-invariant** by construction — it depends only on `h`'s residues mod the cell's primes. The decay is **not** arithmetic; it comes from the normalisation, and the depth-dependence from something the closed form omits (#181) | ⭐ **The frame the refutation hands over**: the error bar's own decay is `N^{-1/2}`-ish, and the measured exponents straddle it — `0.143, 0.215, 0.271, 0.369` are all **below** `1/2` and `0.629` is above. So in every cell but the shallowest the mask **grows relative to its own error bar**, which is A4's per-cell share seen in `N`. Measure the exponent of `dm/se`, not of `dm` || A6 | Novelty of **Theorem D / D′** | The proofs stand and were re-derived at every structural step | Whether the no-go is already in the literature. `LITERATURE.md` records the search this program could do, which is not a specialist's | MathSciNet / zbMATH, or one analytic number theorist reading two pages |

---

## Register B — closures whose premise a later correction moved

`CLOSURE_REAUDIT.md` was itself written (increment 199) to re-audit
every closure against **one** correction, `#30`. Ninety-six corrections
have been recorded since, and **no second pass has been made.**
Increment 282 checked exactly one closure and found the conclusion
stood while its evidence was `1.8 s.e.` rather than the `29 s.e.`
recorded. That is one of thirteen standing closures.

The table below began as a **triage, not a re-audit**. Increments
309–317 worked through it, and **every row is now settled**: Proposition
E and Forge R5 re-derived and strengthened, RV #6 found to have been
settled at 199 already, RV #3's corrected form derived and met with a
power margin, and the kill-tests given the floors and error bars none of
them had. ⚠️ What survives is not the triage's verdicts but its
corrections to itself: four of its rows were wrong as written (#130,
#134, #147), which is the reason a register of this kind has to be
checked against sources rather than trusted.

| Closure | What it rests on | Corrections that moved it | Status |
|---|---|---|---|
| ~~**Proposition E**~~ — circle method has zero margin on `C(N)` | the **trivial bound** `ψ(N) ∼ N` (PNT) and `sup ≥ ‖·‖₂` (Parseval). **Not** the size of `C(N)` — the triage below got that wrong | **none.** Every quantity is an exact sum over `n < N` or a consequence of PNT | ✅ **RE-DERIVED (inc. 309): stands, and is stronger.** `W(N)` cancels between `√(WQ)` and `rms C = √(ρ·𝔄·W)`, leaving `deficit = √(Q/(ρ𝔄)) ≍ √N` — a clean power of `N`, not a log power. Measured `0.888 → 0.961` in units of `√N` across a factor 100 in `N`; Parseval floor `1.52 → 1.92`, never below 1 (#130) |
| ~~**Forge R5**~~ — circle method on `C(N)`, magnitude | the same two lines | **none** | ✅ **Settled with Proposition E** (inc. 309) |
| **C-III #4** (second half) — pointwise budget deficit `x^{1/3}` | a magnitude comparison | `#30` re-opened the first half at increment 199; the second was left standing on structural grounds | ✅ Structural half unaffected by later corrections |
| ~~**RV #3**~~ — SEAM over-normalized by `√P` | a **normalization** | `#36`, `#68`, `#287` | ✅ **RE-DERIVED (inc. 317). The corrected form is met with a power margin.** As written SEAM demanded `|C| ≪ √M(log N)^{−A}`, a log saving over the **square-root** scale — false. The chain consumes at the **trivial** scale (#30), so the demand is `|C| ≪ M(log N)^{−A}`, and `|C|/M ≈ 0.797/√M` beats every fixed `A` by a **power of `M`**; no drift against `log M` (`z = +1.62`). ⚠️ RV #3's own figure was also wrong — `0.717` is the unconditioned estimator; the clean bucket gives `0.7979` within 3 SE (#149) |
| ~~**RV #6**~~ — Theorem E1 (L²) does not match the pipeline (signed L¹) | norm/currency mismatch | none | ✅ **It was settled at increment 199 and this register said otherwise** (#147): the verdict *table* says "in question", the *body* forty lines below says "the mismatch is real as an observation and not fatal as an objection". ⚠️ Its measurement had no spread; supplied at 316 — the loss ratio's **mean does not exist** (`mean/median = 4.7·10¹¹`), its single-band point estimate reaches **38× its own median**, and the `957×` outlier is a 5%-probability tail. ✅ The conclusion survives on the median: `loss ≈ 1.44√K`, constant to 9.3% for `K ≥ 700` (#148) |
| Adjudication route 4 — Perron costs `N^{1−o(1)}` | a magnitude, margin restated as `N^{1/3}` at 199 | none since | ✅ |
| **Construction C1, C2, C2b, C4** | existence of structure — each measured *no* signal | **none.** All four already use a **coin** null: random signs on the real support | ✅ **Hazard 7 is already satisfied** in these four. `e1_constr_c1.py`'s own header reads "Null: 8 draws of random signs on the real support" (#134) |
| **Forge K2, K3, K4** | existence of structure | an **analytic** `z`, no resampling | ⚠️ **Read the distribution, not the null.** An analytic `z` is fine if the reference distribution is right; hazard 7 does not apply, hazard 8 does — the spread was never checked against a resampled one |
| ~~**Forge R1**~~ | existence of structure | a **surrogate** null: random-frequency templates | ✅ **RE-READ (inc. 313): the verdict is much stronger than it claimed.** Null `0.2196 ± 0.0055`, so the pre-registered `2×` threshold sits `+39.9 sd` above it while the measurement sits `−0.80` below. The data exclude a **7.5% enhancement at 3 sd**, thirteen times sharper than the recorded DEAD (#140) |
| ~~**Forge R2**~~ | existence of structure | ratio against a random control | ⚠️ **The criterion could not have returned ALIVE** — it asks for `≥ 2×` a control that came out **negative** (#141). ✅ The verdict stands on the **measurement** instead: `−0.38 sd`. Which of the two a closure rests on had never been distinguished |
| ~~**Forge R4**~~ | existence of structure | `ratio(B) ≤ 0.5 × ratio(1)` | ⚠️ **RE-MEASURED (inc. 315): a 5-sigma test at `B = 8`, no information at `B = 512`.** `ratio(B)` sums `nb = 2048/B` squares, so `nb = 256, 32, 4` and the relative SE runs `9.8%, 35.4%, ~71%` — **degrading fastest exactly where R4's own signature (a fall at large `B`) would appear**. The recorded fall is under 2 sd at every `B` and the replication reverses its sign at `B = 512` (#145). R4b's quoted SEs are underestimates (#146) || ~~**Construction C4**~~ | existence of structure | `defect_real ≤ 0.5 × defect_null` | ⚠️ **RE-MEASURED (inc. 314): the criterion is a `1.3`-sigma test that pure noise satisfies `8.8%` of the time.** Over six levels that is a **42%** chance of a spurious ALIVE, so "0 alive of 6" had probability 58% under noise. ✅ The DEAD direction stands, but on the **measurement** (`+0.10 sd`, never systematically below the null) and not on the count (#143) || **Forge K1, R2, R2b, R4, R4b** | existence of structure | **no null the classifier can see.** R4b's own header shows why: its reference is a `B = 1` baseline ratio, not a resampling — and it re-ran itself for power after a `+3.7σ` band in this program's history regressed to noise | ⚠️ **Read individually.** "No null detected" is a limit of the classifier, not a verdict |

> ⚠️ **This row previously read**: "Forge K1–K4, R1–R4; Construction C1, C2, C4 — all measured no signal against **permutation nulls**, which hazard 7 invalidates." That was written without opening the files and is **false for twelve of thirteen** (#134): exactly one contains a permutation null, four already contain the coin control hazard 7 asks for, and five contain no resampling at all. `code/audit_killtest_nulls.py` classifies them mechanically.

### ⚠️ Every DEAD verdict has a floor, and none of them stated it

Increment 311 measured one, for K2 (`code/reaudit_killtest_power.py`).
The test flags at `|mean|/(sd/√n) ≥ 4` over `n ≤ 600` consecutive `k`,
so its **detection floor is `4/√n ≈ 0.163` standard deviations** of the
field. Shown, not assumed: at twice the floor every one of ten `h`
values flags (`z` from `+7.2` to `+9.8`), so the test **can** return
ALIVE; unmodified it returns DEAD at reduced `N` exactly as the
original did.

**That changes what every kill-test-based closure says.** "DEAD"
was read throughout this program as *no structure*. It means *no
coherent structure above the floor*, and the floor was never quoted.
The two are different statements, and the difference matters exactly
when a design needs a **small** gain — which is the usual case, since
these designs are trying to beat a Cauchy–Schwarz cost of a factor `d`.

| What is now required of every kill-test-based closure | Status |
|---|---|
| its floor measured, in units of the field's own sd | K2 done (`0.163`); twelve others not |
| its ALIVE branch demonstrated, not asserted | K2 done; twelve others not |
| **the design's own requirement stated in the same units** — how large a gain would make it ALIVE? | **K2 done (inc. 312); twelve others not** |

The third row is the one that decides whether these closures stand,
and for K2 it now reads. Design K2 splits the `k`-average mod `d` at a
Cauchy–Schwarz cost of a factor `d`, so a coherent shift must satisfy
`1 + n_r δ²/c > d`, i.e. `δ_design ≈ √(c·d(d−1)/n)` with
`c = Var₀/n_r` **measured** at `0.516`. Against the floor `4/√n`
the two cross near `d = 7`.

⚠️ **But beating the factor `d` is not the route's bar.** The corrected
E1 target (#30) needs a saving of `(log N)^{2A+2}` over trivial, so the
split must reach `d ∼ (log N)^{2A+2}`. A gain of a factor 6 is not a
log power and cannot help the route whatever the design does with it.
At the `d` the target does need:

| `A` | `d` needed | `δ_design` | against the floor |
|---|---|---|---|
| 1 | `8.0·10⁴` | `3.4·10³` | **14,000×** |
| 2 | `2.3·10⁷` | `9.5·10⁵` | **4.1·10⁶×** |

✅ **K2 is the first kill-test closure in this program to stand
quantitatively**, with its one blind region named: `d ≤ 6`, worth at
most a factor 6, which the route cannot spend (#138).
| C-III #2, #3; RV #1, #2; Adjudication 1, 2, 3, 5 | structural — violated premises, absent congruences, illegitimate transforms | none | ✅ Structural closures do not move with a measurement |

**The pattern.** Every closure at risk is a **magnitude, normalisation
or null** closure. Every closure that is safe is **structural**. That
is the same split this program found at increment 300 — *what is proved
survives an audit, what is fitted does not* — applied to the negative
results instead of the positive ones.

---

## Register C — what was derived from a claim that fell

Fourteen corrections withdrew a *positive* claim. The mechanical pass
finds **35 citations of them across the live documents, of which 22 are
history pointers and 13 are withdrawal sentences** — no document states
a withdrawn number as its own. That is the textual half. The logical
half is below: for each withdrawn claim that anything was built on,
what was built, and whether it survives without it.

| Withdrawn | What was built on it | Does it survive? |
|---|---|---|
| **#36** the fitted variance law `Var C ≈ 0.465·𝔖·N·log N` | every statement of the wall's *scale* — that `C` is square-root sized, that `G` has unit variance, the normalisation in Conjecture L | ✅ **Yes, and by a better route.** Proposition V replaced the fit with the **exact** second moment `V(N) = Σ_v μ²(v)Λ(N−v)²`. Everything downstream now goes through `V`, which needs no fit. The withdrawal removed a crutch, not a support |
| **#47** the decay exponent of `Q = ΣC²/V` | Proposition E's *strengthening* at increment 309 — "the method is short by a power of `N`" | ✅ **Yes, and the derivation was rebuilt to avoid it.** The docstring motivating 309 quoted `C ≍ √N(log N)^{0.29}`, but the derivation that ran does not use it: `W(N)` cancels between `√(WQ)` and `√(ρ𝔄W)`, so the result needs only that `𝔄` and `ρ` are bounded, both measured directly |
| **#67, #69, #112, #119** the mask's share and its scaling `N^{1/4}` | the load-bearing claim **"the mask is lower order and does not threaten `C(N) = o(N)`"** | ✅ **Yes, without any exponent.** That claim needs only `\|m(N)\| = O(√V)`, and the mask is measured *in units of* `√V`: the deepest cell reads `−7.09` at `N ≈ 10⁵` falling to `−3.58` at `1.4·10⁷`. Since `√V ≍ √(N log N) = o(N)`, the conclusion follows from a bounded number, not from a fitted exponent. ⚠️ The paper stated it the other way until #127 |
| **#84, #99** the direction of `ρ`'s trend and the rate `b = 2.68` | would be **"the wall is asymptotically exactly square-root"** | ✅ **Nothing to save — it was never claimed.** That statement is Register A1 and is open. Proposition W, which is algebra, is untouched; what fell is every *measurement* that claimed to confront it |
| **#94, #96, #110** the wall carries the zeta zeros, share `0.39%` | would be **"the wall's fluctuation is not phase-random in `log N`"** | ✅ **Survives as a statement about `Λ`, not the wall.** `C = μ*Λ` and `Λ` carries the zeros by the explicit formula. Conjecture L's bulk and tail results do not use the spectrum at all, so nothing else moves |
| **#118** the mask's amplitude is unresolved at large `N` | increment 303's diagnosis that the coin floor was what hid it | ⚠️ **Partly.** The floor is real, but the *reason* the aggregate missed the mask is that `B` weights cells by size while the mask lives in rare cells. Corrected in `LOCATION_MASK.md` at 305 |

**The result is that nothing load-bearing fell with them** — but that
sentence is worth exactly as much as the reading behind it, which is
one pass by the same author who wrote the claims. It is recorded here
so it can be checked rather than assumed.

## What is not in this file

Corrections themselves (`CLOSURE_REAUDIT.md`), the current position on
any topic (the topic's own document), and the state of the program
(`STATUS.md`). This file is only the list of things that are **open**,
and it is written so that a hole cannot be forgotten simply because the
row that opened it scrolled out of view.
