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

| # | The question | What is known | What is not | What would close it |
|---|---|---|---|---|
| A1 | Does `ρ → 1` — is the wall **exactly** square-root? | `ρ = Var C / V` measured `0.857–0.878` at `N ≤ 1.6·10⁷` and **`0.810 ± 0.018` at `N ≈ 10⁸`** (#121). Proposition W makes `ρ − 1` a prime-pair-weighted Chowla correlation, exactly | Whether it tends to 1. Every measurement that claimed to confront it was withdrawn at #300: the de-trended estimator cannot tell `μ` from a coin | ⚠️ **Possibly nothing.** With one realisation of `μ` the quantity may not be estimable at all (#301): the unbiased estimator's spread `0.125` exceeds the deficit `0.094` |
| A2 | The **size** of the cell floor | Its *shape* is the singular series of the shift, `corr(ρ(h), 𝔖₂(h)) = 0.9997–1.0000`; its *cause* is the cell↔divisibility correspondence, confirmed by a placebo key that collapses `B/T` to the iid `(k−1)/n` (#116) | 𝔖₂ overpredicts the closed-form floor by `1.13×` to `1.82×`, **growing with `N`**. The one suspect I could name — a single `a` per band where `a` decays in `h` — was tested and cleared (#119 run) | A second-order term in the pair correlation, or a systematic effect of the finite window. Not identified |
| A3 | The **law** of the mask's decay | The decay is real: it beats a constant by `>2×` in weighted RSS at 5 of 6 depths, with exact per-point errors | Over a factor `160` in `N` the data do **not** separate `A·N^{−a}` from `A·(log N)^{−b}`, and the fitted exponent varies about `4×` **with the cell** under either parameterisation. **There is no single mask exponent** (#119) | A wider range in `N`, or a derivation. `#69`'s `g = −0.489` is withdrawn and not replaced |
| A4 | The mask's **share of the variance** at large `N` | The mask is resolved at every octave to `1.6·10⁷` (`max_c\|z_c\| = 8.4` at the top) and at `N ≈ 10⁸` the two arms differ by `+0.2238 ± 0.0056` (#118, #126) | A *share* — `#67`'s "14.1% at `10⁵`, 1.15% at `1.6·10⁷`" is withdrawn. The aggregate that produced it weights cells by size, and the mask lives in cells of relative size `6.6·10⁻⁵` | A share statistic that is not size-weighted. None is defined yet |
| A5 | Is any part of the wall's spectrum **`μ`'s**? | `C(N)` carries the zeta ordinates, but so does a coin on the same support through the same `Λ`: ratio `1.30×` with 6 of 20 coin draws at or above the real value (#110). The lines are `Λ`'s | Whether `μ` contributes any spectral component at all. `#94`'s `z = +13.9` and `#96`'s share were measured against a permutation null that destroys everything | A statistic that isolates the `μ`–`Λ` correlation. Proposition W names it; A1 says it may be unmeasurable |
| A6 | Novelty of **Theorem D / D′** | The proofs stand and were re-derived at every structural step | Whether the no-go is already in the literature. `LITERATURE.md` records the search this program could do, which is not a specialist's | MathSciNet / zbMATH, or one analytic number theorist reading two pages |

---

## Register B — closures whose premise a later correction moved

`CLOSURE_REAUDIT.md` was itself written (increment 199) to re-audit
every closure against **one** correction, `#30`. Ninety-six corrections
have been recorded since, and **no second pass has been made.**
Increment 282 checked exactly one closure and found the conclusion
stood while its evidence was `1.8 s.e.` rather than the `29 s.e.`
recorded. That is one of thirteen standing closures.

The table below is a **triage, not a re-audit**. It says which
closures rest on a quantity a later correction touched, and therefore
which must be re-derived before they may be quoted. It does not claim
any of them is wrong.

| Closure | What it rests on | Corrections that moved it | Status |
|---|---|---|---|
| ~~**Proposition E**~~ — circle method has zero margin on `C(N)` | the **trivial bound** `ψ(N) ∼ N` (PNT) and `sup ≥ ‖·‖₂` (Parseval). **Not** the size of `C(N)` — the triage below got that wrong | **none.** Every quantity is an exact sum over `n < N` or a consequence of PNT | ✅ **RE-DERIVED (inc. 309): stands, and is stronger.** `W(N)` cancels between `√(WQ)` and `rms C = √(ρ·𝔄·W)`, leaving `deficit = √(Q/(ρ𝔄)) ≍ √N` — a clean power of `N`, not a log power. Measured `0.888 → 0.961` in units of `√N` across a factor 100 in `N`; Parseval floor `1.52 → 1.92`, never below 1 (#130) |
| ~~**Forge R5**~~ — circle method on `C(N)`, magnitude | the same two lines | **none** | ✅ **Settled with Proposition E** (inc. 309) |
| **C-III #4** (second half) — pointwise budget deficit `x^{1/3}` | a magnitude comparison | `#30` re-opened the first half at increment 199; the second was left standing on structural grounds | ✅ Structural half unaffected by later corrections |
| **RV #3** — SEAM over-normalized by `√P` | a **normalization** | `#36`, `#68`, `#287`: three further normalisation corrections, all in the same error class (scale-normalisation drift) | ⚠️ **Re-derive.** Same species as the fault it reports |
| **RV #6** — Theorem E1 (L²) does not match the pipeline (signed L¹) | norm/currency mismatch | `#30` flagged it as the same error class as itself; no later correction touches the norms | ⚠️ Was already "in question" at 199 and has not been resolved since |
| Adjudication route 4 — Perron costs `N^{1−o(1)}` | a magnitude, margin restated as `N^{1/3}` at 199 | none since | ✅ |
| **Construction C1, C2, C2b, C4** | existence of structure — each measured *no* signal | **none.** All four already use a **coin** null: random signs on the real support | ✅ **Hazard 7 is already satisfied** in these four. `e1_constr_c1.py`'s own header reads "Null: 8 draws of random signs on the real support" (#134) |
| **Forge K2, K3, K4** | existence of structure | an **analytic** `z`, no resampling | ⚠️ **Read the distribution, not the null.** An analytic `z` is fine if the reference distribution is right; hazard 7 does not apply, hazard 8 does — the spread was never checked against a resampled one |
| ~~**Forge R1**~~ | existence of structure | a **surrogate** null: random-frequency templates | ✅ **RE-READ (inc. 313): the verdict is much stronger than it claimed.** Null `0.2196 ± 0.0055`, so the pre-registered `2×` threshold sits `+39.9 sd` above it while the measurement sits `−0.80` below. The data exclude a **7.5% enhancement at 3 sd**, thirteen times sharper than the recorded DEAD (#140) |
| ~~**Forge R2**~~ | existence of structure | ratio against a random control | ⚠️ **The criterion could not have returned ALIVE** — it asks for `≥ 2×` a control that came out **negative** (#141). ✅ The verdict stands on the **measurement** instead: `−0.38 sd`. Which of the two a closure rests on had never been distinguished |
| **Forge R4, Construction C4** | existence of structure | ratio against a baseline | ⚠️ **No error bar is printed beside the quantity the threshold judges.** R4b's own header puts the SE at ~25% at `B=64`; at `B=512` the measured `0.601` sits near the `0.5` threshold with an SE never quoted. C4's only clue to its sd is the across-`Q` spread ≈ 0.25, putting the threshold ~2 sd from the null (#142) |
| **Forge K1, R2, R2b, R4, R4b** | existence of structure | **no null the classifier can see.** R4b's own header shows why: its reference is a `B = 1` baseline ratio, not a resampling — and it re-ran itself for power after a `+3.7σ` band in this program's history regressed to noise | ⚠️ **Read individually.** "No null detected" is a limit of the classifier, not a verdict |

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
