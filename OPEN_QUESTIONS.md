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

Register B is the dangerous one. A withdrawn *positive* claim costs
nothing but a claim; a withdrawn *negative* one silently re-opens a
route the program stopped looking at.

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
| **Proposition E** — circle method has zero margin on `C(N)` | the *size* of `C(N)` against `ψ(N) ∼ N`, and `sup ≥ ‖·‖₂` | `#47`: `C(N) ≍ √N (log N)^{0.29}`, i.e. **no sub-square-root cancellation** — larger than a bare `√N`. `#36`: the variance law was missing `√(log N)` | ⚠️ **Re-derive.** The direction is favourable (a larger `C` widens the gap the circle method must close), so the closure most likely strengthens — but "most likely" is not a derivation |
| **Forge R5** — circle method on `C(N)`, magnitude | same as above; "no power slack on that route" | `#47`, `#36` | ⚠️ **Re-derive**, same reason |
| **C-III #4** (second half) — pointwise budget deficit `x^{1/3}` | a magnitude comparison | `#30` re-opened the first half at increment 199; the second was left standing on structural grounds | ✅ Structural half unaffected by later corrections |
| **RV #3** — SEAM over-normalized by `√P` | a **normalization** | `#36`, `#68`, `#287`: three further normalisation corrections, all in the same error class (scale-normalisation drift) | ⚠️ **Re-derive.** Same species as the fault it reports |
| **RV #6** — Theorem E1 (L²) does not match the pipeline (signed L¹) | norm/currency mismatch | `#30` flagged it as the same error class as itself; no later correction touches the norms | ⚠️ Was already "in question" at 199 and has not been resolved since |
| Adjudication route 4 — Perron costs `N^{1−o(1)}` | a magnitude, margin restated as `N^{1/3}` at 199 | none since | ✅ |
| Forge K1–K4, R1–R4; Construction C1, C2, C4 | **existence of structure** — each measured *no* signal | ⚠️ `#110`/`#111` (hazard 7): a null that destroys everything identifies nothing. These kill-tests used permutation nulls | ⚠️ **Re-derive the nulls.** A kill-test that found "no signal" against a destroy-everything null has not shown the signal is absent from `μ` specifically. `#50` audited their *counting* and found it clean; the *null* was not audited |
| C-III #2, #3; RV #1, #2; Adjudication 1, 2, 3, 5 | structural — violated premises, absent congruences, illegitimate transforms | none | ✅ Structural closures do not move with a measurement |

**The pattern.** Every closure at risk is a **magnitude, normalisation
or null** closure. Every closure that is safe is **structural**. That
is the same split this program found at increment 300 — *what is proved
survives an audit, what is fitted does not* — applied to the negative
results instead of the positive ones.

---

## What is not in this file

Corrections themselves (`CLOSURE_REAUDIT.md`), the current position on
any topic (the topic's own document), and the state of the program
(`STATUS.md`). This file is only the list of things that are **open**,
and it is written so that a hole cannot be forgotten simply because the
row that opened it scrolled out of view.
