# Re-audit of every closure against the corrected target (increment 199)

*Correction #30 (increment 198) found that this program had been
carrying a misstated target: the consumable was printed as
Σ_{k∼K}|D(k)|² ≪ (log N)^{−A} Σ M_k — the **square-root** scale —
where the chain actually requires*

> **Σ_{k∼K}|D(k)|² ≪ (log N)^{−2A−2} Σ_{k∼K} M_k²**,

*i.e. a fixed log-power saving over the **trivial** bound. The two
differ by a factor M ≍ N/K, which for K ≤ N^{1/3} is ≥ N^{2/3}.*

**Therefore every death sentence that was passed on a magnitude,
budget, or scale-mismatch argument was measured against the wrong
yardstick and must be re-opened.** This document does that, closure by
closure. It is written to find our own errors, not to defend the
record.

## The decisive quantity

Under the corrected target the proof has a **slack of N/K ≥ N^{2/3}**:
a technique may be lossy by up to that factor and still deliver. Under
the printed target the required bound sat *below* what nature supplies,
so the slack was zero and **any** loss was fatal. That single change
is what re-opens things.

Note the slack is specific to the E1 route. For the scalar C(N)
itself, trivial is ψ(N) ∼ N and the target is N(log N)^{−A}, so the
slack there is only a log power — Proposition E (circle-method zero
margin) is unaffected.

## Verdict table

| Closure | Basis of death | Re-audit |
|---|---|---|
| **C-III #1** — §5 assembly, "scale mismatch x^{2/3}" | **magnitude** | **VOID** — see below |
| C-III #2 — classification hides a type-II region | structural | stands |
| C-III #3 — the dual object is the output of no legitimate transform | structural | stands |
| **C-III #4** — pointwise budget deficit x^{1/3} | **magnitude + structural** | **first half in question**; second half stands |
| RV #1 — gate arithmetic, Q-slot premise violated | structural (premise) | stands |
| RV #2 — role exchange, no pair congruence | structural (hypothesis) | stands |
| RV #3 — SEAM over-normalized by √P | **normalization** | stands, but see "same error class" |
| **RV #4** — h-range collapse, "excess K^{2/3}" | **magnitude** | **in question** |
| **RV #5** — second Vaughan layer multiplies the conductor | **magnitude + structural** | **in question** (first half) |
| **RV #6** — Theorem E1 (L²) does not match the pipeline (signed L¹) | **norm/currency** | **in question** — this is the same error class as #30 itself |
| Adjudication routes 1, 3 | structural (lemma hypotheses) | stand |
| Adjudication route 2 | structural ×2 + "triple-log saving below spec" | stands (triple-log is still not a fixed log power) |
| Adjudication route 4 | structural + "Perron costs N^{1−o(1)}" | stands, but the margin is now N^{1/3}, not "hopeless" |
| Adjudication route 5 | structural (exceptional set) | stands |
| Forge K1–K4, R1–R4 | existence-of-structure kill-tests | stand — each measured *no* signal, not *insufficient* signal |
| Forge R5 (circle method) | magnitude, but on C(N) | stands — no power slack on that route |
| Construction C1, C2, C4 | existence-of-structure | stand |

## C-III #1 is void

The verdict reads, verbatim:

> "E1 is an L²(k) statement at GAUSSIAN (square-root) scale: it demands
> |D(k)|² ≈ M_k. Every leaf of the tree outputs LINEAR-scale pointwise
> log-savings (|D(k)| ≪ (x/k)(log)^{−C}). Squaring gives
> Σ|D|² ≪ (N²/K)(log)^{−2C} — off target by N/K ≥ x^{2/3} at every
> J, K, A."

The premise in the first sentence is exactly the misstatement of
correction #30. The corrected target is

> Σ|D|² ≪ (log N)^{−2A−2}·N²/K,

and the tree's own output, as the verdict itself computes it, is

> Σ|D|² ≪ (N²/K)(log)^{−2C}.

**These are the same statement**, and the leaves meet the target
whenever C ≥ A+1. The "off target by N/K ≥ x^{2/3}" was the distance
to the square-root scale, which is not the target and which nothing in
the chain requires. **Refutation #1 does not hold.**

Consequence for the headline: C3_REVIEW claims "four independent kill
coordinates — two of which stand even if Lemma S is granted for free",
and those two were #1 and #4. #1 is void and #4's arithmetic half is in
question. **The claim that the tree dies even with a free Lemma S is
withdrawn.**

## What this does and does not re-open

**Does not**: C-III is still refuted *as drafted*. Coordinate #3 is
fatal to the specific construction — the draft's hybrid object
(Kloosterman factors, √-phase, 1/c weights and a c-range at once) is
the output of neither legitimate transform — and coordinate #2 shows
the classification leaves an unhandled type-II region. Those are
independent of any target.

**Does**: C-III is no longer *impossible on budget grounds*. It moves
from "dead even with a free Lemma S" back to "an open route requiring
a correct construction", and its remaining obstruction (C-III #4,
second half) is now a *named* target rather than a scale gap: all
saving must come from u′-family sign cancellation, which expands into
shift-averaged binary μμ correlations needed **at log-power strength**
— MRT/Pilatte territory, where the best known savings are
(log)^{1−c} for a small fixed c, i.e. short of a fixed power but in an
actively moving area.

## The same error class, twice

RV #3 killed the SEAM formalization for being "over-normalized by √P
— it demands more than the chain needs". That is the identical species
of mistake as #30: **this program has twice written its own target at
the wrong power of the natural scale.** RV #3's verdict stands (the
SEAM statement as written was indeed too strong), but the corollary is
that **SEAM should be re-derived against the corrected target before
it is called false** — an over-strong conjecture that is false is not
evidence against the correctly-normalized one.

Named as a standing hazard: *scale-normalization drift* — writing a
target at the square-root scale when the chain consumes it at the
trivial scale, or vice versa. Both known instances were caught only by
comparing a stated target against our own measurements of the same
quantity. That comparison is now mandatory before any target is used
to adjudicate anything.

## What must be re-derived (not adjudicated here)

RV #4, #5, #6 and C-III #4's arithmetic half cannot be settled from
the verdict texts alone: they quote powers (K^{2/3}, conductor
inflation, an L²-vs-L¹ mismatch) whose derivations live in documents
that themselves used the printed target. Each needs its budget redone
against Σ M_k². They are recorded here as **open re-audit items**, not
as reinstated or as void.

Honest summary: **one closure void, four in question, the rest stand.**
The negative map is smaller than it was this morning, and the part of
it that shrank is precisely the part that had been measured with the
wrong ruler.
