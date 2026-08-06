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

## Round 2 (increment 200): the four open items, settled

### RV #6 — which norm the chain consumes: **downgraded, not fatal**

RV #6 objected that "the tex's Theorem E1 (L²) does not match the
pipeline (signed L¹)". Settled from scratch, since it is a
norm/currency call and that is the very class of error #30 was.

The chain's actual consumable is the **signed** sum
T_II = Σ_{k∼K} b_k D(k) ≪ N(log N)^{−A}, with b_k the specific
arithmetic weight the Vaughan decomposition produces (1, log k, …),
not an arbitrary bounded sequence. The L² statement enters **only**
through Cauchy–Schwarz, |T_II| ≤ ‖b‖₂(Σ|D|²)^{1/2} — a proof strategy
we chose, not a requirement. So L² is *sufficient* (the implication
runs the safe way) and *strictly stronger than needed*.

Measured price of that step (`code/norm_audit.py`, N = 10⁸):

| K | 10³ | 3·10³ | 6·10³ |
|---|---|---|---|
| (Σ\|D\|²)^{1/2} | 2661 | 1423 | 774 |
| \|Σ_k D(k)\| | 795 | 550 | 14 |
| ‖1‖₂·(Σ\|D\|²)^{1/2} | 46095 | 24654 | 13403 |
| **Cauchy–Schwarz loss** | **58×** | **45×** | 957× |

The loss is ≍ √K, as it must be for a sign-random field (√1000 = 32,
√6000 = 77; the 957 is one band where the signed sum happened to land
near zero). Both norms are cleared by nature with room — the signed
sums sit 5–7 orders below N — so the √K loss is affordable.

**Verdict: the mismatch is real as an observation and not fatal as an
objection.** Downgraded from "refuted (major)". The strategic
corollary matters more than the verdict: **the chain needs only the
signed L¹ bound, and Cauchy–Schwarz discards exactly the sign
structure of b_k** — the structure whose power R4 already exhibited
(for b_k ≡ 1 the *unrestricted* sum is exactly μ(N−1)). The program
has been aiming at something strictly harder than its own target.

### C-III #4 arithmetic half — **void**

The verdict's own words are "**No-margin** arithmetic: pointwise
consumption loses by N_n^{1/2} ≥ x^{1/3−O(1/J)}". The no-margin
premise is the misstated target. Against the corrected slack N/K ≥
N^{2/3}, a loss of x^{1/3} is affordable, with N^{1/3} to spare.

The structural half also softens. "The object is typically exactly
√-sized, so absolute-value summation reaches at best trivial scale" is
no longer a death sentence: **the corrected target is exactly one
log-power below trivial**, so reaching trivial scale leaves a gap in
the currency rather than a gap in the exponent. What remains is the
verdict's own last clause — all saving must come from u′-family sign
cancellation, which expands into **shift-averaged binary μμ
correlations at log-power strength**.

**Verdict: arithmetic half void; structural half downgraded from
"fatal" to "reduces C-III to a named open problem".**

### RV #4 and #5 — **moot**

Their budgets (the "excess K^{2/3}", the conductor inflation) were
computed against the printed target and are probably void by the same
arithmetic — at K = N^{1/3}, K^{2/3} = N^{2/9}, far inside a slack of
N^{2/3}. But re-deriving them would be archaeology without payoff:
they belong to the proof program that RV **#1 and #2 kill
structurally** (a violated Q-slot premise; an absent pair congruence),
and those two are untouched by any target. Re-opening #4 and #5
revives nothing.

**Verdict: budget arithmetic likely void, immaterial.**

## Where this leaves C-III

All four kill coordinates re-audited: **#1 void, #4 void/downgraded,
#2 and #3 fatal to the draft as written but not to the route.** The
headline "the tree dies even with a free Lemma S" is withdrawn in
full.

What C-III now needs, stated exactly:

1. a **legitimate transform** in place of the draft's hybrid object
   (coordinate #3 stands: Kloosterman factors, √-phase, 1/c weights
   and a c-range are not simultaneously the output of truncated
   Voronoi or of the delta method);
2. a **complete classification** covering the type-II region the draft
   's bookkeeping hid (coordinate #2);
3. **quantitative averaged Chowla**: shift-averaged binary μμ
   correlations with a *fixed log-power* saving. Best known is
   (log)^{1−c} for a small fixed c.

Item 3 is worth stating plainly because it changes the character of
the obstruction. The adjudication's central finding was that the
**dilate** average admits no diagonalizing character family, whereas
the **shift** average is precisely MRT's home ground (translation,
Fejér kernel). If a legitimate transform converting the one into the
other exists, the remaining gap is *quantitative strength inside an
active area*, not the absence of any coupling surface. That is a
different kind of open problem from the one this program had recorded.

**None of this is progress toward Goldbach.** It is the withdrawal of
over-claimed refutations. The route is open again and needs exactly
the thing nobody can currently do.

## Round 3 (increment 201): item ① probed directly

C-III's remaining need ① is "a legitimate transform". Stated cleanly,
the geometry is this. In the L² off-diagonal put a = N−mk, b = N−m′k;
then **m′a − mb = (m′−m)N**, so as k runs the pair traverses lattice
points on a fixed line of rational slope m/m′. A shift correlation
Σ_a μ(a)μ(a+h) is the same object on a line of **slope 1**. So

> dilate problem = μμ correlations along lines of arbitrary rational
> slope; shift problem = the slope-1 case.

Item ① asks for a transform carrying the first to the second. Two
things are now measured (`code/c3_transform_probe.py`).

**(i) The obvious shift reading is an exact relabelling, not a
transform.** Re-index by d = m′−m and fix (d,k): the shift h = dk is
then *fixed* and a = N−mk runs over the progression a ≡ N (mod k), so
the inner factor genuinely is a fixed-shift correlation
Σ_{a≡N(k)} μ(a)μ(a−dk). But **the shift dk is a multiple of the
modulus k**, and rescaling the progression sends μ(a)μ(a−dk) straight
back to μ(N−mk)μ(N−(m−d)k). Verified numerically: the off-diagonal
computed directly and via the (d,k) reading agree **exactly** (184.0
vs 184.0, difference 0.000e+00). The shift reading exists and is
circular.

**(ii) The slope family is not tamer than the shift family.**
Pre-registered: LEAD iff mean |T(m,m′)|/√K ≤ 0.6, CLOSED if the two
families agree within 10%.

| K | slope mean | shift mean | ratio |
|---|---|---|---|
| 10³ | 0.4561 | 0.4826 | 0.945 |
| 3·10³ | 0.4714 | 0.4591 | 1.027 |

**CLOSED.** Both sit at half-normal on their own support (0.46/√0.37 ≈
0.76 against 0.798, the mask accounting of Conjecture L), and they
agree with each other. The two families are the same difficulty; the
difference between them is provability, not size.

**What ① therefore requires.** Not a re-indexing (circular) and not a
motivation from the slope family being tamer (it is not). It must be a
genuine *analytic* transform — a summation formula whose dual variable
supplies the shift average. The named candidates for that are exactly
the ones already measured shut: the inverse-domain spectrum (C2,
mask-exact under accurate nulls), manufactured modularity (C4, 0/6
levels), and the dispersion/delta route (REVIEW_VERDICT #1–#2,
structurally refuted). **So ① is open in principle and empty among the
transform classes this program can name** — which is a weaker claim
than "C-III is dead", and a stronger one than "C-III is open".

## Tally after both rounds

| | count |
|---|---|
| void | 2 (C-III #1, C-III #4 arithmetic) |
| downgraded, not fatal | 2 (RV #6, C-III #4 structural) |
| moot (immaterial even if void) | 2 (RV #4, RV #5) |
| standing | the rest |

The negative map is smaller than it was, and everything that shrank
had been measured with the wrong ruler. What survives is structural:
violated premises, absent congruences, illegitimate transforms, and
kill-tests that found *no* signal rather than *insufficient* signal.
