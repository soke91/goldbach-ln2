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
structurally refuted).

## Round 4 (increment 204): ① closed for changes of variable

The centered coordinates settle half of ① by proof rather than by
measurement. Put A = N−a, B = N−b, so A = mk and B = m′k; then

> **m′A − mB = 0 exactly**

(verified over 200,000 triples, zero mismatches,
`code/c3_pencil_check.py`). So the dilate family is the **pencil of
lines through the origin** and the shift family B = A − h is a family
of **parallel lines of slope 1**. A pencil is characterised by its
vertex — here a finite point — and a parallel family is a pencil whose
vertex is at infinity. **Affine maps send finite points to finite
points**, so no affine map carries one family to the other; and since
μ lives on ℤ, the structure-preserving changes of variable available
are exactly the integral affine ones. No shear helps either: it
permutes slopes and fixes the origin.

**So ① is closed for pointwise changes of variable**, which explains
the measured circularity rather than merely recording it. It remains
open for summation formulas — where the obstruction is coordinate 3's,
in general form: Voronoi-type formulas compress only against a
**smooth** outer weight, and ours is μ.

**Status of ①**: closed in one of the two classes a transform can
belong to, unpopulated in the other. This is a classification of
transform types, not a theorem about all conceivable mathematics, and
is stated at that strength. **C-III therefore stands at: ① half
closed, ② unattempted, ③ a named external open problem.**

## Round 5 (increment 205): ② settled — the hidden leaf is the bulk

C-III's requirement ② was "a classification covering the type-II
region". Attempting it settles it, negatively for the draft and more
sharply than the review had stated.

The draft's CO bookkeeping needs the μ-side a to satisfy
a ≤ y^{O(1)} = M^{o(1)}. Measuring the absolute Heath–Brown weight
W(a) = Σ_j C(J,j)·A_j(a)·D_j(M/a) across a-sizes
(`code/c3_hb_mass.py`):

| J | a > M^{0.05} | a > M^{0.10} |
|---|---|---|
| 3 | 0.939 | 0.819 |
| 4 | 0.949 | 0.826 |
| 6 | 0.960 | 0.836 |
| 8 | 0.947 | 0.750 |

and raising M to 10⁶ moves these to 0.961/0.865 (J=4) and
0.969/0.849 (J=8). **So ~95% of the identity's weight — rising with
M, flat in J — lies outside the region the draft's classification
covers.** The excluded pieces carry a rough coefficient α(a), a
product of μ-blocks with no divisor structure, hence no Voronoi
entry: §3's spectral door serves only the remaining few percent.

**Verdict on ②**: the classification can be completed, and completing
it shows the tree's dominant leaf is the bilinear μμ core with no
spectral entry. This is scoped to the one-sided Heath–Brown opening as
drafted, over J = 3…8 and two M; a different opening could reallocate
mass, and nothing here rules that out.

**Can the draft be repaired by reparameterising?** No, and the reason
does not depend on parameters. The identity with cut z needs z^J ≥ x,
and its j-th term has μ-side a = a_1⋯a_j ≤ z^j, which at j = J is x —
for every admissible (z, J). So a ≤ z^{O(1)} is compatible with the
identity only if its high-j terms are dropped, which destroys it. And
those are the terms carrying the weight: measured mass by block count
is 0.824 at j = 3 (J=3), 0.967 at j ∈ {3,4} (J=4), 0.806 at
j ∈ {5,6} (J=6), 0.847 at j ∈ {6,7,8} (J=8).

**What a repair would require.** The excluded region is a bilinear
form Σ_a Σ_b α(a)β(b)μ(N−abk) with α rough and a long. Each piece is
individually trivial — for a near M the progression w ≡ N (mod ak)
holds O(1) terms — so all the content is cancellation *across* a,
which is a type-II estimate for μ(N−·): the wall. Cauchy–Schwarz in b,
the only available move, returns the dilate-pair correlation.

This is not special to the one-sided opening. Every identity
decomposing μ produces a type-II term, and an identity whose every
piece had either a long free variable or a divisor-structured
coefficient on the rough side would dispose of the parity obstruction.
**"Complete the draft" and "break the wall" are the same task.**

**A methodological correction, same species as #30.** The first pass
pre-registered a single threshold, a > M^{1/3}, and got a J-dependent
answer (0.326 → 0.032 as J ran 3 → 8) that pointed the opposite way.
That threshold does not match what the claim requires: "a ≤ y^{O(1)}"
means a = M^{o(1)}, i.e. small θ, not θ = 1/3. **Choosing a threshold
that does not match the claim being tested is the same error as
writing a target at the wrong scale** — and it was caught the same
way, by widening to the whole profile instead of trusting one number.

## Round 6 (increment 227): the hazard-4 sweep

Hazard 4 is new, so the same question #30 raised has to be asked again:
**does any standing closure rest on a null that was a size heuristic?**
Every null in `code/` was classified.

| null construction | scripts | exposed? |
|---|---|---|
| Monte-Carlo draws from the model, scored with the *same* estimator | `e1_constr_c1/c2/c2b/c4`, `e1_spectral_null`, `e1_wishart_null`, `e1_forge_kt3/kt4`, `e1_dilation` (permutation), `e4_wheelnull` | **no** — the null is generated, not guessed |
| exact sampling distribution of the statistic, or SE computed from the data as `std/√n` | `sweep_A` A6–A10, `sweep_B` B1–B5, `sweep_B2/B3`, `hyp_round2b`, `lab_gb_multiplicativity`, `e1_forge_kt2` | **no** |
| verdict is directional ("does the ratio approach 1"), with an a priori `(log N)^{−1}` band printed for scale only | `hl_S1_check`, `hl_assembly` | **no**, but the band is heuristic — the verdicts do not turn on it, and both landed at 0.982 and 1.002 |
| CLT exponent 0.5 as the null, with a hand-chosen band [0.45, 0.60] | `sweep_C` | **low** — the exponent is exact under the null; only the band width is judgement, and the reading was "every variant sits at square root" |
| **size estimate `√(M_p log N)` in place of the second moment** | `lab_prime_factor_split` | **YES** — recomputed |
| **Euler-product and `√(N/w) log N` size estimates** | `lab_grouping_asymmetry` | **YES** — corrections #33, #34, fixed in the same session |

**Result of the recomputation.** Session 6's null was replaced by the
exact second moment `V_p = Σ_{v<N,p|v} μ²(v)Λ(N−v)²`. The constant moved
from 0.83 to 0.87 and `S_abs/S_null` remains flat (0.883 → 0.871 across
a factor 8 in N); `S_abs/triv` is null-free and unchanged. The heuristic
sat 5–7% high — except in the top dyadic range, where it returned a
"typical size" of 1.0968 times the trivial bound, which is impossible;
the exact second moment returns 1.0000 there, correctly reporting that
one-term sums admit no cancellation at all. **The verdict stands: the
prime-factor split has positive margin.**

**Nothing re-opens.** And there is a structural reason it was unlikely
to: hazard 4 bites on *magnitude* nulls, whereas the closures that
survived rounds 1–5 are the ones whose criterion was "is there any
signal at all", which is threshold-free and therefore null-shape-free.
The two exposed scripts are both from this week's invention work, where
the questions are quantitative for the first time.

**Coverage, stated.** This swept nulls, not closures. A closure resting
on a *correct* null but a wrong model is untouched by it — and several
E1 verdicts use Conjecture L's own iid model as the null, which is
circular by construction and already recorded as such in `sweep_B`.

## The four hazards this campaign actually suffers from

Both are about **stating a criterion that does not match the thing
being tested**, and both are caught the same way — by widening from a
single number to the whole profile, or by computing the null first.

1. **Scale-normalisation drift.** Writing a target at the square-root
   scale when the chain consumes it at the trivial scale, or the
   reverse. Two instances: the E1 target (correction #30) and the SEAM
   formalisation (REVIEW_VERDICT #3). Rule: state which scale a target
   is normalised against, in the same sentence, and compare it against
   our own measurement of the same quantity before using it to
   adjudicate anything.

2. **Thresholds below the noise floor.** Choosing a pre-registered
   cutoff without computing what the statistic reads under no
   structure. Five instances in a single session: the a > M^{1/3} cut
   (increment 205), the spread-only test for a mask on C(N) (210), the
   |corr| ≥ 0.10 cross-N test and the top-decade test (212), and a
   verdict string that misreported its own rule (211). Rule:
   **compute the null value of the statistic first and put it in the
   pre-registration next to the threshold.**

3. **Identity checks with unenforced side conditions.** Verifying an
   arithmetic identity while imposing only some of its hypotheses, and
   reading the resulting failures as a defect. Three instances in one
   session: the ladder identity needing p² ∤ m (increment 218, ~10%
   false failures); the μ² insertion, where the paper inserts μ² into
   the product *before* the split and we applied it after (221) — this
   one would have become a **false claim of a second defect in a
   published paper**; and the μ–μ dilation needing p ∤ v *and*
   p ∤ (N′−v) (224, 5.2% false failures). Rule: **when checking an
   identity with side conditions, enumerate them all and enforce every
   one before reading any mismatch.**

4. **Nulls taken from size heuristics instead of from the data.** Rule
   2 says compute the null first; it does not by itself say *how*, and
   an a priori estimate of the null is still an assumption. Two
   instances at increment 226 (corrections #33, #34): an Euler-product
   reading gave `A(W) ~ c/log W` where the truth is PNT-strength, and
   `|E_w| ~ √(N/w) log N` overshot the measured residual by a factor 7
   because it is dominated by the w where the model does not apply.
   Rule: **where the statistic's own second moment is computable from
   the same loop, use it as the null.** Print the heuristic beside it
   if it is informative, never in place of it.

Every closure in this repository predating these rules was re-checked
against them; the ones that survived are the ones whose criteria were
"is there any signal at all", which is threshold-free.

## Supersession record (all of it, in one place)

This repository carries no correction banners elsewhere. Every
document states its current position directly; what was withdrawn to
get there is listed here.

| # | Withdrawn statement | Replaced by |
|---|---|---|
| 25 | "gates pass"; "the remainder is one 1/30-wide seam"; "the rest is transcription" | REVIEW_VERDICT.md items 1–3 |
| 26 | a z = 9 spectral excess (iid-entry Wigner null applied to a Gram matrix) | the correct Wishart null, dead centre |
| 27 | "C-III dies even with a free Lemma S" | C3_REVIEW.md as it now stands: draft refuted, route open |
| 28 | Theorem A's bound is N e^{−c√log N} | ≪_A N(log N)^{−A}; BV's internal Siegel–Walfisz range caps it |
| 29 | in Prop. D″, the top-r piece dominates by a power of log N | the pieces are the same order (r₂/r₁ = 0.771 → 0.810); the closure rests on nonnegativity |
| 30 | E1 target `Σ|D|² ≪ (log N)^{−A} Σ M_k` (square-root scale) | `Σ|D|² ≪ (log N)^{−2A−2} Σ M_k²` (trivial scale) — a factor N/K |
| 31 | "no route verdict changes" after #30 | this document: 2 void, 1 downgraded, 2 moot, 13 standing |
| 32 | ③ is "a named external open problem" with a qualitative gap; MRT's averaged Chowla is qualitative | an **exponent** gap: MRT gives ≈(log H)^{−1+o(1)}, Pilatte a fixed power, and the chain wants ≈(log N)^{−6} — provisional, since the 6 is read off Huang–Li (TRANSFORM_LAB §2) |
| 33 | null for the singular-series average: `A(W) ~ c/log W` from the Euler product 1 − f(q)/q ~ 1 − 1/q | `A(W)` tracks `Σ_{w≤W} μ(w)/w`, PNT-strength — measured ~10⁻³ against 1/log W ≈ 0.09. The Euler reading describes the smooth-truncated product, not the Möbius sum |
| 34 | null for the w-grouping residual: `Σ_w|E_w| ~ √2·N log N`, so the ratio to trivial should be ≈1 and flat | the data-driven second-moment null `Σ_w √(Σ_p t_p²)`; measured ratio 1.18–1.26, flat, and 0.26 of trivial and falling. The heuristic overshot by a factor 7 |
| 35 | the criterion "does `S_abs/N → 0`?" tested as `(log N)^{−1/2}` **against** `c + b/log N`, returning DEAD | the right comparison is free-exponent (→0) against `c + b/log N` (→c). Rejecting the *rate* −1/2 is not evidence for a positive *limit*. Corrected verdict at N ≤ 6.4·10⁶: **UNDECIDED** — the two fit equally well (resid sd 0.00081 vs 0.00083) and separate by 0.0035 at N ≈ 5·10⁸, below the noise. What *is* settled: the exponent is ≈ **0.773**, i.e. **faster** decay than the square-root-per-p prediction |
| 39 | increment 243's position: the location mask has **no closed form**, only an enumerated table | it has a **derivation**. n prime and q\|N force q∤(N−n), so v is confined to the integers coprime to rad(N), where the primes (μ = −1, no sign variation) are over-represented and the Möbius sum carries a bias. With `κ(N) = Π_{q\|N} q/(q−1)` **derived** — numerator ψ(N) ~ N, denominator N·Π(1−1/q) — the prediction `C(N) ≈ κ₀·κ(N)·M_{rad N}(N)` has **one** free constant and reaches weighted R² **+0.569** on the cell means, against +0.219 for a fitted nine-parameter additive model. Not complete: it under-predicts the deepest cells by ~1.7× |
| 38 | Conjecture L's mask applies to the **scale** of a μ-family | it applies to the **location** too, and there it is the larger effect. The mean of `C(N)/√(V(N))` is a deterministic function of *which* small primes divide N — not of 𝔖(N), and not additively: an additive one-coefficient-per-prime model predicts −0.49 for N = 510510 where the truth is −8.99, a factor 18. Full modular enumeration over q ∈ {3,…,23} takes the excess kurtosis from +0.468 to **+0.014**, every tail inside 3 SE, and the extreme from z = +16.8 to **+0.61** |
| 37 | increment 237's reading: the wall's fluctuation has a **heavy tail**, refuting Conjecture L's Gaussian half | a **location** effect, not a tail one. The outliers are primorials — 510510 = 2·3·5·7·11·13·17 and its relatives — and **all twelve have C < 0**. When n is prime, N−n is forced coprime to rad(N), so v is drawn from integers with no small prime factor, ω(v) shifts down, and session 8's alternating balance tilts toward ω = 1 where μ = −1. Centring on an 𝔖(N)-dependent mean removes **96%** of the excess kurtosis (+0.468 → +0.019) and the extreme from z = +16.8 to **+2.4, which passes**. Conjecture L needs a **mean** term — a local mask on the *location*, not only on the scale |
| 36 | the wall's law `C(N) = √(𝔖(N)·N)·G(N)` with `G ∼ N(0,1)` | **`Var C(N) ≈ 0.465·𝔖(N)·N·log N`** — the recorded law is missing a factor `√(log N)`. Measured over every even N ≤ 4·10⁶: sd(G) drifts 2.069 → 2.306, and `sd²/log N` is flat to **0.5%** (fitted exponent 0.895; α = 0 gives CV 7.5%, α = 1 gives 1.0%). The normaliser that needs no fitting is the exact second moment `Σ_v μ²(v)Λ(N−v)²`, against which the measured variance is 1.01 → 0.87. **Scale-normalisation drift, third instance** |

Two of these are the same species — #30 and REVIEW_VERDICT #3 (SEAM's
√P over-normalization) are both **scale-normalisation drift**: writing
a target at the square-root scale when the chain consumes it at the
trivial scale. Both were caught only by comparing a stated target
against our own measurements of the same quantity. That comparison is
now mandatory before any target is used to adjudicate anything.

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
