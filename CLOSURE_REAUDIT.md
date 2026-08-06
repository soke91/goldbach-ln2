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

## The six hazards this campaign actually suffers from

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

5. **Reading a declining trend as convergence to the null, when the
   statistic has no power to distinguish.** Four instances, all on the
   same object and all found only at increments 248–254: `sweep_B`'s B4
   sign balance (z = −1.55 on ~1500 values; **−94.3 and +61.6** on
   1.95·10⁶); B5's `corr(C,𝔖)` falling −0.194 → −0.138 → −0.070 across
   three windows and read as convergence (**z = −193.5** at full power,
   and **−137.4** after the very fix that was supposed to remove it);
   the "mean drift, already insignificant at N ≈ 9·10⁵"; and §10's
   `3 | N` split, where "𝔖-normalisation takes it below 3σ" was n ≈ 500
   and full power gives z = −70 after that same normalisation. Rule:
   **a trend toward the null is not evidence of the null unless the
   last point has the power to reject it.** Quote the power, not the
   trend.

6. **A step that silently does not run, while the record says it did.**
   The most dangerous of the six, because nothing looks wrong. Two
   instances: a pre-registered control prime dropped by an `n ≥ 12`
   guard that printed nothing (increment 243), and a document rewrite
   whose script raised `AssertionError` on a string that did not match
   while the surrounding shell went on to commit a message describing
   the rewrite as done (increment 254, corrected at 255). Rule:
   **chain the edit to the commit so a failure aborts both, and make
   every skipped branch print the fact that it skipped.** A control
   that does not announce its own absence is worse than no control.

   **It then recurred three times immediately after being named.**
   `sed -i 's/^\*Increment 255.../.../'` matches nothing when the file
   says 254, and `sed` exits 0 on no match, so increments 255–257 each
   reported a bump that did not happen and STATUS.md sat at 254 for
   three commits. The rule above is necessary and not sufficient: a
   chained command that *succeeds while doing nothing* defeats it.
   **Every edit must be followed by a read-back that asserts the new
   state**, which is what the fix at 258 does. Naming a hazard does not
   prevent it; only a check that fails loudly does.

   **A third form: a check that cannot fail.** Increment 272's
   identity test printed a "C(N) − Λ(N−1)" line built by rearranging
   the very quantity it was meant to confirm, so it agreed to all
   digits by algebra and verified nothing. It also compared
   `Σ_p log p D_p` (weight 1) against ω-classes built with a `1/log v`
   weight — two different objects, and the mismatch showed up as a
   1.5·10⁴ "failure" of a true identity. Both fixed; with the weights
   matched and a genuinely independent third computation the identity
   holds to 6·10⁻¹⁰. Rule: **a verification must be able to come out
   false.**

   **A fourth form: a shell heredoc collapsing `\\` to `\`.**
   At increment 280 a patch reported success and wrote `\text` into
   STATUS.md as a literal TAB and `\asymp` as a BEL, because the
   heredoc ate one level of escaping before Python saw it. The
   assertions passed — they tested for the *surrounding* text, not the
   escapes — and only reading the file back caught it. Repairing it
   exposed that this very paragraph had been spliced into the middle of
   a sentence by an earlier edit of the same kind, unnoticed since.
   Rule: **assert on the exact bytes that were the point of the edit**,
   and prefer a file-writing tool to a shell heredoc for anything
   containing backslashes.

   **A fifth form, and it broke the rule written one increment
   earlier: `open(path, "w")` truncates before it can fail.** Writing
   `𝔖` (U+1D516, outside the BMP) as a surrogate pair `𝔖`
   raises `UnicodeEncodeError` at the `.write()` — but the file has
   already been emptied by the `open`. This destroyed FINDINGS.md at
   increment 280 (recovered from context; it is deliberately outside
   git, so there was no other copy) and CLOSURE_REAUDIT.md at 281
   (recovered by `git checkout`). Rules: **write to a temporary path
   and rename**, or use the editing tool rather than a shell heredoc —
   which is exactly what the previous paragraph said to do, one
   increment before it was ignored. Naming a hazard does not prevent
   it, and neither does having just named it.

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
| 86 | increment 288 left `ρ → 1?` unsettled and said more computing would not settle it | **identifying `ρ − 1` settles what it is, if not yet what it tends to.** Expanding `C(N)²` with `u = N−p`, `h = p′−p`: **`ρ − 1 = (1/V)Σ_{h≠0} c(h)S(h)`** with `c(h)` a weighted prime-pair count and `S(h)` the **binary Chowla correlation** (**Proposition W**). So the wall's excess over square-root *is* a prime-pair-weighted Chowla correlation, and Chowla-type input gives `ρ → 1` — **the wall is exactly square-root, over-delivering by a power of `log` and no more.** Three pre-registered quantitative tests, all passed: **(A)** `M(h)` sits at **1.051–1.068×** the random-sign floor `√(0.32264(X−h))` — not `√X`, since the sum only sees `n` with *both* `n` and `n+h` squarefree — stable across five decades of shift; **(B)** the reconstructed `ρ−1 = −0.0976` against the measured −0.17…−0.19, a factor 0.54 (bar: 3×; the gap is the dropped prime-power tail and the uniform-`u` approximation); **(C)** the sign is negative. ✅ **And the wall leans on the provable end of Chowla**: shifts `h < 10³`, where the conjecture is hardest and the averaged theorem weakest, carry **1.1%** of the gross mass; 48.9% sits at `h ≈ 10⁵–10⁶`. ⚠️ The ranges cancel 2.0× (net `−2.2·10¹³` against gross `4.4·10¹³`), so shares of the *net* would mislead and the first draft's table — which printed them, reaching 97.6% and −47.5% — was replaced |
| 85 | *(what is **not** settled, stated so it is not overclaimed in the other direction)* does `ρ → 1`, i.e. is the wall exactly square-root? | **not established, and the pre-registered rule says so.** Three two-parameter models on the same points: `ρ = c` (RMS 0.0258), `ρ = c(log N)^{−β}` (0.0060), `ρ = a − b/log N` (0.0049). The best is not below half of each rival, so the run returns **INDETERMINATE** — and every parameter still **walks**: `a` runs 1.174 → 1.088 → 1.063 → 1.070 → 1.054 → **1.039** across the window, `β` runs −0.501 → −0.286, `c` runs 0.782 → 0.810. ⚠️ **An upward drift is not evidence for 1 any more than a downward drift was evidence for 0** — hazard 5 applies symmetrically, and this entry exists because the previous one is the kind of finding that invites exactly that overclaim. What *is* established is the **direction**, and the direction says the extra cancellation is disappearing, not growing. ✅ Unforced cross-check: `sd(C) = √(ρV)` recomposes increment 281's independently fitted `β = 0.5457 ± 0.0032` as `0.5348 + 0.0100 = 0.5448`, difference **−0.0009** — different scripts, different fitted quantities, so this is evidence and not the algebraic identity flagged at #71. The decomposition also locates the exponent: **0.5348 of it is the arithmetic scale `√V` and only 0.0100 is the `ρ` trend** |
| 84 | the wall's cancellation ratio `Var C / Σ_v μ²Λ²` was recorded running **1.006 → 0.873**, read as the wall beating a coin by a growing margin | **the sign of that trend is wrong: it was the mask.** With the location mask removed by the same enumeration used since #67, the ratio **RISES**, 0.760 → 0.837, while the raw ratio falls 1.006 → 0.858. The two move in **opposite directions and converge** — gap 0.246 at `N ≈ 10⁵`, **0.020** at `1.6·10⁷` — exactly as finding #69 requires, the mask decaying like `N^{−1/2}` relative to the fluctuation. So `ρ < 1` is real and large (the off-diagonal is genuinely negative, the wall does beat a coin), but **the excess over square-root is shrinking with `N`, not growing**. The question only became askable at increment 287: until Proposition V the denominator in use was a fitted stand-in (#74, #75, #83), and a ratio to a wrong denominator has no direction worth reading |
| 83 | *(the sharpest correction of the session, and it explains three earlier ones)* this program has used `𝔖(N) = 2C₂∏_{q\|N,q>2}(q−1)/(q−2)` as **the** local factor everywhere, including for the wall's scale | **for the wall's scale it is the wrong function, by a factor of 760.** Increment 283 left `V(N) = Σ_v μ²(v)Λ(N−v)²` exact but opaque; expanding it gives a `(log p)²`-weighted count of **squarefree shifted primes**, whose local density is `1` at `q\|N` (since `q²\|N−p` forces `p=q`) and `1 − 1/(q(q−1))` at `q ∤ N` (a unit class mod `q²`). So the factor is `𝔄(N) = ∏_{q∤N}(1 − 1/(q(q−1)))` — Mirsky 1949, **recalled not claimed** (**Proposition V**). Measured on every even `N ≤ 1.6·10⁷`, with the analytic part removed *exactly* by dividing through `W(N) = Σ_{w<N}Λ(w)²`: `R/𝔄 → 1.000000` with sd falling 0.00154 → 0.000145, and per radical cell `R/𝔄 = 1.00000` to five decimals for all six cells tested. Rescaled to the same mean so only *shape* is judged: residual sd **0.000323** for `𝔄` against **0.245235** for `𝔖`. ⚠️ **The singular series of the problem is not the local factor of the noise.** `𝔖` is Hardy–Littlewood's, correct for the Goldbach *count*, and was carried across to the *variance* of `C(N)`. This one mis-identification explains #74 (93.7% of `(𝔖N)/V`'s variance is cell-explained — `𝔖` is the wrong cell function), #75 (`α` would not converge because a log-power was absorbing a wrong **arithmetic** factor), and why `κ ≈ 0.465` still fitted on average (`𝔖` and `𝔄` correlate) |
| 82 | *(a second error of my own, in the audit written to find the first)* the first draft of `audit_zero_support.py` explained why §6's dispersion figure 0.801 was undamaged by arguing that **"a number sitting at 0.80 cannot have been depressed by a large factor"** | **that assumes the conclusion.** It is the shape of reasoning this table exists to catch, written *inside the tool built to catch it*. Replaced with a measurement: the same family summed over **all `m`** — which is what `dispersion_engine.py` actually does — gives an annihilation rate of **0.0%** and `(a)/(b) = 1.0000` at both bands. §6's number is safe, and now for a checkable reason. Also reclassified: `e1_seam_law.py` is listed unguarded by the scan but is **instrumented** — its line carries `# 0 when killed` and the next statement *prints the annihilation rate*. The scan looks for filters and misses a published diagnostic, so `NO` in that column means "no filter", **not** "nobody noticed" |
| 81 | *(the fault of #80, swept across the corpus)* `abs(Σ)/√max(support,1)` is a **code idiom**, not a one-off, and this program's oldest recurring fault — a mask term read as a measurement — is exactly it | **37 sites in `code/`, 27 with no filter — but syntax is the wrong discriminator, and the right one is measurable.** The same statistic under three conventions, `N ≈ 10⁷`: over **primes `p`** the annihilation rate is **45%**, giving `(a) = 0.42` against `(b) = 0.78` — the fault is worth a factor **0.55**, and `(a)/(b) = 1 − (annihilation rate)` by construction, so **it is exactly as large as the mask is dense**. Over **all `m`** the rate is **0.0%** and the ratio is **1.0000**. So the idiom is harmless in dense summation ranges and severe in sparse ones; the static scan can prioritise but cannot judge. This also places §6's raw "0.348–0.370, stronger than random": it is the `(a)` convention, and §7's zero-accounting (`e1_zero_account.py`, convention `(b)/(c)`, giving 0.97–1.02) had already superseded it — **the correction was made; what was missing is that the documents never said which convention any number used**. `results/audit_zero_support.txt` now does |
| 80 | *(my own error, twice, inside the fix for #78)* the repaired stamp immediately failed on V4 (seam band, `r = 0.401` against half-normal 0.798), and **I diagnosed it wrong twice** | **first**: "the interval was set by analogy" — true, and hazard 4, but not the cause of 0.401. **Second**: "unlike V5, V4 pools the variance-suppressed shared classes" — plausible, documented in §6/§7, and **refuted by splitting them**: the shared class runs at `r = 0.808`, half-normal. **The actual cause**: the old line `abs(vals.sum())/sqrt(max(v,1))` fed **zero-support pairs into the mean as `r = 0`**. Those are M.3's predicted annihilations — pairs about which the law claims nothing — and counting them as *measured zeros* drags the average down. **45% of pairs (118 of 265) have support < 50**, and the arithmetic closes: `0.55 × 0.81 ≈ 0.45` against the old 0.401. This is the campaign's oldest recurring fault — **a mask term read as a measurement** — sitting in the CI stamp. The free/shared split is kept anyway, since matching V5's discipline is right and it is what refuted diagnosis two |
| 79 | V1 of the CI stamp drew its 30 test values as `rng.integers(...) // 6 * 6 + 2` | **every one of them is `N ≡ 2 (mod 6)`, so `3 ∤ N` always.** The location mask is largest exactly at `N` divisible by many small primes — increment 240 measured deep cells at −5 to −7 sd — and those `N` were **structurally excluded from the stamp**, which would therefore have passed unchanged if the mask were catastrophic. A **deep arm** (`30030 \| N`, so `3·5·7·11·13 \| N`) now runs alongside, and the two arms differ visibly: `r = 0.659` shallow against **`0.935` deep**. ⚠️ At 30 `N` per arm the standard error is ≈0.11, so that gap is ~1.8 s.e. — **suggestive, not established**, and it is recorded as the reason the arm exists rather than as a result |
| 78 | *(the fault this program named at increment 272, sitting in its own primary certification)* `verify_all.py` — what STATUS calls "the heavy CI stamp for the measurement corpus" | **it could not come out false.** Zero `assert`, zero `sys.exit`, zero `FAIL`: the criteria ("기준 0.80", "기준 1.0 ± 0.1") were **characters inside output strings**, not code, and the script exited 0 whatever the numbers were. Every V-stamp is now a `(name, text, value, pre-registered interval)` row judged mechanically, printing PASS/FAIL and exiting 1 on any failure, with a **sensitivity block** that pushes each value out of its interval and shows the verdict flips — asserting a check can fail is worth nothing (increment 276). Also fixed: the script died on `UnicodeEncodeError` under a cp949 console, which turns a stamp into an *exception* rather than a *failure* and erases the verdict entirely. ✅ **The one thing the audit expected to find was not there**: V1's `Σ log²p·[μ≠0]` and V5/V6's directly-counted support are **exact** second moments, so increment 283's `𝔖N` stand-in never reached certification |
| 77 | *(a prediction upgraded from "all observed" to exact)* M.3 — `q²\|gcd(k,N)` ⟹ `D(k) = 0` identically — was recorded as "1212 predicted, all observed" | **checked in both directions on 401,000 pairs and it is exact.** Predicted zeros 115,950; observed zeros 115,950; predicted-but-nonzero **0**; unpredicted zeros **0**. The second direction is the one that had never been tested: an unpredicted zero would be a *second* annihilation mechanism, and there are none. ⚠️ Sample note recorded in the output rather than left to be misread: every `N` on this grid is divisible by 4 (both `3·10⁶` and the step 2500 are), so `4\|gcd(k,N)` alone kills a quarter of the `k` and the 28.9% zero rate is a property of **this grid**, not of the field |
| 76 | *(the first time this audit found the campaign right)* Conjecture L's **original** half — `μ`-field = mask × exactly-Gaussian — was stamped at "kurtosis 2.99–3.03" on ~12,000 pairs, i.e. `±0.045`, while increment 283 used `±0.0017` to break the `C(N)` extension | **it survives at 26× the precision.** On **285,050** pairs: excess kurtosis **−0.0034 (z = −0.4)**, skew −0.0075 (z = −1.6), `E\|Z\|/sd = 0.79760` against `√(2/π) = 0.79788` (**z = −0.2**), variance ratio 0.99781. **No class structure** either: `gcd(k,N)` splits give `\|z\| ≤ 2.5` over four classes, variance ratios 0.978–1.020, consistent with sampling. **And the reason this half was safe is instructive**: repeating it with a *band-mean* support in place of the exact per-`k` count gives excess kurtosis **+0.4645 (z = +50.6)** and `E\|X\|/sd` off by −13.1σ, because support varies 40.4% across `k`. §7 counted the zero terms directly and stepped over the trap; the `C(N)` extension used a fitted `𝔖N` scale and fell in. **Same program, both halves, and the difference is whether the exact quantity was available *and used*** |
| 75 | *(a puzzle dissolved rather than solved)* increment 280 could not identify `α` in `Var C = κ𝔖N(log N)^α` — it walked with every fitting window and no reachable `N` would settle it | **because there is no `α` to identify.** `κ𝔖N(log N)^α` is a **fitted stand-in for an exact quantity**: the wall's scale is `V(N) = Σ_v μ²(v)Λ(N−v)²`, which needs no fit at all. The stand-in matches `V` in the *mean* at reachable `N` (`κ log N ≈ 6.5` against the measured `V/(𝔖N) ≈ 6.05`), which is why the fit looked fine, and gets the *fluctuation* wrong, which is why the exponent would not converge. ⚠️ **Increment 238 already said this** — *"the normaliser that needs no fitting is the exact second moment"* — and the conjecture was then written with `𝔖N` anyway, and three later increments spent effort fitting an exponent of a quantity that was never the right object. **A fitted law standing next to an exact one is a choice, and this program made it without noticing** |
| 74 | Conjecture L claims the wall's fluctuation is **"exactly Gaussian … kurtosis 3"**, and the campaign displays `C(N) = m(N) + √(κ𝔖(N)N log N)·G(N)` | **the claim is true and the displayed formula is false; they use different normalisers.** Every even `N ≤ 1.6·10⁷`, cell means removed, band standardised: under the displayed `𝔖N`-based scale the excess kurtosis is **+0.1704 (z = 98.1)** and `E\|X\|/sd` misses `√(2/π)` by −0.0062 (z = −29). Under the **exact** `V(N) = Σ_v μ²(v)Λ(N−v)²` it is **−0.0005 (z = −0.3)**, with `E\|X\|/sd` off by −0.00018 (z = −0.8) — Gaussian to a precision this program has never previously reached, on 6.3·10⁶ values, needing cell **means** only. A depth ladder rules out mask leakage (going from primes ≤ 11 to ≤ 31 removes **−3.3%**, i.e. nothing); per-cell standardisation removes 101%, and the reason is that `(𝔖N)/V` varies with sd 0.049 about a mean of 0.165 with **93.7% of that variance cell-explained**. Increment 240's recorded "+0.014 after masking" was measured under `V` and is right; every document that displays the `𝔖N` form inherits a normaliser under which the conjecture's own central claim fails |
| 73 | the criterion statistic behind that closure, `corr(C,D)`, was reported as **0.60 pooled** across five groups, and the supporting exponent as `\|R\| ~ N^{0.599}` | **both are wrong, and one of them in the direction that flatters the closure.** Pooling across bands of different scale *attenuates* a correlation; per band the de-masked value is **0.80 falling to 0.71**, and pooling every band of the full census reproduces the same artefact (0.7064). The true correlation was never as low as 0.60. Removing the location mask **raises** the correlation by **+0.042** — the mask was *suppressing* the association, the opposite of the natural guess that a shared deterministic term inflates it. Corrected exponent: **`\|R(N)\| ~ N^{0.6458 ± 0.0065}`**. Independently, `β_C = 0.5461 ± 0.0033` here reproduces increment 281's `0.5457 ± 0.0032` by a **different estimator** (`E\|·\|` de-masked vs `sd`), which is the one cross-check in this pair that is not forced |
| 72 | *(a recorded CLOSURE re-opened)* §9 refuted "the ln 2 comet and the final scalar are the same object" on a pre-registered criterion, citing `\|R\| ~ N^{0.599}` against `\|C\| ~ N^{0.503}` | **the closure stands; the evidence recorded for it did not.** Full census of every even `N ≤ 1.6·10⁷`, mask removed: `β_R − β_C = +0.0997 ± 0.0034`, **29 s.e.** — R really does outgrow C. But replicating §9's own design 500 times (five groups of 80 *consecutive* even N, offsets varied) gives that difference as **+0.079 ± 0.054**, so the recorded gap 0.096 was **1.8 s.e. of that design** and **7.0% of replicates return a difference ≤ 0**. A refutation was recorded on evidence that would have come out the other way one time in fourteen. ⚠️ **What actually makes the closure safe is the one thing §9 did not check**: the de-masked correlation **falls**, `−0.0198 ± 0.00035` per unit `log N`, moving *away* from the 0.9 threshold. §9's five noisy groups appeared to show it *rising* (0.457, 0.746, 0.771, 0.831, 0.726); had that been real, the closure would have been premature. **Right conclusion, inadequate evidence, and the adequacy was decided by luck** |
| 71 | *(a check I nearly shipped, caught before the run finished)* increment 281's consistency test compares the amplitude exponent β against the one implied by §12's variance exponent α | **that comparison is near-algebraic and confirms nothing.** β is the slope of `log sd(C)` on `log N`; α is the slope of `log sd(Z)²` on `log log N`; and `sd(C) = sd(Z)·~√(𝔖N)`. They are **the same fitted numbers in two parametrisations**, so `β = 0.5 + (α/2)(Δloglog/Δlog)` holds by construction — and the run duly returned agreement to **6·10⁻⁴** against α = 1.30, which reads like strong corroboration and is worth nothing. **This would have been the fifth check in this program that could not come out false.** The output now says so *above* the numbers. The one informative row is `α = 0`, which the old law predicts and the data refute; a genuinely unforced diagnostic is `E|C|/sd = 0.79145` against the Gaussian `√(2/π) = 0.79788`, a consistent **−0.81%** with the sign expected from the residual excess kurtosis |
| 70 | `\|C(N)\| ~ N^{0.503}` — **"square-root, to three digits"** (§9, and the headline of STATUS §3b) | **the precision is refuted and the exponent is not 0.503.** §9's design is five groups of **80 consecutive even N** (each group spans 160 in N); a scale from 80 samples has sampling CV ≈ 7.9%. Re-running **that exact design** 500 times with varied offsets: **β = 0.516 ± 0.0426**, 5th–95th percentile [0.447, 0.585], and the recorded 0.503 sits at **percentile 39** — an ordinary draw. The replication spread is **43× the 0.001** that "three digits" asserts. Full census of every even `N ≤ 1.6·10⁷` with the mask removed gives **β = 0.5457 ± 0.0032**, so `β = 1/2` is **excluded at 14.5σ**, and the mask alone shifts β by +0.014. `\|C(N)\|` grows **measurably faster than √N** — which is what §12's own variance law requires, and §9 was read as confirming the opposite. **Nothing about the program's position changes**: the requirement was only `o(N)`. What changes is that a number quoted to three digits had three digits of noise |
| 69 | *(a finding, not a correction)* the location mask's **own** scaling in N had never been measured — only its shape at fixed N | **`Var_mask(Z) ∝ N^g` with `g = −0.489 ± 0.005`, walking `−0.443 → −0.489` as the window grows — i.e. `g → −1/2`, so `m(N) ≍ √𝔖(N)·N^{1/4}`.** The value is still drifting but **the sign is not** (~100σ), and the sign is what matters: the mask is **lower order than the fluctuation** and therefore does **not** threaten `C(N) = o(N)`. With `E₃` already known to cancel 82% of it in the Goldbach count, the mask is harmless twice over — a real feature of the wall, not an obstruction. It also explains the bias in #67 mechanically: a term whose variance share falls like `N^{−1/2}` must bias a fitted log-exponent downward by an amount that shrinks with the range, which is why the raw and de-masked fits close from 0.70 apart (3 bands) to 0.29 apart (8) |
| 68 | the wall's variance exponent has been quoted as `log N` to the first power, with `κ ≈ 0.465`, in CONJECTURE_L.md and MEASUREMENTS.md | **the exponent is not identified, in either direction.** Fitted `α` **walks with the window**: 0.83→1.02 raw, **1.60→1.30 de-masked**, neither converged, both drifts ~10× the standard error. ⚠️ The de-masked fit reports `z = +8.58` against `α = 1` — **quoting that as ‘α=1 excluded’ would have been a false ALIVE of exactly the kind this table exists to catch** (hazard 5, turned on the estimate instead of the null). Worse for the form: over this range a pure power `N^ε` fits **as well or better** (residual ratio 0.45 raw, 1.91 de-masked, against a bar of 2 fixed before the run), so the log factor is *consistent with* the data, not *shown by* it. **What survives is `α ≠ 0`** — the variance grows faster than `𝔖N`, which is #36's actual content. And it will not be settled by computing: `log log N` spans 0.33 over a factor 160 in N, so halving the standard error needs `N ~ 10¹⁰` |
| 67 | the variance law was fitted at increments 236–238; **the location mask was found at increment 240** — two increments later — and the law was never re-fitted | **the fit was contaminated by a term this program itself discovered.** `sd(G)` is `g.std()`, which removes one band-wide mean and nothing else, so `m(N)` sat inside the measured variance throughout. Removing it by the same modular enumeration that defines it (cells keyed by which of `{3,…,23}` divide `N`, per-cell means subtracted within each band, exact `n−k` dof correction so the removal cannot shrink the variance by itself): the mask supplies **14.1% of the variance at `N≈10⁵`, 1.15% at `1.6·10⁷`**, and because that share *shrinks*, leaving it in **tilts the exponent by 0.29** at full range and 0.70 on the first three bands. Found by checking commit dates, not mathematics — **a measurement is only as current as the last thing discovered before it** |
| 66 | *(an observation nobody asked for, and it reframes every no-go here)* from what N does Theorem D′ actually **bite**? | **N ≈ 10⁴⁸⁰ at Lichtman's level, and further out below it.** Solving `c√(η log N) = A log log N` (c=1, A=3): `η=0.40` crosses at `N≈10⁴⁸⁰`, `η=0.10` at `10³⁰⁷¹`, `η=0.01` at `10⁵³⁷⁴⁴`. Below those the switch route loses **nothing that matters** and the theorem is silent. Two independent solves agree (required η at `10⁵⁰⁰` is 0.388, just under 0.40; the 0.40 crossover is `10⁴⁸⁰`, just under `10⁵⁰⁰`). This does not weaken any no-go — they are asymptotic statements and are true — but **they constrain METHODS, not any computation anyone will run**, and the documents were easy to read as though they bit at reachable N. Now stated wherever they appear |
| 65 | the prose under Theorem D′ said closing the gap *“would require **θ_E = 1 exactly** — … each progression holds **O(1) terms** and the statement carries no information”* | **wrong on both halves.** The theorem itself is correctly hedged (*for each **fixed** θ_E<1*); the prose then skipped the regime `θ_E = θ_E(N) → 1` **at a rate**. Solving instead of gesturing: `exp(c√(η log N)) ≤ (log N)^A ⇔ η ≤ (A/c)²(log log N)²/log N` (**Proposition D⁗**), so the boundary is `1−θ_E ≪ (log log N)²/log N` — **strictly stronger than EH, strictly weaker than θ_E=1**. And that regime is **not vacuous**: progressions retain `exp(C(log log N)²)` terms, which beats every fixed power of `log N`. ⚠️ This prose stood from increment 196 and was used at increment 235 as **one of three grounds for closing the program's highest-value open item**. The closure survives — Lichtman's `3/5` is a fixed constant and reason 1 (object mismatch) is independent and decisive — but *“vacuous”* says the target is meaningless where *“stronger than EH”* says it is meaningful and out of reach, and only the second leaves a well-posed question |
| 64 | *(a robustness result, not a correction)* Theorem D rests on one unconditional estimate, Huang–Li's Lemma 1 `ρ_n(x) ≪ e^{−c√log x}`, which invites the objection that the no-go **records our ignorance and would dissolve if the estimate improved** | **it runs the other way, and the proof is one line.** The paper already displays `|B_w| ≤ ‖b‖₁·max_d|ρ_{dN}(K/d)|` (`theorem_A.tex` l. 833), so the loss factor is `≥ 1/max_d|ρ|` — **monotone decreasing in any upper bound for ρ**. Every improvement of Lemma 1 strengthens Theorem D; none can weaken it (**Proposition D‴**). Since `1/ζ` enters to the first power, RH gives `ρ_n(y) ≪ y^{−1/2+ε}`, and `K/d ≥ N^{1/2+δ}` turns the conclusion from `exp(c√((1/2+δ)log N))` into **a power of N** — at N=10⁵⁰, 3·10¹² in place of 44. `lab_rho_decay.py` measures ρ's decay directly (discrimination rule fixed before the run): deep n give **H_pow** decisively (RMS 0.081 vs 0.202 at n=2310), shallow n come back INDETERMINATE **because ρ changes sign** there — the signature of a ζ-zero-governed tail, not regime ambiguity. ⚠️ The fit is **not** evidence that β=1/2 is the truth and the deep-n cleanliness is plausibly pre-asymptotic; **the robustness conclusion does not depend on the fit**, only on the monotonicity |
| 63 | LITERATURE.md's placement table said Theorem D's new content was making the classical genre precise *“for this reduction, **under full EH**”* | **wrong, and it misstates the theorem in the direction that flatters it.** Theorem D's ingredient list is **Bombieri–Vinogradov at `Q = N^{1/2−δ}`**; full EH is what **Bombieri's asymptotic sieve** assumes, not what Theorem D assumes. The two were conflated because increment 233 placed Theorem D “in the genre of” Bombieri's sieve and never compared the hypotheses. Laying them side by side shows they **do not overlap**: Bombieri's sieve is positive, at level 1, and needs the weight to obey parity; Theorem D is negative, at every `θ′ < 1`, and quantifies over all weights. Found only because Proposition D‴ required reading both statements at once — **a claim about the literature had gone 45 increments without anyone re-reading the theorem it described** |
| 62 | *(a check, not a correction)* the (18) defect report has been held back on the reasonable ground that **the authors may already know** | **checked, and there is no sign that they do.** arXiv:2005.03811 has exactly two versions, v1 (2020-05-08) and v2 (2022-08-28); v2's revision note is *"title changed, typos corrected"*, not a mathematical correction, and it is the latest. The paper is a Springer book chapter (CANT 2020, doi 10.1007/978-3-030-67996-5_17). No erratum and no citing work noting the issue was found. So the step stands uncorrected roughly four years on. **Caveats**: absence in a web search is weak evidence for a book chapter, whose errata are poorly indexed, and the authors may know privately. **Nothing has been sent** — that is the author's decision, and the check is recorded in the report's own header so it travels with it |
| 61 | *(a structural response, not a correction)* four times in this run a verification was shipped or nearly shipped that **could not come out false** — the silent heredoc patch (259), the algebraically forced identity line (272), the audit's own false alarms (274), a hardcoded `True` inside the verifier (275) | **asserting that a check can fail is worth nothing; showing it is worth something.** `verify_propositions.py` now re-evaluates each identity with one side perturbed by 1 part in 1000 and requires the verdict to flip. All five do, each residual landing ≈1000× its tolerance (M.2, P.1, P.3b: 10.95 against 0.011; P.2, P.3a: 206.5 against 0.207). The four counting checks test `violations == 0`, which a single injected violation flips by construction. **A verifier nobody has run against a known-false input is not known to work** |
| 60 | *(an audit that came back clean, plus a fault in the audit itself)* were the load-bearing figures now standing in the documents transcribed correctly from `results/`? | **yes — 24 of 24, and 3 derived figures recomputed within tolerance.** Transcription is a fault class that nothing else here would catch: a number can be measured right, recorded right, and quoted wrong. ⚠️ The audit's **first run reported 8 false alarms out of 24** because the documents typeset the minus as U+2212 and often quote a rounded form, while the check compared raw ASCII substrings — **the very fault class the tool was built to find, in the tool**. Fixed by normalising the sign and accepting any rounding. Stated limit, in the output: it checks that a figure was not invented or mistyped, **not** that it is the right figure for the claim it supports |
| 59 | *(a derivation, not a correction)* the cancellation nulls used throughout §4 and increments 267–269 model `D_p` as **independently signed** terms | **the signs are forced.** Proposition P.7: with r the least prime not dividing N, for `(N−1)/p < r²` the terms of `D_p` are `v = p` with `μ = −1` and `v = qp` with `μ = +1` — **one negative against a sum of positive**. Verified with 0 genuine violations in ~70 000 checks, the 21 exceptions all being the `N−v = q^k` case P.5 states. So "how far from random" is the wrong question at these p; the right one is the balance `pos/neg`, which reads **0.021, 2.61, 7.63, 19.52** at 1, 2–3, 4–8, ≥9 terms for deep N against 0.558, 2.08 for a typical one. Cancellation lives only in the narrow window where that ratio passes 1 |
| 58 | *(a derivation, not a correction)* does P.5's floor refute P.4? | **No.** Corollary P.6: the floor is `2A(N)N/(√δ·log N)` with `δ = φ(rad N)/rad N` — predicted 0.2577 against a measured 0.2699 for the deep group. That is `Θ(N/log N)` typically and `Θ(N√(log log N)/log N)` at maximal depth, **`o(N)` in every case**. Depth raises the floor by `1/√δ` = 2.28, 2.88, 4.05, 4.96 as `rad N` runs over the primes up to 13, 10², 10⁴, 10⁶ — unbounded but very slowly growing, and never breaking. It also **derives §5's ceiling in general**: §5 got `C(N) ≍ N/log N` from `j = 1` alone, P.5–P.6 get the same order from every `j` at once |
| 57 | *(not a correction — a derivation that closes increments 264–268)* the depth effect on P.4 was measured but unexplained | **Proposition P.5**: `Λ(N−v) ≠ 0` forces `(v, rad N) = 1` (M.1), so with `v = mp` the nonzero terms of `D_p` are exactly those with `(m, rad N) = 1`, giving `(φ(rad N)/rad N)(N/p)` of them. The floor on P.4's demand therefore scales as `√(rad N/φ(rad N))` — **2.284** at `N = k·30030` against ≈1.41 typically. Measured floor ratio **1.69** against a predicted **1.62**. And at deep N every class with `j = ⌊N/p⌋ ≤ 15` has `S/M = 1.0000` exactly, because `m = 1` is the only `m ≤ 15` coprime to 30030: **39% of the demand is single-term and structurally irreducible** |
| 56 | *(a pre-registered prediction of mine, refuted in direction)* the deep-N excess in `Σ log p\|D_p\|` sits at **large p**, since `D_p` is a single forced-sign term there | **it sits entirely at small and middle p.** By band, deep N gives `S/S_sign` = **14.39** at p ∈ [16,32] falling through 10.03, 7.17, 4.01, 2.59, 1.42 to **exactly 1.0000** above p ≈ 3·10⁴; all even N shows no excess anywhere (0.84–1.04). The reasoning failed because a **one-term sum has ratio exactly 1** under a sign-randomised null — the null uses the same magnitude, and sign is irrelevant to a single absolute value. Two consequences: **the deep-N excess IS the mask**, at precisely the p where increment 263 found the removable part; and for `p > N/2` the statistic is **structurally blind**, which is the same fact as §5's ceiling `Σ_{p>N/2}Λ(N−p)` |
| 55 | two specification faults in this program's **own** depth ladder (increments 264–266) | **(a) the null.** `R_null = Σ log p √V_p` is the right random-sign level only where `D_p` has one term; where it has many, `E\|D_p\| = √(2/π)√V_p = 0.798√V_p`. Replaced by a **sign-randomised** null — same support, same Λ weights, μ replaced by random signs — which needs no formula and handles every p. It comes out **close to the old one** (the mass sits at large p, where the one-term case applies), so the correction is only **+5 to +8%**: `R/R_sign` = 0.959, 1.324, 1.442, 1.500, 1.551, **1.584** against the old 0.882 … 1.522. **Every conclusion of #52 and #54 stands**, with the ratios slightly larger. Draw-to-draw sd 0.0001–0.0010. **(b) the labels.** The ladder's rows are *N divisible by the core*, hence **nested supersets** and lower bounds on depth — and the first row is **every even N**, not a depth-0 class. Increments 264–266 wrote them as depth classes |
| 54 | increment 264's depth ladder, whose deepest row had **n = 2** and which left open whether R reaches 1 | repeated at X = 2·10⁷ with **n = 10 per depth**. The two axes separate: **R decays −13.3% on average per factor 12.5 in N at fixed depth**, while the depth steps shrink to **+0.025** by depth 7. Along the primorial sequence a depth step multiplies N by the new prime, so depth 6 → 7 costs **−15.4%** in N against **+3.6%** in depth — **the N-decay outruns the depth growth**, and the measured pair agrees (0.7817 at depth 6, N ≈ 10⁶ → 0.7069 at depth 7, N ≈ 1.9·10⁷). P.4's demand decays even along the hardest sequence. What stays refuted is the **premise**: `R/R_null` is flat in N at fixed depth (1.522 → 1.537), so square-root-per-p is not what carries deep N |
| 53 | *(a pre-registered calibration that failed, recorded because the proxy would otherwise have been used)* the depth ladder restricted to `p ≤ 10⁵` as a cheap stand-in for the full statistic | **it does not earn its use.** The criterion was that the two depth trends agree. Direction and monotonicity do agree, but at depth 6 the proxy reads `R/R_null = 3.073` against the full statistic's **1.522** — an exaggeration of about 2× — and the proxy's own value shrinks with X (R at depth 0 falls 0.176 → 0.065 from X = 1.6·10⁶ to 2·10⁷), so it is not comparable across N either. It may be read for **direction only**; the quantitative question, whether R reaches 1, needs the full statistic |
| 52 | transform P's margin premise: **square-root cancellation per p**, measured as `S_abs/S_null` **flat at 0.87** | flat only **on average over N**. Along the primorial ladder at N ≈ 10⁶ the ratio climbs **0.959, 1.030, 1.234, 1.363, 1.433, 1.488, 1.522** as rad(N) takes 0 to 6 odd primes: at deep N the per-prime sums cancel **worse than random**, because the location mask gives each `D_p` a systematic component. `R = S_abs/Σ log p M_p` climbs 0.391 → 0.782 with `R = 1` being no cancellation at all. The increments shrink, so R may saturate below 1, but the last row has **n = 2** and the decisive measurement — R along the primorial sequence at N ≥ 10⁷ — is not done |
| 51 | Proposition **P.4's margin**, reported as `S_abs/N ≈ 0.30 and falling` | that is an **average over all even N**, and P.4 is a statement about **every** N. Split by depth in the location mask at N ≈ 10⁶: shallow `N = 2q` gives **0.2576**, deep `N = k·30030` gives **0.7773** — three times the reported figure, with the all-N average sitting near the shallow end because deep N are rare. Not a new obstruction (the deep group declines at the same relative rate: 0.8137, 0.7695, 0.7488 across a factor 3 in N), but the margin claimed in TRANSFORM_P §4 is the margin at a typical N and not at the hardest. Also measured and negative: **de-masking does not help** — replacing `Σ log p\|D_p\|` by `Σ log p\|D_p − m_p\|` gives **1.7677** against a raw 0.7773 (permuted floor 1.1114), because the mask lives at `p ≤ 8192` while the demand's mass is at large p |
| 50 | *(an audit that came back clean, recorded because a clean audit is a result)* after #48–#49 it was a live question whether the **negative map's kill-tests** shared the counting fault | **they do not.** `sweep_A` and `sweep_B` count both the mean and the sd-ratio z; the forge kill-tests K2 and K4 report 0 flags with every printed |z| ≤ 1.8, so their two-flags-required rule never had to discard anything; `hyp_round2b` reports max |z| = 1.92 against a null of 0.798 mean |z|. The closures stand. **One caveat**: in K2 the third ALIVE criterion `|ρ₁| ≥ 0.15` is set to `NaN` when fewer than 50 consecutive-k pairs survive, which happens at h = 30 and h = 210, and `abs(nan) >= 0.15` is False — so *could not be evaluated* was scored as *no signal* at two of eight h. The other two criteria did evaluate there and the verdict is not in doubt, but the record should say two of three, not three of three |
| 49 | `sweep_B2`'s verdict, from **3 of 5** flags surviving `÷𝔖` | the same omission as #48 — `sweep_B3` was derived from this helper, so both dropped the three printed mean-z. Counting all eight: **6/8 raw, 4/8 after ÷𝔖**, the added survivor being the `ω≥5` **mean at z = −5.25**. Read together with #48 the pair now says something clean: **√𝔖 is the scale mask** (all three sd ratios → 1.040, 0.998, 1.007), **𝔖 overshoots** (they invert to 1.287, 0.708, 0.820), and **neither touches the location** (mean −6.35 after ÷√𝔖, −5.25 after ÷𝔖). `sweep_A` and `sweep_B` were checked and count both z, so their verdicts are unaffected |
| 48 | `sweep_B3`'s verdict "PARTIAL — 𝔖(N) explains some of the structure", from **1 of 5** flags surviving `÷√𝔖` | the counter saw **5 of the 8 statistics it printed**: `stats()` computed the three mean-z, printed them, and appended only the sd-ratio z to the flag list. Counting all eight gives **6/8 raw and 2/8 after ÷√𝔖**, and the two survivors are `corr(·,𝔖) = −5.34` and the **ω≥5 mean at z = −6.35 — larger than any counted statistic**. Both are location; all three sd ratios collapse to 1.040, 0.998, 1.007. The sharper verdict is **√𝔖 removes the scale structure completely and leaves the location structure untouched** — the fifth independent sighting of the mask in the old sweeps, and the only one where the evidence was printed and simply not counted |
| 47 | *(no prior claim — a new measurement, recorded here because its first fit was wrong)* the ratio `Q = Σ C² / Σ_v μ²(v)Λ(N−v)²` decays like `(log N)^{−0.62}` | **`(log N)^{−0.43}`.** The first fit ran across all eight bands including two with `Q > 1` — the pre-asymptotic regime, where the wall is *less* cancelled than random and no asymptotic model applies. Fitting a power law through a curve still turning over inflated the exponent by 45%. On the `Q < 1` bands the decay still beats a constant (weighted residual **71.4** against **405.3**), so the decay is real, but it is slow: `Q ≈ 0.68` at 10¹² and `≈ 0.27` at 10¹⁰⁰. **No sub-square-root cancellation**: `C(N) ≍ √N·(log N)^{0.29}` |
| 46 | MEASUREMENTS §13: the off-diagonal of the wall's variance identity is **positive and comparable**, "about half the variance of C(N) is genuine shifted-Möbius correlation" | **it is negative.** The reading rested on X = 2000 and 4000. Reproduced there (0.554, 0.498) and then continued: `off/LHS` = **0.077** at 10⁵, **−0.144** at 10⁶, **−0.208** at 4·10⁶, crossing zero near X ≈ 3·10⁵. The shifted-Möbius correlation **removes ~14% of the diagonal** rather than supplying half the total. Separately, the same identity's diagonal measures `diag/(n·N·log N)` = 0.751 → 0.806 → A(N), so §13 **already implied** `Var C ≍ N log N` while §12 asserted `𝔖(N)·N`: the document contained its own refutation for fifty increments (see #36) |
| 45 | MEASUREMENTS §12: sweep_B's five flags are "all what a 𝔖-dependent **scale** predicts", so "the mask is **√𝔖**", with the leftover **mean drift** "already insignificant at N ≈ 9·10⁵" | **four of the five are scale; the fifth is location.** A scale mask cannot correlate C with 𝔖. At full power (n = 1.95·10⁶ against sweep_B's 1500) `corr(C,𝔖) = −0.1386, z = −193.5`, and sweep_B's own fix `C/√𝔖` still reads **z = −137.4**; exact scale removal `C/√V` leaves **−211.8**. Dividing by √V flattens the sd across 𝔖 bands (0.88–0.94) while the mean runs 0.012 → **−1.682** — pure location. **The campaign had the location mask in hand three times at sweep_B** — the 𝔖-correlation, the "mean drift", and the 3\|N mean split of §11 — and read all three as decaying residuals of a scale law. §12 rewritten |
| 44 | the mask should appear as an **observable excess** in the Goldbach count, since Theorem C gives `r̃(N) − 𝔖(N)N ≈ −𝔖(N)C(N)` | **it does not — E₃ absorbs it.** The relation is real (corr `+0.701`, slope `+1.263` against a predicted 1, permuted corr `+0.0001`), but resolved by depth the slope falls 2.15 → 1.67 → 1.40 → 0.91 → 0.43 → **0.18**: exactly where the mask is largest it is least visible. Since `R = P + E₃` by definition, a slope of 0.18 means `cov(E₃,P)/var(P) = −0.82`, i.e. **E₃ cancels 82% of the mask's contribution**. Consistent with `rms E₃/rms 𝔖C` growing 1.02 → 1.36 and `E₃ ∼ N^0.65` against the mask's `N^0.5`. The mask is a real second-order term of C(N) and **not** an observable bias in r̃ |
| 43 | the wall's location mask should have an analogue on the **supply side**, since `q\|k, q∤N` forces `N−mk` coprime to q | **refuted.** No group of k shows a mean shift (all \|z\| < 1.2; prime-by-prime max 1.65 against a permuted control reaching 1.66). The reason is structural: the wall's mask comes from **Λ** forcing v coprime to rad(N), a *primality* constraint, whereas both factors of D(k) are μ, which imposes squarefreeness. The analogue is therefore a **support** mask, not a location shift: `q²\|gcd(k,N)` ⟹ `q²\|N−mk` for every m ⟹ **D(k) = 0 identically**. All 1212 predicted zeros observed. En route, `P(D>0) = 0.346` looked like a z = −12.7 refutation of Conjecture L's Gaussian half; it was the atom at zero, and excluding zeros gives 0.5044, z = +0.46 |
| 42 | increment 248's residual after mask removal (z ≈ ±5) is a **missing piece of the mask** | it is **skewness**, and no mask is missing. Raising the enumeration from q ≤ 13 to q ≤ 37 moves z(3\|N) only 5.51 → 4.18 and leaves z(3∤N) flat at ≈ −5.2, so truncation is refuted. Subtracting a cell **mean** zeroes the mean, not the median, and for mean-zero data `P(X>0) − 1/2 = −γ₁/(6√(2π))` with no fitted constant: measured +0.00317 against +0.00258 for 3\|N (ratio 1.23) and −0.00233 against −0.00103 for 3∤N (2.27). **The enumerated mask is the whole of the deterministic part.** Recorded with it: the Edgeworth term was first written with a **plus**, which inverted every verdict |
| 41 | sweep_B item **B4**: the sign balance of C(N) reads `P(C>0) = 0.4800, z = −1.55`, i.e. **no signal** | a **false negative from underpowering**. B4 ran on ~1500 values, SE 1.3·10⁻², against an effect the mask makes structured rather than uniform. At 1.95·10⁶ values (SE 3.6·10⁻⁴): `3\|N` gives **0.44153, z = −94.3**; `3∤N` gives **0.52703, z = +61.6**; the gap is z = −112.6. By depth, N divisible by five of {3..23} has `P(C>0) = 0.109` — **89% negative**. Removing the enumerated mask returns all of it to within z ≈ 5 of 1/2 |
| 40 | the derived mask `R_A(N) = Π_{q∤N}(q−3)/(q−1)`, from `P(q\|v) = 1/(q−1)` for `q∤N` with independent indicators | **refuted quantitatively, right qualitatively.** Its sharpest prediction — the factor at q = 3 is exactly 0, so `3∤N` kills the mask — fails: the mask does not vanish but **changes sign** (mean R_A = −6.10·10⁻⁴ for 3\|N against **+2.45·10⁻⁴** for 3∤N, ratio 2.49 against a criterion of 5). Per-prime factors measure 12.98, 4.70, 2.91, 2.62, 2.28, 2.16, 2.02 against a predicted 2.00, 1.50, 1.25, 1.20, 1.143, 1.125, 1.10 — decreasing as predicted but toward **2, not 1**, which no per-prime effect can do. One fitted constant gives R² = +0.381 and under-predicts the deep cells 3–4×. The independence of the divisibility indicators, given that N−v is prime, is the identified culprit |
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
