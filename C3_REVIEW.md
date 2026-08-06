# C-III review verdict — the draft is refuted, the route is not

*Fresh-context refutation-mandate review of C3_DRAFT.md (increment
188), line-checked against the Heath–Brown identity, truncated Voronoi
series, μ-BV and the exponential-sum literature; re-audited against
the corrected E1 target (increments 199–201). Two of the original four
kill coordinates were budget calls measured against a target that was
itself misstated, and they do not survive. The other two do, and they
are fatal to the draft as written.*

## The four coordinates, as they stand

### 1. §5 assembly — **does not hold**

The original objection was that E1 is an L²(k) statement "at Gaussian
(square-root) scale, demanding |D(k)|² ≈ M_k", while every leaf of the
tree outputs linear-scale pointwise log-savings
|D(k)| ≪ (x/k)(log)^{−C}; squaring gives Σ|D|² ≪ (N²/K)(log)^{−2C},
declared "off target by N/K ≥ x^{2/3}".

The premise is false. E1 normalises against Σ M_k² ≍ N²/K, the
**trivial** bound, not against the square-root scale:

> Σ_{k∼K}|D(k)|² ≪ (log N)^{−2A−2} Σ_{k∼K} M_k².

The tree's own output, as the objection itself computes it, is
Σ|D|² ≪ (N²/K)(log)^{−2C} — the same form, meeting the target whenever
C ≥ A+1. The distance that was measured was the distance to a scale
the chain does not require.

### 2. §2 classification — **fatal to the draft**

The draft claimed total μ-block mass a ≤ y^{O(1)}; the Heath–Brown
identity actually allows a up to y^J ≈ x. The true classification has
three leaves: P-I′ (some 1-block ≥ x^{1/2+δ} — survives), the CO
corner (a genuinely small — near measure-zero), and an unclassified
**type-II region** (medium μ-mass, no long smooth block) where most
generic mass lives, in which the rough×rough core returns. The draft's
bookkeeping hid that region.

### 3. §3 dual bookkeeping — **fatal to the draft**

CO is a one-free-variable sum, and "Voronoi on the divisor side with w
outside" is not a defined operation. The two legitimate realisations
are: (A) truncated Voronoi + Abel — yields the √-phase but with no
Kloosterman factors, no c-range, and a non-summable d(u′)u′^{−3/4}
weight with forced bandwidth U ≥ N_n(log)^{2A}; and (B) the delta
method — yields Kloosterman factors and c-ranges but a linear w-phase,
whose repair is the already-refuted dispersion route. The draft's
hybrid — Kloosterman factors, √-phase, 1/c weights and a c-range
simultaneously — is the output of neither.

### 4. §4 budget — **arithmetic does not hold; one clause survives**

The objection opened with "**no-margin** arithmetic: pointwise
consumption loses by N_n^{1/2} ≥ x^{1/3−O(1/J)}". The no-margin
premise is the same misstatement as in coordinate 1. Against the
actual margin N/K ≥ N^{2/3}, a loss of x^{1/3} is affordable with
N^{1/3} to spare.

Its second clause also softens. "The object is typically √-sized, so
absolute-value summation reaches at best trivial scale" is not a death
sentence once the target sits exactly one log-power *below* trivial.
What survives is the last clause, and it is the real content: all
saving must then come from u′-family sign cancellation, which expands
into **shift-averaged binary μμ correlations at log-power strength**.
That is a named open problem, not an impossibility. (T-S0b's sampling
also missed the mass-dominant conductor-saturated block, F ≍ length,
wavelength O(1), by orders of magnitude — a measurement caveat that
stands.)

## Verdict

**The draft is refuted** by coordinates 2 and 3, which are structural:
an incomplete classification and an object that is the output of no
legitimate transform. **The route is not refuted.** Nothing above
shows that C-III cannot work; what it shows is that this construction
does not.

## What C-III needs, exactly

1. **A legitimate transform** in place of the draft's hybrid object.
   Its geometry is now pinned down: with a = N−mk and b = N−m′k one
   has m′a − mb = (m′−m)N, so the dilate problem is μμ correlations
   along lines of arbitrary rational slope and the shift problem is
   the slope-1 case. Two things are measured about it
   (`code/c3_transform_probe.py`): the obvious re-indexing by
   d = m′−m produces a genuine fixed-shift correlation but with shift
   a multiple of the modulus, so rescaling returns the original object
   — verified as an exact identity, 184.0 versus 184.0, difference
   0.000e+00 — and the slope family is **not** tamer than the shift
   family (means 0.4561/0.4714 against 0.4826/0.4591, ratios 0.945 and
   1.027, both half-normal on their own support). So the transform
   cannot be a re-indexing, and cannot be motivated by the slope
   family being easier. It must be a genuine analytic transform whose
   dual variable supplies the shift average — and the named candidates
   for that (inverse-domain spectrum, manufactured modularity, the
   dispersion/delta route) are each already closed by measurement or
   by structure.
2. **A classification** covering the type-II region of coordinate 2.
3. **Quantitative averaged Chowla**: shift-averaged binary μμ
   correlations with a *fixed* log-power saving. Best known is
   (log)^{1−c} for a small fixed c.

Requirement 3 is worth stating on its own, because it changes the
character of the obstruction. The adjudication's central finding was
that the **dilate** average admits no diagonalising character family;
the **shift** average is precisely the home ground of the
Matomäki–Radziwiłł–Tao machinery. If a legitimate transform exists,
what remains is quantitative strength inside an active area rather
than the absence of any coupling surface.

## What survives independently

- **P-I′ (corrected)**: pieces with a smooth block ≥ x^{1/2+δ}, closed
  by a k-multiple-moduli variant of Möbius–BV — classical, consistent
  with the surviving m-face.
- The one-sided-opening structural observation, and the T-S0
  measurement discipline (including the re-reading of a HEALTHY signal
  as a budget diagnostic — a method upgrade).
- **A map coordinate**: the conductor-saturated Möbius-vs-√-phase sum
  (F ≍ X) is a genuinely open single-variable problem outside current
  Möbius-orthogonality technology (Green–Tao nil-phases, MRT
  dynamics) — though even its solution would not close CO without
  binary μμ input.

## Terminal sentence

> Beyond the spectral door waits the same wall — the bilinear μμ
> pairing. The door is not proved shut; what is proved is that this
> draft does not open it, and that the wall behind it is the
> shift-averaged binary correlation at log-power strength.

*C3_DRAFT.md is retained with a refutation banner; this file is the
verdict of record. Supersessions and the corrections that produced
them are recorded in one place: CLOSURE_REAUDIT.md.*
