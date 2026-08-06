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

### 2. §2 classification — **fatal to the draft, and worse than stated**

The draft claimed total μ-block mass a ≤ y^{O(1)}, and derived from it
the CO modulus bound q = ak ≤ x^{1/3+O(1/J)}. The Heath–Brown identity
actually allows a up to y^J ≈ x. The true classification has three
leaves: P-I′ (some free block ≥ x^{1/2+δ} — survives), the CO corner
(a genuinely small), and a **type-II region** (medium μ-mass, no long
free block) in which the rough×rough core returns. The draft's
bookkeeping hid that region.

**How much mass it hides — measured** (`code/c3_hb_mass.py`). Since
bookkeeping must bound each piece, the relevant weight is the absolute
one: W(a) = Σ_j C(J,j)·A_j(a)·D_j(M/a), where A_j(a) counts
factorisations of a into j squarefree blocks each ≤ z = M^{1/J} and
D_j(x) = Σ_{b≤x} d_j(b) is the free side. The draft needs
a ≤ y^{O(1)} = M^{o(1)}, so the small-θ columns are the operative
ones:

| J | a > M^{0.05} | a > M^{0.10} | a > M^{1/3} |
|---|---|---|---|
| 3 | **0.939** | 0.819 | 0.326 |
| 4 | **0.949** | 0.826 | 0.257 |
| 6 | **0.960** | 0.836 | 0.173 |
| 8 | **0.947** | 0.750 | 0.032 |

**94–96% of the identity's weight sits at μ-side sizes the draft's
assumption excludes, at every J including the J = 8–12 the draft
proposes, with no decay in J.** The apparent decay in the last column
is an artefact of comparing a fixed cut M^{1/3} against max a = z^J,
which shrinks from 195112 to 65536 as J runs 3 → 8 while M stays
200000; at the thresholds the claim actually requires there is no such
effect. The mass also concentrates in the largest block counts
(j = J−1, J carry the bulk), which is precisely where a can be large.

Raising M moves the figures the wrong way for the draft, which settles
the finite-size question:

| M | J | a > M^{0.05} | a > M^{0.10} |
|---|---|---|---|
| 2·10⁵ | 4 | 0.949 | 0.826 |
| 10⁶ | 4 | **0.961** | **0.865** |
| 2·10⁵ | 8 | 0.947 | 0.750 |
| 10⁶ | 8 | **0.969** | **0.849** |

So the type-II region is not a leak at the edge of the classification;
it is the bulk of it. And the pieces there carry a *rough* coefficient
α(a) — a product of μ-blocks — with no divisor structure, so they have
no Voronoi entry: the spectral door of §3 serves only the corner.

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
   See the next section: this requirement is now closed for one of the
   two classes a transform can belong to, and unpopulated in the
   other.
2. **A classification** covering the type-II region of coordinate 2.
3. **Quantitative averaged Chowla**: shift-averaged binary μμ
   correlations with a *fixed* log-power saving. Best known is
   (log)^{1−c} for a small fixed c.

## Requirement 1 in detail: the transform

**The geometry.** Put A = N − a and B = N − b, so that A = mk and
B = m′k. Then

> **m′A − mB = 0 exactly**,

verified over 200,000 triples with zero mismatches
(`code/c3_pencil_check.py`; the uncentered form m′a − mb = (m′−m)N is
verified alongside it). So in centered coordinates the dilate family is
the **pencil of lines through the origin**, one line per slope m′/m,
traversed as k runs — while the shift family b = a + h is B = A − h, a
family of **parallel lines of slope 1**. Requirement 1 asks for a
transform carrying the first to the second.

A transform can be of exactly two kinds. Both are now accounted for.

### (a) Pointwise change of variable — **closed**

A pencil is characterised by its vertex, which here is a finite point;
a parallel family is a pencil whose vertex lies at infinity. An affine
map sends lines to lines and **finite points to finite points**, so no
affine map carries one family to the other. Since μ lives on ℤ, the
structure-preserving changes of variable available are exactly the
integral affine ones. In particular no shear (A,B) ↦ (A, B − jA) helps:
it permutes slopes and fixes the origin, so the family remains a
pencil.

This settles, by proof rather than by measurement, what the earlier
probe had found empirically: the re-indexing by d = m′ − m produces a
genuine fixed-shift correlation, but with shift a multiple of the
modulus, so rescaling returns the original object — an exact identity,
184.0 versus 184.0, difference 0.000e+00
(`code/c3_transform_probe.py`). Nor can such a transform be motivated
by the slope family being easier: it is not. Measured means
0.4561/0.4714 against the shift family's 0.4826/0.4591, ratios 0.945
and 1.027, both half-normal on their own support.

### (b) Summation formula — **open, and unpopulated**

The remaining possibility is a genuine analytic transform whose dual
variable supplies the shift average. Here the obstruction is the one
coordinate 3 already identified, in general form: a Voronoi-type
formula compresses a sum only against a **smooth** outer weight, and
our outer weight is μ — rough. The two legitimate realisations are
therefore the two that coordinate 3 examined, and both fail for stated
reasons: truncated Voronoi + Abel yields the √-phase but no Kloosterman
factors and a non-summable weight, while the delta method yields
Kloosterman factors but a linear phase whose repair is the refuted
dispersion route. The other named candidates in this class — the
inverse-domain spectrum and manufactured modularity — are closed by
measurement against accurate nulls.

**Status of requirement 1**: closed for changes of variable; open for
summation formulas, with every named candidate in that class already
closed. That is a classification of transform types, not a theorem
about all conceivable mathematics, and it is stated at exactly that
strength.

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
