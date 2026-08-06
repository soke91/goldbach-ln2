> **⚠ ADVERSARIAL REVIEW VERDICT (increment 188): this tree is
> REFUTED at four independent coordinates — assembly scale mismatch
> (x^{2/3}), a Heath-Brown bookkeeping error hiding the type-II leaf,
> an illegitimate Voronoi dual, and a pointwise budget deficit
> (x^{1/3}) sealed by our own T-S0b measurement. Two refutations
> stand even granting Lemma S for free. See C3_REVIEW.md. Survives:
> P-I' (corrected, classical), the one-sided-opening observation, and
> the new map coordinate "beyond the spectral door waits the same
> wall." Retained as a record of the attempt (correction #27).**

# C-III Draft — the spectral-entry reduction tree (increment 186)

*A pure-theory attempt at the only surviving representation class.
Status: DRAFT — no step below is claimed proved unless marked
[classical]; every gate is an open task; the whole document goes to
fresh-context adversarial review (the increment-143 protocol) before
any status upgrade. The author-model's prior gate arithmetic was
refuted once (REVIEW_VERDICT.md); this draft therefore claims no
margins anywhere.*

## 0. Target

E1 (weak): Σ_{k∼K}|D(k)|² ≪ (log N)^{−2A−2} Σ_k M_k² for
D(k) = Σ_{√N<m≤N/k} μ(m)μ(N−mk), K ≤ x^{1/3} dyadic — a fixed
log-power saving over the **trivial** bound Σ M_k² ≍ N²/K.

> This draft was written against the square-root-scale form
> (Σ_k M_k), which is a factor N/K stronger than the chain consumes.
> That misstatement is what invalidated two of the four coordinates in
> the review of this draft; the surviving two are structural. See
> C3_REVIEW.md and CLOSURE_REAUDIT.md.

## 1. One-sided Heath–Brown opening [classical]

Apply the Heath–Brown identity for μ to the factor μ(m) ONLY, with J
blocks at truncation y = x^{1/J} (J to be chosen, e.g. 8–12). Every
resulting piece of D(k) has the form

> (intact rough factor μ(N−mk)) × (product of ≤ J short μ-blocks
> ≤ y) × (product of ≤ J−1 unrestricted smooth 1-blocks),

with m = (blocks product). Key structural point: **the second μ is
never decomposed** — roughness is concentrated in one factor, and the
decomposed side is smooth-heavy (each piece has either a long
1-variable or a multi-block divisor structure).

## 2. Piece classification

- **P-I (one long 1-variable)**: inner sum = linear μ-sum over an AP
  with modulus (short-blocks)·k ≤ y^{O(1)}K. For J moderately large
  this modulus is ≤ x^{1/2−δ}: Möbius–Bombieri–Vinogradov territory
  [classical]. This is the m-face that survived the refutation.
- **P-II (two or more medium 1-blocks, none individually long)**: the
  smooth variables assemble into a divisor-type weight. Collapsing
  them: the piece is a linear combination (short-μ-block coefficients
  α(a), |α| ≤ d_J) of the **Central Object**

  > **CO(a, k): Σ_{w ≡ N (mod ak)} μ(w) · d_I((N−w)/(ak))**,

  the μ-twisted divisor-in-progression sum, modulus q = ak with
  a ≤ y^{O(1)}, hence q ≤ x^{1/3 + O(1/J)}. [reduction bookkeeping:
  TO VERIFY — ranges, block counts, and the exact d_I window]

Note what did NOT appear: no rough×rough core arises from a one-sided
opening — the tree has exactly two leaf classes, P-I (closed,
classical) and CO.

## 3. The spectral entry: Voronoi on the divisor weight

d_I in arithmetic progressions admits Voronoi summation: the weight
d_I((N−w)/q) for w ≡ N (mod q) expands into a main term (computable;
absorbed into the mask/main-term ledger) plus dual terms carrying
Kloosterman-type arithmetic factors and **Bessel oscillations in w**:
phases of the shape e(±2√(u′(N−w))/c) with dual parameters (u′, c),
c related to q. [classical for smooth test functions; the transfer to
our rough μ-weight is by exchanging the roles — Voronoi is applied to
the d-side, and w remains the outer variable — TO VERIFY: the precise
dual-term normalization and the c-range]

After this exchange, CO becomes a sum over dual parameters of

> **Σ_w μ(w) · e(±2√(u′(N−w))/c) · (smooth weights)** —

one-variable Möbius sums against square-root (Bessel) phases.

## 4. The terminal object — Missing Lemma S

> **Lemma S (needed, not known).** Uniformly over the Voronoi-dual
> ranges (u′, c) arising in §3, and summably against their 1/c-type
> weights, Σ_{w∼W} μ(w) e(2√(u′(N−w))/c) ≪ W (log W)^{−A}.

Remarks, honest on both sides:

- **Why this is a genuinely new coordinate.** Lemma S is a
  ONE-variable μ-sum — the binary pairing is gone. Sums
  Σμ(n)e(αn^θ) (fixed θ < 1) are an existing analytic literature
  (Vaughan decomposition + van der Corput estimates give savings in
  various ranges), and the parity obstruction in its combinatorial
  binary form does not obviously apply to a single-μ analytic sum.
  None of the fifteen closed designs died here; this coordinate was
  not reachable before the one-sided opening + Voronoi exchange.
- **Where the wall may re-enter (Checkpoint X — the adversarial
  review's first target).** The dual phase couples w and u′ through
  √(u′(N−w)) — a bilinear coupling in disguise. If the required
  uniformity in u′ forces the estimate to hold for the FAMILY (a
  u′-averaged square, say), the hyperbolic pairing may reconstitute
  exactly the wall this program has mapped, now through the spectral
  door. Whether Lemma S is consumed pointwise-in-u′ (each dual term
  separately, summed trivially against 1/c-weights — in which case
  the existing μ-vs-nonlinear-phase technology is the right
  comparison) or in family (in which case likely circular) is a
  bookkeeping fact of §3 that must be settled FIRST.
- **Budget question (Checkpoint Y).** Even pointwise: (number of
  dual terms) × (per-term bound) must beat the trivial bound of CO.
  Voronoi trades length for dual length; the c- and u′-ranges and
  the 1/c weights decide whether a (log)^{−A} per-term saving
  survives the dual mass. No margin is claimed; this is an open
  computation with the program's full gate-arithmetic discipline
  (and its history of one refuted gate computation) applied.
- **Consistency check available now**: the object of Lemma S is
  measurable. Before any proof attempt, the program can stamp
  whether Σμ(w)e(2√(u′(N−w))/c) actually exhibits square-root
  cancellation uniformly over sampled (u′, c) — if even nature
  refuses, Lemma S is dead on arrival and the tree closes here.

## 5. Assembly shape (conditional, not claimed)

If Lemma S (pointwise version) holds with summable uniformity, then
CO ≪ q-average-compatible (log)^{−A} losses, P-II closes, and with
P-I classical the k-average square in E1 follows by Cauchy–Schwarz
over the piece decomposition — with A degraded by the O(log^{O(J)})
piece count. [assembly bookkeeping: TO VERIFY — nothing here is
claimed beyond the shape]

## 6. Task list (in order)

1. **T-S0 (measurement, cheap)**: stamp Lemma S's object numerically
   — square-root cancellation of μ against Bessel phases over
   sampled (u′, c); kill-test before theory.
2. **T-S1 (bookkeeping)**: §3 dual-term normalization — settle
   Checkpoint X (pointwise vs family consumption).
3. **T-S2 (budget)**: Checkpoint Y dual-mass arithmetic, no-margin
   discipline.
4. **T-S3**: literature adjudication of μ-vs-n^{1/2}-phase sums
   (what is actually proved, in which uniformity).
5. **T-S4**: fresh-context adversarial review of the whole tree.

## T-S0 record

- (increment 186) First run INVALID by self-audit: 9/240 pairs
  survived a miscalibrated oscillation filter — verdict suspended.
- (increment 187) **T-S0b, corrected sampling (240/240 pairs)**:
  mean r = 0.871 ∈ [0.6, 1.0], max 2.026 < 4; all three oscillation
  bands healthy (0.827 / 0.886 / 0.913 — a consistent 9–14% above
  half-normal, no runaway band, no heavy tail). **HEALTHY** — nature
  supplies square-root cancellation to Lemma S's object uniformly
  over the sampled dual ranges. Theory tasks proceed; the mild
  above-half-normal constant is noted for the budget arithmetic
  (T-S2), not assumed away.

*Opened at increment 186. The draft's value even in failure: §2's
one-sided tree and §4's Lemma S give the spectral-side failure
coordinate this program's map lacked.*
