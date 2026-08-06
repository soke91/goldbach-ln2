# Technique Forge — engineering the missing invariance (opened at increment 171)

*The adjudication (AMPLITUDE_ADJUDICATION.md) reduced "no technique
exists" to a precise design specification. This document is the forge:
candidate designs, fast kill-tests against the measurement corpus, and
the escalation rule. Discipline: every candidate gets a numeric
kill-test BEFORE any proof prose; a candidate that fails its kill-test
is recorded and closed (the cascade lesson: one-step R² = 0.68 looked
alive, two-step R²₂ = 0.19 killed it).*

## The design target (from the adjudication, verbatim spec)

Construct a transform/kernel/invariance T acting on the dilate family
{D(k)} such that the pair constraint of T-averaged correlations
becomes LINEAR in some character family (additive or multiplicative),
while T costs at most log-power in the E1 budget. Equivalently: what
translation is to the h-average and dilation-invariance of the log
measure is to Tao's average, T must be to the k-average.

## Assets unique to this program

- **A1. The exact dilation ladder** (identity error 0, verified):
  the p-divisible sub-sum of D(k) maps exactly onto −D^{(p∤)}(pk).
  The only exact invariance the family is known to possess.
- **A2. The determinant observation**: the pair constraint
  mu′ − m′u = N(m − m′) is a 2×2 determinant equation — the natural
  home of GL₂ spectral theory (Kloosterman/DI). Known blockage: our
  parameters violate the premises of the existing dispersion lemmas
  (REVIEW_VERDICT items 1–2); the forge question is whether a NEW
  weight/congruence structure can be manufactured, not borrowed.
- **A3. Conjecture L**: the target field is measured
  mask × unit-Gaussian at every level — any candidate technique may
  assume-and-verify against this exact statistical profile.

## Candidate designs and kill-tests

| # | Design | Kill-test (fast, numeric) | Status |
|---|---|---|---|
| K1 | **Multiplicative Fejér kernel on the ladder orbit**: represent D(k) as an optimally weighted combination of its orbit {D(sk) : s squarefree, s \| 30030} via A1; if the representation is near-exact, orbit-overlap across k linearizes the k-average | Full-orbit least-squares R² of D(k) on {D(sk)}: R² → 1 = alive; R² saturates < 0.9 = dead (cascade already suggests death; this is the definitive orbit-level version) | **DEAD** — full 63-divisor orbit R² = 0.466 (one-layer 0.442); residual energy = 0.499 of the unit-Gaussian budget. Half of D(k)'s energy is invisible to its entire multiplicative orbit: the exact ladder transfers sub-sums, but the p-coprime core it leaves behind is itself half the field at every orbit depth. Consistent with (and stronger than) the cascade death. Closed (`code/e1_forge_kt1.py`) |
| K2 | **Determinant/Kloosterman route with manufactured congruence**: impose an artificial pair congruence k ≡ k′ (mod d) by splitting the k-average into progressions — recovers the conductor-collapse mechanism at cost d; question is whether some d-average wins a log-power | Measure the conductor-collapsed correlation strength as a function of d: does Σ_{k≡k′ (d)} C_{k,k′} gain more than the d-splitting cost? | **DEAD** — 10 h-values (1…210), 0 pre-registered flags: means ≤ 1.8z, variance deviations ≤ 1.6σ, autocorrelations ≤ 0.128 mixed-sign. Congruent pairs are statistically h-blind, exactly as pure Conjecture L predicts; the manufactured congruence buys nothing against the factor-d cost. Closed (`code/e1_forge_kt2.py`) |
| K3 | **Wishart/operator route**: E1 needs only the band ratio; the C-matrix is measured Wishart-clean. Design: prove the ratio bound from a moment-method bound on the Gram spectrum of the row family {μ(N−pk)}_p — rows are deterministic, but row-inner-products are the SAME correlations (circularity check needed at 4th-moment level) | Compute trace-moment ratios tr((RRᵀ)^j)/Wishart-prediction for j = 2,3,4: sub-Wishart 4-cycle cancellation = alive; match = dead (circular) | **DEAD** — tr(M²)/tr(M³)/tr(M⁴) all dead-center on the Wishart null (z = −0.48 / +0.42 / −0.03). No sub-Wishart cancellation to exploit; the moment route consumes 2j-fold μ-correlations with nothing gained. Side product: Conjecture L re-confirmed at 4-cycle moment level. Closed (`code/e1_forge_kt3.py`) |
| K4 | **N-average descent**: the N-averaged theorem is provable (adjudication 5-(2)); design a descent that trades the N-average for a k-average within a single N via the ladder's N-independence | Does the k-averaged dual field T(m, δ) retain, within one N, the shift-structure the N-average exploits (coherence / m-autocorrelation)? | **DEAD** — 0/8 flags (coherence z ≤ 1.22, autocorr ≤ 0.051, m2 ≈ 1 across four δ). The N-average is load-bearing; the k-average supplies no substitute decorrelation. Closed (`code/e1_forge_kt4.py`) |

## Escalation rule

A design survives its kill-test → write the reduction in full,
adversarially review it fresh-context (the increment-143 protocol),
and only then attach it to the chain documents. A design fails → one
paragraph in this table, closed, never re-litigated.

## Ledger

- (increment 171) K1 kill-test launched.
- (increment 172) **K1 DEAD** by pre-registered threshold: full-orbit
  R² = 0.466 < 0.9. The ladder is an exact identity but not a
  spanning one — the multiplicative orbit sees only half the field's
  energy. Design lesson: any viable T must couple to the p-coprime
  core directly, not through divisibility sub-sums. Next candidate on
  resume: K2 (manufactured congruence / determinant route).
- (increment 174) **K2 DEAD**, 0/10 pre-registered flags: the
  congruence field is h-blind. Design lesson: the f-substitution's
  power in the source machine comes from the congruence arising
  *inside* an average over moduli (so the collapse is traded against
  an existing sum), not from the congruence per se; imposing it
  externally only pays the splitting cost. Next: K3 (Wishart moment
  method — pre-register: dead unless 4-cycle trace moments show
  SUB-Wishart cancellation, i.e., exploitable structure beyond pair
  independence; a mere match confirms L but leaves the moment route
  circular, since tr(M^{2j}) consumes 2j-fold μ-correlations).
- (increment 175) **K3 DEAD**: all trace moments dead-center on the
  Wishart null. No sub-Wishart cancellation; the moment route is
  circular. (Side product: L confirmed at 4-cycle level.)
- (increment 176) **K4 DEAD**: the dual field is δ-blind and
  m-uncorrelated; the N-average's linearization has no within-N
  substitute.
- (increment 194) Round 3 opened on the increment-193 asset. New exact
  identity recorded and verified. **R4 DEAD** on its pre-registered
  criterion; a same-direction deficit at 2048 k was settled at 8000 k
  (increment 194b) as block-count SE, with lag-1 autocorrelation
  +0.011 ± 0.011. Fifth constraint on T recorded.

## Round-1 synthesis (increment 177)

Four designs, four clean deaths, four constraints on any viable T:

1. **(from K1)** T cannot act through divisibility sub-sums — the
   multiplicative orbit spans only half the field's energy; T must
   couple to the p-coprime core directly.
2. **(from K2)** T cannot manufacture its congruence externally —
   collapse structure pays only when it arises inside an intrinsic
   average that is already present.
3. **(from K3)** T cannot bootstrap from the field's own moments —
   every moment consumes higher μ-correlations and the field carries
   no sub-Wishart surplus to fund the exchange.
4. **(from K4)** T cannot descend from the N-average — its
   linearization is load-bearing and the k-average holds no shadow of
   it.

**Round-1 conclusion.** The four internal resources of the dilate
family — its multiplicative orbit, its congruence classes, its Gram
spectrum, and its N-family embedding — are each *exactly as flat as
Conjecture L predicts*: the field's internal statistics contain no
lever. Any linearizing invariance T must therefore import
cancellation from OUTSIDE the μμ-correlation universe (zeros of
L-functions consumed in a new way, automorphic spectra attached to
the determinant structure by a construction not yet in the
literature, or an algebraic identity of a kind the ladder is not).
This is the forge's honest round-1 boundary: the technique to be
created is not an assembly of the family's own parts. Round-2
designs must begin from an external cancellation source.

## Round 2 — external cancellation sources (opened at increment 178)

| # | Design | Kill-test | Status |
|---|---|---|---|
| R1 | **Zero-spectrum visibility**: the explicit formula says μ(m) carries the oscillations m^{iγ} (γ = ζ-zero ordinates). If the field D(k) has a component aligned with the zero-templates T_γ(k) = Σ_m m^{iγ}μ(N−mk), the explicit formula gives an external handle on that component | Project the D-field onto span{Re T_γ, Im T_γ} (first 30 zeros) and compare captured energy against random-frequency templates at matched dimension. Pre-registered: ALIVE iff zero-capture ≥ 2× random-capture; else dead | **DEAD** — R²_zeros = 0.2152 vs random 0.2196 ± 0.0055 (ratio 0.98; both at the 60/300 overfitting baseline). The zero oscillations are entirely invisible through the μμ-pairing: the second factor μ(N−mk) scrambles the m^{iγ} coherence completely. The direct explicit-formula channel closes. Closed (`code/e1_forge_r1.py`) |
| R2 | **Determinant/Kloosterman phase**: the pair constraint is the determinant equation; test whether C_{k,k′} correlates with the determinant phases e(N·k̄′/k) (the object DI/Kuznetsov machinery would control if a construction attached it) | Coherent twist gain + per-pair regression vs matched random-phase controls | **DEAD** — the regression test (strictly more general than the single-phase gain) read dead-zero in two independent runs (R² = −0.0001 / +0.0004 vs controls ± 0.0002, 48,000 coprime pairs total): the pair field is completely phase-blind to its own determinant structure. (The per-k gain aggregation hit an uninvestigated nan twice; verdict rests on the regression, which upper-bounds any single-phase gain at noise level.) Closed (`code/e1_forge_r2.py`, `e1_forge_r2b.py`) |
| R3 | **Character transform of the k-average**: Σ_k χ(k)D(k) diagonalizes the inner sum into character-twisted μ-sums over APs mod m | **DEAD BY ANALYSIS** (recorded without a run): for the chain's range K ≤ x^{1/3} the AP moduli are m ~ x/K ≥ x^{2/3} — the transform relocates the difficulty into thin-progression μ-sums (the seam regime, provably void of current technology per the resistance profile); on the mirror side K ≥ x^{2/3} where moduli would be BV-small, the chain does not consume the estimate. Parseval also forbids a statistical gain (Σ_χ \|D̂(χ)\|² is basis-invariant). Closed |

## Round-2 synthesis (increment 180)

The three classical external cancellation sources close in three
distinct ways:

- **Zeros (R1): invisible.** The m^{iγ} coherence that the explicit
  formula provides for μ(m) alone is completely scrambled by the
  pairing with μ(N−mk) — capture indistinguishable from random
  frequencies to the third decimal.
- **Characters (R3): relocating.** The transform that diagonalizes
  the k-average exists, but it deposits every diagonalized component
  into thin-progression μ-sums (moduli ≥ x^{2/3}) — the seam regime,
  outside all current technology.
- **Determinant spectra (R2): phase-blind.** The field carries no
  coherent component along the Kloosterman phases its own determinant
  structure distinguishes.

**Combined with Round 1** (no internal lever: orbit, congruence,
moments, N-descent all exactly flat), the forge had by increment 180
measured, in every direction that has an existing mathematical name,
that the dilate field offers **no coupling surface among the sources
these designs examined**: nothing internal to exploit, nothing external
that can see it. The scope is the forge's own designs and the five
adjudicated routes; the C-III route is separate and is open (see
C3_REVIEW.md).

## Round 3 — the divisor-switch asset (opened at increment 194)

*Round 2 closed with "the technique to be created is a spectral
representation of μ-pairs themselves". Round 3 exists because
increment 193 produced a genuinely new asset that was not available
when rounds 1–2 were designed: the divisor switch of Theorem A, the
one mechanism this program has found that converts a bilinear μμ
object into a linear one.*

**The new asset, stated exactly.** Applied to the dilate field over
its full ranges, the switch gives an exact identity:

> Σ_{k≥1} Σ_{m: mk ≤ N−1} μ(m)μ(N−mk)
>  = Σ_{u<N} μ(N−u) Σ_{m|u} μ(m)
>  = Σ_{u<N} μ(N−u)·[u=1] = **μ(N−1)**.

Perfect cancellation: O(1) for a double sum of ~N log N terms, far
beyond square root. E1 imposes exactly two restrictions that make the
inner divisor sum incomplete — the type-II cut m > √N and the dyadic
band k ∼ K. R4 asks whether any of the perfect cancellation survives
them.

**The mirror observation (derivation, no measurement needed).** Switch
the banded L¹ sum: Σ_{k∼K} D_full(k) = Σ_{u<N} μ(N−u) Σ_{m|u,\ u/2K<m≤u/K} μ(m).
Here m ≍ N/K ≥ N^{2/3} (since the chain needs K ≤ N^{1/3}), so the
surviving Möbius sits on the **long** variable while the complementary
variable k ∼ K ≤ N^{1/3} is short. Theorem A's assignment is exactly
the opposite: there μ landed on the short variable m ≤ N^{1/2−δ} and
Bombieri–Vinogradov closed it. **Same switch, opposite
variable-length assignment** — this is the structural reason the
demand side was free and the supply side is not, and it is a
constraint on any future use of the switch, not an obstruction that
was previously named.

| # | Design | Kill-test | Status |
|---|---|---|---|
| R4 | **L¹ coherence of the restricted field**: does the exact identity's perfect cancellation leave a residue after the two restrictions? If block sums of D(k) over consecutive k are sub-Gaussian, the L¹ identity funds part of the L² bound and the switch becomes a lever | Block sums S_B(j)=Σ_{k∈block}D(k) at B ∈ {1,8,64,512}; statistic ratio(B)=Σ_j\|S_B(j)\|²/Σ_k supp(k), which is B-independent under Conjecture L. Control: same statistic for the full-m field D_full. **Pre-registered: ALIVE iff ratio(B) ≤ 0.5·ratio(1) for ≥2 block sizes B ≥ 8 at BOTH N; DEAD otherwise** | **DEAD** — 0 hits at either N. The identity itself is confirmed exactly (brute force at N = 5000, 20000: the double sum equals μ(N−1) on the nose), but none of its perfect cancellation survives the restrictions. First run (2048 k) showed a same-direction deficit at both N (B=8: 0.87/0.94, B=64: 0.77/0.80) inside the estimator's own SE; settled at full power (8000 k = the entire band where the type-II field is non-empty, since m > √N forces k < √N): B=8 ratios 0.958/1.023, and the sharper diagnostic — lag-1 autocorrelation of D(k)/√supp(k) — reads **+0.0104 / +0.0127 against SE 0.0112**, i.e. dead zero and, if anything, mildly positive rather than the negative coherence a surviving residue requires. The deficit was block-count SE. Closed (`code/e1_forge_r4.py`, `e1_forge_r4b.py`) |

| R5 | **The circle method, applied directly to the wall** C(N) = ∫S_Λ S_μ e(−Nα)dα. Not a lever on the dilate field but the other classical mechanism for the scalar itself; worth pricing exactly now that the switch is closed over its whole design space | Compute the two bills the method must pay — Cauchy–Schwarz ‖S_Λ‖₂‖S_μ‖₂, and pointwise × L¹, sup\|S_μ\|·‖S_Λ‖₁ — against the trivial bound ψ(N) ~ N. **Pre-registered: ROOM iff the margin N/(sup\|S_μ\|·‖S_Λ‖₁) grows at least like a power of log N; NO ROOM if bounded or decaying** (`code/circle_margin.py`) | **DEAD — zero margin, and the cap is identity-level.** (i) Cauchy–Schwarz gives ~(6/π²)^{1/2}N(log N)^{1/2}, *above* the trivial bound by a growing factor (measured 2.30 → 2.80 over N = 2¹⁴…2²⁰). (ii) Any pointwise route is capped by **Parseval**: sup_α\|S_μ\| ≥ ‖S_μ‖₂ = (6/π²)^{1/2}N^{1/2} for free, so sup·‖S_Λ‖₁ ≫ N — trivial again. Measured margin 0.168/0.175/0.158/0.152, decaying. Davenport's uniform S_μ ≪_A N(log N)^{−A} is useless here: it saves against N while the pairing needs scale N^{1/2}, which Parseval forbids. Computation validated by ‖S_μ‖₂/√N = 0.7797 = √(6/π²) exactly. Closed |

## Round-3 synthesis (increment 194, extended at 196)

The program's newest asset closes too, and it closes informatively.

- **The identity is real and is perfect.** Σ_{k}Σ_{m: mk≤N−1}
  μ(m)μ(N−mk) = μ(N−1), verified exactly. On the *unrestricted*
  ranges the bilinear μμ object is not merely square-root cancelling
  but O(1) — the divisor switch is the strongest cancellation
  mechanism this program has found anywhere.
- **None of it localizes.** After E1's two restrictions the field is
  L¹-incoherent to measurement precision: block ratios flat, lag-1
  autocorrelation +0.011 ± 0.011. The perfect cancellation lives
  entirely in the interaction between the type-I window (m ≤ √N) and
  the small-k range, i.e. outside the window the chain consumes.
- **The mirror**, derived above: the same switch puts μ on the short
  variable on the demand side (m ≤ N^{1/2−δ}, hence BV, hence
  Theorem A) and on the **long** variable on the supply side
  (m ≍ N/K ≥ N^{2/3}). The switch is not a technique that happens to
  fail here; it is a technique whose single requirement — Möbius on
  the short variable — is exactly what the type-II cut forbids.

**Fifth constraint on any viable T** (joining the four of round 1):

5. **(from R4)** T cannot be the divisor switch, or any relative of
   it. The switch's cancellation is a property of complete divisor
   sums over the full range; the type-II cut makes every divisor sum
   incomplete and, measured directly, leaves no coherence behind.

**Addendum (increment 196): the two classical mechanisms are now
priced, and both prices are identity-level.**

- **Divisor switching**: closed over its entire weight space by
  Theorem D, and by Theorem D′ the closure survives assuming
  Elliott–Halberstam for Λ at any level θ_E < 1 — the gap is
  N^{1−θ_E}, so closing it needs θ_E = 1, i.e. equidistribution to
  moduli of size N, which is vacuous.
- **The circle method**: zero margin (R5). The cap is Parseval,
  sup ≥ L², which no technique can move.

Both mechanisms fail for reasons that are *identities*, not estimates:
μ ∗ log = Λ on the demand side, the incompleteness of truncated
divisor sums for the switch, and sup ≥ ‖·‖₂ for the circle method.
That is the sharpest form in which this campaign can state the wall:
**every closure it found is at the level of identities, and the
remaining object, C(N) = o(N), is equivalent to the Goldbach asymptotic
itself.**

**The forge's terminal statement** (unchanged by round 3, and now
tested against the program's own newest mechanism as well as against
borrowed ones). In the automorphic world, shifted
convolutions Σa(n)a(n+h) are controlled because the coefficients
a(n) COME WITH a spectral representation (Kuznetsov applies to them
natively); μ has no such representation, and every route above
failed exactly at the point where it tried to borrow one. **The
technique to be created is a spectral (or equivalent
structural) representation of μ-pairs themselves** — not a
projection onto existing spectra. This is where measurement
methodology reaches its limit: kill-tests can close designs, but the
construction of a new representation is a purely theoretical act.
The forge's deliverable is this boundary, drawn with numbers.
