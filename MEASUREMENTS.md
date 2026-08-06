# Measurements — the empirical core of this repository

*This document contains only reproducible measurements and exact,
machine-verified identities. No proof-program claims. The repository's
proof-program attempt (PROOF_SKETCH_E1.md) was refuted in its core
reductions by an independent adversarial review — see
REVIEW_VERDICT.md. Nothing in this document depends on it; the review
explicitly affirmed the measurement corpus.*

## 1. The final-axiom landscape (likely novel numerics)

By Huang–Li (arXiv:2005.03811), Goldbach for large even N follows from
Bombieri–Vinogradov plus the Möbius-twisted hypothesis EH_μ at any
level past x^{1/2}, fixed residue, whose object is c(n) = Λ(n)μ(N−n).
We measured that object's landscape directly:

- **Fixed-residue discrepancy of c(n)** vs random-walk benchmark:
  mean 0.6–1.1 through θ = 0.30 → 0.70 (moduli as x^θ), five values of
  N at two scales (10⁸, 10⁹) — **no visible change at the √x barrier**
  (`code/ehmu_final.py`, `code/ehmu_final_1e9.py`).
- **The innermost sum** T(N) = Σ Λ(n)μ(N−n) (whose o(N) is itself
  unproven): textbook half-normal over 220+120 values of N at two
  scales (means 0.751/0.746) (`code/core_sum_scan.py`).
- **Smooth (well-factorable-type) moduli**: equally healthy (0.64–0.74,
  479 moduli) (`code/smooth_moduli_probe.py`).
- **Plain Möbius in APs**: (0.27·ln Q)×random-walk to θ = 1/2; flat
  3.2–3.8×√(x/q) deep into θ ≤ 0.88 (`code/ehmu_probe2.py`,
  `code/ehmu_beyond.py`).

## 2. The P-profile: sifted primes on the Buchstab curve

R_P(s) := (sifted-prime count)/(independence model) equals
**e^γ·ω(u) to ±2·10⁻⁴** through the sieve-blind zone s ∈ [1.5, 2.5],
with the fine structure (the ω-oscillation, including a −0.07% dip)
reproducing to four decimals across a decade and the residual scaling
as 1/log z (`code/p_profile.py`). Edge deviations at u ≤ 3 match known
Buchstab/HL corrections. Shifted primes sift exactly like random
integers, at every measurable point.

## 3. Exact structure law at roughness x^{1/8}

s(n) = s₀(x) · ∏_{q|n, q>y} (q−1)/(q−2) to 3–4 decimals in controlled
per-q experiments at three scales (y = 19/23/31); sieved factors cancel
exactly (ratio 1.000); fully smooth extremes show no downward
deviation (`code/s_th18_qlaw.py`).

## 4. Exact identities (machine error 0)

- The dilation ladder A_p(k) = −D^{(p∤)}(pk): error exactly 0 on all
  sampled pairs (`code/entropy_ladder.py`).
- Residue decompositions: error exactly 0 (1200+ pairs, three N-values,
  moduli to 2×10⁷) (`code/g1_completion.py` and extensions).

## 5. Thin-progression Möbius (the hardest regime)

At moduli L ~ 10⁵ over range 10⁸ (L > √range — where unconditional,
Selberg-5/8, and GRH-√y technologies are all void): 10,000 classes
give |S|/√n mean **0.738 — below half-normal (0.798)** with the tail
exactly Gaussian (22 exceedances of 3σ vs 27 expected)
(`code/seam_thin.py`, `code/seam_thin_deep.py`).

## 6. Dispersion-type engines

- Off-diagonal bilinear correlations C_{k,k′}: exact half-normal
  (0.801 vs 0.798, 500 form pairs) (`code/dispersion_engine.py`).
- Shifted-Möbius Dirichlet polynomial |G(t)|/√J: mean 0.875 (Rayleigh
  0.886), sup 1.60 over t ∈ [0, 2000] — cleaner than a random-sign
  control (`code/shifted_mu_poly.py`).
- Cross-scale entropy channel: corr(dual(k), dual(pk)) =
  −0.24/−0.40/−0.23/−0.16 for p = 2/3/5/7, controls ≈ 0
  (`code/e1_dilation.py`).
- Corrected-normalization ledger (post-REVIEW_VERDICT): for the
  prime-indexed correlations C_{k,k′} = Σ_{p∼P} μ(N−pk)μ(N−pk′)
  (2 values of N × 2 K-bands × 1500 pairs): mean |C|/√n_p =
  0.348–0.370 (BELOW half-normal 0.798 — stronger than random),
  variance ratio E|C|²/n_p = 0.36–0.40, kurtosis of C = 5.3–5.8
  (heavy-tailed relative to its own variance — mixture structure).
  Against the correctly normalized requirement
  Σ|C| ≪ K²P(log)^{−2A}, nature over-delivers by a factor √P
  (`code/e1_corrected_norm.py`).
- The mixture resolved (6000 pairs, four conditionings): pairs with
  gcd(kk′, N) = 1 are EXACTLY Gaussian (|C|/√n_p mean 0.740 vs
  half-normal 0.798, variance ratio 0.872, kurtosis 3.12); pairs
  sharing factors with N are variance-SUPPRESSED (0.21–0.40) — the
  global heavy tail is an inter-class mixing artifact (tail-membership
  lift < 1 for every arithmetic conditioning). Structure:
  C_{k,k′} = (local factor determined by gcd with N) × Gaussian —
  the θ=1/8 structure law reappearing at dispersion level
  (`code/e1_tail_anatomy.py`).
- Zero-accounting (12,000 pairs × 2 N, per-pair zero terms counted
  directly): the free class closes COMPLETELY — after removing the
  forced μ=0 terms (q | N, q | k forces q | N−pk; q² | N kills
  entirely, observed z = 1.000), gcd-free pairs give
  E|C|²/support = 0.97–1.02 with kurtosis 2.99–3.03: exactly
  Gaussian. Honest correction to the "stronger than random" reading:
  roughly half the suppression is forced-zero bookkeeping (classical,
  computable); the remainder — shared classes keep m2_eff = 0.13–0.67
  on their nonzero support — is a genuine negative correlation among
  surviving terms (`code/e1_zero_account.py`).

## 7. The factorization law (unifying stamp, increments 144–149)

Every family this program probed — dilated shifted-prime correlations
C_{k,k′} = Σ_p μ(N−pk)μ(N−pk′) and thin-progression Möbius sums —
obeys one measured law:

> **μ-field = (deterministic local mask, computable by finite modular
> enumeration) × (exactly Gaussian fluctuation at half-normal scale).**

Evidence chain: (i) every viable exact (v₂, v₃)-cell of the dispersion
field gives variance ratio 0.99–1.06 and kurtosis 2.8–3.1 on its
nonzero support, including the worst joint cell (n = 4240:
1.031/2.94); annihilated cells are predicted deterministically
(v_q(N) = 1 ∧ v_q(k) = 1 ⟹ μ ≡ 0) (`code/e1_exact_cells.py`);
(ii) thin progressions, zero-accounted: all classes 0.957–0.960 with
kurtosis 3.03–3.05, non-squarefree gcd = predicted annihilation; the
raw 0.581 splits exactly as 0.782·√(support 0.607)
(`code/e1_thin_closure.py`).

**Blind verification of the mask** (a priori, no fitting): computing
the mask by exact enumeration over the units mod q² for q ≤ 50 plus
the tail factor ∏_{q>50}(1 − 2/q²) = 0.99228, then testing on 4000
fresh pairs: corr(predicted, observed support) = 1.0000, max error
0.027, bucket agreement exact (predicted 0.835 → observed 0.835;
predicted annihilation → observed annihilation), and the amplitude
prediction 0.798·√s_pred matches the measured |C|/√n_p to 1.5%
(`code/e1_mask_model.py`). The mask is an algorithm, not an
observation.

**Matrix-level closure (with correction #26)**: the full 400×400
C-matrix on a k-band has spectral norm 36.76 vs the correct
factorization null (random-sign rows on the real support pattern,
Gram pipeline) 36.89 ± 0.68 — z = −0.19, dead center; entry mean
0.006 ± 0.005 (no singular-series main term); top eigenvectors
delocalized. An earlier z = 9 "gap" was a null-design category error
(iid-entry Wigner null applied to a Gram matrix; corrected the same
day, before publication) (`code/e1_spectral.py`,
`code/e1_wishart_null.py`). The law therefore holds at every level
measured: pair statistics, cell statistics, and spectral norm.

**Retroactive corrections this law forces**: every "sub-half-normal /
stronger-than-random" reading earlier in this document (the 0.738
thin stamp of §5, the 0.35 dispersion means of §6) is mask
accounting, not super-random cancellation. On its nonzero support
nothing we measured beats a coin; nothing loses to one either. The
proof-relevant content is now maximally distilled: all arithmetic
structure is classical local bookkeeping; the single unproven
statement is square-root cancellation of the Gaussian part.

**Complete historical closure** (increments 155–157): the law also
explains, with no additional channels, (i) the integer-indexed field
— E1's actual object D(k) and its pair version — no mean field
(per-term mean 0.00017 ± 0.00024), no k-class structure,
leave-one-out r1 = 0.800; (ii) the SEAM band, the program's deepest
apparent anomaly: 43.4% of its pairs are deterministically
annihilated, predicting a naive average of 0.798 × 0.566 = 0.452 vs
0.450 measured; the viable pairs read exactly half-normal
(0.795, variance ratio 1.004). Every sub-random reading in this
program's history (0.29–0.74, all families, all bands) is now
explained by the single factorization law
(`code/e1_integer_field.py`, `code/e1_pair_local.py`,
`code/e1_seam_law.py`).

## 8. The ln 2 constant (where this program began)

Dispersion of g/HL under a fixed canonical estimator: fitted limit
0.6931 ± 0.0061 ≡ ln 2 (18 octaves), two independent null
constructions, window-shape invariant; destroyed by gap permutation —
an invariant of ≥3-point additive prime correlations
(`code/e4_dense_sample.py` et al.).

## 9. The wall's own scalar: C(N) ≍ √N, and what the comet measures

Theorem C gives, unconditionally,
r(N) = 𝔖(N)(N − C(N)) + R(N), where
C(N) = Σ_{n<N}Λ(n)μ(N−n) is the scalar the whole chain reduces to and
R(N) = E₃(α) + O_A(N(log N)^{−A}) is the discrepancy. Measured over
400 even N in five octave groups from 1.2·10⁵ to 1.9·10⁶
(`code/h_deficit.py`):

| N₀ | corr(C, D) | mean\|C−D\|/N | sd(C/N) |
|---|---|---|---|
| 1.2·10⁵ | 0.457 | 0.00838 | 0.01071 |
| 2.4·10⁵ | 0.746 | 0.00730 | 0.00612 |
| 4.8·10⁵ | 0.771 | 0.00455 | 0.00495 |
| 9.6·10⁵ | 0.831 | 0.00372 | 0.00310 |
| 1.9·10⁶ | 0.726 | 0.00294 | 0.00271 |

with D(N) := N − r(N)/𝔖(N) the relative Goldbach deficit. Fitting
across the octaves:

> **|C(N)| ~ N^{0.503}** — square-root, to three digits
> **|R(N)| ~ N^{0.599}**

⚠️ **"To three digits" is refuted, and the exponent is not 0.503**
(increment 281, `code/lab_wall_exponent_reaudit.py`). This design uses
five groups of **80 consecutive even N** — each group spans 160 in `N`,
and a scale estimated from 80 samples has sampling CV ≈ 7.9%. Re-running
*this exact design* 500 times with the group offsets varied (fixed seed)
gives

> **β = 0.516 ± 0.0426**, 5th–95th percentile **[0.447, 0.585]**.

The recorded 0.503 sits at **percentile 39** — an entirely ordinary
draw. The replication spread is **43× the 0.001** that "three digits"
asserts. Not an argument that the design is imprecise: a demonstration.

**The real number**, from the full census of every even `N ≤ 1.6·10⁷`
with the location mask removed:

> **β = 0.5457 ± 0.0032** (raw, mask left in: 0.5315 ± 0.0019).
> `β = 1/2` is excluded at **14.5σ**; the mask alone shifts β by +0.014.

So `|C(N)|` grows **measurably faster than `√N`**, which is what §12's
variance law says it must and what this section read as confirming the
opposite. The requirement is still only `o(N)`, so nothing about the
program's position changes — but "square-root to three digits" was
three digits of noise around a number that is not `1/2`.

⚠️ **One consistency check here is not a check.** `β` and §12's `α`
are the *same fitted numbers in two parametrisations* (`sd(C)` is
`sd(Z)` times ~`√(𝔖N)`), so `β = 0.5 + (α/2)(Δloglog/Δlog)` holds by
construction; the agreement to 6·10⁻⁴ against α = 1.30 confirms nothing.
The one informative row is `α = 0`, which the old law predicts and the
data refute. A genuinely unforced diagnostic is `E|C|/sd`, which is
**0.79145** against the Gaussian `√(2/π) = 0.79788` — a consistent
**−0.81%**, the sign expected from the residual excess kurtosis.

Two readings, one positive and one negative.

- **The chain's slack on its own final object is √N.** It needs only
  C(N) = o(N); nature delivers N^{1/2}. That is the same shape as
  every other margin in this program: nature over-delivers by a power,
  and no technique certifies any of it.
- **The comet is not measuring C(N).** Since
  D = C − R/𝔖 and R grows faster than C, the deviation of the Goldbach
  count from 𝔖(N)N is dominated by R, not by 𝔖C. So the ln 2 corpus
  of §8 and the final-boss scalar are **not** the same object — a
  hypothesis raised here and refuted by its own pre-registered
  criterion (corr(C,D) ≥ 0.9 required; 0.60 pooled). The modest
  correlation is the signature of R dominating, not of noise.

At these N the O_A term in R is itself of size (log N)^{−2} ≈ 0.005 of
N, comparable to what is measured, so R cannot be separated into E₃
and the error; the exponent above is for the total discrepancy.

### Re-audit of this closure (increment 282)

`code/lab_comet_closure_reaudit.py`, full census of every even
`N ≤ 1.6·10⁷`, mask removed. **The closure stands. The evidence
recorded for it did not.**

**The conclusion is right, and now established.** `β_R − β_C = +0.0997
± 0.0034` — **29 s.e.** R does grow faster than C, so D is dominated by
R and the comet is not measuring C(N).

**But §9's own design could not have shown that.** Replicating it 500
times (five groups of 80 consecutive even N, offsets varied) gives a
difference of `+0.079 ± 0.054`: the recorded gap 0.096 is **1.8 s.e. of
that design**, and **7.0% of replicates return a difference ≤ 0**. A
refutation was recorded on evidence that would have come out the other
way one time in fourteen.

**The criterion statistic itself was mis-estimated.** §9 pooled five
groups and got `corr(C,D) = 0.60`. Pooling across bands of different
scale *attenuates* a correlation; per band the de-masked value is
**0.80 falling to 0.71**, and pooling every band here reproduces the
same artefact (0.7064). The true correlation was never as low as 0.60.

**And the mask does the opposite of what one would guess.** Removing it
*raises* the correlation by **+0.042** — the mask was suppressing the
association between C and D, not manufacturing it.

**What makes the closure safe is the trend, which was the one thing not
checked.** The de-masked correlation **falls** monotonically,
`−0.0198 ± 0.00035` per unit `log N`. It is moving *away* from the 0.9
threshold, so the criterion will not be met at larger N. Had it been
rising — the pattern §9's own five noisy groups appeared to show
(0.457, 0.746, 0.771, 0.831, 0.726) — the closure would have been
premature.

> **Corrected exponent: |R(N)| ~ N^{0.6458 ± 0.0065}**, not `N^{0.599}`.
> `β_C = 0.5461 ± 0.0033`, independently reproducing increment 281's
> 0.5457 ± 0.0032 by a different estimator.

The caveat above is unchanged and still governs: at these N the `O_A`
term is comparable to what is measured, so `β_R` is the exponent of the
**total** discrepancy and not of `E₃`.

## 10. Three hypotheses on the wall, and one signal chased down

All pre-registered before running (`code/hyp_round.py`,
`code/hyp_h2b.py`).

**H1 — do the actual Vaughan weights buy cancellation the L² can't
see?** The chain consumes the *signed* sum Σ_{k∼K} b_k D(k), and
Cauchy–Schwarz throws away b's structure; only b = 1 had been tested.
One weight is special: over the full ranges Σ_k μ(k)D(k) =
Σ_u μ(N−u)(μ∗μ)(u), a linear μ-sum against a bounded multiplicative
function. Measured |Σ b·D| / (‖b‖₂·rms D) over five band/N
combinations for b ∈ {1, log k, μ(k), μ(k)log k, Λ(k)}: values scatter
from 0.32 to 1.95 with no weight consistently small (b = μ gives 0.43,
0.39 at one N but 1.52 at the next band). **DEAD** — every weight sits
at the random-sign value, as Conjecture L predicts.

**H3 — is the L² mass concentrated on few k?** Participation ratio
P = (Σ|D|²)²/(K·Σ|D|⁴), which is 1/3 for a Gaussian field. Measured
0.293, 0.232, 0.221, 0.274, 0.265. **DEAD** — mildly heavier-tailed
than Gaussian, nowhere near the ≤ 0.15 that isolation would need.

**H2 — is there a mask on C(N)?** The registered criterion was a
spread test, and it returned DEAD (sd ratio 1.049 across v₂ bins).
But the *means* split by 3 | N at ≈ 6.7σ — a signal the criterion had
not been written to catch, the same under-specification as the
single-threshold cut in §9's neighbourhood. Chased down in three
disjoint N-windows:

| window | mean(3\|N) | mean(3∤N) | z | z after ÷𝔖(N) |
|---|---|---|---|---|
| 3·10⁵ | −0.346 | +0.661 | −4.38 | −4.91 |
| 7·10⁵ | −0.248 | +0.471 | −3.17 | −3.17 |
| 1.4·10⁶ | −0.647 | −0.111 | −2.29 | −0.76 |

**Real, and lower-order, and it is the mask.** Measured at full power
over every even N in [10⁵, 4·10⁶] rather than in three windows of ~500
(`code/lab_split_decay.py`):

| band | n | split | z | perm z | /√N | z of C/√𝔖 | z of C/𝔖 |
|---|---|---|---|---|---|---|---|
| 10⁵–2·10⁵ | 50000 | −620.0 | **−52.5** | −0.70 | −1.649 | −50.0 | −45.3 |
| 4·10⁵–8·10⁵ | 200000 | −785.8 | **−66.2** | +0.84 | −1.045 | −60.7 | −53.6 |
| 1.6·10⁶–3.2·10⁶ | 800000 | −959.5 | **−79.4** | +0.26 | −0.638 | −70.4 | −60.4 |

Two of the three original readings stand.

1. **The split is real** — 50 to 80 σ per band against a permutation
   control at |z| ≤ 1.8.
2. **It is lower-order.** `split/√N` falls as `N^{−0.315}` (the three
   windows gave −0.41; a two-parameter power law fitted to points
   carrying 23%, 32% and 44% relative error determines little). The
   absolute split grows like `N^{0.19}`, so it is indeed not a mask on
   `C(N)/√N`.
3. **But it is the location mask**, seen through a shallow cut. `3 | N`
   averages over mostly-shallow N; isolating the deep cells gives
   `C/√N` around −20, falling far more slowly (LOCATION_MASK.md). And
   the claim that 𝔖-normalisation removes it does not survive power:
   at full power `C/√𝔖` still reads z = −70 and `C/𝔖` z = −60 — a
   reduction of 10–15%, not a removal. "Below 3σ at the largest
   window" was n ≈ 500.

## 11. Two more hypotheses, and a methodological correction

**H10 — do nearby N carry correlated copies of the dilate field?**
Moving N to N+h shifts every Möbius argument by h, so a nonzero
corr_k(D_N, D_{N+h}) would be a binary Chowla signal. Over 18 triples
(three N, h ∈ {2,6,30}, two bands, 1200 k each) the standardised
z = corr·√n has **mean |z| = 0.814 against a standard-normal null of
0.798**, max |z| = 1.92. **DEAD**, and the dead result has content:
nearby N carry *independent* copies of the field, which is precisely
why the N-averaged route buys an exceptional set rather than a
fixed-N statement.

**H12 — does some range of n carry the mass of C(N)?**
For equal-log-length ranges, g = |Σ_range f|/√(#terms) is constant
under Brownian behaviour. Measured max/median g = 3.03, 3.31, 2.03,
3.03 across four N — but the **argmax range differs every time**
(7, 3, 5, 6), and max/median ≈ 3.2 is exactly what eight half-normal
draws give. **DEAD**: the walk is Brownian and no range is special.

**The methodological correction, and it is the important part.** Both
hypotheses first returned ALIVE, and both verdicts were artefacts of
thresholds chosen without computing the null:

- H10 used |corr| ≥ 0.10 on 300 points, where the null SE is
  1/√300 = 0.058 — a 1.7σ threshold, below the noise floor. Two hits
  in twelve tests is what chance delivers (expected 1.1).
- H12 used "one decade of n carries ≥ 50%", but the top decade is
  always the last, spanning ~89% of the range; under Brownian
  behaviour its expected share is ~0.95, so the measured 0.53–0.74 was
  *less* than random.

Counting the single-threshold cut of §9's neighbourhood and the
spread-only test of §10, this session produced **five under-specified
pre-registrations**. The pattern is now named and the rule is
explicit: **compute the null value of the statistic before choosing
the threshold, and state it in the pre-registration.** Every test in
§11 carries its null in the same line as its criterion.

## 12. A twenty-hypothesis sweep, and a law for the wall's own scalar

Twenty hypotheses, run in three parallel sweeps, each criterion
carrying its null on the same line (`code/sweep_A.py`, `sweep_B.py`,
`sweep_C.py`, `sweep_B2.py`, `sweep_B3.py`).

**Sweep A — ten hypotheses on the local structure of D(k)**, over
n = 2861 values of k at N = 10⁷: splits by gcd(k,N), ω(k), v₂(k),
smoothness, k mod 3 and mod 4; skewness, kurtosis, a runs test, lags
1–8, and two long lags. **0 flags at |z| ≥ 4 out of ~22 statistics**,
with mean|d| = 0.7890 against the half-normal null 0.798. The dilate
field has no local structure at any coordinate tested.

**Sweep C — five variants of the scalar**, exponent fitted over four
octaves: pure Möbius pair Σμ(n)μ(N−n) gives α = 0.495, primes-only
0.537, twisted 0.535, against C(N)'s own 0.503; the mask-only object
Σ Λ(n)(μ²(N−n) − 6/π²) gives α = 0.999, linear as it must. **No
softer target**: every relative of C(N) that fluctuates is
square-root sized. (Two design faults recorded: the Mertens
calibration line is unusable because M(N) varies too slowly for
averaging over consecutive N, and the (−1)ⁿ twist is nearly a global
sign flip, so C5 duplicates C2 rather than testing anything.)

**Sweep B — seven hypotheses on C(N), and this one fired.** Five flags
at |z| ≥ 4: kurtosis 4.446 (z +11.4), corr(c, 𝔖) −0.208 (z −8.0), sd
ratio 1.408 for ω(N) ≥ 5 (z +6.4), 0.806 for ω(N) ≤ 3 (z −5.2), and
1.261 for N ≡ 0 mod 5 (z +5.7).

**Four of the five are scale, and √𝔖 fixes them.** Dividing by 𝔖
**overshoots** (the ratios invert to 0.708, 1.287, 0.820), and solving
R = (S₁/S₂)^β across the three splits gives β = 0.497, 0.539, 0.460, so
the scale mask is **√𝔖**. Dividing by it:

| after ÷√𝔖 | N₀ = 2·10⁵ | 5·10⁵ | 9·10⁵ |
|---|---|---|---|
| sd ratio, ω≤3 | 0.997 | 1.040 | 1.047 |
| sd ratio, ω≥5 | 1.145 | 0.998 | 0.966 |
| sd ratio, N≡0(5) | 1.083 | 1.007 | 0.961 |
| kurtosis | 3.333 | 3.397 | **2.983** |
| corr(·, 𝔖) | −0.194 | −0.138 | **−0.070** |

**The fifth flag is not scale, and √𝔖 does not remove it.** A scale
mask rescales a fluctuation; it cannot create a correlation between
C(N) and 𝔖(N). The `corr(·, 𝔖)` row above falls from −0.194 to −0.070
across three windows of ~1500 values each, which is a decline but not a
disappearance, and at full power it is not close to zero:

| statistic | corr with 𝔖 | z (n = 1.95·10⁶) |
|---|---|---|
| C | −0.1386 | **−193.5** |
| C/√𝔖 | −0.0984 | **−137.4** |
| C/√V, exact scale removal | −0.1517 | **−211.8** |

Dividing by the exact per-N second moment √V flattens the scale
completely (sd 0.90, 0.89, 0.88, 0.91, 0.94 across 𝔖 bands) and leaves
the mean varying by two orders: 0.012, 0.121, 0.419, −0.002, −0.324,
**−1.682**. That residual is a **location mask**, and it is the largest
single signal in the whole sweep.

> **C(N) = m(N) + √(κ·𝔖(N)·N·log N)·G(N)**, with `κ ≈ 0.465`, `m(N)`
> the location mask of LOCATION_MASK.md, and G Gaussian once **both**
> masks are applied — every tail inside 3 SE and the extreme at
> z = +0.61.

⚠️ **The `log N` and the `κ` are not measured** (increment 280). The
exponent in `(log N)^α` walks with the fitting window — 1.60 → 1.30 with
the mask removed, 0.83 → 1.02 with it left in — and neither has converged;
the drifts are ten times the standard errors. Over this range a pure power
`N^ε` fits **as well or better** (residual ratio 0.45 raw, 1.91 de-masked,
against a bar of 2 fixed before the run), so the log factor is *consistent
with* the data rather than *shown by* it, and `κ` is whatever `α` is
assumed to be. **What is measured is that `α = 0` is excluded**: the
variance grows faster than `𝔖N`. Read the formula as a shape with one
established feature, not as a fitted law.

This extends Conjecture L to the wall's own scalar, which the campaign
had never done, but on **three** masks and not one: support, scale and
location. Two things in the earlier reading were wrong and both are
recorded (corrections #36 and #45). The scale carries a factor
√(log N) that was absorbed into a constant over too narrow a range of
N. And the "mean drift, already insignificant at N ≈ 9·10⁵" was the
location mask seen and dismissed: at 1.95·10⁶ values the 3 | N sign
split reads z = −94.3 against +61.6 off it, and the 3 | N mean split of
§11 is that same mask, not a residual drift of the scale law.

## 13. The wall's variance is a weighted shift-averaged Chowla sum

The law of §12 is a statement about the variance of C(N), and the
variance unfolds exactly (`code/wall_secondmoment.py`):

> **Σ_{N∈W} |C(N)|² = Σ_h r_W(h)·S_W(h)**, where
> r_W(h) = Σ_{n−n′=h}Λ(n)Λ(n′) counts prime pairs at shift h and
> S_W(h) = Σ_v μ(v)μ(v+h) is the binary Chowla correlation.

Verified to machine precision (differences 2·10⁻⁸ and 1.2·10⁻⁷ at
X = 2000 and 4000). So **the variance of the wall's own scalar is a
prime-pair-weighted, shift-averaged Chowla sum** — and shift-averaging
is exactly the regime Matomäki–Radziwiłł–Tao control, unlike the
dilate-averaging of §everything-else.

Two readings, one of which corrects a prediction of ours.

- **The off-diagonal is negative, and small.** Measured over a proper
  range (`code/lab_offdiagonal_reaudit.py`), `off/LHS` reads 0.554 at
  X = 2000 and 0.498 at X = 4000 — reproducing the first reading —
  then **0.077 at 10⁵, −0.144 at 10⁶ and −0.208 at 4·10⁶**. It falls
  through zero at X ≈ 3·10⁵. On upper halves, where the asymptotics
  are cleaner, the same sequence is 0.475, 0.076, −0.097, −0.146. So
  the shifted-Möbius correlation does not supply half the variance; it
  **removes about 14% of the diagonal**, and the earlier reading was
  an artefact of X ≤ 4000.
- **The diagonal alone already sets the scale.** `diag/(n·N·log N)`
  measures 0.751, 0.787, 0.800, 0.806 across the same windows,
  converging on `A(N) = 0.787`. So this identity, correctly evaluated,
  implies `Var C(N) ≍ N log N` on its own — and §12's `𝔖(N)·N` was
  short by a log. **The document carried its own refutation** and the
  contradiction went unnoticed for fifty increments (correction #36).
- **The gap to MRT is a factor X, not a log power.** MRT gives
  Σ_{|h|≤H}|S(h)| = o(HX), an average saving of o(1) per shift against
  the trivial |S(h)| ≤ X; pinning the identity to the precision a
  fixed-N statement needs would require order 1/X per shift. And a
  second moment yields almost-all-N in any case, which is the
  exceptional-set result already on the map.

The value of the identity is as a map coordinate: it is the one place
where the wall's own object connects to shift-averaged rather than
dilate-averaged technology, and the connection is exact. What it does
not do is close anything.

## One-shot reproduction

```bash
python code/verify_all.py   # ~15 min: mu sanity, ladder identity,
                            # dispersion, thin-progression stamps
```

*Every number above regenerates from the scripts named. Questions,
corrections, and literature pointers are explicitly welcome.*
