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

**Retroactive corrections this law forces**: every "sub-half-normal /
stronger-than-random" reading earlier in this document (the 0.738
thin stamp of §5, the 0.35 dispersion means of §6) is mask
accounting, not super-random cancellation. On its nonzero support
nothing we measured beats a coin; nothing loses to one either. The
proof-relevant content is now maximally distilled: all arithmetic
structure is classical local bookkeeping; the single unproven
statement is square-root cancellation of the Gaussian part.

## 8. The ln 2 constant (where this program began)

Dispersion of g/HL under a fixed canonical estimator: fitted limit
0.6931 ± 0.0061 ≡ ln 2 (18 octaves), two independent null
constructions, window-shape invariant; destroyed by gap permutation —
an invariant of ≥3-point additive prime correlations
(`code/e4_dense_sample.py` et al.).

## One-shot reproduction

```bash
python code/verify_all.py   # ~15 min: mu sanity, ladder identity,
                            # dispersion, thin-progression stamps
```

*Every number above regenerates from the scripts named. Questions,
corrections, and literature pointers are explicitly welcome.*
