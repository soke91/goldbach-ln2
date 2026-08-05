# Measurements — the empirical core of this repository

*This document contains only reproducible measurements and exact,
machine-verified identities. No proof-program claims. For those (with
their own honesty flags), see PROOF_SKETCH_E1.md; they are secondary
to what follows.*

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

## 7. The ln 2 constant (where this program began)

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
