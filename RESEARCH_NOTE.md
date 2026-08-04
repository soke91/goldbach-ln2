# An Empirical Sub-Poissonian Constant for Goldbach Representation Fluctuations, with Numerical Evidence That It Equals ln 2

**Status: numerical/exploratory research note — no theorems about the constant are proven here.**
*Produced in a human–AI exploratory session (2026-08-03). All code and logs available; 15 documented self-corrections. Expert scrutiny is explicitly requested; overlap with existing literature is possible and pointers are welcome.*

---

## Abstract

Let g(n) denote the number of unordered representations of an even integer n as a sum of two primes, and let HL(n) be the first-order Hardy–Littlewood prediction. We study the fluctuations of z(n) = g(n)/HL(n) with a fixed ("canonical") estimator: dyadic windows, residue class n ≡ 2 (mod 6), removal of the linear drift in ln n, and the normalized dispersion r = σ(z)·⟨HL⁻¹ᐟ²⟩⁻¹-type ratio calibrated so that an independent (Bernoulli/Cramér) model gives 1. Across 18 dyadic windows up to 5.4×10⁸ we measure a monotone approach of r to a limit; a weighted least-squares fit gives

**r∞ = 0.6931 ± 0.0061 (stat), to be compared with ln 2 = 0.69315.**

Two further independent constructions — a "Hardy–Littlewood null" comparison world (wheel-structured Bernoulli), and the plateau of a wheel-modulus scan with an interaction-complete regression basis — give variance-ratio values 0.468 ± 0.010 and 0.472 ± 0.010, consistent with (ln 2)² = 0.4805. The dispersion is invariant under changing the window shape ratio from 1.3 to 6 (ruling out a window artifact). We conjecture that the canonical suppression constant equals ln 2 exactly.

We further present evidence that the constant is **invisible to two-point statistics**: permuting the sequence of prime gaps (globally, or locally in blocks of 8), which preserves the entire gap multiset, destroys the suppression completely (ratio returns to 1.000 ± 0.007); an elementary calculation with Montgomery-type pair-correlation input yields only O(1/log) corrections. If correct, the constant is an invariant of **additive correlations of order ≥ 3** of the primes — an object beyond GUE pair-correlation phenomenology.

---

## 1. Setup and definitions

- P = set of primes. For even n, g(n) = #{(p,q): p ≤ q, p+q = n, p,q ∈ P}.
- First-order Hardy–Littlewood: HL(n) = C₂·𝔖₀(n)·n/ln²n with C₂ = 0.6601618…, 𝔖₀(n) = ∏_{p|n, p>2} (p−1)/(p−2).
- **Canonical estimator.** For a dyadic window W = [2^k, 2^{k+1}) restricted to n ≡ 2 (mod 6):
  1. z(n) = g(n)/HL(n);
  2. subtract the best linear fit of z against ln n over W (drift removal);
  3. r(W) = std(residual)·[mean(HL⁻¹ᐟ²)]⁻¹ — i.e., the dispersion measured in units of the "Poisson" scale HL⁻¹ᐟ².
- **Null calibration.** In matched independent models (odd sites m carrying "prime" with probability 2/ln m, all structure identical in expectation), the same pipeline returns r ≈ 1 (verified; see §4 for the refined null).

All prime data computed by sieve; g(n) via FFT convolution of prime indicators (exact integer counts, cross-checked against direct summation at 10⁻⁹ relative error).

## 2. Main empirical result

Dyadic-window dispersion values (pure-fluctuation variant; windows 2¹² … 2²⁸, the last four measured on 8000-point subsamples with per-point exact g):

r rises monotonically 0.506 → 0.595 over 17 octaves. Fitting r = a + b/ln n:

| fit range | a (limit) |
|---|---|
| all 18 octaves, weighted | **0.6931 ± 0.0061** |
| tail 10 | 0.6864–0.6903 |
| tail 8 | 0.6898–0.6928 |

ln 2 = 0.693147. The successive fitted limits over growing data ranges were 0.633 → 0.665 → 0.677 → 0.684/0.695 — a monotone march toward ln 2 (consistent with a fit that underestimates in the presence of higher-order corrections).

**Window-shape invariance.** At fixed central scale 4×10⁶, changing the window ratio r_w ∈ {1.3, 1.5, 2, 3, 4, 6} changes the measured dispersion only from 0.5650 to 0.5824 (a scale-mixing drift), while ln r_w varies by a factor 6.8. The constant is not an artifact of octave (ln 2-wide) windows.

**Conjecture 1 (ln 2 conjecture).** For the canonical estimator, r∞ = ln 2 (equivalently, the variance ratio tends to (ln 2)² = 0.4805).

## 3. Two further independent constructions

**(a) Hardy–Littlewood null world.** Comparing the real primes against wheel-structured Bernoulli worlds (sites coprime to 30, density boosted accordingly) under an *identical* projection estimator whose regression basis spans both worlds' arithmetic structure (polynomial × singular-series interactions plus per-prime indicator columns), the Fano-type variance ratio is **0.4680 ± 0.0100**.

**(b) Wheel-modulus scan.** Extending the null family to wheels of modulus z ∈ {3,5,7,11,13} with an interaction-complete basis (products of indicator columns; without them a spurious dip appears at z = 7 and is fully explained by unspanned joint-residue-class structure), the ratio forms a plateau **0.472 ± 0.010** at z = 11–13.

Both agree with (ln 2)² = 0.4805 within ~1σ. We emphasize a caution learned repeatedly: with *incorrectly specified nulls or bases*, the measured constant can land anywhere in ≈0.16–0.52; the three constructions above are the ones whose null/basis specification we can defend, and they agree.

## 4. Structural results (verified; provable or conditional-standard)

1. **Decomposition identity** (exact algebra, machine-verified to 10⁻⁸): with Λ the von Mangoldt function and e = Λ − 1, R = Λ⋆Λ (additive convolution) satisfies R(n) = (n−1) + 2[ψ(n−1) − (n−1)] + (e⋆e)(n). The fluctuation is dominated (99.8% of variance) by the pair term.
2. **Conditional pair singular series** (elementary derivation + numerical verification at the 1% level on totals, 3–15% pointwise): in the band n ≡ 2 (mod 6) the admissible lower primes are p ≡ 1 (mod 6), forcing pair distances d ≡ 0 (mod 6), with conditional pair factor κ(d) = C₅·∏_{q≥5, q|d} (q−1)/(q−2), C₅ = ∏_{q≥5}(1 − (q−1)⁻²) = 0.88022.
3. **Spectral profile = interval-variance (Goldston–Montgomery) function.** The structure-subtracted spectrum of a prime-indicator window, band-averaged off the rationals, follows φ(λ) ≈ 1 − ln λ / ln x across ~4 decades of wavelength λ (checked at multiple scales). This is the GUE/pair-correlation interval-variance profile appearing as the frequency-resolved content of the prime indicator.
4. **Whitening.** After band-restriction on an arithmetic progression (step 6) and low-rank regression, the residual field is spectrally white (flat ratio across all λ; kurtosis 3.14; no autocorrelation): the estimator "bakes in" the constant globally. Consequently the constant is formed where the dense arithmetic (rational spikes and their tails) meets its low-rank *coherent* subtraction — sharp masking of arcs is impossible (they are dense), and eight attempted shortcut routes all converged to this same obstruction (documented).

## 5. Evidence that the constant is beyond pair correlations

- **Gap-permutation surrogates.** Randomly permuting the sequence of prime gaps (globally: preserves the exact gap multiset; or locally within blocks of 8 gaps) yields dispersion ratio 0.999 ± 0.007 — the suppression vanishes although all two-point/interval statistics of the process are (essentially) preserved.
- **Matched-spectrum Gaussian worlds** (random phases with the correct zero-density power spectrum on RH scaling) show *no* suppression (ratio ≈ 1.28–1.39, consistent with the Beta-function geometry constant 2 ln 2 of the additive convolution).
- **Elementary pair-input calculation.** Feeding Montgomery/Goldston-type pair correlation into the variance of Σ_p 1_P(n−p) gives corrections of size O(1/ln) — no O(1) suppression (session-grade derivation; the corresponding numerical check requires full residue-conditioning machinery and is left open).

**Conjecture 2 (beyond-pair).** The constant of Conjecture 1 is determined by connected additive correlations of the primes of order ≥ 3 (plausibly the 4-point function under the constraint p₁+p₂ = p₃+p₄ = n), and is not a functional of the pair correlation alone. Equivalently, it lives in the same regime as minor-arc fourth-moment control of exponential sums over primes.

## 6. Proof program and obstruction map (for the interested analyst)

The full internal roadmap (lemmas established, two remaining computations "Theorem A: closed form of the canonical estimator's spectral measure" and "Theorem B: the GM-profile integral equals (ln 2)²"), together with a list of **15 documented pitfalls** (parity of sampling bands; circular-convolution artifacts; collision-thinning in jitter surrogates; ordered/unordered singular-series factors; aliasing of spiky fields under block averaging; sharp arc masks being (N,Q)-crossings rather than constants; window overlap; density of arcs; etc.) exists as dated internal session logs (in Korean, 30+ increments); available on request via an issue.

A geometric identity we note for interest: [∫_{1/2}^{1} dv/v]² = (ln 2)² — the squared logarithmic measure of the Goldbach partner half-range [n/2, n]; a stationary-phase analysis of the zero double-sum localizes at the midpoint p ≈ q ≈ n/2, but we could not (yet) turn this into a derivation.

## 7. Limitations and honest caveats

- **Nothing about the constant is proven.** Conjectures 1–2 rest on numerical evidence up to 5.4×10⁸ (dispersion series), 1.34×10⁸ (per-point subsamples), auxiliary experiments at 2×10⁶–1.7×10⁷, plus a Goldbach verification sweep to 10¹¹ (no counterexamples; consistent with prior published verifications to 4×10¹⁸).
- **Estimator dependence.** Different estimator frames (different structure models/nulls) yield different constants (observed 0.26–0.52). The ln 2 claim is specifically about the canonical estimator of §1; its distinguished status is argued (HL-normalization only; null-calibrated; window-invariant) but a frame-independent formulation is part of the open program.
- **Extrapolation risk.** The limit is a fitted asymptote of a slowly converging sequence; the fitted value moved from 0.633 to 0.693 as data grew. We report this drift openly; it is consistent with (and was predicted by) a 1/ln n approach, but larger computations would strengthen or refute the claim. Falsification condition: the fitted limit departing from ln 2 with more data.
- **Literature.** We have not performed a systematic literature comparison; variance of Goldbach representations has been studied (e.g., in the ψ-weighted setting, and conditionally under RH), and the constant may exist in some form in prior work. Corrections and references are welcome.
- **Provenance.** The work was carried out in an AI-assisted exploratory session; despite the documented self-correction discipline, undetected errors are possible. All claims are tied to runnable scripts.

## 8. Reproducibility

All experiments are plain Python (numpy only): sieve + FFT convolutions; largest run ≈ 20 minutes on a desktop. Key scripts: series (e4_dense_sample.py, e4_octave27/28.py), window invariance (e4_window_test.py), HL-null and wheel scan (e4_wheelnull.py, e4_zscan3.py), gap surrogates (e4_gapshuffle.py), spectral profile (e4_phi_profile.py), decomposition/whitening (q3d_decompose.py, e4_A_exact.py), κ(d) verification (e4_breakthrough.py), parity ratio & concentration (p2_ratio_theorem.py, p2_concentration.py). Internal session logs (Korean) available on request.

## Addendum A: the parity-barrier ratio equals 1 : ln 2 — a theorem with proof

**Theorem 1 (the 1/(1+ln 2) law).** Let A(x) = {2 < m ≤ x : every prime factor
of m exceeds x^{1/3}}. Then every element of A(x) is either a prime or a
product of exactly two primes (each > x^{1/3}), and as x → ∞ the proportion of
primes in A(x) tends to **1/(1 + ln 2) = 0.59064…** (so the semiprime share
tends to ln 2/(1+ln 2) = 0.40936…).

*Proof.* If m ∈ A(x) had three prime factors (with multiplicity), then
m > (x^{1/3})³ = x, a contradiction; so A(x) consists of primes in
(x^{1/3}, x] and semiprimes qr with x^{1/3} < q ≤ r. By the Prime Number
Theorem the prime count is π(x) − π(x^{1/3}) = (1+o(1))·x/ln x. The semiprime
count is Σ_{x^{1/3} < q ≤ x^{1/2}} (π(x/q) − π(q) + O(1)); since q ≤ x^{1/2}
implies ln(x/q) ≥ ½ ln x, the PNT applies uniformly, giving
(1+o(1))·Σ_q x/(q ln(x/q)) + O(x^{1/2+ε}). Writing q = x^u and using Mertens'
theorem (Σ_{q≤t} 1/q = ln ln t + M + O(1/ln t)) with partial summation,
Σ_{x^{1/3}<q≤x^{1/2}} 1/(q ln(x/q)) = (1+o(1))·(1/ln x)·∫_{1/3}^{1/2}
du/(u(1−u)) = (1+o(1))·ln 2/ln x, since ∫ du/(u(1−u)) = [ln(u/(1−u))] gives
ln 1 − ln(1/2) = ln 2. Hence semiprimes number (1+o(1))·(ln 2)·x/ln x, and the
prime share tends to 1/(1+ln 2). ∎

(The statement is classical in character — counts of almost-primes with
restricted factors are standard — but we include the proof for completeness,
as this constant is the measured "parity-barrier ratio" of our experiments.)
Numerically the error decays faster than 1/ln x: share − limit =
0.0402 / 0.0144 / 0.0061 at x = 10⁶/10⁷/10⁸.

**Conditional corollary (Goldbach partners).** For partner sets m = n − p the
measured share is boosted (+0.0087 at 10⁸) as predicted qualitatively by
prime pair-correlation (Hardy–Littlewood) weighting; a quantitative version is
conditional on the pair-correlation conjecture. Whether this ln 2 (from the
canonical 1/3 dichotomy exponent) is related to the fluctuation constant of
Conjecture 1 remains open.

## Addendum C: the parity-barrier portrait (concentration, closed forms)

Let s(n) be the prime share among n^{1/3}-sieve survivors of the Goldbach
partner set, S(n) the survivor count, and r_sb := std(s)/√(s̄(1−s̄)/S̄) the
sub-binomial ratio. Measurements:

1. **Hyper-concentration (sub-binomial).** std(s) = 0.0043 / 0.0016 / 0.0007
   at n ~ 10⁶ / 10⁷ / 10⁸ — below the binomial prediction at every scale.
2. **The sub-binomial ratio is ≈ 0.85 ± 0.04 (systematics included).**
   Across several sampling designs (windows, detrending, class restriction)
   r_sb ranges 0.80–0.88; design-consistent runs give 0.836(17)/0.884(20)/
   0.845(27) at 3×10⁶/10⁷/3×10⁷. Sub-binomiality (r_sb < 1) is firm at ≥4σ;
   the specific value √(ln 2) = 0.8326 is *compatible but unconfirmed*
   (methodology-sensitive at the ±0.04 level — recorded as an open question,
   not a conjecture).
   **Mechanism of the hyper-concentration (established):** s(n) = g(n)/S(n)
   identically (numerator = the Goldbach count itself), and g, P2, S all carry
   the same multiplicative singular-series profile — measured correlation
   corr(g, P2) = +0.9989 — so the ratio cancels the giant common mode; the
   sub-binomial residual is what remains after this cancellation. A component
   identity r_sb² = (1−s̄)F_g + s̄F_P2 − 2Cov(g,P2)/S̄ is verified numerically
   (0.632 predicted vs 0.637 direct).
3. **Closed-form defense distance.** Since s̄ = 1/(1+ln 2) gives
   s̄/(1−s̄) = 1/ln 2 algebraically, the σ-distance from the all-P2 collapse
   (s = 0, which a Goldbach counterexample would require) is
   D(n) = s̄/std(s) = √(S/ln 2)/r_sb; if Conjecture 3 holds,
   **D(n) = √(S(n))/ln 2**. Verified: predicts 922 at n ~ 10⁸ (measured 917),
   337 at 10⁷ (measured 355). D grows like √(n)/ln²n-scale.
4. **Gaussian tails.** Over 1500 samples the minimum of s sits a routine 2.7σ
   below the mean — no anomalous left tail.

5. **Conservation + parity trade-off (the dynamical structure).** Measuring
   structure-normalized Fano factors of the three counts at 10⁷:
   F̃_g = 0.296 (matching the canonical Goldbach dispersion of Conjecture 1 in
   this scale range — the two phenomena connect), F̃_S = **0.136** (the sieve-
   survivor *total* is hyper-rigid — a new rigidity object), F̃_P2 = 0.440,
   and, after common-mode removal, **corr(ẽ_g, ẽ_P2) = −0.52**: primes and
   near-primes *trade places under a nearly conserved survivor total*.
   The variance budget Var(S) = Var(g)+Var(P2)+2Cov closes to first order
   (0.17 vs 0.14 with a first-order structure model for P2).
6. **The frame's constant table.** F̃_g ≈ 0.27–0.29, F̃_S ≈ 0.13–0.14, and
   F̃_P2 = 0.457(9) at 3×10⁶ rising to 0.490(14) at 10⁷ (7500-sample runs).
   An earlier smaller-sample scan suggested F̃_P2 = (ln 2)² to 0.2%; tighter
   statistics reveal a finite-size upward drift instead (self-correction #17),
   so the correct statement is: F̃_P2 is scale-drifting toward an asymptote in
   the 0.48–0.50 range, *compatible with (ln 2)² but unconfirmed* — settling
   it requires the same multi-decade fitted-limit treatment as the main
   dispersion series. (A recurring lesson of this project: apparent constant
   matches at single scales are unreliable; only fitted limits count.)
7. **First fitted-limit table (7 scales, 10⁶–6.4×10⁷, 12,750 samples).**
   All three Fanos drift upward; weighted 1/ln n fits give limits
   F_g → 0.43, F_S → 0.21, F_P2 → 0.63 (long extrapolations; ±0.05–0.1).
   Suggestive ordering near ((ln 2)², ·, ln 2) — recorded as direction, not
   claim; a comet-grade multi-decade campaign would be needed to settle them.
   Separately, F_S is *window-invariant* at fixed scale (0.137 ± 0.002 across
   window ratios 1.15–2.2 — the same invariance signature as Conjecture 1),
   while a GM interval-variance-at-window-scale explanation is excluded by
   that very invariance (its sharp prediction 0.19→0.06 fails).

Together with Addendum A (s̄ = 1/(1+ln 2)): the parity ratio is pinned at an
elementary constant, fluctuates *less* than independence allows, and the
counterexample-enabling collapse recedes as √S/ln 2·(1/r_sb); dynamically, the
barrier is a *conservation law with anti-correlated parity exchange*.

## Addendum D: an attack map for the sharpened target s(n) > 0 (v2)

Theorem 1 + concentration reduce Goldbach (via the parity barrier) to proving
s(n) > 0 for every large n. A brief empirical/bookkeeping map of routes:

- **Price tag of classical sieves.** The proven population ratio S : P2 =
  (1+ln 2) : ln 2 = 2.443 : 1 means classical upper/lower sieve bounds close
  the gap iff their combined loss factor is < 2.443; standard linear-sieve
  losses (~2× each side) exceed this — a compact numerical statement of what
  the parity barrier costs, and of what any improvement must achieve.
- **Variance route ceiling.** Mean (Theorem 1) + conditional second moments +
  Chebyshev yield only exceptional-set bounds ("almost all n") — the known
  ceiling; individual n is out of reach on this route.
- **Exhaustive envelope is sub-Gaussian-tame.** Over 30,000 consecutive
  candidates at 10⁶, the worst s is 0.5962 — only 3.5σ below the mean (less
  than the ~4.5σ a Gaussian sample of this size would produce), and above the
  asymptotic mean 0.5906.
- **No arithmetic signature: the ratio has already quotiented out structure.**
  A full regression over all 30,000 candidates shows corr(s, 𝔖₀) = +0.004 —
  the singular series explains 0.0% of the variance of s (an earlier
  small-sample "thin-class habitat" reading did not survive this test). This
  is consistent with the hyper-concentration mechanism: g and S share their
  arithmetic mode, so s = g/S cancels it. Residues mod 30 shift class means
  by only ~0.2% relative; the one real effect is that classes with 5 | n are
  ~10% tighter in σ. The remaining fluctuation of s is a single,
  structure-free channel — class-conditioning offers no reduction, and the
  target becomes a tail bound on that universal channel.
- **Sub-binomial variance is a finite-scale effect (protocol v3).** The
  cross-n variance ratio R² = Var(s)/Var_binomial sits below 1 at every
  accessible scale, but a clean measurement required three rounds of
  artifact removal: (i) a window-fixed sieve limit y = X^{1/3} admits
  n-dependent P3 contamination (widening the window at 10⁷ moved R² from
  0.71 to 2.24); (ii) even with per-n roughness, a prime-cube boundary
  q³ inside the window steps s by 1–2σ when q activates; the final protocol
  places windows inside cube gaps (q₁³, q₂³) with drift removal. The clean
  curve rises monotonically, R² = 0.47 → 0.80 over 10⁵–10⁸, with
  (1 − R²)·ln x ≈ 5.2 roughly constant — consistent with **R² → 1 and a
  suppression amplitude decaying like ~5.2/ln x**. A blind extrapolation
  test at 10⁹ (cube-gap window after 997³, 300 samples) predicted
  R² ≈ 0.75 and measured 0.740 ± 0.061 — a direct hit across four decades.
  We therefore do not claim any nontrivial limiting constant; prime
  sampling among survivors appears asymptotically binomial in this
  statistic. (The within-n χ² statistic
  below behaves differently — see its scale analysis.)
- **Roughness-interpolation map and the feasibility profile.** Writing
  y = x^θ and u = 1/θ, the prime share among y-rough survivors is
  1/(u·ω(u)) (Buchstab ω; Theorem 1 is the u = 3 section, since
  ω(3) = (1+ln 2)/3). Against the linear-sieve lower bound f(1/(2θ))
  (level x^{1/2}), the win condition "upper-bound the composite survivors
  within a factor U_req(θ) = f(1/(2θ))·(1 − 1/(u·ω(u)))⁻¹ of truth" has a
  feasibility profile peaking at **θ ≈ 1/8 with U_req ≈ 1.26**; θ ≥ 1/4 is
  an impossibility region (U_req < 1). The parity barrier, in these
  coordinates, is a 26%-margin upper-bound problem at roughness x^{1/8}.
- **Fiber decomposition of the enemy.** Splitting composite survivors at
  θ = 1/3 into fibers F_q (semiprimes divisible by q): a rank-1 common mode
  (the computable singular-series structure) carries 93.8% of fiber
  variance; after removing it, fibers exchange budget (residual total
  variance = 0.44 of the independent-fiber sum) and each fiber is
  sub-Poisson (Fano ≈ 0.59). The enemy total fluctuates at about half the
  σ of independent-Poisson bookkeeping — a first quantification of what
  fiber-independent sieve bounds discard.
- **Exhaustive χ² envelope (strongest individual-n datum).** Define
  χ²(n) = Σ_q (F_q − E_q)²/E_q against the universal fiber profile
  (self-normalized within n). Over all 30,000 candidates at 10⁶:
  χ²/(dof) = 0.515 ± 0.062, **worst case 0.796 — every single n stays
  below the Poisson value 1**, with the deepest excursion (4.55σ) matching
  the Gaussian expectation exactly. At 10⁷ (3000 samples) the typical value
  is 0.539 and the margin to 1 grows from 7.8σ to 10.7σ. Unlike variance
  statements, χ²(n) is defined within a single n, so this is empirically
  the exact *shape* of statement a proof would need. Scale behaviour, on
  standardized fiber ranges (y, 10y] and cube-gap windows, is a slow
  near-linear rise in ln x: typical values 0.515 → 0.530 → 0.546 → 0.560
  at 10⁶/10⁷/10⁸/10⁹ (slope ≈ 0.0068 per e-fold). A pre-registered
  three-model test at 10⁹ measured 0.5603 ± 0.0014: the
  Poisson-convergence path 1 − c/ln x (predicting 0.597) is rejected at
  26σ, a saturating-limit fit (0.550) at 7σ, while the linear-in-ln x
  description hits within 1.2σ. Extrapolating the linear rate, Poisson
  level would not be reached before x ~ 10³⁷; sampled worst cases fall
  with scale (0.80 → 0.73 → 0.69 → 0.65). The within-n fiber-uniformity
  advantage thus persists at every practically accessible scale — exactly
  where the cross-n variance advantage dies out — though its ultimate
  asymptote (saturation below 1 versus astronomically slow approach to 1)
  remains open.
- **Exact structure law at the target roughness (θ = 1/8 zone).** Using a
  segmented-sieve scanner (no full sieve array; feasible up to ~10¹²), the
  prime share s(n) among n^{1/8}-rough survivors at n ~ 10¹⁰ obeys, to
  measurement precision (3–4 decimals in a controlled per-q experiment),
  **s(n) = s₀(x) · ∏_{q | n, q > y} (q−1)/(q−2)** — the twin-type boost
  from each unsieved prime factor of n, while sieved prime factors cancel
  exactly (ratio 1.000; an accidental A/B test with the boundary prime
  confirmed both regimes). After this deterministic factor, residual
  fluctuation returns to binomial level. Strategically this tames the
  target zone: the boost only raises s, so the worst case for s(n) > 0 is
  the smooth class (no medium factors, factor product = 1), whose
  fluctuation is binomial-tame — the individual-n problem at θ = 1/8
  reduces to controlling s₀(x) for one structureless class. The law holds
  unchanged at y = 19/23/31 (x = 10¹⁰/10¹¹/10¹²; e.g. q = 101 ratio
  1.0099 vs predicted 1.0101 at 10¹²), fully y-smooth extremes
  (2^a·5^c·…·19^h) show no downward deviation, and s₀ marches
  0.2544 → 0.2421 → 0.2379 toward the Buchstab value 1/(8ω(8)) = 0.2226.
  Fiber uniformity also transfers to the target zone: χ²/dof =
  0.485/0.528/0.572 at 10¹⁰/10¹¹/10¹² (worst cases all below 0.78) —
  sub-Poisson at θ = 1/8 up to 10¹², rising only at log-crawl pace. The
  three ingredients of the program — the 26% feasibility margin, the exact
  structure law with its identified worst class, and sub-Poisson fiber
  uniformity — now coexist at the same roughness.
- **Anatomy of the chain: the wall localized.** Decomposing the composite
  survivors by smallest prime factor spf = x^α, the upper-bound loss is
  the integral ∫F((1/2−α)/α)·w(α)dα (linear-sieve F, Buchstab weight w),
  ≈ 4.08 classically at u = 8 versus the required 1.26. The mid-zone
  α ∈ [1/8, 1/3] (mass 0.65) has surviving sieve levels and is where the
  measured fiber uniformity lives — the plausibly classical part. The
  core α ∈ [1/3, 1/2] (mass 0.35) consists of semiprimes with both
  factors large (the cofactor is automatically prime); meeting the budget
  there requires a **binary-count upper constant below 2** (best known
  ≈ 3.9), itself a parity-breaking statement, and this requirement is
  θ-invariant. The program thus condenses the entire parity content into
  one classical open problem, with the fluctuation terms measured to be
  negligible (~2% of the margin) and every structural ingredient tamed.
- **Open route (three-link chain).** (i) Prove "χ²(n) ≤ K·dof for all n
  with some K < 1" — anchored near Barban–Davenport–Halberstam /
  Montgomery–Hooley mean-square theory, with the fixed-residue (single
  a = n mod q per modulus) obstacle; (ii) keep the switching loss (dropping
  the cofactor-primality condition) under the 26% margin at θ = 1/8;
  (iii) assemble via Cauchy–Schwarz into an individual-n bound on the
  composite-survivor total.

## Addendum B: first-pass related literature

Averages of Goldbach representations and their connection to zeta zeros are
classical and active: Fujii's formula and refinements
(arXiv:1601.06902, arXiv:1712.00737), oscillation sizes for R(n) − n𝔖(n)
(arXiv:2006.14742), Goldbach in arithmetic progressions vs zeros of Dirichlet
L-functions (arXiv:1704.06103), RH-conditional second-moment results
(Languasco–Perelli and others). Our zero-detection experiments (§ logs)
reproduce this known theory numerically. We did not find, in a first-pass
search, a prior statement of the normalized dispersion constant of Conjecture 1
or the gap-permutation (beyond-pair) observation; pointers welcome.

## Acknowledgment / request

If you are an analytic number theorist: we would be grateful for (i) pointers to prior work containing this constant or refuting it, (ii) an opinion on Conjecture 2, and (iii) whether the two remaining computations of the proof program are tractable under standard conjectures (RH, pair correlation, moment hypotheses for S(α)).
