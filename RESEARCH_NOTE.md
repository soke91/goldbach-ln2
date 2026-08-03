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

## Addendum A: the parity-barrier ratio equals 1 : ln 2 (derived & verified)

The measured ≈60:40 prime:P2 split among n^{1/3}-sieve survivors (§ our earlier
logs) has a closed form. For m ~ x with no prime factor ≤ x^{1/3}, the smaller
factor of a P2 is q = x^u with u ∈ (1/3, 1/2], and the density integral is
∫_{1/3}^{1/2} du/[u(1−u)] = ln 2. Hence asymptotically

**prime share = 1/(1+ln 2) = 0.59064…, P2 share = ln 2/(1+ln 2) = 0.40936…**

Verified: unconditioned random odd m near 10⁸ give 0.5967 (finite-size);
Goldbach partners m = n−p give 0.6054, the excess (+0.0087) being the expected
pair-correlation boost. The derivation is classical sieve density computation
(surely known to specialists); the observation we record is the identification
of the measured parity-barrier ratio with 1 : ln 2 — and the (open) question of
whether this ln 2 (which stems from the canonical 1/3 dichotomy exponent) is
related to the fluctuation constant of Conjecture 1.

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
