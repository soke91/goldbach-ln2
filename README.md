# goldbach-ln2

**An empirical sub-Poissonian constant for Goldbach representation fluctuations,
with numerical evidence that it equals ln 2.**

> Status: **numerical / exploratory — a conjecture, not a theorem.**
> Produced in an AI-assisted exploratory session (2026-08-03), with 15 documented
> self-corrections. Expert scrutiny and literature pointers are explicitly welcome —
> if this constant is already known (or wrong), please open an issue.

## TL;DR

Let g(n) be the number of unordered prime pairs with p+q = n, and HL(n) the
first-order Hardy–Littlewood prediction. Measuring the dispersion of
z = g/HL with a fixed canonical estimator (dyadic windows, n ≡ 2 mod 6,
log-drift removal, Poisson-scale normalization), across 18 dyadic windows
up to 5.4×10⁸:

**fitted limit of the dispersion ratio = 0.6931 ± 0.0061 ≈ ln 2 = 0.69315**

Confirmed by two independent constructions (Hardy–Littlewood-structured null
worlds: variance ratios 0.468 ± 0.010 and 0.472 ± 0.010 vs (ln 2)² = 0.4805)
and invariant under window-shape changes (ratio 1.3 → 6).

Additional finding: permuting prime gaps (which preserves all two-point
statistics) destroys the suppression entirely — evidence the constant is an
invariant of **additive prime correlations of order ≥ 3**, beyond
pair-correlation/GUE phenomenology.

Full write-up: [RESEARCH_NOTE.md](RESEARCH_NOTE.md).
Proof roadmap + 15-pitfall obstruction map: [docs/PROOF_PROGRAM.md](docs/PROOF_PROGRAM.md).
Complete session logs: [docs/](docs/).

## Reproduce (Python + numpy only; minutes on a laptop)

```bash
pip install -r requirements.txt
cd code

# Main dispersion series (high octaves via exact per-point sampling; ~10-30 min)
python e4_dense_sample.py

# Window-shape invariance test
python e4_window_test.py

# Hardy–Littlewood null world comparison / wheel-modulus scan
python e4_wheelnull.py
python e4_zscan3.py

# Gap-permutation surrogates (beyond-pair evidence)
python e4_gapshuffle.py

# Spectral profile phi(lambda) ~ Goldston–Montgomery function
python e4_phi_profile.py

# kappa(d) conditional singular-series verification
python e4_breakthrough.py

# Bonus: detect Riemann & Dirichlet L-zeros inside Goldbach data
python q3p_spectroscopy.py
python q3p_L.py
```

Scripts print all reported numbers to stdout. `code/goldbach/` is the small
shared library (sieve + Hardy–Littlewood utilities). Some scripts use paths or
constants tuned during the session; see docs/E4_PROGRESS.md for the exact
provenance of every reported figure.

## What is claimed / not claimed

- Claimed (facts): the computed values, verifications (Goldbach to 10¹¹ here;
  known literature: 4×10¹⁸), and elementary lemmas (κ(d) formula etc.).
- Conjectured: the canonical suppression constant equals ln 2; it is a ≥3-point
  additive invariant.
- **Not claimed**: any progress on the Goldbach conjecture itself; any theorem
  about the constant.

## Falsification

The claim dies if: (a) the fitted limit departs from ln 2 with larger
computations, or (b) an estimator-theoretic account reproduces the constant
from pair correlations alone, or (c) the constant is shown to be
frame-artifactual. Section 7 of the note lists caveats honestly (estimator
dependence, extrapolation risk, possible literature overlap).

## License

Code: MIT (or your preferred license — set before publishing). Text: CC BY 4.0.
