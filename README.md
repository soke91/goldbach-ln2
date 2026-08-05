# goldbach-ln2

**An empirical program on the Goldbach conjecture: from a sub-Poissonian
fluctuation constant (≈ ln 2) to a complete empirical map of the
Huang–Li conditional chain — including a direct measurement of the final
axiom (EH_μ) at its own coordinate.**

> Status: **numerical / exploratory — measurements and conjectures, not
> theorems.** Produced in AI-assisted exploratory sessions (2026-08-03 → 05)
> with 24 documented self-corrections and three pre-registered blind
> extrapolation tests (all hits). Expert scrutiny is explicitly welcome —
> if anything here is known, wrong, or provable, please open an issue.

## The program in one picture

By a known conditional theorem (Huang–Li 2022, continuing Pan 1982;
arXiv:2005.03811), Goldbach for large even N follows from
Bombieri–Vinogradov **plus one sentence**:

> **EH_μ(N^{θ′}) for some θ′ > 1/2**: the correlation sequence
> c(n) = Λ(n)·μ(N−n) equidistributes in arithmetic progressions with
> moduli past the √N barrier — needed only for the fixed residue class
> N mod q.

This repository measures that landscape (and everything feeding it):

| Result | Measurement |
|---|---|
| **Final axiom, directly** | Fixed-residue discrepancy of c(n): mean 0.6–1.0× random walk through θ = 0.30 → **0.70** (3 values of N at 10⁸) — **no visible change at the √N barrier** |
| Möbius landscape | Discrepancy = (0.27 ln Q) × random walk to θ = 1/2; flat 3.2–3.8×√(x/q) deep into the unproven zone θ ≤ 0.88 |
| **Conjecture P** (sifted primes) | Sifted-prime counts sit on the Buchstab curve e^γω(u) to ±0.0002 through the sieve-blind zone, fine structure reproducing across a decade |
| Minimal target (P_loc) | Two aggregates per n (S, P2) to ~60% accuracy suffice for g > 0 at 10⁸; nature delivers 0.48% worst-case (**factor-125 slack** over 3000 n) while classical lower bounds are blind |
| Structure law at θ = 1/8 | s(n) = s₀(x)·∏_{q|n, q>y}(q−1)/(q−2) to 3–4 decimals at three scales; worst class identified (smooth n), no downward anomalies |
| χ² fiber uniformity | Sub-Poisson (0.49–0.57) from 10⁶ to 10¹²; Poisson-convergence rejected at 26σ in a pre-registered test |
| The ln 2 constant | Canonical dispersion ratio of g/HL → 0.6931 ± 0.0061 ≡ ln 2 (18 octaves, 3 independent constructions); a ≥3-point additive-correlation invariant |
| Loss anatomy | Classical upper-bound loss integral 4.08 vs required 1.26 at u = 8; all parity content condenses into the P2 core / small-s sifted asymptotics |

Full write-up: [RESEARCH_NOTE.md](RESEARCH_NOTE.md) (self-contained, with
methods, corrections, and honest caveats). Current single-page state:
[STATUS.md](STATUS.md). The proof program for the final bound:
[PROOF_SKETCH_E1.md](PROOF_SKETCH_E1.md) →
[paper/e1_proof.tex](paper/e1_proof.tex) (ten certified/adjudicated
rows; the sole remainder is the **SEAM Conjecture** — family-averaged
thin-progression Möbius cancellation at pair-moduli x^0.6–x^(2/3),
stated formally in the paper with its full resistance profile and
three closing routes).

## What is claimed / not claimed

- **Claimed (facts)**: the computed values and verifications; the exact
  identities used (Theorem-1 dichotomy, Buchstab bookkeeping); the
  documented reductions.
- **Conjectured**: the ln 2 constant; Conjecture P; the empirical structure
  laws.
- **Not claimed**: any theorem toward Goldbach. The conditional frame
  (Goldbach ⟸ BV + EH_μ) is Huang–Li's theorem, not ours; our contribution
  on that axis is measurement.

## Reproduce (Python + numpy; laptop-scale)

```bash
pip install -r requirements.txt
cd code

# Final axiom landscape (c(n) = Lambda*mu(N-n) fixed-residue discrepancy)
python ehmu_final.py

# Mobius-in-APs landscape to theta = 0.88
python ehmu_probe2.py && python ehmu_beyond.py

# Conjecture P profile (sifted primes vs Buchstab through the blind zone)
python p_profile.py

# P_loc slack ledger (two aggregates per n vs the winning threshold)
python ploc_scan.py

# Structure law at theta = 1/8 (segmented sieve, controlled per-q experiment)
python s_th18_qlaw.py

# chi2 fiber uniformity ladder
python s1_chi2_envelope.py

# The original ln 2 dispersion series
python e4_dense_sample.py

# One-shot reproducibility suite (mu sanity, ladder identity, dispersion,
# seam band — all four stamps in ~15 min)
python verify_all.py
```

Each script prints its reported numbers to stdout. `code/goldbach/` is the
shared sieve library. ~40 scripts document every measurement in the note.

## Falsification

Each claim carries its own kill condition in the note: fitted limits that
drift, structure that fails at the next scale, or a classical account of any
constant. Three campaigns in this program died exactly this way (R² → ln 2,
thin-class habitat, Λ²-family headroom) and are documented as corrections —
the surviving claims are the ones that passed.

## License

Code: MIT. Text: CC BY 4.0.
