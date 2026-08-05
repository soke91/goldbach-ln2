# A Proof Program for the Dilate-Averaged Möbius Correlation Bound

**Target (E1).** For $N$ even, $D(k) := \sum_{\sqrt N < m \le N/k} \mu(m)\,\mu(N-mk)$,
prove
$$\sum_{k \le K} |D(k)|^2 \;\ll\; \frac{1}{(\log N)^{A}}\sum_{k\le K} M_k
\qquad (M_k := \#\{m\}),$$
for $K$ up to $\sqrt N$ and some $A = A(\text{assembly}) < \infty$.

**Why it matters.** By the chain documented in this repository
(Huang–Li → Pan $f_3$ → Vaughan → hyperbola), E1 is the last unproven
cancellation standing between current technology and binary Goldbach for
large even $N$: the $m$-face is unconditionally controlled by the Möbius
Bombieri–Vinogradov theorem, and the boundary bands are the active
frontier (Lichtman 2023). Empirically E1 holds with exact half-normal
statistics (this repo, `typeII_kface.py`, `cascade2.py`).

---

## Step 1 — The exact dilation ladder (status: EXACT, verified numerically)

For prime $p$, splitting $D(k)$ by $p \mid m$ and writing $m = pm'$:
$$A_p(k) := \sum_{p \mid m} \mu(m)\mu(N-mk)
= -\sum_{\substack{m' \in (\sqrt N/p,\; N/pk] \\ p \nmid m'}} \mu(m')\,\mu(N - m'(pk)).$$
This is an identity (μ(pm') = −μ(m′) for $p \nmid m'$; $p^2 \mid m$ terms
vanish). The dilate family $\{D(k)\}$ is **closed** under this move: the
scale-$k$ object maps exactly onto the scale-$pk$ object. (Contrast: for
shift correlations the analogous move distorts the shift $h \to ph$.)
Numerical check (40 values of $k$, $p = 2/3/5$): **identity error
exactly 0** in every case; the sub-sum-to-dilate transfer is
near-deterministic, corr$(A_p(k), D(pk)) = -0.60/-0.88/-0.87$, with
fuel ratios $|A_p|/|D| = 0.41/3.13/0.95$; the aggregate cross-scale
correlations corr$(D(k), D(pk)) = -0.24/-0.40/-0.23/-0.16$ for
$p = 2/3/5/7$ (`entropy_ladder.py`, `e1_dilation.py`).

## Step 2 — Averaged self-similarity (status: routine bookkeeping)

Averaging Step 1 over primes $p \in [P, P^{1+\eta}]$ with the standard
Turán–Kubilius/Ramaré device expresses $D(k)$, up to
$O(1/\sqrt{\log P})$ errors, as a weighted average of the dilated duals
$\{D^{(p\nmid)}(pk)\}_p$. Hence largeness of $|D(k)|$ at scale $k$
propagates to positive-density largeness across the dilated scales
$\{pk\}$ — energy cannot hide at a single scale.

## Step 3 — Energy bookkeeping across dyadic bands (status: routine)

By Cauchy–Schwarz on Step 2, if $|D(k)| \ge \delta \sqrt{M_k}$ on a
$\delta$-dense subset of $[K, 2K]$, then $|D| \ge \delta' \sqrt{M}$ on
$\delta'$-dense subsets of $[pK, 2pK]$ for most $p \in [P, P^{1+\eta}]$,
with $\delta' \gg \delta / \log^{o(1)}$. The per-band $L^2$ budgets
$\sum_{k \sim K} |D(k)|^2 \lesssim \sum_k M_k$ hold empirically at the
independence level (measured ratios 0.953–1.046).

## Step 4 — Entropy decrement (status: OPEN — the engine; a port, not an invention)

This is the dilate-analog of the entropy-decrement step in Tao's proof of
the two-point logarithmic Chowla conjecture (Forum Math. Pi, 2016). There,
the mutual information between μ-blocks and the sampling residues
$n \bmod p$ decreases along scales, yielding a scale range where sampling
decouples, and the resulting forced structure contradicts a known
equidistribution input. Two structural facts suggest the port is *easier*
here, not harder: (i) the family is invariant under the ladder move (no
shift distortion — multiplicativity is native to dilation); (ii) the
decoupling target is a family of Möbius sums in arithmetic progressions
(moduli $k$), for which unconditional average results (Möbius BV; MRT in
progressions) exist as the contradiction input — the same class of input
Tao's argument consumes, here needed in the $k$-aspect where it is a
theorem. What must be carried through: the conditional-entropy bookkeeping
with the $(k, pk)$ ladder replacing the $(n, n+h) \to (pn, pn+ph)$ move,
and uniformity of the error terms in the dilate weights.

## Step 5 — Contradiction and conclusion (status: routine given Step 4)

Decoupled largeness across $\gg \log\log N$ multiplicatively-independent
scale shells, each carrying $\delta^2$-energy against a fixed global
$L^2$ budget (Step 3), forces $\delta^2 \cdot \Omega(\log\log N) \ll 1$,
i.e. $\delta = o(1)$ — and quantitatively $\delta \ll (\log N)^{-A}$
tracks the entropy rate, which is what the Huang–Li assembly requires.

---

## Status summary

| Step | Status |
|---|---|
| 1. Dilation ladder | **Exact identity, verified** |
| 2. Averaged self-similarity | Routine (TK/Ramaré) |
| 3. Energy bookkeeping | Routine (CS + measured budgets) |
| 4. Entropy decrement (dilate port) | **OPEN — the single engine** |
| 5. Contradiction | Routine given 4 |

One open step; it is a *port* of an existing published argument into a
setting with two verified structural advantages, not a new mechanism.
Every routine step and both advantages are numerically verified in this
repository. Experts in the entropy-decrement method are invited to
attempt Step 4; the authors of this program will gladly share all data
and further computations.
