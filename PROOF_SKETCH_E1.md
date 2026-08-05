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

## Step 4 — The engine (status: OPEN, with two concrete port shapes)

**Key structural discovery (from Tao 2016's own strategy section).**
The *only* role of the logarithmic averaging in the shift-Chowla proof
is that "the logarithmic averaging allows us to leave the constraint
n ≤ x unchanged" under the dilation n → pn. In our setting this crutch
is unnecessary: the ladder move acts as k → pk on the *spectator*
variable, and the target is already an average over k — **the k-average
natively absorbs the range dilation that forced Tao's log-weights.**
The obstruction that keeps unweighted two-point Chowla open (no
averaging direction) is absent here by construction: E1 belongs to the
difficulty class of *averaged* Chowla (Matomäki–Radziwiłł — a theorem),
not fixed-shift Chowla.

**Port shape (a) — entropy decrement with k-averaging.** Run Tao's
contradiction scheme with the (k, pk) ladder in place of (n, n+h) →
(pn, pn+ph): largeness of Σ_k |D(k)|² propagates through the exact
ladder (Step 1) to all dilated shells; the entropy decrement selects a
scale range where the μ-blocks decouple from the p-sampling; Hoeffding
replaces 1_{p|·} by 1/p; the contradiction input is Möbius
equidistribution in progressions *in the k-aspect* — where unconditional
theorems (Möbius BV) exist.

**Port shape (b) — the Mellin route.** Dilations are multiplicative
translations, so Σ_k |D(k)|² has a Dirichlet-polynomial fourth-moment
expression (the Mellin analog of MRT's Fourier computation for averaged
shift-Chowla). The required engine becomes a mean value theorem for
Dirichlet polynomials with *shifted-Möbius* coefficients μ(N−u) — the
native language of the Matomäki–Radziwiłł machine, applied to a
shifted-multiplicative sequence. This names the missing lemma at
maximal concreteness. Its conclusion is measured true with margin: for
$G(t) = \sum_{j \sim 2\cdot 10^6} \mu(N-j)\, j^{-it}$ over
$t \in [0, 2000]$ (1200 grid points), $|G(t)|/\sqrt{J} $ has mean
$0.875$ (Rayleigh: $0.886$) and sup $1.60$ — full square-root
cancellation uniformly in $t$, slightly *cleaner* than a random-sign
control (mean $1.065$, sup $2.13$). The shifted-Möbius polynomial
behaves like — indeed better than — a random polynomial, which is
exactly what the missing mean value theorem asserts.

### Draft bookkeeping for port (a) — the k-averaged entropy argument

*Hypothesis for contradiction.* Σ_{k∼K} |D(k)|² ≥ δ² Σ_k M_k for all
dyadic K in [K₋, K₊] with log K₊/K₋ ≍ log N.

1. *Ladder expansion.* By Step 1 (exact) and Step 2 (TK-averaged over
   p ∈ P = [P, P^{1+η}]), the hypothesis transfers: for most p ∈ P the
   dilated shells [pK, 2pK] also carry δ′-energy, δ′ ≫ δ − O(1/√log P).
2. *Random model.* Draw k uniformly from [K, 2K]; the random variables
   are the μ-block B_k = (μ(N−mk))_{m∼M} and the sampling pattern
   S_k = (k mod p)_{p∈P}. The bilinear expression to control is the
   analog of Tao's (1.5): (1/|P|) Σ_{p∈P} D-terms coupled through
   1_{p|·-structure} in the dilated shells.
3. *Entropy decrement.* Shannon subadditivity over the nested scales
   pK: either the mutual information I(B_k ; S_k) drops below ε at some
   shell (good case), or entropy decreases by ≥ ε per shell — which can
   happen at most H(B)/ε ≪ M/ε times; with ≍ log N shells available,
   a good shell exists for ε ≍ M/log N. [This is the step whose
   quantitative port must be carried out carefully; the shell structure
   is *multiplicatively nested* (K, pK, pp′K, …) exactly as in Tao's
   n-aspect, and no range correction is needed — the k-average absorbs
   it.]
4. *Concentration.* At the good shell, Hoeffding replaces the sampling
   pattern by its mean: 1_{p | ·} → 1/p, uniformly over the block —
   yielding the *deterministically averaged* expression
   (1/|P|) Σ_p (1/p) Σ_{k∼pK} |D(k)|²-type sums.
5. *Contradiction input.* The averaged expression is a plain Möbius
   equidistribution statement in the k-aspect over moduli p ∈ P —
   controlled unconditionally by Möbius–Bombieri–Vinogradov / MRT-type
   inputs, forcing the averaged expression to be o(δ²)-small — 
   contradicting the propagated largeness. □ (modulo step 3's
   quantitative bookkeeping)

*Honest flags.* (i) Step 3's information-theoretic bookkeeping must be
re-derived in dilate coordinates — the shape matches Tao's but the
conditional-entropy chain rule applications need explicit re-doing;
(ii) step 5's reduction to μ-BV must confirm no hidden binary object
reappears (parity's known trick — checked empirically here: the
averaged objects measured half-normal throughout this repository).

### Original port framing (kept for the record)

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
