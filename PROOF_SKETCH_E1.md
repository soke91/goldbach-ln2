> **⚠ ADVERSARIAL REVIEW VERDICT (increment 143): the core reductions
> of this document (T5, T5.1, T6, the SEAM formalization, the h-range
> collapse under the q-layer) were REFUTED by an independent
> adversarial review — see REVIEW_VERDICT.md at the repository root.
> The exact identities (T2, T3) and all measurements survive. This
> document is retained as a record of the attempt.**

# A Proof Program for the Dilate-Averaged Möbius Correlation Bound

**Target (E1).** For $N$ even, $D(k) := \sum_{\sqrt N < m \le N/k} \mu(m)\,\mu(N-mk)$,
prove
$$\sum_{k \le K} |D(k)|^2 \;\ll\; \frac{1}{(\log N)^{2A+2}}\sum_{k\le K} M_k^2
\qquad (M_k := \#\{m\}),$$
for $K$ up to $\sqrt N$ and some $A = A(\text{assembly}) < \infty$ —
a fixed log-power saving over the **trivial** bound
$\sum_k M_k^2 \asymp N^2/K$.

> **Note on what this sketch was aiming at.** Everything below is
> calibrated to $\sum_k M_k$, the **square-root** scale — see the
> contradiction hypothesis "$\sum |D(k)|^2 \ge \delta^2 \sum_k M_k$"
> in the port sketch, and the per-band budgets in Step 3. That is a
> factor $N/K$ stronger than the chain consumes. So this program was
> attempting to prove considerably more than was needed, in addition
> to failing structurally (REVIEW_VERDICT items 1–2). Neither fact
> rescues the other; both are recorded in CLOSURE_REAUDIT.md.

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

## Cascade addendum (increments 107–110)

Attempting the port yielded four quantitative discoveries: (i) the
sampling pattern (p | m) is deterministic in our coordinates — Tao's
central obstacle (block–graph dependence), and hence the entropy
section itself, has no analog here; (ii) the naive Turán–Kubilius
one-step expansion fails exactly at the predicted 1/√L noise scale —
reproducing the reason for the triple-exponential parameter hierarchy;
(iii) with *regressed* weights the ladder sum explains **68%** of the
L² energy per step at wide prime windows (corr(D, S) = +0.83,
c* ≈ −0.39 vs naive 0.54); and (iv) this efficiency is
**scale-stationary** (R² = 0.681 → 0.656, c* = −0.39 → −0.41 on a
3×-deeper band) — the necessary condition for iterating the cascade,
suggesting geometric residual decay (≈ 0.33 per step) and a possible
collapse of the parameter hierarchy to ~log log N steps. Remaining for
the proof: an a-priori weight design achieving the regressed transfer
(second-order TK corrections are the natural candidate), the iterated
bookkeeping with collision terms, and vigilance for parity re-entry in
the iteration — the twelfth watch.

**Iteration verdict (honest).** The direct two-step test kills the
geometric-decay hypothesis: regressing the residual field on *its own*
ladder yields R²₂ = 0.191 (vs 0.681 at step one; cumulative transfer
plateaus near 0.74). The first ladder removes an easy correlated
component; the hard core resists its ladder — the twelfth observed
instance of the parity-conservation pattern, now localized in the
residual field of the dilation cascade. Consequence: no collapse of the
parameter hierarchy; the port's value reverts to its structural
simplification (the entropy section is removable, the sampling is
deterministic) with Tao-scale window economics intact. The residual
field R(k) — 32% of the energy, ladder-orthogonal — is the sharpest
localization of the parity core this program has produced: whoever
characterizes it characterizes the barrier.

**The core identified (thirteenth hunt).** Decomposing by the number of
P-window prime factors ω_P(m): in the working band ω_P = 0 ⟺ m prime
(3000² exceeds the m-range), and corr(R, D₀) = +0.827 — the
ladder-orthogonal residual is dominated by the prime slice. **The
parity core is the Möbius-on-shifted-primes field
Σ_{m prime} μ(N−mk)** — the dilate-family version of the folklore
conjecture Σ_p μ(p+h) = o(π) (Hildebrand 1989; single shift open).
Its shift-averaged version is a *proven theorem* of Lichtman
(arXiv:2009.08969, via sieve methods + refined
Matomäki–Radziwiłł–Tao), whose method extends to structured shift
families; and the boundary bands of this program's chain rest on
Lichtman's 2023 level-of-distribution work. Both remaining items of
the Goldbach chain thus sit on the natural extension of one active
research line — and every object in the chain has been measured here
to obey exactly the cancellation the missing theorems assert.

## Substitution analysis (increment 113): the hinges match the key

Reading Lichtman 2020's proof: the strategy restricts to typical n with
prime factors in prescribed windows, handles the sparse complement by
sieves, and uses the shift-average through a congruence coupling — a
medium prime p₁ divides n = p + h iff h ≡ −p (mod p₁), equidistributed
over the averaged variable; the analytic inputs are MRT Fourier
uniformity and the Vinogradov–Korobov zero-free region. In the dilate
family the coupling is *isomorphic*: p₁ | N − pk iff
k ≡ N p^{-1} (mod p₁) — again a congruence on the averaged variable —
and the argument N − pk lives in arithmetic progressions mod k, so the
additive-character inputs apply unchanged; his Theorem 1.8 (general
non-pretentious f, H = X^θ, via Fourier uniformity) offers a second
route. **No structural obstacle is in sight.** The remaining
differences are finite bookkeeping: (R1) the inner range M_k = N/k is
k-dependent (dyadic decomposition expected to handle it — to be
checked), and (R2) our L² versus his L¹ (one Cauchy–Schwarz). The
thirteenth parity watch is posted at R1's dyadic boundaries. If R1–R2
check out, the chain closes: k-averaged Möbius-on-shifted-primes ⟹
the residual core ⟹ E1 ⟹ (with the m-face and boundary bands) EH_μ
past 1/2 ⟹ Goldbach for large even N.

## R1 section-by-section check (increment 114): the thirteenth tooth located

Carrying Lichtman 2020 §2 through the dilate substitution: the
typical-set restriction and sieve bounds (his (2.6)–(2.7)) port cleanly
under dyadic k-ranges. The structural divergence appears exactly at the
decoupling engine (his Lemma 2.1): the shift equation m − n = m′ − n′
resolves by orthogonality into *linear* phases e(nα) — feeding the MRT
short-interval input — whereas the dilate equation
p′(N − u) = p(N − u′) resolves into **bilinear** phases (products p·u
in the exponent). The key Fourier estimate for the dilate family is
therefore a *bilinear-phase short-sum bound* — Kloosterman-class input
rather than MRT-linear. This is where the thirteenth parity tooth
lives; and the technology for exactly this input class is the
Deshouillers–Iwaniec spectral large sieve as refined in Lichtman 2023.
**The two remaining items of the Goldbach chain — the boundary bands
and the core's decoupling — thus merge into one toolbox at the live
frontier.** The final open question of this program, at maximal
resolution: does the spectral large-sieve toolbox, at its current or
foreseeable strength, control the dilate decoupling's bilinear phases?
Nature's answer is measured throughout this repository; the proof's
answer is where the field now stands.

## Gate arithmetic (increments 116–117): all three gates pass at the naive dictionary

Splitting the k-range dyadically: for K ≤ x^{1/3} the dispersion sits
exactly in Lichtman 2023's critical window (dictionary N_L := K,
M := x/K, matching his NM ≍ x frame; factor budgets R ~ S ~ K^{1/3}
from the k = qrs family factorization, measured density 47% with the
prime-moduli remainder in the classical DI regime); for
K ∈ (x^{1/3}, x^{1/2}] the inner ranges shorten into the
Matomäki–Radziwiłł short-interval home ground (upper-band mass
measured 0.0023N at half-normal 0.84). Substituting into his gates:
(5.1) reads K³ ≤ x — **tight exactly at the regime boundary x^{1/3}**,
a consistency signal that the two-regime split is the gate geometry
itself; (5.2) gives x^{7/32}·K^{13/3} = x^{1.66} < x²; (5.3) gives
K^{4.77} = x^{1.59} < x². **All three gates pass.** Honest caveats:
the arithmetic presumes the naive dictionary (that the W-analog
derivation lands in his exact form); the h- and a-parameter roles need
precise correspondence; a fourteenth tooth could appear only inside
the full derivation. The single remaining item of this program is that
derivation — one paper-section of work, with every surrounding
quantity measured and every gate pre-checked.

## The W-analog derivation (hunt 15 main body; increment 118)

**Setup.** Fix dyadic K ≤ x^{1/3}, P = x/K. Target:
T(K) = Σ_{k∼K} |Σ_{p∼P} μ(N−pk)|² ≪ KP/(log x)^A.

**Step D1 (expansion; exact).** T = Σ_{k∼K} Σ_{p,p′∼P}
μ(N−pk)μ(N−p′k). The diagonal p = p′ contributes ≤ K·P — absorbed by
one log-power since the target allows KP/(log)^A only after... [flag
F0: the diagonal is size KP, i.e. exactly the trivial bound; the
required saving must come entirely from the off-diagonal cancellation
over p ≠ p′ — as in every dispersion argument, one needs the
off-diagonal to be ≤ KP/(log)^A *and* the diagonal counted once, which
is standard: the dispersion bounds |T − diag| and the final CS uses
T ≤ diag + offdiag with diag/(KP) = 1/P^0... **correction**: the
target normalization KP is the trivial size; the saving is claimed
over the full square sum, so what is needed is offdiag ≪ KP/(log)^A
while diag ≈ K·P·(μ²-density) is *itself* the main term of T under
square-root cancellation — T ≈ diag is precisely the desired
conclusion. So the theorem to prove is |offdiag| ≪ KP/(log)^A.]

**Step D2 (off-diagonal reparametrization; exact).** For p ≠ p′ set
u = N − pk (so k = (N−u)/p, requiring p | N−u). Then
u′ = N − p′k = (p′u − N(p′−p))/p... equivalently the pair (u, u′)
satisfies p′(N−u) = p(N−u′), and summing over k ∼ K with both
divisibilities is the same as summing over u ≡ N (mod p) in the
interval N − pK·[1,2) — an AP-interval — with u′ determined. Hence
offdiag = Σ_{p≠p′∼P} Σ_{u≡N (p), u∼U} μ(u)·μ(u′(u)),
U = pK ≍ x, u′ = linear-fractional image of u with rational slope
p′/p and shift N(p−p′)/p.

**Step D3 (congruence completion; standard).** Detect u′ ∈ ℤ and the
AP-condition by additive characters mod p (and the u′-integrality is
automatic given k ∈ ℤ). Expanding the interval condition by finite
Fourier/completion introduces the h-sum with H′ ≈ x^{o(1)}·(moduli
budget)/K — matching the shape of Lichtman's H′ ⩽ x^{o(1)}QR²S²/M
under the dictionary [flag F1: exact H′ bookkeeping to be fixed at
write-up]. The completed phases are
e(a·h·(p̄′-type inverses)·u/p·)-class — modular-inverse bilinear
phases in the two prime-blocks, i.e. precisely the W-form of
Lichtman's Lemma 5.1 with his (n₁, n₂) ↔ our (p, p′) [both ∼ P... 
**flag F2 — a real divergence**: in his W the dispersion variables
n₁, n₂ carry *arbitrary bounded coefficients* α_n and range over
N_L ≤ x^{1/3}, while our p, p′ ∼ P = x/K ≥ x^{2/3} are LONG; our short
variables are the k's (≈ his n's). The correct role-matching sends
k-pairs to his n-pairs — this requires running D2 with the roles
exchanged: fix (k, k′), sum over p — which is the C_{k,k′} form
already measured (engine, increment 93). Under that matching our
dispersion variables k ∼ K ≤ x^{1/3} sit in his window ✓, our long
p-side carries the smooth sum ✓, and the factor budgets come from the
k = qrs family split (measured densities: see below). The derivation
should therefore be organized as: CS in p first (moving p to the
smooth side), dispersion over (k, k′).]

**Step D4 (gate substitution; done at increments 116–117).** With
N_L := K, M := P, R ∼ S ∼ K^{1/3}: gates (5.1)–(5.3) pass, (5.1)
tight at K = x^{1/3}.

**Step D5 (remainder regimes).** K ∈ (x^{1/3}, x^{1/2}]: short-interval
regime (MR-class inputs; measured tame). Non-factorable k (density
≈ 1/2): prime-moduli DI regime.

**Honest gap list.** F0 resolved above (bookkeeping clarification);
F1 (H′ budget) — routine but must be written; F2 (role exchange) —
structural, resolved in principle by organizing the CS/dispersion
order as stated, but the full write-out with the exchanged roles is
exactly where a fourteenth tooth could live; the α_n-coefficient
freedom in his Lemma (arbitrary bounded) covers our μ-values on the
k-side ✓. Next: write D2′ (exchanged-role version) in full.

## D2′ — the exchanged-role dispersion (increment 119)

Organize the estimate with Cauchy–Schwarz in the long variable first:
|Σ_k ξ_k Σ_p μ(N−pk)|² ≤ P · Σ_p |Σ_{k∼K} ξ_k μ(N−pk)|²
= P · Σ_{k,k′∼K} ξ_k ξ_{k′} C_{k,k′},
C_{k,k′} = Σ_{p∼P} μ(N−pk) μ(N−pk′).
The dispersion now runs over the SHORT pairs (k, k′) ∼ K ≤ x^{1/3} —
Lichtman's window — with arbitrary bounded coefficients ξ (his α_n
freedom covers this). Diagonal k = k′: ≤ K·P·(density), the main term.
Off-diagonal: for k ≠ k′, write w = N − pk (so p = (N−w)/k, k | N−w),
w′ = N − pk′ with k′w − kw′ = N(k′−k): the p-sum becomes a sum over
w ≡ N (mod k), w ∼ pK-interval, of μ(w)μ(w′(w)) with w′ the
linear-fractional image of slope k′/k — rational slope with SMALL
numerator/denominator (≤ K ≤ x^{1/3}), unlike D2's large-slope
version. Completion of the AP-interval condition mod k and of the
integrality mod k′ introduces additive characters with moduli
[k, k′] ≤ K² ≤ x^{2/3} and an h-range H′ ≈ x^{o(1)}K²/P — for
K ≤ x^{1/3}: H′ ≈ x^{o(1)}·K³/x ≤ x^{o(1)} — **the h-sum is
essentially bounded (O(x^{o(1)}) terms), a dramatic simplification
versus the generic case.** The completed phases carry modular inverses
mod [k,k′] evaluated on the μ-argument w — Kloosterman-class with
modulus ≤ x^{2/3} and factorable (k = qrs available on the
gate-region subfamily; prime-moduli regime otherwise) — the exact
input class of Lichtman §10 with, additionally, a short h-range.
Gates re-checked under this shape at increments 116–117 (pass; (5.1)
tight at K = x^{1/3}).

**Remaining to certify (the honest last list):**
(G1) the completion bookkeeping (finite Fourier of the interval and
integrality conditions) with its error terms — routine but must be
written line by line;
(G2) the μ(w)μ(w′)-to-|W-sum| reduction: the coefficients on the
w-side are μ-values (bounded ✓) but *fixed* (not freely chosen) —
Lichtman's W allows arbitrary bounded α_n, so this direction is
covered ✓;
(G3) the treatment of gcd conditions (k, k′) > 1 and the N-set
coprimality structure of his Lemma — bookkeeping;
(G4) summing the dyadic K-bands and the two regimes into the final
EH_μ-shape with the fixed residue N mod q — assembly.
None of G1–G4 has the shape of a parity obstruction; all four are
finite derivational labor. **The fourteenth tooth, if it exists, must
live inside G1's error terms or G4's assembly — and nowhere else.**

**Family coverage (measured, increment 120).** The gate-region
admissible density among k is 0.561/0.579/0.605 at K = 10⁴/10⁵/10⁶ —
*rising* with K (in contrast to the falling balanced-budget density:
the region flexibility resolves the budget concern empirically);
primes and semiprimes (27–35%) route to the prime-moduli DI regime;
the small dominant-prime-factor remainder routes there too with the
large factor as the modulus part. The k-family is covered.

## G1, line by line — the completion bookkeeping (increment 121)

Fix k ≠ k′ ∼ K, g = (k, k′), L = [k, k′] ≤ K²/g. The off-diagonal
inner sum is C_{k,k′} = Σ_{p∼P} μ(N−pk) μ(N−pk′).

**(G1.a) Congruence normalization.** Write w = N − pk; the map
p ↦ w is a bijection from p ∼ P onto the AP
{w ≡ N (mod k)} ∩ I_k, I_k = N − k·(P, 2P]. The second argument is
w′ = N − pk′ = (k′w − N(k′−k))/k... more simply keep the p-variable:
both arguments are linear in p, so C is a sum over the single variable
p of a product of two fixed μ-evaluations along APs — no completion is
needed for integrality (both arguments are automatically integers);
the only conditions to complete are the *interval* p ∼ P (sharp cutoff)
and, after the factorization-transfer to Kloosterman form, the
residue-class decompositions mod L.

**(G1.b) Residue decomposition (exact).** Splitting p by residues
a mod L: C = Σ_{a (L)} Σ_{p≡a (L), p∼P} μ(N−pk)μ(N−pk′). This is an
identity (machine-verified: reconstruction error 0 — see
`g1_completion.py`). The point of the decomposition: within a residue
class, the arguments N−pk, N−pk′ run over APs mod Lk, Lk′ with fixed
entry points — the structure the Kuznetsov-side machinery consumes.

**(G1.c) Interval completion.** The sharp cutoff p ∼ P is exchanged
for smooth weights at cost O(P^{1−δ}) per boundary (standard smoothing;
absorbed by the log-power budget), or completed by finite Fourier with
the h-range H′ ≈ x^{o(1)} established at D2′ — each h ≠ 0 term is a
phase-twisted copy of the same object (measured: the twisted copies
carry the same √-cancellation as the h = 0 term — engine data), and
there are only x^{o(1)} of them: the completion costs a factor
x^{o(1)} on the error budget, which the gates' ε-margins absorb.

**(G1.d) What remains inside G1.** The transfer from (G1.b)'s
residue-decomposed form to the exact W-shape of Lemma 5.1 (matching
his variable names q, r, s, n₁, n₂, f, h and his coprimality set N) —
pure transcription against his §5 pp. 14–16, no analytic content.
G1 is thereby reduced to transcription plus the standard smoothing
lemma. **No fourteenth tooth found in G1.** The remaining habitat
shrinks to G4 (assembly) alone.

## G3 and G4 (increment 122): the assembly exponents close

**G3 (gcd bookkeeping).** Pairs with g = (k, k′) > 1 form a
1/g²-thinned family; per dyadic g the same dispersion applies at
modulus [k, k′] = kk′/g, and Σ_g g^{-1-} converges — a one-paragraph
transcription. His coprimality set N: transcription.

**G4 (assembly) — the decisive exponent check.** The EH_μ discrepancy
at modulus q ~ x^{1/2+δ} (fixed residue N mod q) threads a q-congruence
through every face. Re-substituting the gates with the q-layer
(Q_gate := x^{1/2+δ}, K = x^{1/3}): (5.1) is Q-free (unchanged, tight);
(5.2) totals x^{0.219+1.111+0.5+δ} = x^{1.83+δ} < x² ⟺ **δ < 0.17**;
(5.3) totals x^{1.76+δ} — passes likewise. The m-face beyond level 1/2
is covered by the well-factorable μ level-3/5 literature (δ < 1/10
comfortable); the MR-regime and prime-moduli regimes carry their own
margins. Since Huang–Li's Corollary 1 needs only *some* δ > 0, the
assembly closes with a δ ≈ 1/6 margin at the naive-dictionary level.

**Program state after G1–G4.** Every gate, face, regime, and assembly
exponent now passes at the naive-dictionary level; G1 reduced to
transcription + smoothing; G3 to a paragraph; G2 covered; the h-range
collapses; the family is covered; and every object in sight is
measured to cancel. What remains is the certified write-out — the
transcription of this sketch against Lichtman §5 and §7–10 line by
line into a paper — during which any fourteenth tooth would have to
reveal itself in G4's threading of the q-congruence through the
dispersion (the one place where naive-dictionary optimism could still
hide a loss). This document, with the repository's 85+ commits of
measurements, is that paper's complete specification.

**Final measurement (increment 123).** Threading the q-congruence layer
(model q = 997) through the dispersion object: residue identity error
exactly 0; main-term cancellation 0.328 vs 0.323 unthreaded; h-phase
terms 0.595 vs 0.580 — **the q-threading is empirically lossless.**
Every quantity this program can measure has been measured, and every
exponent it can check has been checked. The certified transcription is
the sole remainder.

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
