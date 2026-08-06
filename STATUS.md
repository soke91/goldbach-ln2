# STATUS — single-page state of the program

*Increment 201 (2026-08-06). 4 days, 31 recorded corrections,
18 recorded closures (13 standing), 15 documented "teeth".*

**Consolidated working paper: `paper/negative_map.tex`.**
**Correction and supersession record: `CLOSURE_REAUDIT.md`** — every
statement this program has withdrawn, and what replaced it, in one
place. The other documents carry no correction banners; they state the
current position directly.

## Where things stand

By Huang–Li (arXiv:2005.03811) plus Bombieri–Vinogradov, binary
Goldbach for large even N follows from one hypothesis,
$EH_\mu(N^{\theta'})$ for a single $\theta' > 1/2$. This program mapped
that hypothesis from both sides. The demand side is closed at the level
of identities. The supply side is densely mapped but not closed: one
route is open.

### 1. The demand side — closed by theorems

Huang–Li consume $EH_\mu$ at exactly two places, $E_3(\alpha)$ and
$E_4(\alpha)$, distinguished only by the weight $w_k$ attached to
$\mu(k)$.

- **Theorem A**: for $w_k = 1$ the signed, fixed-class sum is
  $\ll_A N(\log N)^{-A}$ for every $A$, **unconditionally**.
  Ingredients all classical (BV at level
  $N^{1/2-\delta}(\log N)^{4A+8}$, Huang–Li's own Lemma 1, a density
  identity verified in exact rational arithmetic). First derivation of
  this campaign to survive the adversarial-review protocol, after three
  were refuted. **Corollary**: their $E_4$ / Lemma 4 consumption
  becomes unconditional and the whole $EH_\mu$ demand collapses to the
  single scalar $E_3$.
- **Theorem C**: for $w_k = \log k$ the same divisor switch returns the
  binary Goldbach sum itself, because $\mu * \log = \Lambda$. The
  unconditional identity $E_3(\alpha) = \sum\Lambda\Lambda -
  \mathfrak S(N)(N - \sum\Lambda\mu) + O_A(N(\log N)^{-A})$ — their own
  (22) — makes the weakest sufficient form of the demand **equivalent**
  to the conclusion.
- **Theorem D / D′**: the interior of the weight space is empty.
  Writing $b = \mu*w$, extraction needs $b$ to have mass at the
  truncation point $K = N^{\theta'}$, while a level of distribution
  $\theta_E$ admits $b$ only below $N^{\theta_E-(1-\theta')}$; the gap
  is $N^{1-\theta_E}$ and the loss $\exp(c\sqrt{(1-\theta_E)\log N})$
  exceeds every power of $\log$. **No weight extracts $C(N)$ by divisor
  switching — even granting the full Elliott–Halberstam conjecture.**
- **Proposition D″**: the smooth-weight family closes too. For
  $w_k = f(\log k)$ one has $b = \mu*\log^D = \Lambda_D$, and for a
  monomial every term of the complete part is nonnegative, so it is
  $\asymp N(\log N)^{D-1}$ with fixed sign; cancelling across monomials
  would have to be tuned against the Goldbach sum's own asymptotics.
- **Proposition E**: the circle method has zero margin on $C(N)$ —
  Cauchy–Schwarz lands above the trivial bound by a growing factor, and
  any pointwise route is capped by Parseval, $\sup \ge \|\cdot\|_2$.
- **Net progress toward Goldbach: zero.** Theorem A removes only the
  half of the demand that carries no Goldbach content.
- By-product: a genuine defect in the published equation (18) — a
  dropped $n$-dependent truncation — with its repair
  (`paper/defect_report_18.md`, not sent).

### 2. The supply side — mapped, not closed

The consumable is the signed type-II sum $\sum_{k\sim K} b_k D(k) \ll
N(\log N)^{-A}$; via Cauchy–Schwarz it suffices to prove
$\sum_{k\sim K}|D(k)|^2 \ll (\log N)^{-2A-2}\sum_{k\sim K}M_k^2$ — a
fixed log-power saving over the **trivial** bound. Nature supplies
square-root cancellation, far more than is asked; nothing certifies any
of it.

- **Conjecture L**: every $\mu$-family probed factorizes as
  (deterministic local mask, computable by finite modular enumeration)
  × (exactly Gaussian fluctuation on the surviving support). Mask
  blind-verified (corr 1.0000); Gaussian half holds at pair, cell,
  matrix and E1-ratio level across $10^8$–$2\cdot10^9$ and four
  $N$-structures. Every "sub-random" reading in this program's history
  is mask accounting.
- **Five routes adjudicated** against the source papers' verbatim lemma
  hypotheses: all blocked structurally. The common obstruction is that
  the bilinear pair constraint of $\mu(m)\mu(N-mk)$ is diagonalized by
  no character family, additive or multiplicative, while the
  $k$-average supplies no linearizing invariance.
- **Technique Forge**, three rounds, nine designs, nine deaths: no
  internal lever (orbit, manufactured congruence, Gram moments,
  $N$-average descent), no external coupling surface (zeros invisible
  through the pairing, characters merely relocating the difficulty into
  thin progressions, determinant phases blind), and not even the
  program's own divisor switch — which gives an *exact* identity on the
  full ranges, $\sum_k\sum_m\mu(m)\mu(N-mk) = \mu(N-1)$, none of which
  localizes into the type-II window.
- **The Construction**: C-I abelian, C-II inverse-domain and C-IV
  manufactured modularity are each closed by measurement against
  accurate nulls.
- **C-III is open.** Its draft is refuted — an incomplete
  classification, and a dual object that is the output of no legitimate
  transform — but the route is not. It needs a legitimate transform, a
  classification covering the type-II region, and quantitative averaged
  Chowla at fixed log-power strength (best known $(\log)^{1-c}$).
  What is measured about the transform: the dilate problem is $\mu\mu$
  correlations along lines of arbitrary rational slope, the shift
  problem is the slope-1 case, the obvious re-indexing is an exact
  relabelling (circular), and the slope family is not tamer than the
  shift family. So the transform must be genuinely analytic — and the
  named candidates for that are already closed.

### 3. The refuted proof-program, retained as a record

PROOF_SKETCH_E1.md, `paper/e1_proof.tex`, `paper/e1_transcription.md`,
with REVIEW_VERDICT.md as the verdict of record: the gate arithmetic is
vacuous (Q-slot premise violated), the role exchange into the
dispersion lemma is unfounded (no pair congruence), and the SEAM
formalization is over-normalized. **Read the verdict before the
sketch.** The exact identities and every measurement survive.

### 4. The measurement corpus — the primary artifact

MEASUREMENTS.md, ~65 scripts. The ln 2 constant, Conjecture P /
Buchstab profile, the structure law at $\theta = 1/8$, the $\chi^2$
ladders, the final-axiom landscape through the $\sqrt N$ barrier, the
thin-progression stamps, the factorization-law stamps. One-shot
reproduction: `python code/verify_all.py`.

## What is and is not claimed

- **Claimed**: Theorems A, C, D, D′ and Proposition E, with proofs;
  the exact identities; every measurement (reproducible, affirmed by
  the adversarial review); the closure map, at the strength stated for
  each entry.
- **Conjectured**: Conjecture L; the ln 2 constant; Conjecture P; the
  empirical structure laws.
- **Not claimed**: any theorem toward Goldbach; and no blanket claim
  that the supply side admits no coupling surface — one route is open.

## The map, in one sentence

The Huang–Li hypothesis has exactly two couplings of its two Möbius
factors: **divisibility-coupled** (the demand side, $k \mid N-n$),
where the divisor switch is free but hands back the Goldbach sum itself
and no weight evades it; and **difference-coupled** (the supply side,
$\mu(m)\mu(N-mk)$), where every coupling surface examined is absent and
one — a genuine analytic transform carrying arbitrary rational slope to
slope 1 — remains unexamined because nobody knows how to build it.

## How to continue

- **Read first**: `paper/negative_map.tex`; then `CLOSURE_REAUDIT.md`
  for what has been withdrawn and why.
- **Entry point**: MEASUREMENTS.md, then `python code/verify_all.py`.
- **The theorems**: THEOREM_A.md → `paper/theorem_A.tex`.
- **The open problem**: CONJECTURE_L.md (amplitude half), with
  AMPLITUDE_ADJUDICATION.md and TECHNIQUE_FORGE.md for why each route
  fails, and C3_REVIEW.md for the one that is open.

*The wall, at final resolution, is a binary correlation that nature
cancels everywhere we can measure and no current technique — including
the routes this program invented and then refuted — can certify. The
dossier of measurements, the unconditional theorems on the demand side,
and the failure coordinates are the deliverable. — the program*
