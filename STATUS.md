# STATUS — single-page state of the program

*Last updated: increment 194 (2026-08-06). 4 days, 28 recorded
corrections, 16 pre-registered closures, 15 documented "teeth".*

**Consolidated working paper: `paper/negative_map.tex`** — the whole
campaign in one document (the two couplings, Theorem A and the
equivalence, Conjecture L, all sixteen closures with their blocking
coordinates, the five constraints on any future technique, and the
methodology that produced the negative results).

## Where things stand

By Huang–Li (arXiv:2005.03811) plus Bombieri–Vinogradov, binary
Goldbach for large even N follows from a single hypothesis,
$EH_\mu(N^{\theta'})$ for one $\theta' > 1/2$. This program mapped that
hypothesis from both sides. Both sides are now closed, for different
reasons, and the closure of each is the deliverable.

### 1. The demand side — closed by a theorem (increments 189–193)

Huang–Li consume $EH_\mu$ at exactly two places, their $E_3(\alpha)$
and $E_4(\alpha)$, distinguished only by the weight $w_k$.

- **Theorem A** (THEOREM_A.md, full proof in `paper/theorem_A.tex`):
  for $w_k = 1$ the signed, fixed-class sum is
  $\ll_A N(\log N)^{-A}$ for every $A$, **unconditionally**. Every
  ingredient is classical (BV at level $N^{1/2-\delta}(\log N)^{4A+8}$,
  Huang–Li's own Lemma 1, a density identity verified in exact rational
  arithmetic). It survived the program's adversarial-review protocol —
  the first derivation of this campaign to do so, after three were
  refuted.
- **Corollary**: their $E_4$ / Lemma 4 consumption becomes
  unconditional, and the entire $EH_\mu$ demand collapses to the single
  scalar $E_3$.
- **Theorem C (the permanent closure)**: for $w_k = \log k$ the same
  divisor switch returns the binary Goldbach sum itself, because
  $\mu * \log = \Lambda$. The unconditional identity
  $E_3(\alpha) = \sum\Lambda\Lambda - \mathfrak S(N)(N - \sum\Lambda\mu)
  + O_A(N(\log N)^{-A})$ — Huang–Li's own (22) — makes the weakest
  sufficient form of the demand **equivalent** to the conclusion. There
  was never a weakening to find, and no choice of $\theta'$, truncation
  or smoothing can evade an identity.
- **Net progress toward Goldbach: zero.** Theorem A removes only the
  half of the demand that carries no Goldbach content.
- By-product: a genuine defect in the published equation (18) — a
  dropped $n$-dependent truncation — with its repair
  (`paper/defect_report_18.md`, not sent).

### 2. The supply side — closed by measurement (increments 144–188)

- **Conjecture L** (CONJECTURE_L.md): every $\mu$-family this program
  probed factorizes as (deterministic local mask, computable by finite
  modular enumeration) × (exactly Gaussian fluctuation on the surviving
  support). The mask is blind-verified (corr 1.0000; annihilation
  fractions predicted exactly at fresh $N$); the Gaussian half holds at
  pair, cell, matrix and E1-ratio level across $10^8$–$2\cdot10^9$ and
  four $N$-structures. **Every "sub-random" reading in this program's
  history is mask accounting.**
- **The amplitude half is what the chain needs, and it is unreachable.**
  A fresh-context adjudication (AMPLITUDE_ADJUDICATION.md) found all
  five candidate routes blocked at named-lemma level; the common
  obstruction is that the bilinear pair constraint of
  $\mu(m)\mu(N-mk)$ is diagonalized by no character family, additive or
  multiplicative, while the $k$-average supplies no linearizing
  invariance.
- **Technique Forge** (TECHNIQUE_FORGE.md), three rounds, eight
  designs, eight pre-registered deaths: no internal lever
  (multiplicative orbit, manufactured congruence, Gram moments,
  $N$-average descent all exactly as flat as Conjecture L predicts),
  no external coupling surface (zeros invisible through the pairing,
  characters merely relocate the difficulty into thin progressions,
  determinant phases blind), and — round 3, increment 194 — **not even
  the program's own newest mechanism**: the divisor switch of Theorem A
  gives an *exact* identity on the full ranges
  ($\sum_k\sum_m \mu(m)\mu(N-mk) = \mu(N-1)$, perfect cancellation,
  verified), but none of it localizes into the type-II window (lag-1
  autocorrelation $+0.011 \pm 0.011$ over 8000 $k$ at two $N$). The
  mirror explains the whole asymmetry: the same switch puts Möbius on
  the **short** variable on the demand side and on the **long**
  variable on the supply side.
- **The Construction** (CONSTRUCTION.md), four classes, three closed by
  measurement (C1 abelian, C2 inverse-domain, C4 manufactured
  modularity) and the fourth (C-III, spectral/Voronoi) refuted by
  adversarial review at four independent coordinates (C3_REVIEW.md).

### 3. The refuted proof-program, retained as a record

PROOF_SKETCH_E1.md, `paper/e1_proof.tex`, `paper/e1_transcription.md`.
An independent adversarial review (REVIEW_VERDICT.md) found the gate
arithmetic vacuous, the role-exchange into the dispersion lemma
unfounded, and the SEAM formalization over-normalized and falsified by
our own data. **Read REVIEW_VERDICT.md before the sketch.** The exact
identities and every measurement survive.

### 4. The measurement corpus — the primary artifact

MEASUREMENTS.md, ~60 scripts. The ln 2 constant, Conjecture P /
Buchstab profile, the structure law at $\theta = 1/8$, the $\chi^2$
ladders, the final-axiom landscape through the $\sqrt N$ barrier, the
thin-progression stamps, the factorization-law stamps. One-shot
reproduction: `python code/verify_all.py`.

## What is and is not claimed

- **Claimed**: Theorem A and its corollaries (proved, reviewed,
  numerically confirmed); the exact identities; every measurement
  (reproducible, affirmed by the adversarial review); the failure map.
- **Conjectured**: Conjecture L; the ln 2 constant; Conjecture P; the
  empirical structure laws.
- **Not claimed**: any theorem toward Goldbach. Theorem A is
  unconditional but Goldbach-neutral, and the chain's remaining
  consumable — the amplitude half of Conjecture L — is untouched.
- **Withdrawn** (increment 143): "gates pass", "remainder = one
  1/30-wide seam", "the rest is transcription". **Withdrawn**
  (increment 193): the claim that Theorem A's true bound is
  $Ne^{-c\sqrt{\log N}}$; BV limits it to $N(\log N)^{-A}$.

## The map, in one sentence

The Huang–Li hypothesis has exactly two couplings of the two Möbius
factors: **divisibility-coupled** (the demand side, $k \mid N-n$),
where the divisor switch is free but hands back the Goldbach sum
itself; and **difference-coupled** (the supply side, $\mu(m)\mu(N-mk)$),
where the field offers no coupling surface, internal or external, in
any direction that has a mathematical name.

## How to continue

- **Read first**: `paper/negative_map.tex` — the consolidated working
  paper; everything below is its source material.
- **Entry point**: MEASUREMENTS.md, then `python code/verify_all.py`.
- **The theorem**: THEOREM_A.md → `paper/theorem_A.tex`.
- **The verdict**: REVIEW_VERDICT.md — read before the sketch.
- **The open problem**: CONJECTURE_L.md (amplitude half) with the
  reasons every route fails in AMPLITUDE_ADJUDICATION.md and
  TECHNIQUE_FORGE.md.
- **History**: 193 dated increments with 28 corrections and 15
  documented teeth — the map of every route that died, so no one
  repeats them.

*The wall, at final resolution, is a binary correlation that nature
cancels everywhere we can measure and no current technique — including
the routes this program invented and then refuted — can certify. The
dossier of measurements, the one unconditional theorem, and the failure
coordinates are the deliverable. — the program*
