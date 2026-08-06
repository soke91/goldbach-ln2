# STATUS — single-page state of the program

*Last updated: increment 200 (2026-08-06). 4 days, 31 recorded
corrections, 18 recorded closures (2 void, 1 downgraded, 2 moot, 13 standing), 15 documented "teeth".*

**Consolidated working paper: `paper/negative_map.tex`** — the whole
campaign in one document (the two couplings, Theorem A and the
equivalence, Conjecture L, all eighteen closures with their blocking
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
- **Theorem D (increment 195): the interior of the weight space is
  empty.** Theorems A and C are the two endpoints of a design space —
  the weight $w_k$ attached to $\mu(k)$. Writing $b=\mu*w$, extraction
  requires $b$ to have mass *at* the truncation point $K=N^{\theta'}$
  (Huang–Li's Lemma 1 damps everything below it by
  $e^{-c\sqrt{\log(K/d)}}$), while Bombieri–Vinogradov admits $b$ only
  below $N^{\theta'-1/2}$. The two thresholds are separated by exactly
  $N^{1/2}$, and the resulting loss $\exp(c\sqrt{\tfrac12\log N})$
  exceeds every power of $\log$. **No weight extracts
  $C(N)=\sum\Lambda(n)\mu(N-n)$ by divisor switching plus BV** — a
  no-go over the entire weight space of that method (not a claim about
  other methods).
- **Proposition D″ (increment 197): the smooth-weight family closes
  too.** Theorem D assumed $b$ low-supported; the natural family it
  excluded is $w_k = f(\log k)$, where $b = \mu*\log^D = \Lambda_D$ and
  the complete part splits by $\omega(u)$ into a Goldbach-type piece
  and Chen-type pieces. For a monomial every term is nonnegative
  ($\Lambda\ge0$, $\Lambda_D\ge0$), so it is $\asymp N(\log N)^{D-1}$
  with fixed sign; cancelling across monomials would have to be tuned
  against the asymptotics of the Goldbach sum itself. The one
  analytically canonical tuning (mean-zero $b$) moves it by 5%.
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
- **⚠️ Correction #30 (increment 198): the target itself was
  misstated.** The adjudication printed the consumable as
  $\sum_{k\sim K}|D(k)|^2 \ll (\log N)^{-A}\sum M_k$ — with $M_k$ where
  $M_k^2$ belongs, a factor $N/K$. As printed it demands *better than
  square-root cancellation*, and this repository's own data refute it:
  $\sum|D|^2/\sum M_k = 0.305/0.310/0.319$ at $N=10^8$ against a demand
  of $8.8\cdot10^{-6}$. The correct requirement, from Cauchy–Schwarz on
  $T_{II}=\sum_k b_kD(k)$, is
  $\sum_k|D(k)|^2 \ll (\log N)^{-2A-2}\sum_k M_k^2$ — a fixed log-power
  saving over the **trivial** bound, which square-root cancellation
  clears with margin $(N/K)(\log N)^{-2A-2}\to\infty$ (invisible at
  accessible $N$). The five *adjudication routes* are unaffected —
  they are blocked structurally, at named-lemma level, not on margin —
  but several **other** closures were budget calls and had to be
  re-audited; see the next item. (An earlier version of this line said
  "no route verdict changes", full stop. That was too quick — the
  re-audit voided two closures, downgraded one and mooted two.)
- **The amplitude half is what the chain needs, and it is unreachable.**
  A fresh-context adjudication (AMPLITUDE_ADJUDICATION.md) found all
  five candidate routes blocked at named-lemma level; the common
  obstruction is that the bilinear pair constraint of
  $\mu(m)\mu(N-mk)$ is diagonalized by no character family, additive or
  multiplicative, while the $k$-average supplies no linearizing
  invariance.
- **Technique Forge** (TECHNIQUE_FORGE.md), three rounds, eight
  designs, nine pre-registered deaths: no internal lever
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
  adversarial review — **but see the re-audit**.
- **⚠️ Re-audit of every closure (increment 199, CLOSURE_REAUDIT.md).**
  Correction #30 changed the yardstick, so every death sentence passed
  on a magnitude or scale-mismatch basis was re-examined.
  **One is void**: C-III's kill coordinate #1 explicitly assumed "E1
  demands |D(k)|² ≈ M_k", the square-root scale; the tree's own output
  $\sum|D|^2 \ll (N^2/K)(\log)^{-2C}$ is exactly the corrected target's
  form. **Four were opened for re-audit** (REVIEW_VERDICT #4, #5, #6
  and C-III #4's arithmetic half) and round 2 settled all four.
  **The rest stand**, being structural: violated lemma premises,
  absent congruences, and kill-tests that measured *no* signal rather
  than *insufficient* signal.
  **Round 2 (increment 200) settled the four**: C-III #4's arithmetic
  half is also **void** (its own premise is "no-margin arithmetic");
  **RV #6 is downgraded, not fatal** — the chain consumes the *signed*
  sum, L² enters only via Cauchy–Schwarz at a measured price of ≍ √K
  (58×, 45×), so L² is sufficient but strictly stronger than needed;
  **RV #4 and #5 are moot**, since #1 and #2 kill their program
  structurally. Net: **C-III is re-opened as a route** (its #2 and #3
  remain fatal to the draft as written) and what it needs is now
  stated exactly — a legitimate transform, a complete classification,
  and quantitative averaged Chowla at fixed log-power strength. The
  character of the obstruction changes with it: the **dilate** average
  admits no diagonalizing character family, but the **shift** average
  is MRT's home ground, and C-III's value is that it moves the
  difficulty from the first to the second.
  **Strategic corollary**: Cauchy–Schwarz discards the sign structure
  of $b_k$, the structure whose power R4 exhibited (for $b_k\equiv1$
  the unrestricted sum is exactly $\mu(N-1)$) — the program has been
  aiming at something strictly harder than its own target.

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
