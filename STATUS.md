# STATUS — single-page state of the program

*Increment 332 (2026-08-07). 6 days, 177 recorded corrections,
18 recorded closures (13 standing), 15 documented "teeth",
9 named hazards.*

**Hazard 9** (inc. 330): *two summaries of one object are not comparable
until each one's weight is stated.* Increments 326–329 spent four passes
on a contradiction between `4.39` and `0.90` in which **both numbers
were correct** — one implicitly weighted by `|Λ̂(j/q)|² ≈ (X/φ(q))²`, the
other an unweighted mean over 255 frequencies of which only 22 mattered.
It cost three invented mechanisms, and the first of them was right and
discarded (#174). Distinct from #30: there a quantity met a **target**
at the wrong scale; here two measurements of the **same object** met
each other under different aggregations. ✅ **First swept at 331 (#175)** at the quantity most claims rest on: `ρ = Var C / V` names **three** quantities — a ratio of means, a mean of ratios, and the half-normal arm `(π/2)mean(|C|/√V)²`. They agree to `0.75%` at `N ≈ 1.4·10⁷` and differ by `10.3%` at `10⁵`, ordered `A > B > C`. **Proposition W wants the ratio of totals, so #121's `ρ = 0.810` is `0.841` in its units.** ✅ **Swept across every quoted ratio at 332 (#177), and the two that carry the most weight are exactly unambiguous**: `𝔄(N) = V/W` — Proposition V's local factor and increment 309's `0.8106` — gives `0.00%` between the two summaries in every band, and `Q(N)/N` likewise. Flagged: `ρ` at `3.10%` and the Cauchy–Schwarz deficit `C²/(W·Q)` at `9.10%`, both at the smallest band and both falling below `1%` at the largest. The deficit's ambiguity does not threaten 309, whose conclusion is that `deficit/√N` is **bounded** rather than equal to a constant.

**What is trusted, and what is not: `DEPENDENCY_AUDIT.md`.** The
theorems proved here do not rest on Huang–Li; the strategic frame does.
That frame has now been **independently re-derived here at every
structural step** — inversion, split, $\mu^2$ insertion, switch, the
$S_1$ evaluation, and the assembly giving $\mathfrak S(1-A(N))$ — with
the one defect we found (the dropped truncation at (18)) repaired.
No second defect was found. What is still taken on their authority is
the exponent-level analysis of the error terms.

**Consolidated working paper: `paper/negative_map.tex`** — now carries §"The wall's own law": Proposition V (the exact scale in closed form, Mirsky), the corrected law `C(N) = m(N) + √V·G`, Proposition W (the excess is a Chowla correlation), and the four measured facts — Gaussian bulk, lower-order mask, Gaussian tail, spectral component. Fifteen increments of results had lived only in markdown. ⚠️ Writing it exposed that the paper **had not compiled**: two propositions used an environment that was never defined (#103).
**CI stamp: `code/verify_all.py`** — repaired at increment 285. It had **no assertions and no failure path**: its criteria were characters inside output strings and it exited 0 regardless. Every stamp is now judged against a pre-registered interval, prints PASS/FAIL, exits 1 on failure, and carries a sensitivity block showing the verdict can flip. A **deep arm** (`30030 | N`) was added because the old sample forced `N ≡ 2 (mod 6)` and so excluded every `N` divisible by 3 — the cells where the location mask lives. ⚠️ Increment 306 pointed **hazard 8** at it and found a third defect: the intervals had never been calibrated against the rows' own spread, and under 40 seeds the stamp **failed a quarter of the time** — the recorded green was luck (#120). A sensitivity block shows a test *can* fail; it says nothing about whether the interval is the right *width*. Widths now come from the measured spread (`centre ± 4σ`), centres stay on theory where theory survives, and the re-audit gives **0% failure at 3.4–4.1 σ margins** (`code/audit_stamp_calibration.py`). En route it also gave **`ρ = 0.810 ± 0.018` at `N ≈ 10⁸`**, since V1's old centre was the `ρ = 1` half-normal that increments 288–289 had already refuted (#121).
**Gate linter: `code/lint_gates.py`** (inc. 307) — a file whose output **announces a verdict** must be able to return the bad one. It flags any script printing a verification banner with no `sys.exit(≠0)`, `raise` or `assert`. Built because increment 285's repair of `verify_all.py` was never swept to its siblings: **`verify_propositions.py` — the file `CLOSURE_REAUDIT` #61 holds up as the *answer* to "a check that cannot fail", and which this document cites as a reproduction command — printed `SOMETHING FAILED` and exited 0** for twenty-two increments (#123). Three files flagged, three repaired, corpus clean. Its own self-test caught it counting `sys.exit(0)` as a failure path before it ran on anything.
**Deep-N stamp: `code/verify_deep.py`** (rewritten inc. 307) — until then a verbatim pre-285 copy of the CI stamp, and a file called *verify **deep*** whose sample `// 6 * 6 + 2` excluded every $N$ divisible by 3. It now measures the half-normal ratio at $n = 300$ per arm — ten times `verify_all`'s power — and gates **the gap between them**, the one row in this program that stakes a sign as well as a size. At $N \approx 10^8$: deep $0.9476 \pm 0.0293$, shallow $0.7238 \pm 0.0245$, **gap $+0.2238 \pm 0.0056$ — forty standard errors** (#126). The shallow arm independently reproduces #121's $\rho = 0.810$, since $\sqrt{2\rho/\pi} = 0.718$. Intervals calibrated from the measured spread; 0% failure across 40 seeds at 3.9–4.0 σ.
**Open-questions register: `OPEN_QUESTIONS.md`** (inc. 308) — what the corrections **left open**, which nothing recorded. Two registers: six holes with nothing in them, and — the dangerous one — a triage of **closures whose premise a later correction moved**. `CLOSURE_REAUDIT.md` was written to re-audit every closure against correction #30; **ninety-six corrections have been recorded since and no second pass was made** (#129). The triage splits cleanly: every closure at risk is a **magnitude, normalisation or null** closure, every safe one is **structural**. ⚠️ And increment 311 attached the first **number** to a kill-test verdict: K2's detection floor is `4/√n ≈ 0.163` standard deviations of its field, so **"DEAD" has always meant "no coherent structure above a floor nobody quoted"** (#137). Its ALIVE branch is now shown rather than asserted — at twice the floor all ten `h` flag at `z = +7.2` to `+9.8`. ✅ Increment 312 then measured that quantity for K2 and **the closure stands quantitatively for the first time**: design K2 needs δ ≈ √(c·d(d−1)/n) against a floor of 4/√n, and at the `d ∼ (log N)^{2A+2}` the corrected E1 target actually requires, that is **14,000× the floor at A = 1** and `4.1·10⁶×` at A = 2. Its one blind region is `d ≤ 6`, worth at most a factor 6 — not a log power, so the route cannot spend it (#138). ⚠️ Increment 313 then found the remaining twelve do not all have K2's shape: their criteria are **ratios** (`2×`, `0.5×`), chosen as effect sizes and **never compared against the spread of the thing they threshold**. Restated in standard errors of their own nulls: **R1's threshold sits `+39.9 sd` above its null while its measurement sits `−0.80` below** — the data exclude a `7.5%` enhancement at 3 sd, thirteen times sharper than the recorded verdict claims (#140). **R2's ALIVE branch is unreachable as written**, asking for `≥ 2×` a control that came out negative; its verdict stands on the measurement (`−0.38 sd`) and not on the criterion (#141). **R4 and C4 print no error bar beside the quantity their threshold judges** (#142). ✅ Increment 314 supplied C4's, with 40 null draws per level instead of six, and **the criterion does not survive**: the `0.5×` threshold sits only `−1.29 sd` below its own null, so **pure noise satisfies it 8.8% of the time** and the recorded "0 alive of 6" had probability 58% under noise — unremarkable. The DEAD direction still stands (a trigger-happy test that did not trigger is evidence of absence) but **the closure must rest on the measurement, `+0.10 sd` on average and never systematically below the null, not on the count of flagged levels** (#143). ✅ Increment 315 supplied R4's, and the answer depends entirely on the block size: `ratio(B)` is a sum of `nb = 2048/B` squares, so **R4 is a 5-sigma test at `B = 8` and carries no information at `B = 512`, where four blocks are averaged** — its precision degrades as `√(2/nb)`, fastest exactly where its own signature would appear (#145). R4b's worry about "a mild deficit in the same direction at both `N`" is answered: the fall is **under 2 sd at every `B`**, and the replication reverses its sign at `B = 512`. ⚠️ R4b's quoted SEs are themselves underestimates — bootstrapped, `35.4%` against the stated `25%` at `B = 64` — so every R4 test is **weaker** than the recorded numbers say (#146). It now carries a third register: **claims that still stand but were DERIVED from one that fell** — the opening sentence used to say a withdrawn positive claim "costs nothing but a claim", which is wrong (#135). And its kill-test row, which asserted a property of eleven files without opening them, was **false for twelve of thirteen**: one permutation null, four already carrying the coin control hazard 7 asks for (#134, `code/audit_killtest_nulls.py`). Enforced in the other direction by **`code/audit_withdrawn_forms.py`** — no live document may still *assert* a withdrawn form. It found the working paper carrying two (#127).
**Notation register: `NOTATION.md`** (inc. 321) — every symbol that has collided or could. It exists because **`ρ` denoted two unrelated quantities in this document**: Huang–Li's Lemma 1 function in Theorem D (`0.0137` at `N = 10⁸`) and the cancellation ratio `Var C / V` (`0.810`). Read together, "the loss factor is `≥ 1/max_d|ρ|`" and "`ρ = 0.810`" give `1.23` where the truth is `≥ 73` — the same species as the 𝔖 / 𝔄 collision of #74, #75, #83 (#159). Huang–Li's is now `ρ_HL` everywhere and `lint_docs.py` check **(H)** enforces it. ✅ Asked whether quantum uncertainty or entanglement could help, increment 321 answered in numbers: **uncertainty is Parseval** — `sup|S_μ| ≥ ‖S_μ‖₂ = √(6N/π²)` is Proposition E's route (ii), the inequality that *closes* the pointwise route, and its slack is a bounded `3.5`–`3.9` with fitted slope `+0.0005`; saturating it exactly buys `3.79×` and leaves the best conceivable pointwise product at `1.92 ×` trivial, saving nothing against a target needing `(log N)^A`. **Entanglement** is `ρ − 1`, which Proposition W identifies as a Chowla correlation — and which increment 320 showed cancels out of the comparison entirely. **Neither is a resource** (#160).
**Verdict linter: `code/lint_verdicts.py`** — a conclusion must be *computed*, not *composed*. It flags any `print` of a plain string literal carrying verdict vocabulary outside a branch, since such a line asserts an outcome the run cannot contradict; f-strings with a computed value, prints inside a conditional, and lines marked `# verdict-ok: <reason>` pass. Built after a script printed the result it expected and the run rejected it (#100). First pass over 226 files: 10 candidates, 7 criterion-or-structural and marked with their reason, **3 real and all mine**.
**Document linter: `code/lint_docs.py`** — six mechanical checks over every `.md`/`.tex` file, exit 1 on failure, with a self-test that shows each one detecting a synthetic fault. It exists because escape collapse through shell heredocs had corrupted tracked files five times and been answered five times with a note; it found a sixth, live, on its first run. Two of its own invariants failed that self-test before the third passed.
**Correction and supersession record: `CLOSURE_REAUDIT.md`** — every
statement this program has withdrawn, and what replaced it, in one
place. **The other documents state the current position directly and
do not argue with themselves**: where a statement is superseded its body
is rewritten, not bannered. Increments 280–290 broke that rule — sixteen
correction banners accumulated outside this file, and CONJECTURE_L.md
reached the point of displaying a formula and then denying it three
times below. Repaired at increment 291 (#89): every body that
denied itself has been rewritten to state the current position, with the
superseded form named once and its history pointed here. Six markers
remain outside this file — three in MEASUREMENTS.md, two here, one in
LITERATURE.md — and all six are **forward-looking cautions about live
figures** (which convention a number used; that a cross-check is
near-algebraic; that the no-go theorems are asymptotic), not denials of
the text above them.

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
- **Proposition V** (inc. 287): the wall's exact scale has a closed
  form. $V(N)=\sum_v\mu^2(v)\Lambda(N-v)^2$ is a $(\log p)^2$-weighted
  count of **squarefree shifted primes**, so $V(N)=W(N)\,\mathfrak A(N)(1+o(1))$
  with $W(N)=\sum_{w<N}\Lambda(w)^2$ and
  $\mathfrak A(N)=\prod_{q\nmid N}(1-1/(q(q-1)))$ (Mirsky 1949).
  Verified to $1.000000\pm0.000145$. ⚠️ **The local factor is
  $\mathfrak A$, not $\mathfrak S$** — residual sd 0.000323 against
  0.245235, a factor of 760. The singular series of the *count* is not
  the local factor of the *noise*, and substituting one for the other
  is the single error behind corrections #74, #75 and #83.
- **Theorem D / D′**: the interior of the weight space is empty.
  Writing $b = \mu*w$, extraction needs $b$ to have mass at the
  truncation point $K = N^{\theta'}$, while a level of distribution
  $\theta_E$ admits $b$ only below $N^{\theta_E-(1-\theta')}$; the gap
  is $N^{1-\theta_E}$ and the loss $\exp(c\sqrt{(1-\theta_E)\log N})$
  exceeds every power of $\log$. **No weight extracts $C(N)$ by divisor
  switching — even granting the full Elliott–Halberstam conjecture.**
- **Proposition D‴** (inc. 278): the loss factor is
  $\ge 1/\max_d|\rho_{HL}|$, hence **monotone** in the bound on
  $\rho_{HL}$ — Huang–Li's Lemma 1 function, **not** the cancellation
  ratio $\rho = \mathrm{Var}\,C/V$; the two differ by two orders of
  magnitude and shared a symbol until increment 321 (`NOTATION.md`,
  #159) — every improvement of
  Huang–Li's Lemma 1 *strengthens* the no-go, none can weaken it. Under
  RH it becomes a power of $N$. The route is not blocked by what we
  cannot prove; it is blocked harder the more we know.
- **Proposition D⁗** (inc. 279): the boundary is **not** $\theta_E = 1$.
  Solving gives $1-\theta_E \ll (\log\log N)^2/\log N$ — strictly
  stronger than EH, strictly weaker than $\theta_E = 1$ — and that
  regime is **not vacuous** (progressions retain $\exp(C(\log\log N)^2)$
  terms, beating every fixed power of $\log N$). Corrects prose that had
  stood since increment 196.
- ⚠️ **All of these are asymptotic and bite late**: Theorem D′ acquires
  content around $N \approx 10^{480}$ at $\eta = 0.40$, $10^{3071}$ at
  $\eta = 0.10$. They constrain **methods**, not any computation anyone
  will run (`results/lab_level_threshold.txt`).
- **Proposition D″**: the smooth-weight family closes too. For
  $w_k = f(\log k)$ one has $b = \mu*\log^D = \Lambda_D$, and for a
  monomial every term of the complete part is nonnegative, so it is
  $\asymp N(\log N)^{D-1}$ with fixed sign; cancelling across monomials
  would have to be tuned against the Goldbach sum's own asymptotics.
- **Proposition E**: the circle method has zero margin on $C(N)$ —
  Cauchy–Schwarz lands above the trivial bound by a growing factor, and
  any pointwise route is capped by Parseval, $\sup \ge \|\cdot\|_2$.
  **Re-derived at increment 309 and strengthened.** Proposition V makes
  $W(N)$ cancel between $\sqrt{W Q}$ and
  $\mathrm{rms}\,C=\sqrt{\rho\,\mathfrak A\,W}$, so the method is short
  of $C(N)$ by $\sqrt{Q/(\rho\mathfrak A)}\asymp\sqrt N$ — **a clean
  power of $N$, with the log powers cancelling identically**. Measured
  $0.888\to0.961$ in units of $\sqrt N$ over a factor 100; route (ii)'s
  Parseval floor $1.52\to1.92$, never below 1. "Zero margin"
  understates it. Nothing in the proof depends on any recorded
  correction, so `OPEN_QUESTIONS.md`'s triage mis-assigned it (#130).
- **Net progress toward Goldbach: zero.** Theorem A removes only the
  half of the demand that carries no Goldbach content.
- **Where the wall is, in numbers** (inc. 319). Measured truth
  $C(N)\asymp\sqrt N(\log N)^{0.29}$ against trivial $\psi(N)\sim N$: the
  margin at $N=10^8$ is $N^{0.454}$, **a full power of $N$**, and measured
  $\max|C|/N$ falls $0.056\to0.0082$ across the octaves. **The gap is not
  the size of the extreme.** What is missing is a *pointwise method
  losing less than a power of $N$*: Proposition E puts the circle
  method's loss at $N^{0.498}$ — the margin's own size, which is why
  "zero margin" is exact rather than rhetorical — while Theorem D's is
  only $N^{0.233}$ and it still fails, its obstruction being the absence
  of mass at the truncation point rather than a size deficit.
  ⚠️ **Corrected at 320 (#157): that comparison used the wrong
  reference.** A proof needs *bound ≤ target*, not *bound ≈ truth*. In
  the target's own units, trivial sits $(\log N)^A$ above the target,
  Cauchy–Schwarz $(\log N)^{A+1/2}$ above it, and the truth a power of
  $N$ **below** it. **The entire difficulty is a log power** — the parity
  problem in its own units — and the two routes are then the *same*
  failure in one currency: both must supply $(\log N)^A$ over trivial,
  the circle method instead loses $\sqrt{\log N}$ and the divisor switch
  loses $\exp(c\sqrt{\log N})$. $\rho$ and $\mathfrak A$ **cancel out of
  the comparison entirely**, which is why sharpening the Möbius input
  cannot help. One new exact identity came
  out of the attempt, verified to $1.5\cdot10^{-16}$:
  $\sum_N C(N)^2=\sum_h M(h)P(h)$ with $M$ the Möbius autocorrelation
  and $P$ the prime-pair count — the global form of Proposition W. It
  buys almost-all $N$ by Chebyshev, which is 1938, and no more (#155).
- By-product: a genuine defect in the published equation (18) — a
  dropped $n$-dependent truncation — with its repair
  (`paper/defect_report_18.md`, not sent).

### 2. The supply side — mapped, not closed

**The tail, which is the only thing the requirement constrains**
(inc. 290). A bulk can be Gaussian to five decimals while the tail is
heavy, and $C(N)=o(N)$ constrains **every** $N$, not typical ones.
With the mask removed and $V$ exact, $\max|Z|$ tracks the Gumbel law:
mean deviation from $E[\max]$ is $+0.54\pm0.45$ over eight bands, and
aggregate tail counts against the Gaussian expectation give ratios
**0.999** ($t=3$), **0.997** ($t=4$), 0.878 ($t=5$). The extremes sit at
**generic** $N$ ($2\cdot8317$, $2\cdot138917$, ...), not deep radicals, so
the mask removal is not leaking into the tail. **The margin at the
extreme** — the figure the program should quote, and never has — is
$\sqrt N/(a_n\sqrt{\mathfrak A\log N})$: $10^{4.4}$ at $N=10^{12}$,
$10^{22.8}$ at $10^{50}$. **The requirement is not remotely tight; the
whole difficulty is in proving it** (#87, #88).

**And what that excess IS** (inc. 289, **Proposition W**): expanding
$C(N)^2$ gives $\rho-1=(1/V)\sum_{h\ne0}c(h)S(h)$ with $c(h)$ a
weighted prime-pair count and $S(h)$ the **binary Chowla correlation**.
So the wall's excess over square-root is a prime-pair-weighted Chowla
correlation, and Chowla-type input forces $\rho\to1$: **the wall is
exactly square-root, over-delivering by a power of $\log$ and no more.**
Measured: $M(h)$ sits at $1.051$–$1.068$ times the random-sign floor
$\sqrt{0.32264(X-h)}$ across five decades of shift; the reconstruction
gives $\rho-1=-0.0976$ against a measured $-0.18$; the sign is negative.
**The wall leans on the provable end of Chowla** — shifts $h<10^3$ carry
1.1% of the gross mass, 48.9% sits at $h\approx10^5$–$10^6$.

**How much cancellation nature actually delivers, measured at last**
(inc. 288). Proposition V makes $\rho(N):=\operatorname{Var}C(N)/V(N)$
well posed — it is exactly 1 under random signs on the squarefree
support. With the location mask removed, $\rho$ **rises** 0.760 → 0.837
over $N \le 1.6\cdot10^7$, while the *raw* ratio falls 1.006 → 0.858:
opposite directions, converging (gap 0.246 → 0.020). **The recorded
downward trend was mask contamination and its sign was wrong.** So the
wall does beat a coin ($\rho<1$) but by a margin that is **shrinking**,
not growing. ⚠️ **Withdrawn at increment 300.** The centred estimator behind every
number in this paragraph is **biased**, and the bias accounts for the whole
effect: replacing $\mu$ by a random $\pm1$ on the same support leaves $V$
identical and every other step byte-identical, and the coin reproduces the
real curve with $z$ between $-0.5$ and $+0.4$ in every band. **The centred
estimator cannot tell $\mu$ from a coin.** Nearby $N$ share the same
$\mu(v)\Lambda(N-v)$ terms, so the $C(N)$ are positively correlated across a
band and subtracting the band mean removes real variance. Withdrawn with it:
#84's direction, #99's $b=2.68$, and the quantitative half of #86.
Proposition W's identity is algebra and is untouched; what is gone is every
measurement that claimed to confront it. Under the estimator the coin does
validate --- the uncentred $E[C^2]/V$, exactly 1 for a coin --- the real
ratio runs $1.008$ down to $0.858$. **That is not established either**
(inc. 301): with 40 coin draws, 0 of 8 bands pass $-5\sigma$, 12 of 40 coin
trends are more negative than the real one, and pooled the real $0.906$ sits
at $z=-0.77$ against coin $1.002\pm0.125$ with 9 of 40 draws below it. The
$z=-4.2$ of #107 divided by the standard error of the coin *mean* rather
than by its spread (#108). ⚠️ **The low-variance estimator is biased and the
unbiased one is too noisy** — pooled deficit $0.094$ against a coin scatter
of $0.125$, because the $C(N)$ across a band share most of their terms. So
Proposition W's $\rho-1$ has **never been measured**, and with one realisation
of $\mu$ it is not measurable by either route. How much the wall beats a coin
is currently **unknown** (#106--#109).

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
  $N$-structures. For the **wall** $C(N)$ the Gaussian half holds
  only under the **exact** second moment $V(N)=\sum_v\mu^2(v)\Lambda(N-v)^2$
  (excess kurtosis $-0.0005$, $z=-0.3$); under the $\mathfrak S N$-based
  scale this program displays it **fails at $z=98$** (inc. 283). Every "sub-random" reading in this program's history
  is mask accounting. The **original** half (the $\mu$-families
  $D(k)$) was re-tested at increment 284 on 285,050 pairs and **holds**:
  excess kurtosis $-0.0034$ ($z=-0.4$), $E|Z|/\mathrm{sd} = 0.79760$
  against $0.79788$, no class structure. M.3 is **exact in both
  directions** on 401,000 pairs (115,950 predicted zeros, 115,950
  observed, none unpredicted).
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
- **C-III**: its draft is refuted, and two of its three outstanding
  requirements are now settled against it. **①** is closed for
  pointwise changes of variable (see below) and unpopulated for
  summation formulas. **②** is settled by measurement: completing the
  classification shows the draft's assumption $a \le y^{O(1)}$ excludes
  **~95% of the Heath–Brown weight** (0.939–0.969 above $M^{0.05}$,
  flat in $J$ over 3–8, rising with $M$), and those pieces carry a
  rough coefficient with no divisor structure, hence no Voronoi entry —
  so the spectral door serves only a few percent of the mass. **③**,
  quantitative averaged Chowla at fixed log-power strength (best known
  $(\log)^{1-c}$), is a named external open problem. Scope: ② is for
  the one-sided Heath–Brown opening as drafted; a different opening
  could reallocate mass, and nothing here rules that out.
  On the transform: in centered coordinates $A = N-a$, $B = N-b$ one
  has $m'A - mB = 0$ exactly, so the dilate family is the **pencil of
  lines through the origin** and the shift family is a **parallel**
  family of slope 1. A pencil's vertex is finite, a parallel family's
  is at infinity, and affine maps send finite points to finite points —
  so **no integral change of variable can do it**, which proves what
  the earlier probe measured as circularity. What remains is a
  summation formula, blocked in general by the roughness of the outer
  weight $\mu$ (Voronoi compresses only against smooth weights), with
  every named candidate in that class already closed.

### 3. The refuted proof-program, retained as a record

PROOF_SKETCH_E1.md, `paper/e1_proof.tex`, `paper/e1_transcription.md`,
with REVIEW_VERDICT.md as the verdict of record: the gate arithmetic is
vacuous (Q-slot premise violated), the role exchange into the
dispersion lemma is unfounded (no pair congruence), and the SEAM
formalization is over-normalized. **Read the verdict before the
sketch.** The exact identities and every measurement survive.

### 3b. The wall's own scalar, measured

The chain reduces to $C(N) = \sum_{n<N}\Lambda(n)\mu(N-n) = o(N)$,
which is equivalent to the Goldbach asymptotic. Measured over 400 even
$N$ in five octave groups (MEASUREMENTS §9, `code/h_deficit.py`):
**$|C(N)| \sim N^{0.5457\pm0.0032}$** — measurably above square root,
which is what the variance law requires, against a requirement of merely
$o(N)$. (§9 recorded $N^{0.503}$, "square-root to three digits", from
five groups of 80 consecutive even $N$; replicating that design 500
times gives $0.516\pm0.043$, so all three digits were noise. #70.)

The discrepancy $R(N) = r(N) - \mathfrak S(N)(N - C(N))$ splits
exactly as $R = E + \mathfrak S C$ with $E(N) = r(N) - \mathfrak S(N)N$ the
**classical Goldbach error term**, and $E$ **dominates**: 60% of $|R|$ at
$N\approx10^5$, rising to 70% at $1.6\cdot10^7$ (inc. 292). So the Goldbach
count's deviation from $\mathfrak S(N)N$ is carried by the classical error,
not by the wall — **the ln 2 comet corpus and the final-boss scalar are
not the same object**, a hypothesis raised and refuted here, and now
refuted by a cleaner mechanism than the exponent comparison §9 used.
**And $E$ really is the zeta-zero error, by location and not by size**
(inc. 293). If $E(N)\approx-2\sum_\rho N^\rho/\rho$ then $E(N)/\sqrt N$
carries $e^{i\gamma\log N}$ for every zero, so its spectrum in $\log N$ must
have lines at the ordinates. Binned on a uniform $\log N$ grid and tested
against a **local** null: **7 of 10** of the first ordinates clear their
local 99th percentile ($\gamma_2$ at $2.85\times$), while a phase-randomised
surrogate gives **0 of 10**. First detection of the zeta zeros anywhere
in this program's data (#91, #92). Recalibrated at increment 294 by a
**permutation test** — the local-percentile null of 293 could not be
estimated from a window holding only ~6.5 independent frequencies — and
it survives: $E/\sqrt N$ gives $z=+12.4$, $p\le0.005$ (#95).

**But the lines in the wall are $\Lambda$'s, not the wall's** (inc. 302).
$C(N)$ does carry the ordinates — but replacing $\mu$ by a random $\pm1$ on
the same support, through the same $\Lambda$ and the identical pipeline, gives a
coin mean of $2.99\cdot10^{-3}$ against the real $3.90\cdot10^{-3}$, with **6 of
20 coin draws at or above it**. <!-- withdrawn-ok: the 1566x is quoted
in order to withdraw it --> The ratio is $1.30\times$, not the $1566\times$
reported against a *permutation* null. $C=\mu*\Lambda$ and $\Lambda$ carries the
zeros by the explicit formula, so this was never evidence about $\mu$ (#110).
**Hazard 7, which is behind this and #106 both: a null that destroys
everything identifies nothing.** When the object is built from several
arithmetic inputs, the null must vary exactly one and leave the rest
byte-identical (#111). ✅ What survives across both controls: **every claim
that $\mu$ behaves like a coin stands; every claim that $\mu$ is special has
fallen.**

Sizes, with the log powers solved for rather than left inside an
effective exponent: $|E|\sim\sqrt N(\log N)^{1.51}$,
$|\mathfrak S C|\sim\sqrt N(\log N)^{0.20}$.

> **Re-audited at increment 282: the closure stands, the evidence for it
> did not, and the exponent was wrong.** Full census to $1.6\cdot10^7$:
> $\beta_R-\beta_C = +0.0997\pm0.0034$ (**29 s.e.**), so the conclusion
> is now established — but replicating §9's own design 500 times gives
> that difference as $+0.079\pm0.054$, making the recorded gap **1.8
> s.e.**, with **7% of replicates returning the opposite sign**. The
> criterion statistic was also mis-estimated: $\mathrm{corr}(C,D)$ is
> **0.80 falling to 0.71** per band, not the 0.60 that pooling across
> scales produced, and removing the mask *raises* it by 0.042. What
> makes the closure safe is that the correlation **falls** with $N$
> ($-0.0198\pm0.00035$ per unit $\log N$), moving away from the 0.9
> threshold — the one thing §9 did not check. Corrected exponent:
> $|R(N)| \sim N^{0.6458\pm0.0065}$.

A twenty-hypothesis sweep (MEASUREMENTS §12) then pinned the scalar's
law. Ten hypotheses on the local structure of $D(k)$ returned **0
flags out of ~22 statistics**; five variants of $C(N)$ are all
square-root sized, so there is no softer target; and seven hypotheses
on $C(N)$ itself fired five flags which resolve into one statement:

**And a law for the process, not only the marginal** (inc. 324). A
covariance function *is* a process specification, and increment 304
measured this one in closed form: $\rho(h) = a\,\mathfrak S_2(h) + b$ at
$\mathrm{corr} = 0.9997$–$1.0000$. Since $\mathfrak S_2$ has the
Ramanujan expansion $\sum_q (\mu^2(q)/\varphi^2(q))c_q(h)$, the spectral
measure must be **atoms at the rationals $j/q$ with Hardy–Littlewood
weights**. Tested on $3.87\cdot10^6$ even $N$, chosen a multiple of
$30030$ so every such frequency lands exactly on a bin: the atomic bins
carry **$525\times$** an equal number of others, and on a coin the
per-frequency mass tracks $\mu^2(q)/\varphi^2(q)$ at
**$\mathrm{corr}=+0.9864$** across sixteen moduli. ⚠️ The coin
reproduces it, so **the dynamics are $\Lambda$'s through the shift, not
$\mu$'s** — the fourth structure in the wall to turn out that way. ⚠️ And
it is a law, not a bound: the decomposition it suggests reduces to a
pointwise Siegel–Walfisz for $\mu*\Lambda$, which is the parity problem
again (#162). ⚠️ **Qualified at 325 (#164): "purely atomic" was a claim
about DENSITY per bin**, which holds at $525\times$; the sixteen exact
moduli carry only **11.3%** of the coin's mass and **7.0%** of the
real's. That is what the Ramanujan measure predicts — it has infinitely
many atoms and $\sum_q \mu^2(q)/\varphi(q)$ diverges — and the
cumulative mass tracks the cumulative weight at increment correlation
$+0.9110$, stable to $4.7\%$. Of the eight largest peaks outside the
exact set, three are DC residue, three are leakage skirts of atoms
already counted, and two are **genuine atoms at $q = 17$ and $q = 23$**,
the moduli the set excludes. Nothing unexplained.
⚠️ **And at 326 the coin control came out backwards** (#166). The mask
and the atoms are exactly nested — the cell projection sits inside the
atomic span at containment $1.0000$ — and the mask does *not* exhaust
them: the atomic projection carries $1.52\times$ the cell projection's
energy. But **the coin's atom-to-cell ratio is $16.87$ against the
real's $1.52$**, and in absolute share the coin carries $14.78\%$ of its
energy in the atoms against the real's $7.00\%$. The real's cell energy
is $5.25\times$ the coin's, which is the mask behaving as it should;
its *non-mask* atomic energy is several-fold **below**. So $\mu$
**suppresses** the periodic covariance $\Lambda$ supplies, except where
M.1 creates one — the first place in this program where $\mu$ has looked
*unlike* a coin in a direction that is not an estimator defect. ⚠️ The
mechanism I proposed (Davenport smallness at rationals) was tested and
**refuted**: $|\sum_v\mu(v)e(vj/q)|^2/Q = 0.9492$ against the coin's
$0.9975$, a factor $1.05$ (#167). ⚠️ **Suspended at 327 (#169).** The first suspect — that the
$1/\sqrt{V}$ rescaling created it — is cleared: the **raw** $C(N)$
already shows $5.155\times$ and the rescaling contributes about $1\%$
(#168). But two of my own measurements of the same object now disagree
by a factor 5. The atoms sit at $j/(2q)$, not $j/q$ — the periodogram
runs over even $N$ — and #167 measured at the wrong frequencies; **at
the corrected ones the ratio is $0.91$, still nothing like $5.155$.**
Since the spectrum of $C$ is exactly $|\hat\mu|^2|\hat\Lambda|^2$, a
share ratio of 5 with a per-frequency ratio of 1 requires
$|\hat\mu|^2$ to differ elsewhere in the spectrum, and where is not
identified. **#166 is suspended until they are reconciled — neither
withdrawn nor standing.** ⚠️ **Localised at 328 (#170):** split
absolutely, the **numerator** ratio is $4.390$ and the **denominator**
ratio $1.140$, so the factor lives in the atomic bins and not in the
window's total. ⚠️ And a third mechanism — the even-$N$ subsampling,
which aliases $f$ with $f+1/2$ and, since $\Lambda$ lives on odd
numbers, should leave $2\hat\Lambda(f)\sum_{v\ \mathrm{odd}}\mu(v)e(vf)$
— was **refuted too**, at $0.90\times$ (#171). Three mechanisms
proposed and three refuted. **The rule that follows: stop proposing
mechanisms.** When every direct measurement of $|\hat\mu|^2$ says
$\approx 1$ and an aggregate says $4.39$, the likeliest explanation is
that one of the two computations is wrong, and the next step is to
verify the periodogram against the exponential sum **at a single bin**.
✅ **Resolved at 329 (#172, #173), and it reverses the reading.** The two
summaries carried different weights: the periodogram sums
$|\hat\mu|^2|\hat\Lambda|^2$ and $|\hat\Lambda(j/q)|^2\approx(X/\varphi(q))^2$,
so $q=3$ outweighs $q=143$ by $3600$, while every exponential-sum check
reported an **unweighted mean** over 255 frequencies. Per modulus both
objects move together: the coin exceeds the real by $8.40\times$ at
$q=3$, $15.16\times$ at $q=5$, $4.97\times$ at $q=7$, $3.60\times$ at
$q=11$, falling to $\approx1$ from $q=13$ up — and only 22 of the 255
frequencies sit at $q\le11$, so the mean drowned them. **So #167's
refutation was itself the artefact**: $\mu$'s exponential sums at
small-denominator rationals really are far smaller than a coin's, which
is the classical major-arc cancellation, and **#166 is restored with
that as its mechanism**. ⚠️ It is not a new fact about $\mu$ — it is this
program measuring a classical one for the first time, after three wrong
refutations — and the reconciliation is qualitative: $13.0$ against
$8.4$ at $q=3$, so **#166 stands as a direction, not as a number**.

> **$C(N) = m(N) + \sqrt{V(N)}\cdot G(N)$**, with
> $V(N)=\sum_v\mu^2(v)\Lambda(N-v)^2$ exact, $m(N)$ the location mask,
> and $G$ Gaussian in the bulk and in the tail.

The mask $m(N)$ is real (M.1 is a theorem) and, measured against an
**exactly computed** floor, is **resolved in every band to
$1.6\cdot10^7$**: $\max_c|z_c|$ runs 11.1, 10.9, 12.9, 12.6, 11.2,
10.3, 9.4, **8.4**, and the deepest cell $3\cdot5\cdot7\cdot11\cdot13\mid N$
sits at $-3.58\sqrt V$ in the top band (inc. 305). Increment 303 had
called the large-$N$ amplitude *unresolved*; that was a statement about
the **statistic**, not the mask. The aggregate
$B=\sum_c (n_c/n)(m_c-\overline m)^2$ weights each cell by its **size**
while the mask lives in the **rare** deep cells — the largest-$|z|$ cell
holds $6.6\cdot10^{-5}$ of the top band — so with the floor removed
exactly $B_{\text{mask}}$ still clears twice its error in only **2 of 8**
bands. **#69's statistic could not have measured the mask at large $N$
whatever the floor** (#118). Its **scaling remains withdrawn**: a decay
beats a constant at 5 of 6 depths, but over a factor 160 in $N$ the data
do **not** separate $N^{-a}$ from $(\log N)^{-b}$, and the fitted
exponent varies about $4\times$ **with the cell** under either
parameterisation ($a = 0.14, 0.22, 0.27, 0.37$ at depths 5→2), so **there
is no single mask exponent** to quote (#119). It
does not threaten $C(N)=o(N)$. The form this section carried until
increment 283 — $\sqrt{\mathfrak S(N)N}\,G(N)$ with $G$ of unit variance
and a "mean drift decaying with $N$" — is superseded; that drift was the
mask, and the scale was wrong. History in `CLOSURE_REAUDIT.md`
#36, #67–#69, #74, #83.

**What that coin floor is** (inc. 304). It is not estimation noise and
not a property of "the estimator": it is the **singular series of the
shift**. For a coin $\operatorname{Var}Z(N)=V/V=1$ exactly, so the
covariance *is* the correlation and needs no simulation —
$\rho(h)=(\mu^2*g_h)(N)/\sqrt{V(N)V(N+h)}$ with
$g_h(w)=\Lambda(w)\Lambda(w+h)$, whose numerator is a **prime-pair count
at shift $h$** and therefore carries $\mathfrak S_2(h)$. Measured:
$\rho(h)=a\,\mathfrak S_2(h)+b$ with
$\operatorname{corr}(\rho,\mathfrak S_2)=0.9997$–$1.0000$ in all eight
bands. Two $N$ in one cell agree on which small $q$ divide them, so
$h=N'-N$ is divisible by small primes more often — $1/(q-1)$ against
$1/q$ — exactly where $\mathfrak S_2$ is larger. **Permuting the cell
labels across $N$**, which preserves every cell size and leaves $Z$
byte-identical, collapses $B/T$ from $0.053$ to $0.0034$ at the bottom
band and from $0.038$ to $0.00015$ at the top: the $(k-1)/n$ of an
independent sample. So the floor's flatness in $n$ is explained — **the
excess is a property of *pairs*, and no sample size removes it** (#116).
Its **size** is not: against the closed form $E[B]/E[T]$ the singular
series overpredicts by $1.13\times$ to $1.82\times$, growing with $N$,
and the one suspect I could name was tested and cleared.

Two by-products. The cell floor of *any* coin statistic is now
available **exactly**, from
$E[(\sum_{N\in c}Z)^2]=\sum_v\mu^2(v)u_c(v)^2$ with
$u_c(v)=\sum_{N\in c}\Lambda(N-v)/\sqrt{V(N)}$ — one FFT per cell,
no draws. And **hazard 8**: a pre-registered tolerance means nothing
until the target's own spread is measured. A "$1.5\times$ in 7 of 8
bands" pass was noise from a statistic with **86.9%** per-draw spread
at $R=8$; two runs of the same band differ by $2.4\times$ (#115).

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
- **Entry point**: MEASUREMENTS.md, then `python code/verify_all.py`
  for the measurement corpus, and `python code/verify_propositions.py`
  for the proved statements — M.1–M.3 and P.1–P.7, nine exact checks
  in seconds, each shown -- not merely asserted -- to flip to FAIL
  under a 1e-3 perturbation of one side.
- **The theorems**: THEOREM_A.md → `paper/theorem_A.tex`. Includes
  **Proposition D‴**: Theorem D's loss factor is `≥ 1/max_d|ρ_HL|`,
  hence monotone in the bound on `ρ_HL`, so *every* improvement of
  Huang–Li's Lemma 1 strengthens the no-go and none can weaken it.
  Under RH the conclusion becomes a power of N. The route is not
  blocked by what we cannot prove about `ρ_HL` — it is blocked harder
  the more we know.
- **The wall's deterministic term**: LOCATION_MASK.md — C(N) is not
  mean-zero. Propositions M.1–M.3, the derivation and where it stops,
  and the seven corrections that came out of finding it.
- **The live construction**: TRANSFORM_P.md — the prime-factor split,
  the first transform here with positive margin (Props. P.1–P.4), with
  TRANSFORM_LAB.md as the notebook it came out of.
- **Where our results sit in the record**: LITERATURE.md — provisional
  placement, and the resolution (negative) of whether current
  large-moduli technology already covers what the chain consumes. It
  does not: the object is a correlation Λ(n)μ(N−n), not a single
  function, and the well-factorable machinery has no entry point.
- **The open problem**: CONJECTURE_L.md (amplitude half), with
  AMPLITUDE_ADJUDICATION.md and TECHNIQUE_FORGE.md for why each route
  fails, and C3_REVIEW.md for the one that is open.

*The wall, at final resolution, is a binary correlation that nature
cancels everywhere we can measure and no current technique — including
the routes this program invented and then refuted — can certify. The
dossier of measurements, the unconditional theorems on the demand side,
and the failure coordinates are the deliverable. — the program*
