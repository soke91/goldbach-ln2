# Literature pass (increment 233)

*What this is: eight web searches and two paper fetches, run to place
this program's results against the published record. What it is not: a
literature search. No MathSciNet, no zbMATH, no reference-chasing, no
specialist consulted. Everything below is a **provisional** placement,
and the novelty column should be read as "not found in eight searches",
which is weak evidence.*

## The finding that matters most, and it is not a novelty verdict

**Lichtman, *Primes in arithmetic progressions to large moduli, and
Goldbach beyond the square-root barrier*, arXiv:2309.08522 (2023).**

> "We show the primes have level of distribution 66/107 ≈ 0.617 using
> triply well-factorable weights. … For the Goldbach problem, this is
> the first use of a level of distribution beyond the square-root
> barrier, and leads to the greatest improvement on the problem since
> Bombieri–Davenport from 1966."

Two things follow, in opposite directions.

**It does not collide with us.** The result is for **primes**, gives
**upper bounds** for Goldbach representations, and does not touch
`Σ_{n<N} Λ(n)μ(N−n)`. The Huang–Li route aims at the asymptotic, not an
upper bound. No overlap.

**But its companion result is the closest live technology to what this
chain needs.** The same line of work establishes **level of distribution
3/5 for the Möbius function with triply well-factorable weights**
(Lichtman, *Primes in APs to large moduli II/III*; Maynard's 3/5 for
primes is the antecedent). Huang–Li need `EH_μ(N^{θ′})` for one
`θ′ > 1/2`, and the campaign's own Theorem C shows the weight that
carries the Goldbach content is `w_k = log k`. So there is a concrete,
checkable question sitting between the two:

> **Open, and checkable: does the weight that Huang–Li's consumption
> actually requires fall inside the triply-well-factorable class for
> which level 3/5 is known?**

Caveats, stated so the question is not overclaimed:

- Huang–Li's `EH_μ` is stated with `max_a |·|`, an absolute-value form.
  Well-factorable results are for **signed weighted** sums `Σ_k λ_k
  E(x;k,a)`, which is strictly weaker. The comparison must be against
  what the chain *consumes* (a signed sum — see REVIEW_VERDICT #6, which
  established exactly this), not against the max-form as written.
- Theorem D blocks extraction *by divisor switching*. Lichtman's
  technology is the Deshouillers–Iwaniec spectral large sieve, a
  different mechanism, so Theorem D does not close it.
- `log k` is not obviously well-factorable in any of the senses used;
  the honest state is that nobody here has checked.

**This was the highest-value open item in the program** — higher than
anything in TRANSFORM_LAB, because it asked whether an existing theorem
already covers the gap rather than asking us to invent one.

### Resolved, negatively (increment 235)

Three independent reasons, any one of which settles it.

**1. Object mismatch, and it is decisive.** Well-factorable results —
Lichtman II (arXiv:2006.07088) is for **primes**, level `x^{3/5−ε}`,
extending Bombieri–Friedlander–Iwaniec's `x^{4/7−ε}` — are distribution
statements for a *single* arithmetic function `f(n)` in progressions.
Huang–Li's `EH_μ` is a distribution statement for

> `f(n) = Λ(n)·μ(N−n)`,

a **correlation**, whose second factor's argument is the additive shift
`N−n`. The machinery behind every well-factorable result
(Vaughan/Heath-Brown decomposition → Type I/II sums → Kloosterman sums →
the Deshouillers–Iwaniec spectral large sieve) requires convolution
structure **in n**. `Λ(n)μ(N−n)` has multiplicative structure in two
variables tied by an additive constraint — this program's "two
couplings" — and the machinery has no entry point. Whether the weight is
triply well-factorable is not reached; the question is upstream of it.

**2. Raising the level does not help, and our own Theorem D′ already
priced this in.** For `Λ` of level `θ_E` the loss factor is
`exp(c₁√((1−θ_E)log N))`. At `θ_E = 1/2` that is `exp(c₁√(0.5 log N))`;
at Lichtman's `3/5` it is `exp(c₁√(0.4 log N))`. Both exceed every power
of `log N`. Only `θ_E = 1` exactly would close the gap, and there each
progression holds `O(1)` terms and the statement is vacuous. **A better
level was never the missing ingredient.**

**3. The author's own conclusion, as a sanity check.** Lichtman's
Goldbach application gives **upper bounds** — "the greatest improvement
on the problem since Bombieri–Davenport from 1966", and
Bombieri–Davenport 1966 is an upper-bound constant. If the method
reached `EH_μ` it would give the *asymptotic*, which by our Theorem C is
equivalent to binary Goldbach for large even N. He does not claim that,
and the claim would not be modest. The method demonstrably does not
reach the object.

**Nothing is salvageable elsewhere either.** Transform P's demand is the
same object type (μ against Λ across an additive constraint). Theorem A
already covers every `θ′ > 1/2`; a level of `3/5` would extend it to
`θ′ > 2/5`, which is more than the chain needs and therefore no gain.

**One thing is gained, and it is not small.** The reason the
state-of-the-art technology cannot touch this is *exactly* the
obstruction this program identified independently and calls the two
couplings. That is external corroboration of the negative map, arrived
at from the other direction — and it is the first such corroboration the
campaign has had.

*Cost of resolving it: eight searches and three fetches. Cost of not
checking: months.*

## Provisional placement of our results

| result | genre it belongs to | verdict |
|---|---|---|
| **Theorem A** | Bombieri–Vinogradov + divisor switch moving the work to the short variable — a standard mechanism | **routine as mathematics, plausibly new as a statement.** Its content is specific to a 2020 paper's step; no hit in search. Should be presented as *a lemma about Huang–Li's E₄*, never as a theorem about BV |
| **Theorem C** | `μ ∗ log = Λ`, classical | **bookkeeping.** Its value is directional (the E₃ route is circular), not technical |
| **Theorem D / D′** | **Bombieri's asymptotic sieve** and the parity problem — "no choice of sieve weights extracts primes" | **the strongest thing here, and still expected.** The genre is classical; what is new is making it precise *for this reduction, with BV as the only input, at every level `θ′ < 1`*. **Corrected at increment 278**: this row read "under full EH", which misstates the theorem — Theorem D's ingredient list is Bombieri–Vinogradov at `Q = N^{1/2−δ}`, and full EH is what *Bombieri's asymptotic sieve* assumes, not what Theorem D assumes. An expert would likely call the conclusion unsurprising and the proof necessary |
| **Proposition D″** | monomial weights, `μ ∗ log^D = Λ_D` | routine |
| **Proposition E** | Parseval; `sup ≥ ‖·‖₂` | **an observation.** Most experts would assume it without proof |
| **Transform P, P.1** | **Ramaré's identity** / peel-a-prime-factor: Bombieri's asymptotic sieve, Heath-Brown's identity, Tao's log-weighted device | **presumed known.** Ramaré corrects the multiple counting with `ω_{(P,Q)}(m)+1`; our `log p / log v` normalises instead. A weighted Ramaré split |
| **P.3** (grouping asymmetry) | sieve folklore | no hit; **unverified, plausibly folklore** |
| **ω-alternation** | `μ = (−1)^ω`, Erdős–Kac | **elementary.** Restates the parity problem; not adjudicable numerically |
| **Conjecture L** | Chowla-type randomness, local-mask factorisation | the *mask* half being blind-verified is the real content; the conjecture itself is cheap to state |
| **the (18) defect** | — | **an erratum to a published paper.** Not a novelty question; a factual one, and it stands |

## What this changes in how the program should speak

1. **Theorem A is a lemma, not a theorem.** Rename it in any paper. Its
   ingredients are all classical and its mechanism is standard; what is
   ours is the application and the four corrections that make it hold.
2. **Theorem D should be presented against Bombieri's asymptotic
   sieve**, explicitly, as "the same phenomenon, made precise for the
   Huang–Li divisor switch under EH". Presenting it without that frame
   would read as ignorance of the classical result.
3. **Transform P's P.1 must be attributed** to the Ramaré family. It is
   still the case that the *arrangement* (P.2–P.4, the margin, the
   ceiling) is ours; the identity is not.
4. **Lichtman must be cited** wherever this program says "beyond the
   square-root barrier". It is the state of the art and it postdates
   Huang–Li.

## Theorem D against Bombieri's asymptotic sieve (increment 278)

Increment 233 placed Theorem D "in the genre of Bombieri's asymptotic
sieve" and left it there. Re-reading the two side by side, **they do not
overlap, and the reason is worth recording** — it is the sharpest thing
this program can say about where its own no-go sits.

| | Bombieri's asymptotic sieve | **Theorem D** |
|---|---|---|
| direction | **positive** — computes an asymptotic | **negative** — no weight works |
| level assumed | **1** (full EH) | **θ′ < 1** |
| condition on the weight | must **obey parity** (equal weight to odd and even numbers of prime factors) | *any* `w`, i.e. any `b = μ∗w` |
| what fails, and by how much | parity-violating weights leave an **unknown scalar**; sieve majorants are off by a **factor of 2** | loss is `exp(c√log N)`, later `N^{1/4}` — **beats every power of log** |

Three consequences.

**1. Bombieri's sieve cannot evaluate `C(N)` even at level 1.** `μ` is
the maximally parity-violating weight, so full EH does not deliver
`Σ Λ(n)μ(N−n) = o(N)` through that machine. This is why every thread of
this campaign closed on parity, and it is a *reason*, not a coincidence.

**2. Theorem D is not a restatement of the parity problem.** The parity
obstruction costs a **constant** (Tao: "any upper bounds must be off
from the truth by a factor of 2 or more"). Theorem D costs
`exp(c√log N)`, and by Proposition D‴ a power of N under RH. Different
mechanism — Theorem D's loss comes from Siegel–Walfisz decay of ρ *at
the truncation point*, not from a majorant's blindness to `Ω(n) mod 2`.
The two close the same route for unrelated reasons.

**3. It explains why Huang–Li's hypothesis is not circular.** `EH_μ` is
a statement about **μ** in progressions — parity-sensitive input, and
therefore outside the sieve axioms by construction. That is exactly the
kind of ingredient Tao's survey lists as required to break parity
(bilinear/Type-II information, Vaughan's identity, exceptional
characters). So the chain is not "assume a sieve fact, deduce a
parity-hard conclusion"; the assumption is doing legitimate work. **This
is a point in Huang–Li's favour and is recorded as such** — the campaign
has an interest in the opposite verdict and should say so when the
evidence runs this way.

⚠️ **Limits.** This is a reading of two survey expositions, not of
Bombieri's paper. The "factor of 2" and the level-1 assumption are
quoted from Tao's notes; the comparison table is our construction and
nobody in the area has checked it. It sharpens the *placement* of
Theorem D — it does **not** establish novelty, which remains open.

## What a real check would require

- MathSciNet / zbMATH searches on the actual statements, not keywords.
- Reading Lichtman II and III (arXiv:2006.07088 and the sequel) far
  enough to state the exact weight class for the μ level-3/5 result.
- Reading Ramaré's identity in its original form to confirm the
  relation to P.1 precisely rather than by description.
- Someone who works in the area. Eight searches by a non-specialist
  cannot settle novelty and this document does not claim to.

Sources:
[Lichtman, arXiv:2309.08522](https://arxiv.org/abs/2309.08522) ·
[Primes in APs to large moduli II, arXiv:2006.07088](https://arxiv.org/abs/2006.07088) ·
[Tao, Notes on the Bombieri asymptotic sieve](https://terrytao.wordpress.com/2016/07/17/notes-on-the-bombieri-asymptotic-sieve/) ·
[Tao, Möbius function notes (Ramaré identity)](https://terrytao.wordpress.com/tag/mobius-function/) ·
[Huang–Li, arXiv:2005.03811](https://arxiv.org/abs/2005.03811) ·
[Tao, Open question: the parity problem in sieve theory](https://terrytao.wordpress.com/2007/06/05/open-question-the-parity-problem-in-sieve-theory/)
