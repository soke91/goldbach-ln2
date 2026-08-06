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

**This is now the highest-value open item in the program** — higher than
anything in TRANSFORM_LAB, because it asks whether an existing theorem
already covers the gap rather than asking us to invent one.

## Provisional placement of our results

| result | genre it belongs to | verdict |
|---|---|---|
| **Theorem A** | Bombieri–Vinogradov + divisor switch moving the work to the short variable — a standard mechanism | **routine as mathematics, plausibly new as a statement.** Its content is specific to a 2020 paper's step; no hit in search. Should be presented as *a lemma about Huang–Li's E₄*, never as a theorem about BV |
| **Theorem C** | `μ ∗ log = Λ`, classical | **bookkeeping.** Its value is directional (the E₃ route is circular), not technical |
| **Theorem D / D′** | **Bombieri's asymptotic sieve** and the parity problem — "no choice of sieve weights extracts primes" | **the strongest thing here, and still expected.** The genre is classical; what is new is making it precise *for this reduction, under full EH*. An expert would likely call the conclusion unsurprising and the proof necessary |
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
[Huang–Li, arXiv:2005.03811](https://arxiv.org/abs/2005.03811)
