# Transform Lab — inventing requirement ①

*Opened at increment 218. Requirement ① of C-III is a legitimate
analytic transform carrying the dilate-averaged μμ correlation to a
shift-averaged one. Closed for pointwise changes of variable (the
pencil-versus-parallel argument); open and unpopulated for summation
formulas. This document is the working notebook for building one.*

*Discipline unchanged: nulls before thresholds, verify rather than
assert, and every failed prediction stays in the record.*

## Session 1 — Helfgott–Radziwiłł, and where our ladder is not theirs

### What their machine is

Helfgott–Radziwiłł (`papers/pdf/arxiv_2103.06853.pdf`) prove that the
graph Γ on V = (N, 2N] with edges

> n ~ n ± p for p ∈ P with **p | n**

is a strong local expander almost everywhere: the eigenvalues of
A = Ad_Γ − Ad_{Γ′} restricted to a density-1 subset are O(√(KL)),
where Γ′ carries every edge n ~ n ± p with weight 1/p regardless of
divisibility. Their Corollary 1.1 is the usable form — the
divisibility-restricted sum equals the unrestricted 1/p-weighted sum,
with a saving.

The engine underneath is **Tao's reduction**: two-point Chowla reduces
to Σ_n λ(n)·(Σ_{p|n, p∈P} λ(n+p)) = o(xL), and this is available
because λ is multiplicative — for p | n, λ(n) = −λ(n/p) — so the
correlation at scale N is tied to the correlation at scale N/p.
Expansion then forces the correlation to equal its own average, which
is zero.

Note precisely what makes the graph interesting: **the move is
additive (n ↦ n ± p) while the index is multiplicative (p | n)**. Tao
himself could not prove expansion for it, calling the edges "neither
random enough nor structured enough"; that mixture is exactly the
difficulty, and exactly the value.

### Our ladder is not that

We have the analogous self-similarity, the exact ladder A1:

> p | m ⟹ μ(m)μ(N − mk) = −μ(m/p)·μ(N − (m/p)(pk)),

and it is exact — verified over 83,592 draws with zero mismatches
(`code/hr_correspondence.py`).

**But two things collapse it.**

1. **The second factor never moves.** N − (m/p)(pk) = N − mk
   identically. So the "ladder" is nothing but μ(m) = −μ(m/p) with the
   companion factor held fixed; the pair value is unchanged apart from
   a sign. (Recorded honestly: the first version of the check reported
   a 10% failure rate. That was the test not enforcing the hypothesis
   p² ∤ m — the failures are exactly the density of p²|m, p³∤m with
   both Möbius factors nonzero. A test fault, not a mathematical one.)
2. **The move is multiplicative in both coordinates.** It sends
   (m, k) ↦ (m/p, pk) and therefore **preserves u = mk**. So the graph
   our ladder generates is not HR's graph on integers with additive
   edges; it is the **divisor lattice of a fixed u**, with edges
   "move one prime factor from m to k". Summing over the whole lattice
   is the complete divisor sum Σ_{m|u} μ(m) = [u = 1] — verified — and
   that is precisely the identity this program already had. Its
   "expansion" yields nothing we do not know.

### The correspondence, stated

| | Helfgott–Radziwiłł | ours |
|---|---|---|
| target | Σ λ(n)λ(n+1) | Σ_m μ(m)μ(N−mk); C(N) = Σ Λ(n)μ(N−n) |
| self-similarity | p\|n ⟹ λ(n) = −λ(n/p) | p\|m ⟹ same summand at (m/p, pk), sign flipped |
| the move | n ↦ n ± p — **additive** | (m,k) ↦ (m/p, pk) — **multiplicative**, fixes mk |
| the index | p \| n — multiplicative | p \| m — multiplicative |
| resulting graph | integers, additive edges, divisibility-indexed | divisor lattice of u |
| what expansion buys | the correlation equals its own average | the complete divisor sum, already known |
| external input | MRT short-interval averages | — |

**The finding of this session**: HR's power comes from the *mixture* —
additive move, multiplicative index. Ours is multiplicative on both
sides, so its graph is trivial and exhausted. **A transplant needs an
additive move, and our ladder does not supply one.**

### Where this points

There is a coordinate in which our object already lives on HR's graph.
Writing u = mk, the switch gives

> Σ_{k∼K} b_k D(k) = Σ_{u<N} μ(N−u)·G_b(u),
> G_b(u) = Σ_{k|u, k∼K} b_k μ(u/k).

The vertex set here is the integers u, and HR's edges u ~ u ± p for
p | u are available on it. What is missing is **an analogue of Tao's
reduction (1.3)**: something that converts the shift-0 correlation
Σ_u μ(N−u)G_b(u) into a p-indexed shifted one. Tao's step used the
multiplicativity of λ at p; our G_b is a *truncated* divisor sum and is
not multiplicative, which is exactly where the analogue has to be
built or shown impossible.

**Open sub-question 1 (the first well-posed question of this phase).**
Is there a reduction of Σ_u μ(N−u)G_b(u) to a sum indexed by p | u,
with the shift acting additively on u? If yes, HR's expansion theorem
applies to our object as it stands. If no, the reason will be the next
map coordinate.

Two things to establish before attacking it, both cheap:
- how far G_b is from multiplicative, measured rather than assumed;
- whether the truncation is what breaks multiplicativity, by comparing
  G_b against its untruncated counterpart, which is μ(u)b-transformed
  and *is* well behaved.

## Session 2 — sub-question 1 answered, and two corrections to our record

### Sub-question 1: no, and the reason repeats the wall

Tao's reduction works because λ is multiplicative in the summation
variable: for p | n, λ(n) = λ(p)λ(n/p), so the shift-1 correlation at
scale N equals the shift-p correlation at scale pN restricted to
multiples of p. Try the same on ours.

Introduce the p-index the standard way, using that Σ_{p∈P}[p|u] has
mean L:

> Σ_u μ(N−u)G(u) ≈ (1/L)·Σ_u Σ_{p|u} μ(N−u)G(u).

Now to produce a *shift* we must move one factor by p. Two attempts,
both forced:

- **Index by p | u.** Then the natural move is u ↦ u/p, and the second
  factor becomes μ(N − u/p) — the argument is *dilated*, not shifted.
- **Index by p | (N−u).** Writing v = N−u and using μ(v) = −μ(v/p),
  the sum becomes −(1/L)Σ_v Σ_{p|v} μ(v/p)·G(N−v): again the Möbius
  argument is dilated while the companion sits at N−v.

**Either way the p-trick produces a dilation, never a shift**, because
the multiplicative structure lives on u and the additive structure on
N−u — they are on different factors. This is the campaign's wall
reappearing at the reduction step rather than at the estimate.

The other half of the obstruction is measured rather than assumed.
Tao's step needs the shifted function multiplicative; the untruncated
G is (b∗μ), multiplicative by construction, so any failure is caused
by the truncation k ∼ K alone. Measured over 4000 coprime pairs
(`code/lab_gb_multiplicativity.py`): G(uv) = G(u)G(v) holds for only
38.2% of pairs at K = 60 and 28.0% at K = 300 — and most of those are
cases where both sides vanish — while corr(G(uv), G(u)G(v)) is −0.140
at K = 60 and +0.035 at K = 300, flipping sign. **The truncation
destroys multiplicativity outright.**

### What is still live in the HR route

Two facts keep it from being closed.

1. **Their Corollary 1.1 needs no multiplicativity.** It bounds
   Σ_n Σ_{p|n} f(n)g(n+σp) against its 1/p-weighted counterpart for
   *any* f, g with |f|₂,|g|₂ ≤ 1 and |f|₄,|g|₄ ≤ e^{CL}. So the
   barrier is not the roughness of our factors; it is purely getting
   the object into the form Σ_n Σ_{p|n} f(n)g(n ± p).
2. **The self-similarity Tao uses is present on our Möbius factor.**
   Index by p | (N−u) and move u ↦ u ± p: then N−u and N−u∓p are both
   divisible by p, and with w = (N−u)/p one gets
   **μ(N−u)·μ(N−u∓p) = μ(w)·μ(w∓1)** — the two-point Chowla
   correlation at a smaller scale, exactly Tao's mechanism. What we
   lack is that our companion is G(u), a divisor-type function of u,
   rather than a second Möbius at a shifted point.

**Sub-question 2.** Is there a form of the wall in which *both*
factors are Möbius at additively related points? The second moment of
C(N) is one (§13 of MEASUREMENTS: Σ_h r(h)S(h)), but it costs the
passage to almost-all-N. Is there another?

### Correction #32 — two things our record had wrong about ③

**(a) MRT's averaged Chowla is not qualitative.** Our documents said
the shift-averaged form is stated as o(HX) with the quantitative
strength unchecked. The paper's own abstract says otherwise: *"Our
arguments in fact give quantitative bounds on the decay rate (roughly
on the order of log log H / log H)"*. So the shift-averaged form
already carries a saving of about (log H)^{−1+o(1)}. Pilatte's
(log x)^{1−c} likewise is a saving of (log x)^{−c} for an absolute
c > 0 — a *fixed power of log*, for a fixed shift, logarithmically
averaged.

**(b) The chain needs a fixed log power, not every A.** Huang–Li's
Corollary 1 proof takes **A = 2 + ε/2** and needs only
1 − A(N) ≫ (log N)^{−2}. So binary Goldbach for large N consumes
EH_μ at one fixed exponent, not for all A. Carrying that through the
corrected E1 target gives an exponent of about 2A+2 ≈ 6.

**What this changes.** ③ was recorded as "a named external open
problem" with an implicitly qualitative gap. It is better described
as an **exponent gap**: the technology delivers fixed log powers of
order 1, and the chain wants a fixed log power of order 6. Whether
those compose — the passage from averaged Chowla through E1 to EH_μ
has losses of its own — is a separate question and is not claimed
here. But "qualitative versus quantitative" was the wrong description
and is withdrawn.

**Caveat on (b), and it is not small.** The exponent 2 + ε/2 is read
off Huang–Li's own Corollary 1 proof, so it is only as good as that
proof — and we have already found one defect in the same paper. See
DEPENDENCY_AUDIT.md: our theorems do not rest on Huang–Li, but every
quantitative statement about what the *chain* needs does, and none of
those has been independently derived here. Treat "③ needs about
(log N)^{−6}" as provisional.

## Session 3 — sub-question 2, and why only the N-average has free shifts

**Sub-question 2 was: is there a form of the wall in which both factors
are Möbius at additively related points?** There is, and grouping by
the shift makes it explicit. In the L² off-diagonal put d = m′ − m and
h = dk; then N − m′k = (N − mk) − h, so with v = N − mk,

> Σ_{k∼K}|D(k)|² − (diagonal)
>  = Σ_{k} Σ_{d≠0} Σ_m μ(m)μ(m+d)·**μ(v)μ(v−h)**, v = N−mk, h = dk.

So both factors are Möbius at points differing by exactly h — a genuine
shifted correlation. **But the shift is h = dk, a multiple of k, while
v runs over the progression N mod k.** Rescaling the progression by k
therefore returns μ(N−mk)μ(N−(m−d)k), the dilate object we started
from. This is the increment-201 circularity again, now visible as a
property of the grouping rather than of one parametrisation.

Contrast the second moment over N (MEASUREMENTS §13):
Σ_{N∈W}|C(N)|² = Σ_h r(h)S(h). There the shift h = n − n′ runs over
prime differences with **no tie to any modulus**, and v runs over a
full interval — which is why that one is genuinely shift-averaged, and
why it is the only coordinate where MRT-type technology applies. Its
price is the passage to almost-all-N.

> **Answer to sub-question 2.** Yes at fixed N, but every fixed-N form
> has its shift locked to the modulus of the progression the other
> variable runs over, so rescaling undoes it. Free shifts appear only
> when the average over N is taken, and that average is exactly what
> costs the exceptional set.

That is a sharper statement of the exceptional-set barrier than
"N-averaging is provable but not consumable": it says *why* — the
N-average is the only place the shift decouples from the modulus.

**Sub-question 3, raised and answered in the same session.** Is there
an averaging cheaper than the full N-average that still decouples the
shift from the modulus?

The natural candidate is a **short window**. Take W of length H and
unfold: S_W(h) = Σ_{N∈W} μ(N−n)μ(N−n′) with n − n′ = h, and as N runs
over W the variable v = N − n runs over an interval of length H. So

> **the second moment over a short window unfolds into shifted μμ
> correlations over short intervals** — and short intervals are exactly
> Matomäki–Radziwiłł's home ground. The shift decouples from any
> modulus, at a cost of H rather than N.

That is the right territory. It does not help, and the reason is
arithmetic rather than technological.

- Actual size: Σ_{N∈W}|C(N)|² ≈ Σ_{N∈W} 𝔖(N)N ≈ 𝔖·N·H.
- Trivial bound: Σ_h r_W(h)|S_W(h)| ≤ (max_h |S_W(h)|)·Σ_h r_W(h)
  = H · N².
- Shortfall = HN²/(𝔖NH) = **N/𝔖, independent of H.**

Both the true size and the trivial bound scale linearly in the window
length, so **no choice of window changes the shortfall**; it is the
same factor N found for the full second moment. Shrinking it would
need per-shift cancellation of size H/N — which is Chowla itself.

> **Answer to sub-question 3.** A short window does decouple the shift
> from the modulus, and it lands the object squarely in short-interval
> territory. But the factor-N shortfall is window-independent, so the
> decoupling buys nothing on its own. Cheap averaging is available;
> what is not available is anything that improves the ratio.

The residual value is a coordinate: **if the shift-decoupling is ever
to pay, it must be combined with something that improves the per-shift
bound**, not with a cleverer choice of average.

## Session 4 — the summation-formula entry to ①, closed for the same reason as K2

Session 2 left one thing alive: HR's Corollary 1.1 needs no
multiplicativity of f and g, only L² and L⁴ bounds, so the barrier is
purely getting our object into the shape Σ_n Σ_{p|n} f(n)g(n ± p).
Session 1's p-trick produced dilations. This session asks whether some
*other* entry works.

It does not, and the reason is one this program has already named.

**What their corollary actually does.** It says a sum *already
carrying* a p | n restriction equals its 1/p-weighted unrestricted
counterpart, with a saving. It is a tool for **removing** a
divisibility restriction that the problem hands you — which is what
happens in Tao's reduction, where the restriction arises from the
multiplicativity step and is not put in by hand.

**Our object carries no such restriction.** Σ_u μ(N−u)G_b(u) has no
p | u anywhere. One can manufacture one, using that Σ_{p∈P}[p|u] has
mean L:

> Σ_u μ(N−u)G(u) ≈ (1/L)·Σ_u Σ_{p|u} μ(N−u)G(u),

and then apply the corollary to remove it again — returning exactly
what we started with. **Manufacture-then-remove is a round trip.**

This is K2's lesson verbatim, from the Technique Forge's first round:
*"collapse structure pays only when it arises inside an intrinsic
average that is already present; imposing it externally only pays the
splitting cost."* K2 established it for a manufactured congruence; it
holds identically for a manufactured divisibility index.

> **Status of ①, complete.** Pointwise changes of variable: closed by
> the pencil-versus-parallel argument. Summation formulas: the named
> candidates were already closed by measurement, and the one live
> entry — HR's Corollary 1.1 — is closed because our object has no
> divisibility restriction for it to remove and manufacturing one is
> circular. **① is now closed in both classes, each for a stated
> structural reason rather than a budget.**

What that leaves is not a gap in the enumeration but the enumeration
itself: a transform outside both classes — neither a pointwise change
of variable nor a summation formula in the usual sense. Inventing one
is the whole task, and this lab now has the boundary drawn tightly
enough that a candidate can be tested against it in one step: **does
it move a finite pencil vertex to infinity, or does it remove a
divisibility restriction the problem supplies rather than one we
inserted?** Anything that does neither is new.

## Session 5 — the first stake: why the machinery cannot see our pairing

With both transform classes closed, the useful question becomes *what
exactly* the working machinery requires that our object lacks. It has
a one-line answer, and it is sharper than "the coupling is bilinear".

**Tao and Helfgott–Radziwiłł both run on a simultaneous dilation.**

> λ(pn)·λ(pn+p) = λ(p)²·λ(n)λ(n+1) = λ(n)λ(n+1)

— the correlation at scale pn *is* the correlation at scale n, because
**both factors transform under the same dilation**. Verified exactly,
200,000 draws, zero mismatches (`code/lab_dilation_pairing.py`).

**Our wall pairs μ with Λ, and only μ dilates.** Λ is not completely
multiplicative: Λ(pm) is log p when m is a power of p and 0 otherwise.
Measured: among 23,124 draws with Λ(m) ≠ 0, Λ(pm) ≠ 0 in **74**
(0.0032). The second factor simply does not survive the dilation, so
no simultaneous dilation exists for C(N) = Σ_v μ(v)Λ(N−v).

**Replacing Λ by μ restores the dilation but moves N.** For the pure
pair, μ(pv)μ(p(N′−v)) = μ(v)μ(N′−v) whenever p divides neither v nor
N′−v — exact, 133,173 draws, zero mismatches. But read what it
relates: the correlation at **N = pN′** restricted to p | v, against
the correlation at **N′**. The dilation moves the additive constraint,
so it acts *across* the N-family rather than within one N — which is
the N-average, and we know its price.

> **The first stake of the invention phase.** A transform that works
> must handle **two factors with different transformation behaviour**.
> Every mechanism that currently breaks a parity-type correlation
> requires the two factors to move together; ours cannot, and forcing
> them to move together (by replacing Λ with μ) converts the problem
> into an N-family statement.

This adds a third clause to the candidate test. A candidate transform
is new if it does **none** of: move a finite pencil vertex to infinity;
remove a divisibility restriction we inserted ourselves; or require the
two factors to dilate together.

*(Recorded: the first run of test (3) enforced only p ∤ v and failed on
5.2% of draws — exactly the density of p | (N′−v), the side condition
the second factor needs. Third time this session a verification script
under-enforced an identity's hypotheses; the rule is now explicit —
when checking an identity with side conditions, enforce all of them.)*

## Session 6 — the prime-factor split: the first candidate with margin

The three-clause test is a filter, and the point of a filter is to run
candidates through it. Here is one that passes.

**The transform.** μ vanishes off the squarefree numbers, so *on the
support of μ* the identity `log v = Σ_{p|v} log p` is exact — there are
no higher prime powers to account for. Feeding it into the wall:

> **Transform P.**  `C_log(N) := Σ_{v<N} μ(v)Λ(N−v) log v = Σ_p log p · D_p(N)`,
> where `D_p(N) := Σ_{v<N, p|v} μ(v)Λ(N−v)`.

An exact identity, no error term. Equivalently, with `v = pw`
(`p ∤ w` automatic on the support of μ), `D_p(N) = −Σ_w μ(w)Λ(N−pw)`:
the Möbius variable is dilated and the additive form becomes `N − pw`,
i.e. **Λ restricted to the progression N mod p**.

**It passes all three clauses.** Not a divisor switch, so no pencil
vertex moves. The restriction `p | v` is intrinsic to v, not one we
inserted — the problem supplies the prime factors. And only μ is
dilated; Λ is *redirected to a progression* rather than asked to
transform. Session 5's stake said a viable transform must handle two
factors with different transformation behaviour: this one does, by not
transforming the second factor at all.

**It is lossless.** `Σ_p log p M_p = Σ_v μ²(v)Λ(N−v) log v` exactly,
where `M_p := Σ_{v<N, p|v} μ²(v)Λ(N−v)` is the trivial mass of `D_p`
(verified to 1.5·10⁻¹⁶). The split neither gains nor loses at the
trivial scale — it redistributes the same mass over p. Both closed
classes lost here; this one does not.

### The measurement: does it leave room?

Lossless is necessary, not sufficient. The decisive quantity is the
**absolute aggregate** `S_abs(N) = Σ_p log p |D_p(N)|`, which uses *no
cancellation across p whatsoever*. Beside it, the random-sign null:
`D_p` runs over ≈ `M_p/log N` prime-power terms of size ≈ log N, so
`null_p = √(M_p log N)` and `S_null = Σ_p log p √(M_p log N)`.

| N | S_abs/triv | S_null/triv | **S_abs/S_null** | S_abs/(N log N) |
|---|---|---|---|---|
| 5·10⁴ | 0.4459 | 0.5380 | **0.8287** | 0.3162 |
| 10⁵ | 0.4200 | 0.5097 | **0.8239** | 0.3029 |
| 2·10⁵ | 0.4033 | 0.4872 | **0.8277** | 0.2916 |
| 4·10⁵ | 0.3828 | 0.4615 | **0.8294** | 0.2784 |

Two things, and the second is the one that matters.

1. `S_abs/triv` **decays** — the fraction of the trivial bound that
   survives without any p-cancellation is falling.
2. `S_abs/S_null` is **flat to three digits across a factor 8 in N**
   (0.8287 → 0.8294). The measured aggregate sits at a constant 0.83 of
   its own random-sign prediction. The per-p cancellation is exactly
   square-root strength, neither better nor worse, and stably so.

Asymptotically `M_p ≈ cN/p` gives `S_null ≈ √(cN log N)·Σ_{p<N} log p/√p
≈ 2√c·N√(log N)`, so `S_null/(N log N) ≈ const/√(log N) → 0`.

> **Square-root cancellation per p suffices, with √(log N) to spare,
> and no cancellation across p is needed at all.**

That is a positive margin. Theorem D found the divisor switch tight to
zero slack; Proposition E found the circle method with zero margin.
This is the first transform in the campaign where the budget closes
with room left over.

*Stated with its limit.* Over N ∈ [5·10⁴, 4·10⁵] the factor √(log N)
moves by 9%, so the fitted exponents (−0.851 measured, −0.862 null)
separate "decaying" from "flat" and nothing finer; they are not a
confirmation of −1/2. The finding that carries weight is the flatness
of S_abs/S_null, which needs no extrapolation.

### Where the difficulty actually sits (dyadic profile, N = 4·10⁵)

| p range | mass frac | abs frac | mean ρ_p | null ρ_p |
|---|---|---|---|---|
| 2–4 | 0.0370 | **0.0006** | 0.0158 | 0.0101 |
| 32–64 | 0.0482 | 0.0007 | 0.0147 | 0.0438 |
| 1024–2048 | 0.0587 | 0.0079 | 0.1350 | 0.2402 |
| 32768–65536 | 0.0615 | 0.0488 | 0.7937 | 0.9223 |
| 131072–262144 | 0.0630 | 0.0625 | 0.9929 | 1.0359 |
| 262144–4·10⁵ | 0.0638 | **0.0638** | **1.0000** | 1.0968 |

The trivial mass is spread almost perfectly evenly over dyadic ranges
(≈0.06 each, as `Σ log p/p` demands). The *surviving* mass is not: it
climbs monotonically from 0.0006 to 0.0638. **All the difficulty is at
large p**, and for the top range `ρ_p = 1.0000` exactly — those `D_p`
are one-term sums, where cancellation is not merely hard but
unavailable. That region carries mass ≈ 1/log N and vanishes, slowly.

### What this hands the next session

For `p > √N` the sum is `D_p = Σ_{w<N/p} μ(w)Λ(N−pw)` with `w` short
and `p` long. Reorganised by the short variable, the hard part of the
wall is

> `Σ_{w<√N} μ(w)·G_w(N)`,  `G_w(N) = Σ_p log p · Λ(N − wp)`,

a **μ-average of binary problems**: each `G_w` counts `N = wp + n` in
primes and is positive with main term `≈ 𝔖_w(N)·N/w`. So the large-p
part of the wall is not a cancellation-in-Λ question at all — it asks
whether **`Σ_w μ(w)𝔖_w(N)/w` cancels**, an Euler product over a
multiplicative-in-w singular series. That is computable, and it is the
next measurement.

### Not yet read

MRT's *averaged* form of Chowla (`arxiv_1503.05121`) is stated
qualitatively as o(HX). Requirement ③ needs a fixed log-power saving in
the shift-averaged form, and Pilatte's (log x)^{1−c} is for a **fixed**
shift, logarithmically averaged — a different object. Checking the
quantitative strength actually available for the shift-averaged form is
the other open task, and it is independent of ①.
