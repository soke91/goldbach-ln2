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
cancellation across p whatsoever*. Beside it, the random-sign null
taken from the data rather than from a size estimate: `D_p` is a signed
sum of `t_v = μ(v)Λ(N−v)` over `p | v`, so its scale is the square root
of its own second moment, `V_p = Σ_{v<N, p|v} μ²(v)Λ(N−v)²`, computed
in the same loop.

| N | S_abs/triv | S_null/triv | **S_abs/S_null** | S_abs/(N log N) |
|---|---|---|---|---|
| 5·10⁴ | 0.4459 | 0.5048 | **0.8833** | 0.3162 |
| 10⁵ | 0.4200 | 0.4817 | **0.8718** | 0.3029 |
| 2·10⁵ | 0.4033 | 0.4620 | **0.8729** | 0.2916 |
| 4·10⁵ | 0.3828 | 0.4394 | **0.8712** | 0.2784 |

Two things, and the second is the one that matters.

1. `S_abs/triv` **decays** — the fraction of the trivial bound that
   survives without any p-cancellation is falling.
2. `S_abs/S_null` is **flat across a factor 8 in N** (0.883 → 0.871,
   a 1.4% drift). The measured aggregate sits at a constant ≈0.87 of
   its own random-sign prediction. The per-p cancellation is exactly
   square-root strength, neither better nor worse, and stably so.

Asymptotically `V_p ≈ cN log N/p` gives `S_null ≈ √(cN log N)·Σ_{p<N}
log p/√p ≈ 2√c·N√(log N)`, so `S_null/(N log N) ≈ const/√(log N) → 0`.

*(The size heuristic `√(M_p log N)` — `M_p/log N` terms of size `log N`
— lands 5–7% above the exact second moment here, so it would not have
changed the verdict. It is printed beside the null rather than in place
of it; see hazard 4 in CLOSURE_REAUDIT.md.)*

> **Square-root cancellation per p suffices, with √(log N) to spare,
> and no cancellation across p is needed at all.**

That is a positive margin. Theorem D found the divisor switch tight to
zero slack; Proposition E found the circle method with zero margin.
This is the first transform in the campaign where the budget closes
with room left over.

*Stated with its limit.* Over N ∈ [5·10⁴, 4·10⁵] the factor √(log N)
moves by 9%, so the fitted exponents (−0.851 measured, −0.781 null)
separate "decaying" from "flat" and nothing finer; they are not a
confirmation of −1/2. The finding that carries weight is the flatness
of S_abs/S_null, which needs no extrapolation.

### Where the difficulty actually sits (dyadic profile, N = 4·10⁵)

| p range | mass frac | abs frac | mean ρ_p | null ρ_p |
|---|---|---|---|---|
| 2–4 | 0.0370 | **0.0006** | 0.0158 | 0.0097 |
| 32–64 | 0.0482 | 0.0007 | 0.0147 | 0.0421 |
| 1024–2048 | 0.0587 | 0.0079 | 0.1350 | 0.2308 |
| 32768–65536 | 0.0615 | 0.0488 | 0.7937 | 0.8911 |
| 131072–262144 | 0.0630 | 0.0625 | 0.9929 | 0.9974 |
| 262144–4·10⁵ | 0.0638 | **0.0638** | **1.0000** | **1.0000** |

The trivial mass is spread almost perfectly evenly over dyadic ranges
(≈0.06 each, as `Σ log p/p` demands). The *surviving* mass is not: it
climbs monotonically from 0.0006 to 0.0638. **All the difficulty is at
large p**, and for the top range `ρ_p = 1.0000` exactly — those `D_p`
are one-term sums, where cancellation is not merely hard but
unavailable. The null agrees there (`1.0000`, since `√(t²) = |t|`),
which is the check that the exact second moment gets right and the
size heuristic did not: the heuristic returned a null of 1.0968, a
"typical size" larger than the trivial bound, which is impossible.
That region carries mass ≈ 1/log N and vanishes, slowly.

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

## Session 7 — the same identity, grouped the other way

Transform P is a double sum over (p, w) with v = pw. Session 6 grouped
it by p. Grouping by w instead costs nothing and changes everything:

> `C_log(N) = Σ_p log p · D_p(N) = −Σ_w μ(w) · G_w(N)`,
> `G_w(N) = Σ_{p<N/w, p∤w} log p · Λ(N−pw)`

verified exactly (rel. diff ≤ 1.9·10⁻¹³). And `G_w(N)` is a familiar
object: it counts `N = wp + n` in primes. **The wall is a Möbius
average of binary Goldbach problems**, one for each w, with `G_1(N)`
the Goldbach count itself.

### The asymmetry, and it is exact rather than measured

`G_w ≥ 0` termwise. Therefore

> `Σ_w |G_w| = Σ_w G_w = Σ_v μ²(v)Λ(N−v) log v` = the trivial bound,
> **identically**.

| N | w-grouping Σ\|·\|/triv | p-grouping Σ\|·\|/triv |
|---|---|---|
| 5·10⁴ | 1.00000000 | 0.4459 |
| 10⁵ | 1.00000000 | 0.4200 |
| 2·10⁵ | 1.00000000 | 0.4033 |

One identity, two groupings, opposite verdicts. The w-grouping has
**provably** zero margin — every scrap of its cancellation must come
from μ(w) across groups. The p-grouping's `D_p` contains the Möbius,
cancels internally, and that alone suffices.

### What the w-grouping would need

**Its main terms cancel, and completely.** With `𝔖_w(N) = 𝔖(N)f(w)`,
`f(w) = Π_{q|w}(q−1)/(q−2)`, the main-term sum is `𝔖(N)N·A(W)`,
`A(W) = Σ_{w≤W} μ(w)f(w)/w`. Since f is within O(1/q) of 1, A(W) is a
perturbation of `Σ_{w≤W} μ(w)/w` — PNT-strength, not `1/log W`:

| N | A(W) | Σ μ(w)/w | 1/log W | MT_tot | \|MT\|/triv |
|---|---|---|---|---|---|
| 2.5·10⁴ | −9.8·10⁻⁴ | 9.9·10⁻⁴ | 0.106 | −42.9 | 2.4·10⁻⁴ |
| 10⁵ | −2.9·10⁻⁴ | 5.3·10⁻⁴ | 0.092 | −50.7 | 6.1·10⁻⁵ |
| 2·10⁵ | −4.4·10⁻⁴ | −4.9·10⁻⁴ | 0.087 | −155.5 | 8.8·10⁻⁵ |

The main terms cancel from ~N log N down to O(100) — better than
(log N)². **The singular-series average is not the obstruction.**

**And the residual is at its own square-root scale.** With
`E_w = G_w − 𝔖_w N/w`, measured against the second moment of `G_w`
computed in the same loop (`null_w = √(Σ_p t_p²)`):

| N | Σ\|E_w\| | Σ null_w | **ratio** | Σ\|E_w\|/triv |
|---|---|---|---|---|
| 2.5·10⁴ | 52379 | 44236 | **1.184** | 0.2888 |
| 5·10⁴ | 112398 | 89426 | **1.257** | 0.2929 |
| 10⁵ | 220735 | 185640 | **1.189** | 0.2658 |
| 2·10⁵ | 463475 | 377073 | **1.229** | 0.2625 |

Flat at ≈1.2, exactly as the p-grouping was flat at 0.83 against its
own null. Both groupings, once normalised against the right scale, sit
at square-root strength — the difference is entirely in what has to be
supplied from outside to get there.

> **The design rule.** Group so that **μ is inside the group**. Then
> the group cancels internally and the main-term subtraction comes for
> free; group the other way and the group is nonnegative, has zero
> margin as it stands, and only becomes comparable after main terms are
> put in by hand.

### Two nulls of mine were wrong, and both were caught by the data

Recorded as corrections #33 and #34. I predicted `A(W) ~ c/log W` from
the Euler product — the truth is PNT-strength, three orders below. And
I predicted `Σ|E_w| ~ √2·N log N` with a flat ratio ≈1 to trivial — the
measurement is 7× smaller and *decaying*. Both errors have one cause:
**the null came from an a priori size heuristic instead of from the
data**, when the statistic's own second moment was computable in the
same loop. That is now hazard 4 in CLOSURE_REAUDIT.md. Rule 2 said
compute the null first; it did not say the null must not itself be a
guess. Now it does.

## Session 8 — the wall as an alternating series in ω(v), and why numerics stop here

Session 7 left a puzzle. In the p-grouping the surviving mass reaches
100% at large p, and the reason is plain once said: **for v prime,
μ(v) = −1 with no sign variation at all.** Yet those v alone contribute
`−Σ_p log p·Λ(N−p) ≈ −𝔖(N)N`, which dwarfs the measured `C_log`. So the
wall's smallness is not cancellation *inside* any class of v.

On the support of μ, `μ(v) = (−1)^{ω(v)}`, so

> `C_log(N) = Σ_{j≥1} (−1)^j T_j(N)`,
> `T_j(N) = Σ_{v<N, μ²(v)=1, ω(v)=j} Λ(N−v) log v ≥ 0`

— exact (rel. diff ≤ 1.2·10⁻¹³), every term **nonnegative**, and
`Σ_j T_j` is exactly the trivial bound. The wall is an **alternating
series in the number of prime factors**, and `T_j` is the count of
Goldbach representations `N = v + n` with v a squarefree j-almost-prime.
`T_1` is the Goldbach count itself.

### H_POISSON: the shape is refuted, the size would have sufficed

If ω were Poisson the series would collapse to `e^{−2λ}`. Fitted
properly (λ solves `λ/(1−e^{−λ}) = m`, not λ = m):

| N | R = \|alt\|/ΣT | R_pois | λ fit | P(1) Poisson | P(1) measured |
|---|---|---|---|---|---|
| 5·10⁴ | 0.013627 | 0.1733 | 1.7526 | 0.3674 | 0.2240 |
| 10⁵ | 0.001090 | 0.1609 | 1.8271 | 0.3503 | 0.2124 |
| 2·10⁵ | 0.000697 | 0.1516 | 1.8865 | 0.3371 | 0.2002 |
| 4·10⁵ | 0.009532 | 0.1416 | 1.9550 | 0.3224 | 0.1851 |

The profile is **not** Poisson — measured P(1) is 0.185 against 0.322 —
and the true series alternates 10–200× better than a Poisson would.
But note what R_pois is doing: `R_pois ≈ e^{−λ}/(1−e^{−λ}) ≈ C/log N`,
which is `o(1)`, which is **exactly what the wall needs**. So the
Poisson heuristic is not too weak; it is simply not the truth.

### Where the difficulty actually goes

With `T_j⁰ = Σ_{ω(v)=j} μ²(v) log v` (drop the arithmetic weight) and
`r_j = T_j/T_j⁰` the mean level of Λ(N−v) on the class,

> `Σ_j (−1)^j T_j = r̄·Σ_j (−1)^j T_j⁰ + Σ_j (−1)^j (r_j − r̄) T_j⁰`

The first piece is `r̄` times the **plain Möbius sum** `Σ_v μ(v) log v`,
PNT-small and free (measured 265, −530, 4.4, 172 — a few parts in 10⁵
of the trivial bound). Everything the wall asks sits in the second.

| N | r̄ | r_1 | r_2 | r_3 | r_4 | free/triv | drift/triv |
|---|---|---|---|---|---|---|---|
| 10⁵ | 1.2993 | 1.7692 | 1.5296 | 1.1355 | 0.5542 | −0.00083 | −0.00026 |
| 4·10⁵ | 1.2972 | 1.7393 | 1.5626 | 1.1810 | 0.6766 | 0.00006 | 0.00947 |

**The prime level depends strongly on ω** — a factor 2.6 from j=1 to
j=4, much of it the elementary fact that v must be odd and larger j
makes 2 | v likelier. So the wall becomes: *does the ω-dependence of
the prime density near N cancel in alternation?*

### Why that should work — and why we cannot check it here

An alternating sum is a repeated difference operator, so it annihilates
smooth functions of j: for any polynomial f, `Σ_j (−1)^j f(j)λ^j/j!e^{−λ}`
is `e^{−2λ}` times a polynomial in λ. The drift is therefore harmless
unless `r_j` carries a component oscillating like `(−1)^j` — and such a
component *is* a μ-correlation by definition. Replacing `r_j` by its
weighted least-squares fit:

| N | alt | alt(deg-1 fit) | residual | /triv |
|---|---|---|---|---|
| 2·10⁵ | 1231 | −8474 | 9705 | 0.0007 / 0.0048 / 0.0055 |
| 4·10⁵ | 35777 | −14757 | 50533 | 0.0095 / 0.0039 / 0.0135 |

Both pieces are ~0.4–1.4% of the trivial bound and **comparable to each
other and to the wall**. Neither dominates. The annihilation mechanism
needs many ω-classes to bite, and at N = 4·10⁵ there are five.

> **This route is not numerically adjudicable at reachable N.** Its
> controlling parameter is `log log N`, which moves from 2.38 to 2.56
> across the entire factor-8 range we can compute. To double the number
> of ω-classes we would need N of order `10^{10⁵}`. No amount of
> compute fixes this; it is the same species of limit already recorded
> for Theorem C (not testable below ~10¹³), only far more severe.

### Verdict

The ω-reformulation is **exact, clean, and worth having as a
formulation** — it turns the wall into an alternating series of
explicit nonnegative Goldbach-type counts, with the Möbius part free
and the difficulty isolated in the ω-dependence of the prime level.
It is **not a new resource**: the residual is precisely a `(−1)^{ω}`
component of that level, which is the parity obstruction itself, met
head-on rather than circumvented. Recorded as a formulation, not a
route, and explicitly flagged as beyond numerical reach.

## Session 10 — range grouping: the right shape, a constant-factor gain, and a floor

Proposition P.4 takes absolute values one prime at a time. Grouping the
same identity by dyadic ranges `R = [P, 2P)` gives a strictly weaker
sufficient condition,

> `C(N) − Λ(N−1) = Σ_R C_R(N)`,
> `C_R(N) = Σ_{p∈R} log p·D_p(N) = Σ_v μ(v)Λ(N−v)L_R(v)/log v`,
> `L_R(v) = Σ_{p|v, p∈R} log p`,

and for large P both `p` and `w = v/p` run over long prescribed ranges,
so `C_R` is a genuine bilinear form in (p, w) — the shape Type-II
technology is built for, with μ(w) still inside the group as P.3
demands. Two questions: how much does the weaker condition buy, and
does it reach the part that P.4 cannot?

**It buys a constant factor, not an asymptotic one.**

| N | Σ_R\|C_R\|/triv | S_abs/triv | **coh/abs** | Σ_R\|C_R\|/N |
|---|---|---|---|---|
| 5·10⁴ | 0.3046 | 0.4423 | **0.6887** | 0.2379 |
| 10⁵ | 0.2820 | 0.4170 | **0.6764** | 0.2228 |
| 2·10⁵ | 0.2718 | 0.3996 | **0.6804** | 0.2142 |
| 4·10⁵ | 0.2519 | 0.3797 | **0.6633** | 0.1987 |

The gain is flat at ≈0.67 across a factor 8 in N. Range grouping is a
better sufficient condition by a factor 1.5 and by nothing more.

**And it does not reach the hard part.** For `p > N/2` the only v with
`p | v` is `v = p` itself, so

> `C_R = −Σ_{p∈R} Λ(N−p)` — every term negative.

μ(p) = −1 with no sign variation, so no grouping recovers anything
there. Measured `|C_R|/A_R = 1.0000` exactly in the top range (the
pre-registered bug check, passed), and the size matches the predicted
`𝔖(N)N/(2 log N)` to 8–10%:

| N | top \|C_R\| | prediction | share of triv |
|---|---|---|---|
| 10⁵ | 5377 | 5734 | 0.0680 |
| 4·10⁵ | 18853 | 20471 | 0.0598 |

> **The ceiling.** Any route through transform P that takes absolute
> values — P.4's or the range-grouped version — is bounded below by the
> Goldbach count of the top half, `≍ 𝔖(N)N/(2 log N)`. So transform P
> can deliver `C(N) ≍ N/log N` and **not one power of log better**.

That is not fatal: `N/log N = o(N)`, which is exactly what the
Huang–Li equivalence asks for the Goldbach asymptotic. But it settles
what the transform is for. It is a route to the asymptotic formula,
never to a fixed log power, and any document claiming otherwise for it
would be overstating.

**What this closes and what it leaves.** The measurement question about
the large-p side is answered: regrouping does not make the surviving
mass smaller in order, and the residual there is one-signed by
construction. The remaining question is not "is it small" — it is
**"is the range-grouped bilinear form provable"**, which is analysis,
not measurement, and is where the next work has to go.

## Session 11 — the dispersion test: plain Cauchy–Schwarz dies, the diagonal does not

Session 10 handed the question to analysis. The range-grouped form
`C_R = −Σ_w μ(w)H_{R,w}` with `H_{R,w} = Σ_{p∈R} log p·Λ(N−pw)/log(pw)`
is bilinear with μ on one side and Λ on the other, and the standard
manoeuvre is Cauchy–Schwarz in the rough variable w, then expanding the
square — the dispersion method. Whether that manoeuvre leaves margin is
arithmetic, and arithmetic can be computed.

**Pre-registered**, from the shape alone: `Σ_w Λ(N−pw) ≈ N/p`, so the
square is `≈ c N²/log²N`, i.e. `|C_R| ≲ √c·N/log N` — and with ~log N
ranges the per-range demand is exactly `o(N/log N)`. The prediction was
that Cauchy–Schwarz lands *at* the demand with no margin, the same
verdict Proposition E returned for the circle method.

**It lands on the nose.**

| N | Σ_R\|C_R\| | Σ_R CS_R | **CS/N** | CS/triv | CS/Σ\|C_R\| |
|---|---|---|---|---|---|
| 5·10⁴ | 11896 | 71462 | **1.429** | 1.830 | 6.01 |
| 10⁵ | 22283 | 145031 | **1.450** | 1.836 | 6.51 |
| 2·10⁵ | 42831 | 289105 | **1.446** | 1.835 | 6.75 |

Per range `CS_R/N ≈ 0.08` almost uniformly across all 17 ranges, and
`1/log N = 0.082` at N = 2·10⁵ — the bound equals the demand exactly.
In aggregate the manoeuvre costs **1.45 N**, which is 1.83× the
*trivial* bound: plain Cauchy–Schwarz is worse than doing nothing. The
loss over the true `Σ|C_R|` is a factor 6.0 → 6.8, growing slowly, the
same `√K` species measured in REVIEW_VERDICT #6.

**But the diagonal is comfortable, and that is the useful half.**

| N | Σ_R diag | **/N** | diagonal share of CS |
|---|---|---|---|
| 5·10⁴ | 30550 | **0.611** | 0.4275 |
| 10⁵ | 59487 | **0.595** | 0.4102 |
| 2·10⁵ | 114397 | **0.572** | 0.3957 |

The diagonal — the part of the expansion that no cancellation can
remove — sits at 0.57 N and is **falling**. So the damage is not
structural: it is specifically the **positive off-diagonal main terms**
`Σ_w Λ(N−pw)Λ(N−p'w) ≈ 𝔖(p,p')W`, which Cauchy–Schwarz throws into the
bound with a plus sign after the sign information that made `C_R` small
has already been discarded.

> **What a proof must therefore do.** Not bound `|Σ|²` — *compute* it.
> The off-diagonal main terms are explicit singular series, so the
> dispersion method proper applies: subtract the main terms and bound
> only the variance. The diagonal, which survives that subtraction
> untouched, is at 0.57 N and decaying, so there is room for exactly
> that argument and for no cruder one.

This is the first time this campaign has had a **specification** for a
live route rather than a refutation of one. It is not a proof and it is
not a claim that the variance is boundable — bounding it is a real
analytic problem, and the off-diagonal error terms are correlations of
primes along two linear forms on average over the pair, which is where
the difficulty now honestly sits.

*Zero margin, three times.* Proposition E (circle method on C(N)),
REVIEW_VERDICT #6 (Cauchy–Schwarz in the E1 pipeline), and now this.
The pattern is not that absolute values always fail — Proposition P.4
takes absolute values and keeps a margin of 0.38 and falling. It is
that absolute values taken **after squaring** destroy the sign
information the wall runs on, every time, and land back at the trivial
scale. Take them on already-cancelled quantities and they are affordable.

### Not yet read

MRT's *averaged* form of Chowla (`arxiv_1503.05121`) is stated
qualitatively as o(HX). Requirement ③ needs a fixed log-power saving in
the shift-averaged form, and Pilatte's (log x)^{1−c} is for a **fixed**
shift, logarithmically averaged — a different object. Checking the
quantitative strength actually available for the shift-averaged form is
the other open task, and it is independent of ①.
