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

**Sub-question 3, which this raises.** Is there an averaging cheaper
than the full N-average that still decouples the shift from the
modulus? An average over a short window of N, over a sub-progression
of N, or over an auxiliary parameter not yet in the problem. K4 killed
the N-average *descent* (trading the N-average for a k-average within
one N); this is a different question — not descent, but a smaller
average that still buys decoupling.

### Not yet read

MRT's *averaged* form of Chowla (`arxiv_1503.05121`) is stated
qualitatively as o(HX). Requirement ③ needs a fixed log-power saving in
the shift-averaged form, and Pilatte's (log x)^{1−c} is for a **fixed**
shift, logarithmically averaged — a different object. Checking the
quantitative strength actually available for the shift-averaged form is
the other open task, and it is independent of ①.
