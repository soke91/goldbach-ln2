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

### Not yet read

MRT's *averaged* form of Chowla (`arxiv_1503.05121`) is stated
qualitatively as o(HX). Requirement ③ needs a fixed log-power saving in
the shift-averaged form, and Pilatte's (log x)^{1−c} is for a **fixed**
shift, logarithmically averaged — a different object. Checking the
quantitative strength actually available for the shift-averaged form is
the other open task, and it is independent of ①.
