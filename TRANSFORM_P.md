# Transform P — the prime-factor split of the Goldbach wall

*Formalisation of the candidate found at increment 225 and measured
through increment 230. Everything in §1–§3 is proved; §4 and §5 are
measurement, labelled as such; §5 also states the ceiling the transform
cannot cross, and §6 what is not established.*

Throughout, N is a large even integer,

> **C(N) := Σ_{n<N} Λ(n) μ(N−n) = Σ_{1≤v<N} μ(v) Λ(N−v)**

is the wall (`C(N) = o(N)` is what the Huang–Li reduction asks of us;
see DEPENDENCY_AUDIT.md for what that reduction is and is not verified
to be), and p, q always denote primes.

## 1. The split

**Proposition P.1.** *For every even N ≥ 4,*

> **C(N) − Λ(N−1) = Σ_{p<N} log p · D_p(N)**,
> **D_p(N) := Σ_{v<N, p | v} μ(v) Λ(N−v) / log v**,

*and the identity is exact — there is no error term.*

*Proof.* μ vanishes off the squarefree integers, and a squarefree
v ≥ 2 satisfies `Σ_{p|v} log p = log v` exactly, there being no higher
prime powers to account for. Hence for every v with 2 ≤ v < N and
μ(v) ≠ 0,

> `Σ_{p | v} (log p / log v) = 1`.

Multiply by `μ(v)Λ(N−v)`, sum over 2 ≤ v < N, and interchange the two
finite sums. The v = 1 term of C(N) is `μ(1)Λ(N−1) = Λ(N−1)`, which no
p divides and which therefore stands outside the double sum. ∎

Writing `v = pw` — where `p ∤ w` is automatic on the support of μ — the
same statement reads

> **D_p(N) = − Σ_{w < N/p, p ∤ w} μ(w) Λ(N − pw) / log(pw)**,

so the Möbius variable is dilated while the additive form becomes
`N − pw`, i.e. **Λ restricted to the progression N mod p**. The second
factor is not transformed at all.

*Remark (why weight 1 and not log v).* The identity is more natural
with the weight log v, giving `C_log(N) = Σ_p log p Σ_{v: p|v} μ(v)
Λ(N−v)`. But recovering C(N) from C_log(N) by partial summation costs
`O(N log log N / log N)` from the trivial bound on partial sums, which
caps any conclusion at `C(N) = O(N log log N / log N)` — enough for the
qualitative wall, not for a fixed log power. Dividing the identity by
log v removes the detour at the price of a smooth factor, and P.1 is
the form the chain actually consumes.

## 2. The split is lossless

**Proposition P.2.** *With* `M_p(N) := Σ_{v<N, p|v} μ²(v)Λ(N−v)/log v`,

> **Σ_{p<N} log p · M_p(N) = Σ_{2≤v<N} μ²(v) Λ(N−v)**,

*which is the trivial bound for C(N). Hence the split redistributes the
trivial bound over p and neither enlarges nor shrinks it.*

*Proof.* Identical to P.1 with μ replaced by μ². ∎

The right-hand side is `A(N)·N·(1+o(1))` with `A(N) = Π_{p∤N, p>2}
(1 − 1/(p(p−1)))`, the Λ-weighted squarefree density — the same
constant Huang–Li's Theorem 1 carries. Measured: 0.7809, 0.7901,
0.7878, 0.7888 at N = 5·10⁴ … 4·10⁵ against A(N) = 0.7873.

Losslessness is not automatic and is the first thing this campaign's
earlier transforms failed: the divisor switch loses (Theorem D) and the
circle method has no margin at all on C(N) (Proposition E).

## 3. The grouping asymmetry

P.1 is a double sum over (p, w). Grouping by w instead of by p costs
nothing and gives a **provably** different object.

**Proposition P.3.** *With* `G_w(N) := Σ_{p<N/w, p∤w} log p · Λ(N−pw) /
log(pw) ≥ 0`,

> **C(N) − Λ(N−1) = − Σ_w μ(w) G_w(N)**   *and*
> **Σ_w μ²(w) G_w(N) = Σ_{2≤v<N} μ²(v)Λ(N−v)**,

*so the absolute aggregate of the w-grouping equals the trivial bound
identically. The w-grouping has zero margin.*

*Proof.* The first identity is P.1 rewritten through `v = pw`. For the
second, every `G_w` is a sum of nonnegative terms, so absolute values
change nothing; substituting `v = pw` again, each squarefree v ≥ 2
arises once for each prime factor of v, and `Σ_{p|v} log p/log v = 1`. ∎

`G_w(N)` counts the representations `N = wp + n` in primes, so P.3 says
the wall is a **Möbius average of binary Goldbach problems** — with
`G_1(N)` the Goldbach count itself. Every scrap of the w-grouping's
cancellation must come from μ(w) *across* groups; the p-grouping's
`D_p` contains the Möbius and cancels *inside* each group.

> **Design rule.** Group so that μ is inside the group. Then the group
> cancels internally and the main-term subtraction comes for free.
> Grouped the other way the summand is nonnegative, and P.3 says
> exactly how much that costs: everything.

That the w-grouping is nevertheless not hopeless is worth recording:
its main terms cancel completely — `Σ_w μ(w)𝔖_w(N)N/w = 𝔖(N)N·A(W)`
with `A(W)` tracking `Σ_{w≤W} μ(w)/w`, PNT-strength, measured at ~10⁻³
against `1/log W ≈ 0.09`. The singular-series average is not the
obstruction. But that cancellation has to be supplied from outside,
whereas in the p-grouping it is already inside `D_p`.

## 4. The margin (measured, not proved)

**Corollary P.4.** *If* `Σ_{p<N} log p · |D_p(N)| = o(N)` *then*
`C(N) = o(N)`.

*Proof.* Immediate from P.1, since `Λ(N−1) = O(log N)`. ∎

The point of P.4 is that its hypothesis uses **no cancellation across p
whatsoever**. Writing `S_abs(N) = Σ_p log p |D_p(N)|` and, as the
random-sign null, `S_null(N) = Σ_p log p √(V_p)` with
`V_p = Σ_{v<N, p|v} μ²(v)(Λ(N−v)/log v)²` — the second moment of `D_p`,
computed from the same data, not estimated (`code/lab_transformP_weight1.py`):

| N | S_abs/triv | S_null/triv | **S_abs/S_null** | S_abs/N |
|---|---|---|---|---|
| 5·10⁴ | 0.4423 | 0.4997 | **0.8852** | 0.3455 |
| 10⁵ | 0.4170 | 0.4769 | **0.8743** | 0.3295 |
| 2·10⁵ | 0.3996 | 0.4574 | **0.8735** | 0.3148 |
| 4·10⁵ | 0.3797 | 0.4352 | **0.8725** | 0.2995 |

Two readings, and the second needs no extrapolation.

1. `S_abs/triv` **decays** — the fraction of the trivial bound surviving
   without any p-cancellation is falling.
2. `S_abs/S_null` is **flat** (0.885 → 0.873 across a factor 8 in N).
   The per-p cancellation is exactly square-root strength, stably.

And square-root strength suffices with room. Since `V_p ≈ cN/(p log N)`,

> `S_null ≈ √(cN/log N)·Σ_{p<N} log p/√p ≈ 2√c·N/√(log N) = o(N)`,

so **square-root cancellation per p delivers the wall with `√(log N)` to
spare, using no cancellation across p at all.** Theorem D found the
divisor switch with zero slack and Proposition E found the circle
method with zero margin; this is the first transform in the campaign
whose budget closes with room left over.

**Where the difficulty sits.** The trivial mass spreads almost evenly
over dyadic p-ranges (≈0.06 each, as `Σ log p/p` demands), but the
surviving mass climbs monotonically from 0.0006 at p ∈ [2,4) to 0.0598
at the top. All the difficulty is at large p, and in the top range
`|D_p| = M_p` exactly — those are one-term sums, where cancellation is
unavailable rather than merely hard. That region carries mass ≍ 1/log N
and vanishes, slowly.

## 5. The ceiling: what transform P can and cannot deliver

Grouping by dyadic ranges `R = [P, 2P)` rather than by single primes
gives a strictly weaker sufficient condition than P.4 — `|C_R| ≤
Σ_{p∈R} log p|D_p|` termwise, with `C_R = Σ_v μ(v)Λ(N−v)L_R(v)/log v`
and `L_R(v) = Σ_{p|v, p∈R} log p`. Measured, it buys a **constant
factor 0.67, flat across a factor 8 in N**, and nothing asymptotic.

More importantly, both routes share a floor. For `p > N/2` the only
v < N divisible by p is `v = p`, so on that range

> `C_R = Σ_{p∈R} log p·D_p = −Σ_{p∈R} Λ(N−p)`,

every term negative, and no grouping can recover cancellation from it —
μ(p) = −1 carries no sign variation. Hence

> **Σ_p log p·|D_p(N)| ≥ Σ_{N/2<p<N} Λ(N−p)**,

which is the number of Goldbach representations of N with the prime
exceeding N/2, expected to be `≍ 𝔖(N)N/(2 log N)` and measured within
8–10% of it (18853 against 20471 at N = 4·10⁵).

> **Consequence.** Any argument through transform P that passes to
> absolute values can conclude `C(N) ≍ N/log N` and **no fixed power of
> log better**. That is enough for the Huang–Li equivalence, which asks
> only `C(N) = o(N)` for the Goldbach asymptotic. It is not enough for
> anything that wants `C(N) ≪ N(log N)^{−A}`.

The floor is stated conditionally in one direction: it is `≍ N/log N`
provided N has close to its expected number of Goldbach
representations. Unconditionally the inequality above still holds and
still blocks a fixed log power for any N where Goldbach is not nearly
failing.

## 6. What is not established

- **P.1–P.3 are identities, not estimates.** They prove nothing about
  the wall. What they do is relocate it, losslessly, and P.3 shows that
  the relocation is not neutral: one grouping has margin and the other
  provably has none.
- **The margin is a budget, not an attainment.** "Square-root per p
  suffices" is not "square-root per p is provable". For p = 2 the
  hypothesis of P.4 restricted to that one prime is already
  Goldbach-strength. What has changed is that such p carry an
  asymptotically vanishing share (`O(1)/log N`) of the mass.
- **The measurements are at N ≤ 4·10⁵**, where log N ≈ 12.9. Over that
  range √(log N) moves by 9%, so the decay of `S_abs/triv` separates
  "decaying" from "flat" and nothing finer. The flatness of
  `S_abs/S_null` is the claim that needs no extrapolation.
- **Novelty is unverified.** Splitting a Möbius sum by the prime
  factors of its argument is a standard device; whether this particular
  arrangement, and Proposition P.3 in particular, is in the literature
  has not been checked. See the standing note in STATUS.md.
- **The Huang–Li frame is used, not owned.** That `C(N) = o(N)` is the
  wall is their Theorem 1's equivalence clause. DEPENDENCY_AUDIT.md
  records what has and has not been re-derived here.

## 7. Why this candidate and not another

Sessions 1–5 of TRANSFORM_LAB.md closed two whole classes of transform
and left a three-clause test. A candidate is new if it does **none** of:

| clause | Transform P |
|---|---|
| move a finite pencil vertex to infinity | not a divisor switch — no pencil |
| remove a divisibility restriction *we* inserted | `p \| v` is intrinsic to v; the problem supplies it |
| require the two factors to dilate together | only μ is dilated; Λ is redirected to a progression |

The third clause is the one session 5 paid for: Tao's and
Helfgott–Radziwiłł's machinery both run on `λ(pn)λ(pn+p) =
λ(n)λ(n+1)`, a *simultaneous* dilation, and our μ–Λ pairing admits
none — Λ(pm) is nonzero in 0.3% of the draws where Λ(m) is. Transform P
sidesteps that by never asking Λ to transform.
