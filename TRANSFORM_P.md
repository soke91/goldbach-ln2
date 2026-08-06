# Transform P — the prime-factor split of the Goldbach wall

*Formalisation of the candidate found at increment 225 and measured
through increment 231. §1–§3 are proved; §4 and §5 are measurement,
labelled as such, §5 stating the ceiling the transform cannot cross;
§6 is the analytic status and what a proof would have to do; §7 is what
is not established.*

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

**⚠ The premise fails at deep N, and these are averages over all
even N.** The margin above rests on "square-root cancellation per p suffices", i.e. on `S_abs ≈ S_null`. Measured along the primorial ladder at N ≈ 10⁶ (`code/lab_transformP_depth.py`), with `R = S_abs/Σ log p·M_p` — so `R = 1` is exactly no cancellation:

| rad(N) ⊇ | odd primes | R | R_null | **R/R_null** |
|---|---|---|---|---|
| 2 | 0 | 0.3912 | 0.4078 | **0.959** |
| 2·3·5 | 2 | 0.5508 | 0.4466 | 1.234 |
| 2·3·5·7·11 | 4 | 0.6966 | 0.4862 | 1.433 |
| 2·…·17 | 6 | 0.7817 | 0.5134 | **1.522** |

`R/R_null` climbs from 0.96 to **1.52**: at deep N the per-prime sums cancel **worse than random**, because the location mask gives each `D_p` a systematic component. So the premise of §4's margin — square-root per p — is false exactly where the demand is largest.

**Repeated at X = 2·10⁷ with proper samples** (n = 10 per depth, `code/lab_transformP_depth7.py`; the table above had n = 2 at depth 6), which separates the two axes:

| depth | 0 | 2 | 4 | 6 | 7 |
|---|---|---|---|---|---|
| R at N ≈ 10⁶ | 0.3912 | 0.5508 | 0.6966 | 0.7817 | — |
| **R at N ≈ 1.5·10⁷** | 0.3381 | 0.4738 | 0.5965 | **0.6821** | 0.7069* |
| R/R_null there | 0.954 | 1.220 | 1.422 | **1.537** | 1.576* |

*(The `R_null` column uses `√V_p`, which is the random-sign level only where `D_p` has a single term. A **sign-randomised** null — μ replaced by random signs on the same support — is exact for every p and comes out 5–8% stricter: the depth-6 ratio is **1.584** rather than 1.522, and all even N sits at **0.959**. Correction #55; the readings below are unchanged in direction and slightly larger in size.)*

*(\* n = 1)*. Two facts follow, and together they are the answer.

1. **R decays with N at every fixed depth** — −13.6, −8.6, −14.0, −15.9, −14.4, −14.2, −12.7 percent across a factor 12.5 in N, a mean of **−13.3%**.
2. **R climbs with depth at fixed N**, but the steps are shrinking: +0.055, +0.081, +0.063, +0.060, +0.048, +0.038, **+0.025**.

Along the primorial sequence itself — depth k requires `N ≥ p_k#`, so the two move together — a step from depth 6 to depth 7 multiplies N by 19, which costs **−15.4%**, against a depth gain of **+3.6%**. **The N-decay outruns the depth growth**, and the measured pair agrees: R = 0.7817 at (N ≈ 10⁶, depth 6) falls to 0.7069 at (N ≈ 1.9·10⁷, depth 7). So P.4's demand does decay even along the hardest sequence.

What remains false is the **premise**, not the conclusion: `R/R_null` is flat in N at fixed depth (1.522 → 1.537 at depth 6), so "square-root cancellation per p" is not what carries deep N — something else does, and §4's derivation of the margin from `S_null ≈ 2√c·N/√(log N)` does not apply there.

**And these are averages over all even N, while P.4 is a statement about
every N.** Split by how deep N sits in the location mask
(`code/lab_transformP_demask.py`, N ≈ 10⁶):

| group | S_abs/N |
|---|---|
| shallow, `N = 2q` | 0.2576 |
| **all even N** | **≈0.30** |
| deep, `N = k·30030` | **0.7773** |

The demand for deep N is **three times** the figure §4 reports, and the
all-N average sits near the shallow end because deep N are rare. It is
not a new obstruction — within the deep group the demand declines at
the same relative rate (0.8137, 0.7695, 0.7488 across a factor 3 in N,
against the all-N law `(log N)^{−0.77}`) — but the margin claimed above
is the margin at a typical N, not at the hardest ones, and every
measurement of P.4 before increment 263 averaged over all N.

**De-masking does not help.** Since the location mask is computable, a
weaker sufficient condition would replace `Σ_p log p|D_p|` by
`|Σ_p log p·m_p| + Σ_p log p|D_p − m_p|`, the first term being the mask
and therefore free. Measured with a split estimate (mean from one half
of the group, applied to the other) it makes the demand **worse**:
1.7677 against a raw 0.7773, with a permuted floor of 1.1114. By band
the picture is sharp — the reduction is real at `p ≤ 8192` (0.42 to
0.68) and negative above it (down to −0.66), while the mass sits at
large p. **The mask and the difficulty are in different places.**

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

## 6. Analytic status: what a proof would have to do

The range-grouped form `C_R = −Σ_w μ(w)H_{R,w}`, with
`H_{R,w} = Σ_{p∈R, p∤w} log p·Λ(N−pw)/log(pw)`, is bilinear with μ on
one side and Λ on the other, and the one manoeuvre such a form admits
is Cauchy–Schwarz in the rough variable w followed by expanding the
square — the dispersion method.

**Plain Cauchy–Schwarz has no margin, and is worse than trivial.**
Measured (`code/lab_dispersion_margin.py`), the bound
`CS_R = (Σ_{w≤W} μ²(w))^{1/2}(Σ_{w≤W} H_{R,w}²)^{1/2}` gives

| N | Σ_R CS_R / N | Σ_R CS_R / triv | CS/Σ\|C_R\| |
|---|---|---|---|
| 5·10⁴ | 1.429 | 1.830 | 6.01 |
| 10⁵ | 1.450 | 1.836 | 6.51 |
| 2·10⁵ | 1.446 | 1.835 | 6.75 |

Per range `CS_R/N ≈ 0.08` uniformly, against a per-range demand of
`1/log N = 0.082`: the bound equals the demand exactly, which is the
same zero-margin verdict Proposition E returned for the circle method
on C(N).

**The diagonal is comfortable.** The part of the expansion that no
cancellation can remove measures 0.611, 0.595, 0.572 of N — below the
budget and falling. So the damage is not structural. It is the
**positive off-diagonal main terms** `Σ_w Λ(N−pw)Λ(N−p'w) ≈ 𝔖(p,p')W`,
which Cauchy–Schwarz admits with a plus sign after the sign information
that made `C_R` small has been discarded.

> **Specification.** A proof must not bound `|Σ|²` but *compute* it:
> the off-diagonal main terms are explicit singular series, so the
> dispersion method proper applies — subtract the main terms, bound only
> the variance.

**That was carried out, and it is not enough.** With `h_w =
κ_R·𝔖_w(N)/log(Pw)` — the shape derived from the binary equation
`N = wp + n`, only the constant fitted — writing `C_R = −Σ_w μ(w)h_w −
Σ_w μ(w)(H−h)` and bounding only the second piece gives

| N | Σ_R disp_R / N | plain CS / N | Σ_R \|first piece\| / N |
|---|---|---|---|
| 5·10⁴ | **0.941** | 1.429 | 0.182 |
| 10⁵ | **0.950** | 1.450 | 0.174 |
| 2·10⁵ | **0.940** | 1.446 | 0.165 |

The main term is genuine — it explains 78% of the second moment, and the
first piece is PNT-small and falling, as Proposition P.3's remark
predicted. But the bound is 0.94 N and **flat**, against a demand of
`o(N)`. The variance here is *computed exactly*, not estimated, so no
theorem about it can improve this: the shape of the argument is what
fails.

**The obstruction.** Cauchy–Schwarz over w costs `√(N/P)`, affordable
only for large P; cancellation inside `D_p` needs many v per p, which
happens only for small P. The requirements are opposite and the mass
sits at the crossover. Concretely, the main-term model explains 0.10 of
the second moment at p ∈ [2,4) and 0.83 by p ≍ 10³, while the number of
w available to fit it falls from 6·10⁴ to 1; only p ≍ 10³–4·10³ has
both, and that is two dyadic ranges out of seventeen.

**Net position.** The true value `Σ_p log p|D_p|` is 0.30 N and falling
(§4). The best available argument bounds it by 0.94 N and flat. The gap
between those two numbers is the size of what is missing, and closing it
is not a matter of executing a known method more carefully.

## 7. What is not established

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
- **P.1 should be presumed known.** A first literature pass
  (increment 232, four searches) places it squarely in the
  peel-a-prime-factor family — **Ramaré's identity**, Bombieri's
  asymptotic sieve, Heath-Brown's identity, Tao's log-weighted device.
  Ramaré's version carries a correction factor `ω_{(P,Q)}(m)+1` for the
  multiple counting when v has several prime factors in range; the
  weight `log p / log v` performs the same bookkeeping by normalising
  instead of correcting. **Transform P is a weighted Ramaré split.**
  Proposition P.3 and the margin/ceiling analysis returned no hit, but
  four searches is weak evidence and the honest label is *unverified,
  plausibly folklore*.
- **The Huang–Li frame is used, not owned.** That `C(N) = o(N)` is the
  wall is their Theorem 1's equivalence clause. DEPENDENCY_AUDIT.md
  records what has and has not been re-derived here.

## 8. Why this candidate and not another

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
