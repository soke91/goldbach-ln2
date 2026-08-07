# Notation — one symbol, one meaning

*Increment 321.* This file exists because a symbol in this program
denoted two different quantities, in the same document, with no marker
between them. That is the same species as the `𝔖` / `𝔄` collision of
corrections #74, #75 and #83, where the singular series of the *count*
was substituted for the local factor of the *noise*.

## ⚠️ The collision that prompted this file

| symbol | meaning | typical size at `N = 10⁸` | where |
|---|---|---|---|
| `ρ_HL(x)`, Huang–Li's Lemma 1 | the character-sum quantity bounded by `≪ e^{−c√log x}` | **`0.0137`** | Theorem D, D‴, `THEOREM_A.md`, `paper/theorem_A.tex` |
| `ρ`, the cancellation ratio | `Var C(N) / V(N)`, how far the wall falls short of its exact second moment | **`0.810 ± 0.018`** | Proposition W, Conjecture L, #121, #309, #312, #318 |

**They differ by two orders of magnitude and are unrelated.** Read
together, `STATUS.md`'s "the loss factor is `≥ 1/max_d|ρ|`" and its
"`ρ = 0.810`" give a loss factor of `1.23`, which is nonsense; the true
figure is `≥ 73` at `N = 10⁸`, matching the `exp(√log N) = 73.11`
recorded independently at increment 320.

**Rule.** Huang–Li's is always written `ρ_HL` or with its argument,
`ρ_n(x)`; a bare `ρ` always means `Var C / V`. `code/lint_docs.py`
check (H) enforces it.

## The rest of the alphabet

| symbol | meaning |
|---|---|
| `C(N)` | the wall, `Σ_{n<N} Λ(n) μ(N−n)`. `C(N) = o(N)` is equivalent to binary Goldbach (Huang–Li, Theorem C) |
| `V(N)` | the wall's **exact** second moment, `Σ_v μ²(v) Λ(N−v)²`. Not a fitted stand-in (Proposition V) |
| `W(N)` | `Σ_{w<N} Λ(w)²`, `∼ N(log N − 1)` |
| `Q(N)` | `Σ_{m<N} μ²(m)`, `∼ 6N/π²` |
| `𝔄(N)` | `∏_{q∤N}(1 − 1/(q(q−1)))`, the local factor of the **noise**. `V = W·𝔄(1+o(1))` |
| `𝔖(N)` | the singular series of the **count**. ⚠️ **Not** `𝔄`; substituting one for the other is corrections #74, #75, #83 |
| `𝔖₂(h)` | the twin-prime singular series at shift `h`; the shape of the coin's cell correlation (#116) |
| `m(N)` | the location mask, the deterministic part of `C(N)` (M.1) |
| `G(N)` | `(C − m)/√V`, Gaussian in bulk and tail (Conjecture L) |
| `D(k)` | `Σ_{√N < m ≤ N/k} μ(m) μ(N−mk)`, the type-II dilate field (E1) |
| `M(h)`, `P(h)` | the Möbius autocorrelation and the prime-pair count in `Σ_N C² = Σ_h M(h)P(h)` (#155) |
| `A` | the log-power exponent in the target `N(log N)^{−A}`. ⚠️ Not `𝔄` |
| `θ_E`, `θ′` | the level of distribution and the truncation exponent in `EH_μ(N^{θ′})` |
| `d` | the modulus of the progression split in design K2, costing a factor `d` by Cauchy–Schwarz |
| `B` | in `verify_all`/`lab_*`, the **between-cell variance**; in Forge R4, the **block size**. Two meanings, disjoint contexts |

## What this file is not

It is not a definition list for the mathematics — each object is
defined where it is used. It is a register of **symbols that have
collided or could**, kept so that a statement about one is never read
as a statement about the other.
