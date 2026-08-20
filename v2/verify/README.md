# v2_verify — independent re-verification of v2

Same intent as `v1_verify/`: a claim reproduced by a second
implementation is worth more than one reproduced by rerunning the first.
Each pass is written from the **statement**, not from the script behind
it; where it reaches a different number, the disagreement is the finding,
and neither side is assumed right until it is resolved.

```
passN/code/      re-implementations written from the statements
passN/results/   one output file per script
```

Passes 1–3 checked `conj:L`'s stamps and the ladder exponents. Passes 4
and 5 checked the five papers in `deploy/` — one for the mathematics, one
for whether the code computes the quantity the papers define. Pass 6, the
projection audit, has not run.

Each pass's own record — what it looked at first and why, what it did not
look at, what it proposes the gate should check — is a working document
and is not distributed. It is available on request.

---

## What was found

Confirmed verdicts only, and what each changed. `code/` and `results/`
here are the evidence.

### Passes 1–3 — the stamps of `conj:L`

| Verdict | What it changed | Evidence |
|---|---|---|
| The blind mask prediction reproduces | one of the conjecture's four stamps is now independently verified | `pass2/results/verify_conjL_mask_zeros.txt` |
| Ladder exponents agree with an independent implementation, to the integer | the ladder figures stand | `pass1/results/verify_rung_exponents.txt` |
| The exact $(v_2,v_3)$ cells were not reached | the conjecture's evidence is graded in the text, rather than claimed uniformly | `pass3/results/verify_conjL_exact_cells.txt` |

### Pass 4 — the mathematics

| Verdict | What it changed | Evidence |
|---|---|---|
| **The extraction identity is not an identity** — a term of size $N$ had been dropped and its place written as $O(N^{o(1)})$ | the term is restored, and the no-go is restated so that it grants the hard bound and fails anyway | `pass4/results/a4_extract_identity.txt` |
| The polynomial-weight case attributes a binary-Goldbach lower bound to classical sieve theory; the next clause of the same statement denies it | the claim is weakened to the upper bound, which is what the argument needs | internal comparison |
| Two ratio lists were printed with different numerators while the argument requires one | the numerator is held fixed; the other quantity is named separately | `pass4/results/a1_bh_vs_b.txt` |
| A definition omitted the squarefree restriction its own figures use | the restriction is written into the definition | `pass4/results/a1_bh_vs_b.txt` |
| Two margin figures come from a law the same paper says it does not assert | the section states its scope | internal comparison |
| Cell-floor endpoints and one range factor disagree with measurement | corrected | `pass4/results/a6_cellfloor_top.txt` |
| **The defect reported in \[HL\] eq. (18) is correct** — checked against both arXiv versions, and the missing term reproduced to $8\cdot10^{-16}$ at twelve points | no correction needed; the restatement gains a $\mu^2$, and the truncation budget is replaced by the authors' own | `pass4/results/a7_hl18_delta.txt` |

### Pass 5 — the code

| Verdict | What it changed | Evidence |
|---|---|---|
| **No closed route was closed by a wrong calculation.** Every computation behind a negative verdict matches its statement, with margins of $13\times$ to $10^4\times$ | — | `pass5/results/c01`–`c04` |
| A flat sum was indexed from $k\ge2$ where the statement indexes from $k\ge1$; since $H(N;1)=C(N)$, exactly $C(N)$ was missing | eight ratios and one fitted exponent corrected | `pass5/results/c02_flatsum_k1.txt` |
| The five $N$ behind one demand ratio are all $2^a5^b$, so the threshold is constant across them | the sentence now says what field it was measured on | `pass5/results/c03_demand_field.txt` |
| One figure was quoted from a rounding the same sentence rules out | replaced by the admissible one | `results/audit_hb_weight.txt` |
| One identity check is tautologically true and adjudicates nothing | the statement stands on the exhaustive comparison in the same script, and says so | `results/audit_switch_identity.txt` |
| Sieve layer checked exhaustively against trial division across 22 entry points — one harmless mismatch; FFT precision supports every printed digit; nulls preserve the support | — | `pass5/results/c01`, `c04` |
| Two exponents that looked inconsistent belong to different quantities | printed values stand | `pass5/results/c05_exponent_convention.txt` |

### Pass 6 — the projection

Has not run. It asks whether reducing `v2/paper/` to `deploy/papers/`
left behind any measurement that a surviving sentence rests on.
