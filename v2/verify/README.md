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

It asked whether reducing `v2/paper/` to `deploy/papers/` left behind
any measurement that a surviving sentence rests on. The answer is
mostly no, and the defect it found instead is a different one.

| Verdict | What it changed | Evidence |
|---|---|---|
| **Nothing was dropped that a surviving statement rests on.** All 38 numbered statements crossed, none with an altered statement; of 251 remarks, 201 belong to branches the projection does not raise, and the 11 the projection could rest on were read one by one | — | `pass6/results/coverage.txt`, `pass6/results/stmtdiff.txt` |
| **No withdrawn figure stands as current.** All 271 printed decimals traced to their source context; the 18 whose context carries a withdrawal word were read individually | — | `pass6/results/numbers.txt` |
| **A finding reaches one tree and not the other.** Three confirmed defects had reached neither tree; seven had reached one | all ten applied to both; the gate now checks the signature of each in both trees (G79) | `pass6/results/pass4_landed.txt` |
| One table printed two sampling resolutions in one row, so its claim held at one depth fewer than stated | the declared field's figures stand alone; the ten-times-harder resample is reported separately | `results/lab_cell_singular.txt` |
| An adjective was quoted from the cell that carries no effect | the figure is given at the cell that carries it, and at the other | `results/lab_mask_placebo.txt` |
| A sampling error called "the same throughout" is eight times smaller at one depth | stated by depth | `results/lab_cell_singular.txt` |
| A measurement's printed figures came from a script that does not produce them | the script that does is named and shipped | `pass4/results/a1_bh_vs_b.txt` |
| Three measurements whose unfavourable half was left behind — the threshold's arithmetic dependence, what re-verifying a stamp cost, how large the combined modulus must be | all three carried into the projection | `results/audit_threshold_arithmetic.txt`, `pass2/`, `results/lab_combined_modulus.txt` |
| The gate's first signature for one of the ten matched the corrected form as well as the broken one, so it could never pass | anchored to the position that distinguishes them | `gate/gate.py` |

The pass retracted four of its own findings before reporting, three of
them because it looked for the verification code in the log tree and
concluded from one path that it did not exist.
