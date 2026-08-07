# v1_verify — independent re-verification of v1

Same structure as `v1/`, deliberately. This tree is for checking `v1`'s
claims **without reading `v1`'s code**: a claim reproduced by a second
implementation is worth more than one reproduced by rerunning the first.

```
paper/                notes on what was checked and what was found
code/verify/          re-implementations of the reproduction stamps
code/demand/          re-checks of Theorems 1-6, Propositions 7-8
code/supply/          re-checks of E1, Conjecture 10, the kill-tests
code/wall/            re-checks of Propositions 11-22 and the lemmas
results/              one output file per script, same subdivision
```

## The rule

A re-verification script must be written from the **statement** in
`v1/paper/wall_v1.tex`, not from the script in `v1/code/`. Where it
reaches a different number, the disagreement is the finding, and neither
side is assumed right until it is resolved.

Record each check as: statement, method used here, number obtained,
`v1`'s number, and verdict. A check that merely reruns `v1`'s script is
not a re-verification and does not belong in this tree.

## The corrected paper

`paper/wall_v1_corrected.tex` is `v1/paper/wall_v1.tex` with every
standing finding applied. It is a paper, not an erratum: the affected
passages are **rewritten**, with no banners and no retraction history —
that history lives in `paper/ADVERSARIAL_FINDINGS.md` and in
`v1_log/docs/CLOSURE_REAUDIT.md`, and nowhere else. Nine passages
changed:

| Where | Change |
|---|---|
| Abstract | seventeen closures, nine kill-tested; Proposition 15 described by what it does and does not buy |
| Lemma 13 (`lem:MP`) | restated for the truncated convolution over its full support, with a remark on why the `N <= X` form fails |
| Proposition 15 (`prop:W`) | split into the exact identity and the decoupling; the amplification factor `Gamma ~ N/(A log N)` added; the Chowla inference replaced by what it actually needs |
| §`sec:coin`, "What is measured" | the withdrawn comparison against a measured `-0.18` removed; both normalisations of `S(h)` reported |
| §`sec:coin`, "Where the wall leans" | shift shares tagged with the normalisation they assume |
| §`sec:coin`, "Two controls" | "Proposition W is an identity" narrowed to the exact display |
| `conj:wall` item 3 | tail counts given, with the Poisson bar the `t=5` ratio carries |
| `conj:wall` item 4 | rewritten: the lines are present, the aggregate statistic cannot tell `mu` from a coin, the local one suggests it can, and the attribution is open |
| Proposition 20 (`prop:coh`) | derivation restated on the squared difference, so the three terms no longer cancel it |
| Summary | seventeen closures |

`code/verify/lint_corrected_paper.py` checks it mechanically:
environments and braces balance, every `\ref` resolves, no astral
characters, none of the nine refuted forms survives in the body, and
`v1/paper/wall_v1.tex` is still read-only. It exits nonzero on any of
them.

## Status

Six rounds done. Findings are in `paper/ADVERSARIAL_FINDINGS.md`;
scripts and outputs under `code/` and `results/`. Round 6 cross-checked
every finding against `v1_log/`, the program's internal record.

**Fourteen standing findings, five serious.** Four of my own predictions
were refuted and two findings withdrawn; all recorded. The demand side
survived intact.

| # | Statement | Verdict |
|---|---|---|
| 1 | Lemma 13 (`lem:MP`) | false as stated (LHS/RHS ~ 0.60). `STATUS.md:175` has it right, with no upper limit — **the paper introduced the truncation** |
| 2 | "under that input `rho -> 1`" | invalid; `S(h) = o(1)` is short by a factor `N/log N`. Not noticed anywhere in `v1_log` |
| 3 | the reconstruction `-0.0976` | divides by `W`, not `V`; correct value `-0.12413` |
| 4 | "against a measured `-0.18`, a factor 0.54" | **withdrawn at correction `#106`** and re-asserted in four live documents. A coin reproduces the target — confirmed here, &#124;z&#124; <= 0.94 in every band |
| 5 | "Proposition W is an identity" | overstated; `#86` itself calls the step "the uniform-`u` approximation" |
| 6 | Proposition 20 (`prop:coh`) | derivation cancels itself; conclusion confirmed independently |
| 7 | `conj:wall` item 4 | **withdrawn at correction `#110`** and re-asserted verbatim, including the withdrawn `0.39%`. The guard built for exactly this (`audit_withdrawn_forms.py`) scans the file, has the form registered, and misses it because LaTeX writes `0.39\%` |
| 10 | "ten kill-tested technique designs" | the paper's own table has nine |
| 11 | the `t = 5` tail ratio `0.878` | `#87` records the counts `4 vs 4.6`; the Poisson bar is ±0.43 |
| 12 | statement numbers in `v1/PROVENANCE.md` | off by one from `conj:L` onward; seven citations, all one file, all one direction. `v2/README.md` and this tree inherited it and are fixed |
| 13 | §`sec:R4`'s block ratios | the statistic is not defined; one reading falls by 2.6x (the signature R4 hunts), the other is flat, and neither reaches 2 sigma against its own null. R4's DEAD stands on the lag-1 statistic, not on the blocks |
| 16 | R2's two thresholds | the regression bar is "2x a control whose mean is zero" and decides nothing; the coherent-gain bar sits +5.38 sigma of its own null while the measurement sits +2.78 sigma. DEAD stands on the measurement (z = -0.46); "determinant phases blind" does not |
| 14 | §`sec:c3`(2)'s Heath–Brown table | `A_j`, `D_j`, `z` not defined. The definition is recovered and the conclusion holds to within 1%, but `z`'s rounding alone moves the `J=8` entry by 0.017, so the third decimal is a convention |

Withdrawn by me in round 6: the §3.1 margin point (`v1_log` states it
better) and the K1 `R^2` point (the "half" rests on a separate
residual-energy figure).

**Open, and reopened by round 6:** whether the zeta lines belong to
`mu` or to `Lambda`. `#110`'s aggregate statistic cannot separate them
(reproduced here: ratio 1.42x, 4 of 20 coins at or above real); the
local-background statistic does (real 6/10 against coins 0–3), but
eight coin draws give only `p ~ 0.11`.

Confirmed with no defect found: **Theorem 1**, tested directly for the
first time (`sup_t|T_1(t)|` computed exactly, falls like `N^{-0.346}`
and tracks `e^{-c sqrt(log N)}`); Theorems 3, 5, 6 and Corollary 2 of
`theorem_A.tex` by line-by-line scan; **Proposition 7** (table
reproduced and shown grid-converged); **Proposition 8** (every digit);
**Proposition 11**; **Conjecture 14 item 1**; the `conj:L` band-ratio
and kurtosis stamps.

## Running it

```
python v1_verify/code/verify/verify_all.py            # the gate, 12 rows
python v1_verify/code/verify/lint_numbering.py        # statement numbers
python v1_verify/code/verify/lint_corrected_paper.py  # the corrected paper
```

`verify_all.py` re-derives every finding's load-bearing number at
reduced size, judges each against a pre-registered interval, and exits
nonzero on any failure. `PROVENANCE.md` maps each check to its script
and its result file. The full-size runs are the individual `audit_*.py`.

Independently re-verified on the supply side: K1 (finding 15), K3, R2
(finding 16), R4 (finding 13), C-I, C-III(2) (finding 14). K2 and K4
checked at the level of their thresholds and counts, both calibrated.

Not re-verified: the route adjudications of §7.1, kill-tests R1 and
R3, the representation classes C-II and C-IV, and the reproduction stamps' own
pre-registered intervals. Open items are in
`../v1_verify_log/docs/OPEN_AFTER_ROUND5.md`.
