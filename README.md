# goldbach-ln2

A computational and analytic map of the Elliott–Halberstam route to the
binary Goldbach conjecture, and of the single scalar that route reduces
the problem to:

```
C(N) = sum_{n<N} Lambda(n) mu(N-n) = o(N)
```

Huang and Li proved that binary Goldbach for large even `N` follows from
`EH_mu(N^theta')` for a single `theta' > 1/2`. This repository asks what
that one sentence costs, from both of its sides, and then studies the
wall itself.

**Net progress toward Goldbach: zero.** What is here is a map, a law for
the wall, one unconditional theorem that removes a Goldbach-neutral half
of the demand, one conjecture any future construction must reproduce,
and one defect report for the source paper.

## Layout

| Directory | Contents | In git |
|---|---|---|
| `v1/` | **Version 1 — frozen.** The paper, and the code and results it cites. | yes |
| `v1_verify/` | Independent re-verification of `v1`, mirroring its structure. Seventeen findings, and a corrected paper. | yes |
| `v1_verify2/` | The second pass. Blind recall against the first, plus the corrections the first pass is the only witness to. | yes |
| `v2/` | The continuation. Open questions `v1` states but does not answer. | yes |
| `lib/goldbach/` | Shared helper package, imported by scripts in every tree. | yes |
| `v1_log/` | The program's own record for `v1`: process documents, exploratory code, uncited results. | **no** |
| `v1_verify_log/` | The re-verification's own record, including disagreements while they are unresolved. | **no** |
| `v1_verify2_log/` | The second pass's own record. | **no** |
| `v2_log/` | The same, for `v2`. | **no** |

The three `_log` trees are a lab notebook. They are on disk and they are not
distributed: `.gitignore` excludes them. Nothing in `v1/` depends on
them, so a clone is complete without them.

`v1/` and `v1_verify/` share one structure:

```
paper/                the document and its companion notes
code/verify/          the reproduction stamps
code/demand/          Theorems A, C, D, D' and Propositions E, D''
code/supply/          E1, Conjecture L, the kill-tests, the C-classes
code/wall/            Propositions V, W, the lemmas, the location mask
results/{...}/        one output file per script, same subdivision
```

## Reproducing

```
python v1/code/verify/verify_all.py     # the core corpus, minutes
python v1/code/verify/verify_deep.py    # the deep-N arm at N ~ 1e8
```

Every script is standalone: `python <path>` runs it, and each prints its
own pre-registered pass/fail criteria and exits nonzero on failure.
Python and numpy on a laptop is the whole requirement.

`v1/PROVENANCE.md` maps every numbered statement in the paper to the
code that verifies it and the result file the figure was read from.

## The rules this program ran under

Stated because they are the reason the negative results are believable,
and because each was adopted after the corresponding mistake was made.

1. **Pre-registration.** Every design's decision rule, including the one
   that would refute the hypothesis under test, was written before its
   computation ran.
2. **Adversarial review in fresh context**, against source papers'
   verbatim lemma hypotheses rather than against summaries of them.
3. **Power before belief.** A threshold means nothing until the spread
   of the quantity it judges has been measured.
4. **Weights before comparisons.** Two summaries of one object are not
   comparable until each one's weight is stated.
5. **A count is not an error bar.** The uncertainty of a mean of
   correlated summands does not fall like `1/sqrt(n)`; see
   Proposition 15 of the paper.
