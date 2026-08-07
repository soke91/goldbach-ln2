# v1_verify — independent re-verification of v1

Same structure as `v1/`, deliberately. This tree is for checking `v1`'s
claims **without reading `v1`'s code**: a claim reproduced by a second
implementation is worth more than one reproduced by rerunning the first.

```
paper/                notes on what was checked and what was found
code/verify/          re-implementations of the reproduction stamps
code/demand/          re-checks of Theorems 1-6, Propositions 7-8
code/supply/          re-checks of E1, Conjecture 11, the kill-tests
code/wall/            re-checks of Propositions 12-22 and the lemmas
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

## Status

Empty. Nothing has been re-verified yet.
