# v2 — the continuation

The open questions `v1` states but does not answer, and what has come of
them. This directory is self-contained: no file in it reads from
`v2_log/`, so a clone without the lab notebook is complete.

```
paper/theorem_A.md    the demand side: Theorems A and C, the no-go, the (18) defect
paper/wall_v3.md      the wall: its exact second moment, and the negative map
PROVENANCE.md         statement -> code -> result file, one row each
code/                 one script per measurement
results/              one output file per script
gate/gate.py          the consistency lint the manuscripts are held to
verify/               independent re-verification passes, and what they found
```

## What v2 adds to v1

The demand side is closed at the level of identities. One half is an
unconditional theorem — the Möbius-weighted, fixed-class correlation sum
with weight `w_k = 1` is `<<_A N(log N)^-A` — and the other half is the
same sum with weight `log k`, which an unconditional identity shows to be
equivalent to Huang–Li's equation (22). The interior of the weight space
between those two is empty, and stays empty under the full
Elliott–Halberstam conjecture.

Everything then rests on one scalar,
`C(N) = sum_{n<N} Lambda(n) mu(N-n)`. `wall_v3.md` gives what is exactly
computable about it: its second moment in closed form, an exact identity
for the aggregate second moment, and the reason Chowla's conjecture does
not control the excess — the coefficient amplifying it grows like
`N/log N` and is nonnegative. It also records seventeen pre-registered
routes into the field and what closed each one.

Net progress toward the Goldbach conjecture is zero, and both documents
say so in their own words.

## The gate

`gate/gate.py` is a lint over the manuscripts, not over the mathematics.
It checks that every numbered statement names evidence that exists, that
every printed figure appears in the result file it came from, that one
symbol does not carry two meanings, and that no verdict is printed
without a path by which it could have failed.

```
python gate/gate.py > gate/gate.txt 2>&1; echo $?
```

Each check was added after the corresponding mistake had been made. The
script is the authoritative list — a list copied into prose drifts from
the implementation, and a drifted list is worse than none, because it
looks like it is being kept.

## Reproducing

```
python code/<name>.py > results/<name>.txt 2>&1; echo $?
```

Every script is standalone, prints its own pre-registered pass/fail
criteria, and exits nonzero on failure. Python and numpy on a laptop is
the whole requirement. `PROVENANCE.md` maps each numbered statement to
the script and the result file behind it.

## Verification

`verify/` holds the re-verification passes. Each is written from the
statements rather than from the scripts, so that a second implementation
reaching a different number is a finding rather than a copy.
`verify/README.md` records what was found and what each finding changed.
