# v1_verify2 — a second, independent pass

`v1_verify` ran six rounds and returned seventeen findings. This tree
is the second pass. It exists for three things that the first pass
cannot do for itself.

```
paper/                what was checked and what was found
code/verify/          stamps and lints
code/demand/          Theorems 1-6, Propositions 7-8
code/supply/          E1, Conjecture 10, the kill-tests, the classes
code/wall/            Propositions 11-22 and the lemmas
results/              one output file per script, same subdivision
```

## The three jobs, in order

### Phase 1 — blind recall

Re-verify `v1/paper/wall_v1.tex` from the **statements**, exactly as
`v1_verify/README.md` prescribes, **without opening
`v1_verify/paper/ADVERSARIAL_FINDINGS.md` or any script under
`v1_verify/code/`.** Write the findings to
`paper/FINDINGS.md`. Only then open the first pass and score.

This is the one measurement nobody has: **how much of what is wrong
does one pass of this kind actually catch?** If the second pass
reproduces the first, the process has high recall and the paper is
close to clean. If it returns findings the first pass missed, the miss
rate is the number that matters, and it applies to everything the
first pass declared confirmed.

The first pass's outputs are sealed by hash. Record these in
`paper/FINDINGS.md` before starting, and verify them after, so the
scoring is against a fixed target:

```
5d359db38b00938e8b8de5e1e2537f4e80cdbb825165843a9621339e4db1ff0c
    v1_verify/paper/ADVERSARIAL_FINDINGS.md
454e92e034814525d7e5f3c85013c031888821e5b8aed356686ddadaa034964a
    v1_verify/paper/wall_v1_corrected.tex
```

Do not go looking for seventeen. The count is not a target and knowing
it is already a leak; it is written here because a reader of this file
will see the first pass's commit message anyway.

### Phase 2 — the corrections themselves

`v1_verify/paper/wall_v1_corrected.tex` rewrites nine passages of the
paper and changes numbers in them. **Nobody has checked those.** The
first pass is the only witness to its own repairs, which is exactly
the position `v1` was in before `v1_verify` existed.

Every changed figure needs an independent recomputation. In
particular:

| Where | What was put in its place |
|---|---|
| Lemma 13 | the truncated-convolution form, claimed exact |
| Proposition 15 | `Gamma ~ N/(A(N) log N)`, and the claim that `S(h) = o(1)` is short by `N/log N` |
| §`sec:coin` | `-0.124` and `-0.052` under the two normalisations of `S(h)` |
| `conj:wall` item 3 | tail counts `21441/21463`, `502/503.6`, `4/4.6` |
| `conj:wall` item 4 | the coin and rotation nulls, `6 of 10` local, `0-3` for coins |
| Proposition 20 | `Var/(Q_cc/n_c^2)` between `2` and `60`, and `Var*log N` constant to 1% |
| §`sec:R4` | Definition of the block ratio, and `z = -1.03 .. -0.86` |
| §`sec:c3` | Definition of the Heath-Brown weight, and `0.933 .. 0.977` |
| K1's row | `13 of 63` live columns, and `0.904` on the untruncated field |
| R2's row | `G_1 = 1.513` against `0.994 +- 0.187` |
| R1's row | the six-draw versus forty-draw spread |

A correction that is itself wrong is worse than the error it replaced,
because it now carries a verification stamp.

### Phase 3 — what neither pass has covered

Stated so the gap stays visible:

- kill-test R3 (analytic: does the character transform really deposit
  every component into moduli `>= N^{2/3}`, and is Parseval really the
  obstruction?)
- representation classes C-II and C-IV. The first pass read their
  known defects out of `v1_log` rather than measuring them.
- three of the five route adjudications of §7.1: Lichtman rerun in
  dilate coordinates, the Dirichlet-polynomial fourth moment / Perron
  row, and the partial slices row. The MRT/Lichtman and Tao 2016 rows
  were checked against the sources and hold.
- the reproduction stamps' own pre-registered intervals in
  `v1/code/verify/`. `v1_log/code/audit_stamp_calibration.py` is the
  place to start.
- `v1/paper/theorem_A.tex` was scanned line by line but its numerics
  were never re-run.

## The rules

The first pass's rules carry over and one is added.

1. **Write from the statement, not from the code.** A check that
   reruns someone else's script is not a re-verification.
2. **Pre-register.** Every script states its decision rule and the
   prediction it was written to test, in its docstring, before it
   runs. Four of the first pass's predictions were refuted; that is
   the point of writing them down.
3. **Null before threshold.** A bar quoted as an effect size is not a
   bar. Compute the null's own spread and quote the threshold in it.
4. **Define the statistic.** Findings 13, 14 and 16 of the first pass
   were all one disease: a number quoted with no definition, where two
   defensible readings gave opposite conclusions.
5. **NEW: check whether the design was run.** Finding 15 was a
   regression whose 63 columns were 50 zeros, because a range
   condition emptied them. Before believing any measurement, count how
   much of its design was non-degenerate.

## What already exists and can be reused

Reusing a tool is not reusing a conclusion; these are infrastructure,
and reading them does not leak Phase 1.

```
v1_verify/code/verify/verify_all.py            12-row gate
v1_verify/code/verify/lint_numbering.py        the shared counter
v1_verify/code/verify/lint_corrected_paper.py  structural lint
```

`lint_numbering.py` matters: `v1/PROVENANCE.md` numbers statements one
too high from `conj:L` onward, and both `v2/README.md` and the first
pass inherited it. Resolve numbers from source before quoting them.

## Traps, from the first pass

- The paper's statements share **one counter**, and remarks advance
  it. Cite by `\ref`, never by hand.
- `v1` is frozen and read-only. Record disproofs here; fix in `v2`.
- The kill-tests at `p ~ N/2K` carry no `sqrt N` cut; the type-II
  field does, and it empties `D(v)` for `v >= sqrt N`. That is what
  broke K1.
- `mu` at `2e8` needs a memory-lean sieve. `int8` plus one `bool`
  array, and delete the scratch before the field is built.
- LaTeX escapes defeat plain regexes: the withdrawn-form auditor
  misses `0.39\%` because it looks for `0.39%`.

## Reproducing the environment

Python and numpy only. For the paper, a TeX engine is now available
without a system install:

```
~/.local/tectonic/tectonic.exe -o <outdir> <file.tex>
```

`wall_v1.tex` and `wall_v1_corrected.tex` both compile clean; the
corrected one is 22 pages with no undefined references.

## Status

Empty. Nothing has been re-verified yet. Phase 1 has not started, and
starting it means not reading `v1_verify/paper/` first.
