# v2_verify — the third pass

`v2/paper/wall_v2.tex` is the submission candidate. It has never been
adversarially verified as a document. Its *contents* survived two
passes, but that is not the same claim, and the difference is where the
risk is.

```
paper/                what was checked and what was found
code/verify/          stamps and lints
code/demand/          Theorems A, C, D, D', Corollary B, Propositions E, D''
code/supply/          E1, Conjecture L, the kill-tests, the classes
code/wall/            Propositions V, W, coh, scaleinv; Lemmas MP, cellmom,
                      coin, placebo
results/              one output file per script, same subdivision
```

## Start here: Phase 1 is a blind pass

Re-verify `v2/paper/wall_v2.tex` and its companion
`v2/paper/theorem_A.tex` **from the statements**, with

- `v1_verify/paper/` and `v1_verify/code/`, and
- `v1_verify2/paper/` and `v1_verify2/code/`

**unopened.** Write findings to `paper/FINDINGS.md`, then open both
prior passes and score.

Record these hashes in `paper/FINDINGS.md` before starting and verify
them after, so the scoring is against a fixed target.

**The documents under test:**

```
ada170ba17038ddfdea9f77359b823a42ba469f64fc9d39e71a3f8d325b36d73
    v2/paper/wall_v2.tex
b2a399baeff5564181b81e6233efd4cc7b4bb9d7bc93c96a12d4a2a376931769
    v2/paper/theorem_A.tex
```

**The prior passes, sealed and not to be opened during Phase 1:**

```
5d359db38b00938e8b8de5e1e2537f4e80cdbb825165843a9621339e4db1ff0c
    v1_verify/paper/ADVERSARIAL_FINDINGS.md
454e92e034814525d7e5f3c85013c031888821e5b8aed356686ddadaa034964a
    v1_verify/paper/wall_v1_corrected.tex
0197a4742d2e18dd46b181495869f9401a8b6ea44279d6a9daa1b333640e7b46
    v1_verify2/paper/FINDINGS.md
97f8fec0003f98b28a374104fba21fe34bdfe9fef5195e6352b1f52749f1aab5
    v1_verify2/paper/RECALL.md
b6f1be1ba8f4993c47646fedfa009d9fd9fc85e5716a5427f7f6621b7d876327
    v1_verify2/paper/CORRECTIONS_CHECK.md
```

`v1/` is frozen and is not a target. `v1_log/` is the program's own
record and reading it is a leak during Phase 1, as it was for the second
pass.

## Where the risk actually is

The second pass measured that one adversarial pass catches about
**half** of what is there (`v1_verify2/paper/RECALL.md`). Two passes have
now run, so the inherited material is better covered than anything has
been in this program. **That is not true of the parts the second pass
wrote itself.** Those have exactly one witness, which is the position
`v1` was in before `v1_verify` existed and the position
`wall_v1_corrected.tex` was in before `v1_verify2`.

Ranked, highest risk first:

1. **Passages the second pass authored.** Each is new text, checked by
   nobody:
   - `prop:W`, the amplification `Gamma ~ N/(A log N)` and the claim
     that no bound on `|S(h)|` suffices. The numbers came from the
     *first* pass and were reproduced by the second, but the
     **proposition as stated** is new.
   - `thm:A`'s three-step mechanism, rewritten in `wall_v2.tex` and in
     both places it appears in `theorem_A.tex`.
   - `lem:MP`'s statement in the truncated form, and its remark.
   - `prop:coh`'s rewritten derivation and the "all three terms are
     needed" paragraph.
   - `lem:placebo`'s rewritten statement.
   - `sec:R4`'s block-ratio definition.
   - §"The mask exists" and §"What this version does not claim", both
     written from scratch.
2. **One citation added and never checked against its source**:
   `prop:E` attributes `||S_Lambda||_1 >> N^{1/2}` to Vaughan (1988).
   Verify the reference and the exact form, or the sentence must be
   weakened. This is flagged in `v2/README.md` and is the single
   likeliest defect in the paper.
3. **Three single-witness corrections carried forward**: K1's `13 of 63`
   and `0.904` (which is why K1 is stated as *open*), R2's coherent-gain
   figure, R1's six-draw spread. The paper is written so that nothing
   rests on their numbers, but the *verdicts* do rest on them.
4. **What no pass has covered**: Table `tab:L`'s dilate-pair rows
   (`C_{k,k'}` pair statistics, the `(v_2,v_3)` cells, the blind mask
   stamp, the Wishart spectrum) and four of the five route
   adjudications of `sec:closures`. `conj:L` rests on the first;
   the closure table's top half rests on the second.
5. `v1/paper/e1_proof.tex`, never read by any pass.

## The rules

The first two passes' rules carry over, and two are added.

1. **Write from the statement, not from the code.** A check that reruns
   someone else's script is not a re-verification.
2. **Pre-register.** Every script states its decision rule and the
   prediction it was written to test, in its docstring, before it runs.
   The second pass had six predictions refuted, two of them in the
   paper's favour; that is the point of writing them down.
3. **Null before threshold**, and the null must preserve the structure
   of the field it is a null for. Shuffling a correlated field gives a
   null that is orders of magnitude too tight — that is hazard 11, and
   it is how `conj:wall` item 4 survived one pass.
4. **Define the statistic, and the field.** All three of the second
   pass's contradictions of the first turned on this: what a "cell" is,
   what range a figure was computed on, whether an aggregate was pooled
   or normalised per octave. Both passes silently chose a reading, and
   chose differently. When a statement does not fix these, **compute
   every reading and print them all** rather than picking one.
5. **NEW: check the author's own additions first.** The section above
   ranks them. A paper assembled from verified parts is not a verified
   paper, because the assembly is new work.
6. **NEW: verify every citation that supports a step.** The one added
   citation in this paper has not been checked, and a wrong attribution
   is indistinguishable from a wrong lemma to a reader.

## Scoring

When Phase 1 closes, score against **both** prior passes, not one, and
report:

- how much of `v1_verify` + `v1_verify2` this pass would have found had
  it run on `wall_v1.tex` (the recall series);
- how much of the *new* material it found, which is the number nobody
  has yet;
- and its own retractions. The second pass withdrew two of its own
  findings on comparison and recorded them; do the same.

The recall series so far: **one pass catches ~50%** of what is there.
If the third pass finds little in the inherited material and a lot in
the authored material, that is the expected shape and should be stated
as a confirmation rather than a surprise.

## Reproducing the environment

Python and numpy only. A TeX engine without a system install:

```
~/.local/tectonic/tectonic.exe -o <outdir> <file.tex>
```

Both documents compile clean. Two gates already exist and should be
kept green throughout:

```
python v2/code/lint_paper_refs.py v2/paper/wall_v2.tex
python v2/code/lint_paper_refs.py v2/paper/theorem_A.tex
python v1_verify2/code/verify/lint_corrected_vs_findings.py \
       v2/paper/wall_v2.tex
```

The last one encodes the second pass's seventeen findings and reports
`0` against the current manuscript. It earned its place: the defective
one-line mechanism of `thm:A` had three homes and the first repair
missed the third. **Add this pass's findings to it as they are
confirmed**, so the next manuscript cannot silently reintroduce them.

## Traps, carried forward

- Shell heredocs eat backslashes. Anything containing `.tex`, a regex,
  or an f-string goes through the Write/Edit tools. The second pass hit
  this three times in one session.
- The statements share **one counter**, and remarks advance it. Cite by
  `\ref`, never by hand. `v2/code/lint_paper_refs.py` enforces it, and
  knows that "Huang--Li's Lemma 1" is the source paper's numbering and
  not this document's.
- `mu` at `2e8` needs a memory-lean sieve: `int8` plus one `bool` array,
  and delete the scratch before the field is built.
- The field the wall's figures are computed on is `10^5 < N <= 1.6e7`,
  not "every even `N`". Small `N` dominates several residuals.
- `A(N)` is a product over primes **not** dividing `N`. At
  `N = 4e6 = 2^8 * 5^6` both the `q=2` and `q=5` factors come out of the
  Artin constant; using the even-`N` value alone is wrong by 5%.

## Status

Empty. Nothing has been re-verified yet. Phase 1 has not started, and
starting it means not reading `v1_verify/paper/` or `v1_verify2/paper/`
first.
