# v2 — the continuation

`v1` is frozen. Everything that follows lands here.

```
paper/      the v2 document
code/       new computations
results/    their outputs
```

## The v2 paper

`paper/wall_v2.tex` — **the submission candidate.** It is not `v1` with
patches; it is a narrower paper built only from what survived two
independent adversarial passes (`v1_verify`, `v1_verify2`).

**What it keeps.** The demand side entire (Theorem A and its corollary,
the equivalence, the weight-space no-go, the circle-method margin, the
polynomial-weight closure); the Huang--Li equation (18) defect report,
which was confirmed against the source text; and the wall's *exact*
results — the second moment `prop:V`, the aggregate identity `lem:MP`
in its repaired form, the exact cell floor `lem:cellmom`, the coherence
law `prop:coh`, the scale-invariance obstruction `prop:scaleinv`, and
the two controls. It adds one result the earlier draft did not state:
the amplification `Gamma ~ N/(A log N)` (`prop:W`), which shows Chowla
does not give `Var C ~ V` and that no bound on `|S(h)|` can.

**What it drops, and where that went.** Conjecture `conj:wall` and all
four of its measurements, `sec:floor`'s decay table, the `rho`
reconstruction, the major-arc factors, the autocorrelation excess, and
the margin exponent. §"What this version does not claim" says so in the
paper itself rather than leaving the reader to notice. Those questions
are the ones listed below, and the mask's decay law is added to them.

**The bundle is two files.** `paper/wall_v2.tex` cites
`paper/theorem_A.tex` as its companion note, so the note ships with it
and carries its own copy here (`v1/` is frozen and cannot be edited).
The note's defective one-line mechanism is fixed in **both** places it
appeared — the overview of \S"Notation and the mechanism" and the remark
after the theorem statement. The proof itself was already correct; only
the two summaries of it were not.

That defect turned out to have three homes: `wall_v1.tex`,
`theorem_A.tex`'s overview, and `theorem_A.tex`'s remark. The first fix
missed the third, and the findings lint caught it. Keep the lint in the
gate for that reason.

**Gates.** Two lints, both exiting nonzero on failure, run over both
files:

```
python v2/code/lint_paper_refs.py v2/paper/wall_v2.tex        # 0
python v2/code/lint_paper_refs.py v2/paper/theorem_A.tex      # 0
python v1_verify2/code/verify/lint_corrected_vs_findings.py \
       v2/paper/wall_v2.tex                                   # 0 of 17
python v1_verify2/code/verify/lint_corrected_vs_findings.py \
       v2/paper/theorem_A.tex                                 # 0 of 17
```

For contrast the same findings lint reports 17 of 17 against
`v1/paper/wall_v1.tex`, 16 of 17 against
`v1_verify/paper/wall_v1_corrected.tex`, and 1 of 17 against the frozen
`v1/paper/theorem_A.tex`.

**Before submission.** Four items, none of them mathematical:

1. `\author{[author]}` is a placeholder.
2. The defect report is against arXiv:2005.03811**v2** (2022). Check
   whether a later version already repairs equation (18), and tell the
   authors before publishing a defect report about them.
3. Three of the first pass's corrections are still single-witness — K1's
   reopening, R2's coherent-gain figure, R1's six-draw spread. The paper
   states K1 as *open* and withdraws R1's and R2's quoted precisions, so
   nothing rests on them; do not reintroduce their numbers.
4. One citation was **added** by the repair and has not been checked
   against its source: `prop:E` now attributes
   `\|S_\Lambda\|_1 \gg N^{1/2}` to Vaughan (1988). Verify the reference
   and the exact form of the statement, or weaken the sentence to "an
   unconditional lower bound of this order is known". The old text
   attributed the step to Parseval, which does not cover it.

## What v2 inherits

The four questions `v1` states and does not answer, in the order of how
much each would change the picture.

### 1. Does `rho -> 1`?

`rho = Var C / V` is `1` under random signs, and Proposition 15 makes
`rho - 1` a prime-pair-weighted `mu`-correlation. It does **not**
follow that Chowla's conjecture forces `rho -> 1`, and the arithmetic
of why is the first thing v2 has to carry.

`S(h)` is an average, so `|S(h)| <= 1`, and the coefficient it is
multiplied by is not `O(1)`:

```
Gamma(N) = (sum_{h!=0} c(h)) / V(N) ~ N / (A(N) log N)
```

measured `1.55e3`, `1.85e4`, `3.58e5` at `N = 1e4, 1.6e5, 4e6`, with
`Gamma log N / N -> 1/A = 1.270`. Since `c(h) >= 0`, a hypothesis
`|S(h)| <= eps` yields nothing better than `|rho-1| <= eps Gamma`. So
`S(h) = o(1)` — what Chowla asserts — **is short by a factor
`N/log N`**; the strength required is `S(h) = o(log N / N)`.

Nor does any bound on `|S(h)|` suffice. The absolute budget
`sum_h c(h)|S(h)|/V`, with `S` measured, runs `13.5 -> 30.7` over
`X = 2e4 .. 1.6e5` and grows. Whatever sends `rho` to `1` has to be
**signed cancellation across `h`**, which is a different object from
Chowla's smallness and from the averaged absolute bound of MRT. Naming
what supplies that cancellation is the question.

**And the measurement side is blocked for a second, independent
reason.** Lemma 17 (the coin control) kills the natural estimator:
replacing `mu` by random signs on the same support leaves `V` and every
other input byte-identical, and the centred estimator gives the same
output on both — measured, the real `rho` sits within `0.94` standard
errors of the coin's in every octave. So the often-quoted `rho ≈ 0.810`
is not a measurement of `mu`, and cannot calibrate anything. Any v2
attempt must first exhibit an estimator that *distinguishes* `mu` from
a coin, and show that it does, before quoting a level or a rate. With
one realisation of `mu` this may not be possible at all, and saying so
would itself be a result.

Verification for the two paragraphs above:
`v1_verify/code/wall/audit_propW_chowla_gap.py` and
`v1_verify/code/wall/audit_coin_control_v1claims.py`.

### 2. The mask's decay law

Two things are open and they are different.

- **Which form.** `N^-a` against `(log N)^-b` is not separated over the
  factor `160` in `N` reachable so far. More range, or a derivation.
- **Why the exponent depends on depth.** The exponent rises
  monotonically as fewer small primes divide `N`, at 5 to 30 standard
  errors, and Proposition 22 rules out the mechanism that explains the
  mask's *size*: the singular-series excess is scale-invariant by
  construction and predicts an exponent of zero at every depth. So the
  decay is carried by something not yet identified. The frame handed
  over by that refutation is to fit `dm/se` rather than `dm`, with the
  error bar's own law taken from Proposition 20. Note that Proposition
  20's law is a statement about the **squared difference** of two cell
  averages; `Q_cc/n_c^2` on its own overstates the error bar by a
  factor 2 to 60 depending on the cell, the other two terms of
  Lemma 19 cancelling it to that order.

### 3. Is any part of the wall's spectrum `mu`'s?

The spectral measure is atomic on the rationals `j/q` with
Hardy-Littlewood weights `mu^2(q)/phi^2(q)` — which is `Lambda`'s
structure. `mu`'s own contribution is measured only through the
major-arc deficit (`8.40×` at `q=3`, `15.16×` at `q=5`). Whether
anything beyond that is `mu`'s is open.

The zeta ordinates are the sharp case, and they are open in a specific
way. The aggregate statistic — `R^2` of `G` on
`cos(gamma log N), sin(gamma log N)` for the first ten ordinates —
**cannot** decide it: against Lemma 17's coin the coin mean is
`3.3e-3` against a real `4.6e-3`, and 4 of 20 coin draws reach the
real value. A local comparison, judging each ordinate against control
frequencies from its own neighbourhood, does separate them — 6 of 10
ordinates clear their local 99th percentile against 0 to 3 for each of
eight coins — but eight draws give only `p ~ 0.11`. **The first job
here is cheap: run the local comparison against 200 coins.** If it
holds, it is the first quantity in this program that is `mu`'s and not
`Lambda`'s. Any null that scrambles `G`'s serial structure (a value
permutation) is white and answers a different question; the regressors
sit in the lowest few dozen Fourier bins, where the field keeps its
power.

Verification: `v1_verify/code/wall/audit_zeta_regression_null.py`,
`audit_zeta_local_background.py`, `audit_zeta_coin_local.py`.

### 4. Literature novelty of Theorems 5-6

The no-go over the whole weight space, and its survival under `EH`.
This needs a specialist or a literature search, not a computation, and
until it is done the theorems are stated as ours without a claim of
priority.

## The rules carry over

Pre-registration, adversarial review in fresh context, nulls before
thresholds, weights before comparisons, and a count is never an error
bar. See the root `README.md`.

One more, learned from `v1`. **The paper's statements share a single
counter, and remarks advance it.** `wall_v1.tex` cites only by `\ref`
and is safe; the supporting documents write the numbers by hand and
drifted — `v1/PROVENANCE.md` is off by one for every statement from
`conj:L` onward, and this file inherited three of those. Run
`v1_verify/code/verify/lint_numbering.py`, which resolves the counter
from the source and flags any citation that disagrees, before quoting a
number in a document. The resolved table is in
`v1_verify/results/verify/lint_numbering.txt`.
