# v1 — frozen

The paper and everything it cites. This directory is self-contained: no
file in it reads from `v1_log/`, so a clone without the lab notebook is
complete.

```
paper/wall_v1.tex     the paper
paper/theorem_A.tex   companion note: the full proof of Theorem 1
paper/e1_proof.tex    companion note: the E1 consumable
PROVENANCE.md         statement -> code -> result file, one row each
code/verify/          reproduction stamps
code/demand/          Theorems 1-6, Propositions 7-8
code/supply/          E1, Conjecture 11, the kill-tests, the C-classes
code/wall/            Propositions 12-22, the lemmas, the location mask
results/              one output file per script, same subdivision
```

## What v1 settles

- **The demand side is closed.** One half is an unconditional theorem
  (Theorem 1); the other half is *equivalent* to binary Goldbach
  (Theorem 3). The interior of the weight space between them is empty
  (Theorems 5-6), and the circle method has zero margin
  (Proposition 7).
- **The supply side is closed to measurement.** Eighteen pre-registered
  closures say which directions were checked and why each fails, and
  five structural constraints say what any future technique must
  satisfy.
- **The wall has a law.** Its second moment is exact
  (Proposition 12, local factor `A(N)` and not the singular series);
  its aggregate second moment is an exact identity in two shifted
  correlations (Lemma 13); its bulk and tail are Gaussian under that
  scale and no other; its excess over random signs is a
  prime-pair-weighted Chowla correlation (Proposition 16).
- **The mask is characterised except for its decay.** Shape, cause,
  size and rarity are accounted for; the decay is not, and
  Proposition 22 shows why the mechanism that explains the size cannot
  explain the decay.

## What v1 does not settle

Carried to `v2/`:

1. Does `rho -> 1` — is the wall exactly square-root? With one
   realisation of `mu` it may not be answerable (Lemma 18).
2. The mask's decay law: `N^-a` against `(log N)^-b` is not separated
   over a factor `160` in `N`, and no mechanism is identified for the
   depth-dependence of the exponent.
3. Whether any part of the wall's spectrum is `mu`'s rather than
   `Lambda`'s.
4. The literature novelty of Theorems 5-6 — needs a specialist, not a
   computation.

## Reproducing

```
python v1/code/verify/verify_all.py
python v1/code/verify/verify_deep.py
```

Each script is standalone, prints its own pre-registered criteria, and
exits nonzero on failure.
