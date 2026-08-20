# -*- coding: utf-8 -*-
r"""
pass3, blind: conj:L's "exact cells" stamp

TARGET

OPEN.md's wall item 4 lists three stamps of Conjecture [conj:L] still
carried by a single witness.  pass2 took the blind mask prediction.
This pass takes the second, which v1's Table 'Measured support for
Conjecture L' states in one line:

    Exact cells | every viable (v_2,v_3) cell 0.99--1.06

The quantity is not in doubt.  The conjecture says G is, on the
surviving support, a fluctuation with no class structure, and the
repository's own audit of the same conjecture measures that as the
variance ratio E[D^2] / support -- it prints 0.99781 on exact support
and a table of the same ratio by gcd class.  A ratio near 1 is exact
square-root cancellation on the terms that survive; a cell that
departs from 1 is class structure, which is what the conjecture
forbids.

**What is in doubt is the cell.**  "(v_2,v_3)" does not say of what.
The mask is a function of the v_q-data of (N,k), so at least three
readings are natural:

    A   the valuations of k          (v_2(k), v_3(k))
    B   the valuations of the gcd    (v_2(gcd(N,k)), v_3(gcd(N,k)))
    C   the valuations of N          (v_2(N), v_3(N))

Nothing published picks one.  pass2 found the same shape of defect in
the neighbouring stamp -- its number was right and could not be
reconstructed from what it printed -- and pass2's FINDINGS said the
rule about it would not be written from one instance.  This is the
second instance, so this pass is designed to count it either way.

Nothing is read from v1/code/ and nothing is imported from v2/code/.

THE FIELD

D(N,k) = sum over odd primes p <= 2000 of mu(N - p k), on
N = 3*10^6 + 2500 j for j = 0..400 and k = 1..1000.  That is the grid
pass2 reconstructed for the neighbouring stamp; it is used again here
so the two passes sit on the same field, and it is this pass's own
declaration, not a claim about which grid the stamp used.  Pairs with
empty support are excluded -- those are the mask's annihilations, not
cells of G -- and a cell is called viable when at least 100 surviving
pairs land in it.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  Under at least one of the three readings, every viable cell's
      ratio E[D^2]/support lies in [0.99, 1.06], reproducing the
      stamp.
  C2  Under at least one reading it does NOT, so the stamp's range
      depends on a choice the stamp does not publish.
  C3  The readings are not equivalent: at least two of the three
      disagree on how many cells are viable.
  C4  The stamp's band is one-sided about 1 -- 0.99 to 1.06 is 0.01
      below and 0.06 above -- and that asymmetry is a property of the
      quantity, not of the choice: under every reading, the mean
      ratio over viable cells is above 1.

REFUTATION RULE (fixed before the run)

  C1  REFUTED if no reading reproduces the range.  Then either the
      grid matters as much as the cell does, or the stamp is not
      about this quantity at all, and this pass reports that it could
      not find the stamp's measurement rather than that the stamp is
      wrong.  **This is the outcome to expect if the reading space is
      larger than three**, and three is a guess.
  C2  REFUTED if all three readings reproduce the range.  Then the
      choice does not matter, the stamp is reconstructible after all,
      and pass2's defect does not repeat here -- which would mean the
      rule pass2 deferred should stay deferred.
  C3  REFUTED if all three give the same viable-cell count.  Bookkeeping
      only; it gates nothing.
  C4  REFUTED if any reading's mean ratio over viable cells is at or
      below 1.  This one can fail for a reason that is not
      interesting -- a single small cell dragging a mean -- so the
      per-cell tables are printed in full and the median is printed
      beside the mean, and if they disagree the reading recorded is
      "the asymmetry is not resolved by this field", not a verdict on
      the stamp.

  C1 and C2 are exhaustive between them only in the sense that one of
  them must hold; they are written as two predictions because the
  finding is different in each case, and writing only the one I
  expect would let the other pass unremarked.

  NO NULL IS RUN.  A null is not what this needs: the claim is that a
  ratio sits in a published interval, and the interval is the
  standard.  What could make the reading wrong is too few pairs in a
  cell, and that is handled by the viability threshold and by
  printing every cell's count.
"""

import io
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.abspath(os.path.dirname(__file__))
OUT = os.path.join(HERE, "..", "results", "verify_conjL_exact_cells.txt")

NN, NK, PMAX = 401, 1000, 2000
VIABLE = 100
LO, HI = 0.99, 1.06          # published in v1/paper/wall_v1.tex


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def mobius_upto(n):
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes_upto(int(n ** 0.5) + 1):
        p = int(p)
        mu[p::p] *= -1
        mu[p * p::p * p] = 0
    rest = np.ones(n + 1, dtype=np.int64)
    for p in primes_upto(int(n ** 0.5) + 1):
        p = int(p)
        rest[p::p] *= p
        q = p * p
        while q <= n:
            rest[q::q] *= p
            q *= p
    big = np.arange(n + 1, dtype=np.int64) // np.maximum(rest, 1)
    mu[(big > 1) & (mu != 0)] *= -1
    return mu


def vq(n, q):
    n = np.array(n, dtype=np.int64, copy=True)
    v = np.zeros(n.shape, dtype=np.int64)
    while True:
        d = (n % q == 0) & (n != 0)
        if not d.any():
            break
        v[d] += 1
        n[d] //= q
    return v


def build(Ns, ks, ps, mu):
    """D and its support, both (len(Ns), len(ks))"""
    D = np.zeros((len(Ns), len(ks)), dtype=np.int32)
    S = np.zeros((len(Ns), len(ks)), dtype=np.int32)
    for i, N in enumerate(Ns):
        idx = int(N) - np.outer(ps, ks)
        ok = idx >= 1
        w = np.where(ok, idx, 0)
        m = mu[w].astype(np.int32) * ok
        D[i] = m.sum(axis=0)
        S[i] = np.abs(m).sum(axis=0)
    return D, S


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    Ns = np.arange(0, NN, dtype=np.int64) * 2500 + 3_000_000
    ks = np.arange(1, NK + 1, dtype=np.int64)
    ps = primes_upto(PMAX)
    ps = ps[ps > 2]
    mu = mobius_upto(int(Ns.max()))
    D, S = build(Ns, ks, ps, mu)
    live = S > 0
    say("field D(N,k) = sum over %d odd primes p <= %d of mu(N-pk)"
        % (len(ps), PMAX))
    say("  N = 3e6 + 2500 j, j = 0..%d and k = 1..%d" % (NN - 1, NK))
    say("  %d pairs, %d with nonempty support, %d annihilated"
        % (D.size, int(live.sum()), int((~live).sum())))
    say("  pooled ratio E[D^2]/support over the survivors %.5f"
        % (float((D[live].astype(np.float64) ** 2).sum())
           / float(S[live].sum())))
    say("PRINTBOUND verify_conjL_exact_cells 4 0.00005")

    v2N, v3N = vq(Ns, 2), vq(Ns, 3)
    v2k, v3k = vq(ks, 2), vq(ks, 3)
    g2 = np.minimum(v2N[:, None], v2k[None, :])
    g3 = np.minimum(v3N[:, None], v3k[None, :])
    readings = [
        ("A  the valuations of k",
         np.broadcast_to(v2k[None, :], D.shape),
         np.broadcast_to(v3k[None, :], D.shape)),
        ("B  the valuations of gcd(N,k)", g2, g3),
        ("C  the valuations of N",
         np.broadcast_to(v2N[:, None], D.shape),
         np.broadcast_to(v3N[:, None], D.shape)),
    ]

    d2 = D.astype(np.float64) ** 2
    inside, counts, means, meds = {}, {}, {}, {}
    for name, a, b in readings:
        key = name.split()[0]
        say()
        say("%s" % name)
        say("      v_2  v_3      pairs   E[D^2]/support")
        cells = sorted(set(zip(a[live].tolist(), b[live].tolist())))
        rs, ok_all, nv = [], True, 0
        for (x, y) in cells:
            m = live & (a == x) & (b == y)
            n = int(m.sum())
            if n < VIABLE:
                continue
            nv += 1
            r = float(d2[m].sum()) / float(S[m].sum())
            rs.append(r)
            good = LO <= r <= HI
            ok_all &= good
            say("      %3d  %3d   %8d   %.4f%s"
                % (x, y, n, r, "" if good else "   OUTSIDE"))
        say("  viable cells %d of %d occupied; %d outside [%.2f, %.2f]"
            % (nv, len(cells), sum(1 for r in rs
                                   if not (LO <= r <= HI)), LO, HI))
        say("  mean %.4f, median %.4f, range %.4f to %.4f"
            % (np.mean(rs), np.median(rs), min(rs), max(rs)))
        say("SCATTER cells_%s %.4f" % (key, float(np.std(rs))))
        say("SCALES 1")
        inside[key], counts[key] = ok_all, nv
        means[key], meds[key] = float(np.mean(rs)), float(np.median(rs))

    # -------------------------------------------------------------- C1
    say()
    say("C1  does any reading reproduce the stamp's range?")
    hit = [k for k in inside if inside[k]]
    say("  readings whose every viable cell lies in [%.2f, %.2f]: %s"
        % (LO, HI, ", ".join(hit) if hit else "none"))
    c1 = len(hit) > 0
    say("  C1 %s   (cap: at least one)" % ("hold" if c1 else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  does the choice of cell change the answer?")
    miss = [k for k in inside if not inside[k]]
    say("  readings that do not reproduce it: %s"
        % (", ".join(miss) if miss else "none"))
    c2 = len(miss) > 0
    say("  C2 %s   (cap: at least one fails)"
        % ("hold" if c2 else "REFUTED"))

    # -------------------------------------------------------------- C3
    say()
    say("C3  are the readings distinguishable at all?")
    say("  viable-cell counts %s"
        % ", ".join("%s = %d" % (k, counts[k]) for k in "ABC"))
    c3 = len(set(counts.values())) > 1
    say("  C3 %s   (cap: not all equal)"
        % ("hold" if c3 else "REFUTED"))

    # -------------------------------------------------------------- C4
    say()
    say("C4  is the band's asymmetry about 1 a property of the "
        "quantity?")
    for k in "ABC":
        say("  %s mean %.4f, median %.4f" % (k, means[k], meds[k]))
    c4 = all(means[k] > 1.0 for k in "ABC")
    agree = all((means[k] > 1.0) == (meds[k] > 1.0) for k in "ABC")
    say("  C4 %s   (cap: every reading's mean above 1)"
        % ("hold" if c4 else "REFUTED"))
    if not agree:
        say("  mean and median disagree in sign about 1 somewhere, so "
            "the")
        say("  asymmetry is not resolved by this field and no verdict "
            "on the")
        say("  stamp is read from C4")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  C4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, c3, c4)))
    say()
    if c1 and c2:
        say("the stamp is reproducible and not reconstructible: a "
            "reading exists")
        say("that gives its range, and a reading exists that does "
            "not, and the")
        say("stamp does not say which it used. That is the second "
            "instance of")
        say("the defect pass2 found in the neighbouring stamp.")
    elif c1 and not c2:
        say("the choice does not matter -- every reading gives the "
            "stamp's range")
        say("-- so this stamp is reconstructible and pass2's defect "
            "does not")
        say("repeat here.")
    else:
        say("no reading tried here reproduces the range. The finding "
            "is that")
        say("this pass could not locate the stamp's measurement, "
            "not that the")
        say("stamp is wrong: three readings is a guess at a space "
            "whose size")
        say("the stamp does not publish.")

    head = [
        "STATISTIC: the variance ratio E[D^2]/support pooled",
        "           inside each (v_2,v_3) cell, which is 1 under",
        "           exact square-root cancellation on the",
        "           surviving terms, under three readings of what",
        "           the cell indexes.",
        "FIELD: N = 3*10^6 + 2500 j for j = 0..400 and",
        "       k = 1..1000, 401000 pairs, with the field",
        "       D(N,k) = sum of mu(N-pk) over the 302 odd",
        "       primes p <= 2000. The grid is reconstructed",
        "       from the two facts the stamp publishes about",
        "       its own and is this pass's declaration, not",
        "       a claim about which grid the stamp used.",
        "       Annihilated pairs are excluded and a cell is",
        "       viable at 100 surviving pairs.",
        "STATEMENT: Conjecture conj:L's 'exact cells' stamp, v1's",
        "           table: every viable (v_2,v_3) cell 0.99 to 1.06.",
        "           The quantity is the variance ratio",
        "           E[D^2]/support, which is 1 under exact",
        "           square-root cancellation on the surviving terms.",
        "METHOD HERE: the field rebuilt over odd primes p <= %d on"
        % PMAX,
        "           401 N and 1000 k, annihilated pairs excluded, and",
        "           the ratio pooled inside each cell under three",
        "           readings of what (v_2,v_3) indexes: the",
        "           valuations of k, of gcd(N,k), and of N. A cell is",
        "           viable at %d surviving pairs." % VIABLE,
        "REPOSITORY'S NUMBER: the interval 0.99 to 1.06, with no cell",
        "           index, no grid, and no per-cell table published.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
