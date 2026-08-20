# -*- coding: utf-8 -*-
r"""
pass3, second step: is the prime range the parameter the stamp omits?

WHAT HAPPENED

verify_conjL_exact_cells.py asked whether conj:L's "exact cells" stamp
-- every viable (v_2,v_3) cell 0.99 to 1.06 -- reproduces under any of
three readings of what (v_2,v_3) indexes.  **C1 was REFUTED: none of
the three does.**  Every reading lands mostly below the interval, with
means 0.9844, 0.9750 and 0.9674, and the pooled ratio over all
survivors is 0.98811.

That last number is the one to chase.  The repository's own audit of
the same conjecture prints, for the same quantity on the same kind of
field, **0.99781**, and this pass got 0.98811 on the grid pass2
reconstructed.  A pooled ratio is not a cell choice: it does not
depend on how (v_2,v_3) is read.  So something upstream of the cells
differs, and the field's only unpublished degree of freedom left is
the range the primes p run over.

pass2 found the neighbouring stamp's grid offsets unpublished.  This
step asks whether the prime range is unpublished in the same way --
whether the stamp's number is a function of a parameter it does not
print.

The scan is bounded by arithmetic, not by budget: N - p k must stay
positive, and with N around 3*10^6 and k up to 1000 that caps p below
about 3000.  So the scan runs p <= 250, 500, 1000, 1500, 2000, 2500
and 2900 on the same grid.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The pooled ratio is not constant in the prime range: it moves by
      more than 0.005 across the scan.
  E2  It is monotone in the range.
  E3  The repository's 0.99781 lies inside the interval the scan
      sweeps, so a prime range exists that reproduces it.
  E4  And the stamp's floor 0.99 is reached inside the scan, so the
      stamp's interval is attainable on this grid once the range is
      chosen.

REFUTATION RULE (fixed before the run)

  E1  REFUTED below 0.005 of movement.  Then the prime range is not
      the missing parameter and the disagreement is elsewhere -- most
      likely in the grid, which pass2 already showed is unpublished,
      and this pass would then have found nothing new.
  E2  REFUTED by any reversal.  A non-monotone ratio would mean the
      range interacts with something else and that a single number
      cannot be recovered by naming it.  This is a real possibility
      at the small end, where the support per pair is thin.
  E3  REFUTED if 0.99781 is outside the swept interval.  Then no
      prime range on this grid gives the repository's pooled number
      and the difference is not the range after all -- and the honest
      statement is that this pass has narrowed the missing parameter
      to "not the range", which is weaker but still a result.
  E4  REFUTED if the ratio never reaches 0.99.  Same reading as E3,
      one step milder.

  **What none of these can show** is that the repository used any
  particular range.  A parameter that reproduces a number is not a
  parameter that was used.  The finding available here is about what
  the stamp does or does not determine, and every statement is
  written at that strength.

  NO NULL IS RUN and none applies: the quantity is a ratio compared
  against a published interval.
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
OUT = os.path.join(HERE, "..", "results", "verify_conjL_cells_range.txt")

NN, NK = 401, 1000
SCAN = (250, 500, 1000, 1500, 2000, 2500, 2900)
PUB_POOLED = 0.99781          # published in the repository's own audit
LO, HI = 0.99, 1.06           # published in v1/paper/wall_v1.tex
MOVE = 0.005
VIABLE = 100


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
    D = np.zeros((len(Ns), len(ks)), dtype=np.int32)
    S = np.zeros((len(Ns), len(ks)), dtype=np.int32)
    for i, N in enumerate(Ns):
        idx = int(N) - np.outer(ps, ks)
        ok = idx >= 1
        m = mu[np.where(ok, idx, 0)].astype(np.int32) * ok
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
    mu = mobius_upto(int(Ns.max()))
    v2k, v3k = vq(ks, 2), vq(ks, 3)
    a = np.broadcast_to(v2k[None, :], (NN, NK))
    b = np.broadcast_to(v3k[None, :], (NN, NK))
    cells = sorted(set(zip(a.ravel().tolist(), b.ravel().tolist())))

    say("the same grid, the prime range varied")
    say("  N = 3e6 + 2500 j, j = 0..%d and k = 1..%d" % (NN - 1, NK))
    say("  the published pooled ratio is %.5f and the stamp's "
        "interval is" % PUB_POOLED)
    say("  [%.2f, %.2f]" % (LO, HI))
    say("PRINTBOUND verify_conjL_cells_range 5 0.000005")
    say("  a cell is viable at %d surviving pairs" % VIABLE)
    say()
    say("     p <=   primes   mean support   pooled ratio   cells in "
        "[%.2f,%.2f]" % (LO, HI))
    pooled = []
    for pmax in SCAN:
        ps = primes_upto(pmax)
        ps = ps[ps > 2]
        D, S = build(Ns, ks, ps, mu)
        live = S > 0
        r = (float((D[live].astype(np.float64) ** 2).sum())
             / float(S[live].sum()))
        pooled.append(r)
        good = tot = 0
        for (x, y) in cells:
            m = live & (a == x) & (b == y)
            n = int(m.sum())
            if n < VIABLE:
                continue
            tot += 1
            rc = (float((D[m].astype(np.float64) ** 2).sum())
                  / float(S[m].sum()))
            good += int(LO <= rc <= HI)
        say("  %7d   %6d   %12.2f   %12.5f   %8d of %d"
            % (pmax, len(ps), float(S[live].mean()), r, good, tot))
        say("POINT cellsrange_%d %.6f" % (pmax, r))
    say("SCALES %d" % len(SCAN))

    lo, hi = min(pooled), max(pooled)

    # -------------------------------------------------------------- E1
    say()
    say("E1  does the prime range move the ratio at all?")
    say("  pooled ratio from %.5f to %.5f, movement %.5f"
        % (lo, hi, hi - lo))
    e1 = (hi - lo) > MOVE
    say("  E1 %s   (cap: %.3f)" % ("hold" if e1 else "REFUTED", MOVE))
    say("SPREAD cellsrange_pooled %.5f" % (hi - lo))

    # -------------------------------------------------------------- E2
    say()
    say("E2  is it monotone?")
    d = np.diff(pooled)
    up = int((d > 0).sum())
    e2 = up == len(d) or up == 0
    say("  steps %s" % ", ".join("%+.5f" % x for x in d))
    say("  %d of %d rise" % (up, len(d)))
    say("  E2 %s   (cap: no reversal)" % ("hold" if e2 else "REFUTED"))
    say("SIGNRUN cellsrange_steps %d %d" % (max(up, len(d) - up),
                                            len(d)))

    # -------------------------------------------------------------- E3
    say()
    say("E3  is the repository's pooled number inside the sweep?")
    e3 = lo <= PUB_POOLED <= hi
    say("  %.5f against the swept [%.5f, %.5f]" % (PUB_POOLED, lo, hi))
    say("  E3 %s   (cap: inside)" % ("hold" if e3 else "REFUTED"))

    # -------------------------------------------------------------- E4
    say()
    say("E4  is the stamp's floor reached?")
    e4 = hi >= LO
    say("  the largest pooled ratio in the scan is %.5f against the "
        "floor %.2f" % (hi, LO))
    say("  E4 %s   (cap: reached)" % ("hold" if e4 else "REFUTED"))

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4)))
    say()
    if e1 and e3:
        say("the stamp's number is a function of a range it does not "
            "print. A")
        say("prime range on this grid reproduces the repository's "
            "pooled ratio,")
        say("which does not show that range was used -- a parameter "
            "that")
        say("reproduces a number is not a parameter that was used -- "
            "but it does")
        say("show the stamp does not determine its own value.")
    elif e1 and not e3:
        say("the range moves the ratio and cannot reach the "
            "repository's number")
        say("on this grid, so the missing parameter is narrowed to "
            "'not the")
        say("range'. That is weaker than locating it and is what the "
            "scan can")
        say("support.")
    else:
        say("the range does not move the ratio, so it is not the "
            "missing")
        say("parameter and this step found nothing the grid did not "
            "already")
        say("explain.")

    head = [
        "STATISTIC: the same variance ratio E[D^2]/support, pooled",
        "           over all surviving pairs and per cell under",
        "           reading A, as the prime range varies.",
        "FIELD: N = 3*10^6 + 2500 j for j = 0..400 and k = 1..1000,",
        "       401000 pairs, with the field D(N,k) = sum of",
        "       mu(N-pk) over the odd primes p <= 250, 500, 1000,",
        "       1500, 2000, 2500 and 2900 in turn. The upper end is",
        "       set by N - p k staying positive, not by budget.",
        "STATEMENT: conj:L's 'exact cells' stamp, v1's table: every",
        "           viable (v_2,v_3) cell 0.99 to 1.06; and the same",
        "           conjecture's audit, pooled variance ratio",
        "           E[D^2]/support = %.5f." % PUB_POOLED,
        "METHOD HERE: the same grid as",
        "           verify_conjL_exact_cells.py with the prime range",
        "           varied over %s, the ratio pooled" % (SCAN,),
        "           over all surviving pairs and, for the cell count,",
        "           over reading A (the valuations of k).",
        "REPOSITORY'S NUMBER: %.5f pooled, interval %.2f to %.2f per"
        % (PUB_POOLED, LO, HI),
        "           cell, with no prime range published.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
