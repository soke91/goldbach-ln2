# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Section {#sec:c3}, item (2) -- the Heath-Brown weight
of the identity of level J with cut z.

WHAT IS UNDER TEST

The section says the natural bookkeeping needs the Mobius-side variable
a to satisfy a <= M^{o(1)}, and that the weight of the identity refuses
to sit there.  It prints three figures:

  (a) "the rounding of z alone moves the J=8 entry by 0.017";
  (b) "the top-j share is about 0.82 at J=3";
  (c) "the j in {6,7,8} share about 0.85 at J=8";

and adds the parameter-independent statement that the identity with cut
z needs z^J >= x, while its j-th term has a <= z^j, which at j = J is x
for every admissible (z,J).  No script for any of it exists here, and
the section itself declines to quote the fractions to three decimals
because "its value depends on conventions the source statement leaves
open".  Those conventions are fixed here so the sensitivity can be
measured rather than asserted.

CONVENTIONS FIXED HERE

The absolute weight of the j-th term is its L^1 mass,

    W_j = C(J,j) * #{ (m_1..m_j, n_1..n_j) :
                      m_i <= z, mu(m_i) != 0, prod m_i prod n_i <= x },

so W_j = C(J,j) * sum_t A_j(t) T_j(x/t), where A_j(t) counts ordered
j-tuples of squarefree m_i <= z with product t and T_j(y) = sum_{u<=y}
tau_j(u).  The share of term j is W_j / sum_i W_i.  x = M throughout,
and z is M^{1/J} rounded down, to nearest, or up.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  At J=3, M=1e6, the top-j share lies in [0.78, 0.86] -- "about
      0.82".
  E2  At J=8, M=1e6, the j in {6,7,8} share lies in [0.81, 0.89] --
      "about 0.85".
  E3  The three roundings of z move the J=8 entry by between 0.005 and
      0.05, i.e. the sensitivity is real and of the order of the quoted
      0.017.
  E4  The J=8 entry increases with M over M = 1e4, 1e5, 1e6 -- "the
      fraction increases with M".
  E5  z^J >= x fails for at least one of the three roundings at at
      least one M, so the rounding is not a free convention: it decides
      whether the identity holds at all.

REFUTATION RULE (fixed before the run)

  E1  REFUTED if the share falls outside [0.78, 0.86].
  E2  REFUTED if the share falls outside [0.81, 0.89].
  E3  REFUTED if the spread across the three roundings is below 0.005
      or above 0.05.
  E4  REFUTED if the entry is not non-decreasing across the three M.
  E5  REFUTED if z^J >= x holds for all three roundings at all three M.

  E1, E2 and E5 gate.  E3 and E4 are reported: E3 is the claim the
  paragraph makes about its own precision, and a mismatch there is a
  statement about conventions rather than about the mathematics.

CITED BY: {#rem:hbround} in paper/.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_hb_weight.txt")

MS = [10_000, 100_000, 1_000_000]
PUB_J3 = 0.82
PUB_J8 = 0.85
PUB_SENS = 0.017


def squarefree_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[0] = False
    p = 2
    while p * p <= n:
        s[p * p::p * p] = False
        p += 1
    return s


def tau_prefix(x, J):
    """T_j(y) = sum_{u<=y} tau_j(u) as arrays, j = 1..J."""
    out = {}
    t = np.zeros(x + 1, dtype=np.float64)
    t[1:] = 1.0                                    # tau_1
    out[1] = np.cumsum(t)
    for j in range(2, J + 1):
        nxt = np.zeros(x + 1, dtype=np.float64)
        for d in range(1, x + 1):
            nxt[d::d] += t[1:x // d + 1]
        t = nxt
        out[j] = np.cumsum(t)
    return out


def A_arrays(x, z, J, sqf):
    """A_j(t), ordered j-tuples of squarefree m_i <= z with product t."""
    ds = [d for d in range(1, z + 1) if sqf[d]]
    A = np.zeros(x + 1, dtype=np.float64)
    A[1] = 1.0
    out = {}
    for j in range(1, J + 1):
        nxt = np.zeros(x + 1, dtype=np.float64)
        for d in ds:
            nxt[d::d] += A[1:x // d + 1]
        A = nxt
        out[j] = A.copy()
    return out


def shares(x, J, rounding, sqf, T):
    r = x ** (1.0 / J)
    z = {"floor": int(math.floor(r)), "round": int(round(r)),
         "ceil": int(math.ceil(r))}[rounding]
    z = max(z, 1)
    A = A_arrays(x, z, J, sqf)
    W = []
    for j in range(1, J + 1):
        a = A[j]
        idx = np.flatnonzero(a)
        if idx.size == 0:
            W.append(0.0)
            continue
        W.append(float(math.comb(J, j)
                       * (a[idx] * T[j][x // idx]).sum()))
    tot = sum(W)
    return z, [w / tot for w in W], z ** J >= x


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    xmax = max(MS)
    say("sieving to %d and building tau prefixes to J = 8 ..." % xmax)
    sqf = squarefree_upto(xmax)
    T = tau_prefix(xmax, 8)

    say()
    say("J = 3")
    say("  M          rounding  z      z^J >= x   top-j share")
    j3 = {}
    adm_fail = False
    for M in MS:
        for rd in ("floor", "round", "ceil"):
            z, sh, adm = shares(M, 3, rd, sqf, T)
            adm_fail = adm_fail or (not adm)
            if M == MS[-1]:
                j3[rd] = sh[-1]
            say("  %-10d %-9s %-6d %-10s %.6f"
                % (M, rd, z, "yes" if adm else "NO", sh[-1]))

    say()
    say("J = 8")
    say("  M          rounding  z      z^J >= x   j in {6,7,8} share"
        "   per-j shares")
    j8 = {}
    trend = {}
    for M in MS:
        for rd in ("floor", "round", "ceil"):
            z, sh, adm = shares(M, 8, rd, sqf, T)
            adm_fail = adm_fail or (not adm)
            top3 = sh[5] + sh[6] + sh[7]
            if M == MS[-1]:
                j8[rd] = top3
            if rd == "round":
                trend[M] = top3
            say("  %-10d %-9s %-6d %-10s %-19.6f %s"
                % (M, rd, z, "yes" if adm else "NO", top3,
                   " ".join("%.3f" % v for v in sh)))

    say()
    e1 = 0.78 <= j3["round"] <= 0.86
    say("E1  J=3 top-j share at M=1e6 = %.6f   published 'about 0.82'"
        "   %s" % (j3["round"], "hold" if e1 else "REFUTED"))
    e2 = 0.81 <= j8["round"] <= 0.89
    say("E2  J=8 {6,7,8} share at M=1e6 = %.6f   published 'about 0.85'"
        "   %s" % (j8["round"], "hold" if e2 else "REFUTED"))
    spread = max(j8.values()) - min(j8.values())
    e3 = 0.005 <= spread <= 0.05
    say("E3  spread of the J=8 entry across the three roundings = %.6f"
        % spread)
    say("    published 'the rounding of z alone moves it by 0.017'   %s"
        % ("hold" if e3 else "REFUTED"))
    vals = [trend[M] for M in MS]
    e4 = all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))
    say("E4  J=8 entry by M (round): %s   %s"
        % (", ".join("%.6f" % v for v in vals),
           "hold" if e4 else "REFUTED"))
    e5 = adm_fail
    say("E5  some rounding makes z^J < x, so the identity fails there: "
        "%s   %s" % (adm_fail, "hold" if e5 else "REFUTED"))

    say()
    say("=" * 70)
    ok = e1 and e2 and e5
    say("E1 %s  E2 %s  E3 %s  E4 %s  E5 %s"
        % tuple("hold" if v else "REFUTED" for v in (e1, e2, e3, e4, e5)))
    say("Section {#sec:c3} item (2) reproduces at the stated convention"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the share of the j-th term in the absolute (L^1)",
        "           Heath-Brown weight at level J with cut z, i.e.",
        "           W_j = C(J,j) #{(m_1..m_j,n_1..n_j) : m_i <= z",
        "           squarefree, prod m prod n <= x} normalised by sum_i",
        "           W_i; reported for j = J at J = 3 and for j in {6,7,8}",
        "           at J = 8; the spread of that share across the three",
        "           roundings of z = M^{1/J}; and whether z^J >= x holds",
        "           for each rounding.",
        "FIELD: x = M = 1e4, 1e5, 1e6; J = 3 and 8; z = M^{1/J} rounded",
        "       down, to nearest and up; tau_j prefix sums built by",
        "       repeated divisor convolution to 1e6; A_j by repeated",
        "       convolution with the squarefree indicator on [1,z].",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
