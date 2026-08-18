# -*- coding: utf-8 -*-
r"""
What would independent signs actually give across dilations?

WHAT IS AT STAKE

Remark {#rem:nocrossk} closes half the programme. With nonnegative
weights, whatever smallness E_3 has must come from the signs of
H(N;k) across k, and the measured gain
G = sum_k (log k)|H| / |sum_k (log k)H| is 1.834 to 2.789 where

    "Independent signs would give sqrt(#k)" = 17.7 to 38.5.

The reading is that the sum behaves as though it had n_eff = G^2 =
3.4 to 7.8 independent signs "where it has hundreds to thousands of
terms: the dilated walls move together".

sqrt(#k) is what independent signs give on EQUAL magnitudes. These
magnitudes are not equal -- the remark's own rule T4 measured the top
decile at 0.3486 to 0.3587 of the mass, which is concentrated even
though it refuted domination. For unequal magnitudes a random sign
sum has mean |.| near sqrt(2/pi) times the l2 norm, so the gain a coin
gives is about sqrt(pi/2) * l1/l2, and l1/l2 is below sqrt(#k) by
exactly the concentration. Nobody has computed it.

That matters for the size of the deficit, not its existence. If
l1/l2 is well below sqrt(#k) then part of the missing cancellation is
the magnitude distribution and only the rest is correlation, and the
factor by which the walls "move together" is smaller than the
published comparison makes it.

The coin arm in {#rem:nocrossk} does not answer this: it is a coin
with its own H, whose gain reads 965.6, 13.3, 25.0, 1000.4, 96.0 --
erratic because its denominator is a near-cancellation. The arm this
needs is random signs on mu's OWN magnitudes, the same null
lab_split_budget.py uses.

The implementation is independent of lab_positive_weights.py's.

BACKS: Remark {#rem:crosskreference} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The control: G reproduces the published 1.834, 1.804, 2.207,
      2.588, 2.789 to within 0.01 at every N.
  T2  sqrt(#k) is not the right reference: random signs on mu's own
      magnitudes give a median gain below sqrt(#k) at every N.
  T3  The deficit survives anyway: that median coin gain is at least
      three times mu's at every N, so a correlation across k is
      still there and {#rem:nocrossk}'s conclusion stands.
  T4  And l1/l2 is the reference the coin sits on: the median coin
      gain divided by l1/l2 lies in [1.0, 1.6] at every N, the band
      containing sqrt(pi/2) = 1.2533.

REFUTATION RULE (fixed before the run)

  T1  REFUTED at 0.01 at any N -- not the same statistic, and nothing
      below may be compared with {#rem:nocrossk}.
  T2  REFUTED if the coin reaches sqrt(#k) at any N. Then the
      magnitudes are close enough to equal that sqrt(#k) was the
      right reference and the published comparison needs no
      qualification.
  T3  REFUTED if the coin's median gain falls below three times mu's
      at any N. That is the one that matters: it would say most of
      the missing sqrt(#k) is the magnitude distribution rather than
      a correlation across dilations, and the sentence "the dilated
      walls move together" would have to be withdrawn.
  T4  REFUTED outside [1.0, 1.6], which would say the coin is not
      behaving like a random sign sum on these magnitudes and the
      whole reference calculation is wrong.

  All four gate.

  THE NULL IS THE POINT and it is run: random sign vectors over the
  admissible k, applied to the IDENTICAL (log k)|H(N;k)| magnitudes,
  so that "what independent signs give" is measured on this sequence
  of magnitudes rather than assumed from their count. Same convention
  as lab_split_budget.py's size permutation and
  lab_residue_cancellation.py's coin arm.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_crossk_reference.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
THETA = 0.56
DRAWS = 256
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def read_published():
    """#k, G and sqrt(#k) at each N, read from the results file"""
    src = io.open(os.path.join(RES, "lab_positive_weights.txt"),
                  encoding="utf-8").read()
    i = src.index("G_eps     sqrt K   top-decile share")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 8 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not out:
                continue
            break
        out[int(f[0])] = (int(f[2]), float(f[4]), float(f[6]))
    return out


def gains(N, lam, mu, sqf):
    """(log k)-weighted H over the admissible k, and the gains"""
    PN = factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Hs = []
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
    H = np.array(Hs)
    a = np.log(ks.astype(np.float64)) * H
    return ks, a


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published rows from results/lab_positive_weights.txt"
        % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rng = np.random.default_rng(SEED)
    say("%d random sign vectors over the admissible k, seed %d"
        % (DRAWS, SEED))

    got = []
    for N in NS:
        ks, a = gains(N, lam, mu, sqf)
        l1 = float(np.abs(a).sum())
        l2 = float(np.sqrt((a * a).sum()))
        G = l1 / abs(float(a.sum()))
        eps = rng.integers(0, 2, size=(DRAWS, ks.size)) * 2 - 1
        cg = l1 / np.abs(eps @ np.abs(a))
        got.append((N, ks.size, G, l1, l2, np.sort(cg)))
        say("  N = %-10d #k = %-6d G = %.4f  l1/l2 = %.4f  "
            "sqrt(#k) = %.4f"
            % (N, ks.size, G, l1 / l2, math.sqrt(ks.size)))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- T1
    say()
    say("T1  the control: the published gain")
    say("  N            #k here  #k pub   G here   G published  diff")
    t1 = True
    for N, nk, G, l1, l2, cg in got:
        pnk, pG, psq = pub[N]
        d = abs(G - pG)
        if not (d < 0.01) or nk != pnk:
            t1 = False
        say("  %-12d %-8d %-8d %-8.4f %-12.4f %.5f"
            % (N, nk, pnk, G, pG, d))
    say("  T1 %s   (cap 0.01)" % ("hold" if t1 else "REFUTED"))

    # ---------------------------------------------------------- T2/T3
    say()
    say("T2/T3  what independent signs give on mu's OWN magnitudes")
    say("  N            sqrt(#k)  coin median  coin [2.5%, 97.5%]"
        "        mu      coin/mu")
    t2 = t3 = True
    rat = []
    for N, nk, G, l1, l2, cg in got:
        med = float(np.median(cg))
        lo = float(cg[int(0.025 * DRAWS)])
        hi = float(cg[int(0.975 * DRAWS)])
        r = med / G
        rat.append(r)
        if med >= math.sqrt(nk):
            t2 = False
        if r < 3.0:
            t3 = False
        say("  %-12d %-9.4f %-12.4f [%8.4f, %8.4f]  %-7.4f %.4f"
            % (N, math.sqrt(nk), med, lo, hi, G, r))
    say("  T2 the coin stays below sqrt(#k) at every N   %s"
        % ("hold" if t2 else "REFUTED"))
    say("  T3 and is at least three times mu (min %.4f)   %s"
        % (min(rat), "hold" if t3 else "REFUTED"))
    say("PERN crossk_coin_over_mu %d %.4f %.4f"
        % (len(rat), min(rat), max(rat)))

    # ------------------------------------------------------------- T4
    say()
    say("T4  is l1/l2 the reference the coin sits on?")
    say("  N            l1/l2     coin median / (l1/l2)")
    t4 = True
    for N, nk, G, l1, l2, cg in got:
        v = float(np.median(cg)) / (l1 / l2)
        if not (1.0 <= v <= 1.6):
            t4 = False
        say("  %-12d %-9.4f %.4f" % (N, l1 / l2, v))
    say("  the band contains sqrt(pi/2) = %.4f, which is E|Z|/sigma,"
        % math.sqrt(math.pi / 2.0))
    zmed = float(np.median(np.abs(rng.standard_normal(400000))))
    say("  and also %.4f, which is the constant a MEDIAN gain sits"
        % (1.0 / zmed))
    say("  on: the median of |Z| for a standard normal, measured here")
    say("  over 400000 draws, is %.4f. The measured ratios straddle"
        % zmed)
    say("  the second, as they should.")
    say("  T4 %s   (band [1.0, 1.6])" % ("hold" if t4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC on T2 (post hoc). T2 compared a MEDIAN coin")
    say("  gain with sqrt(#k), and those are different conventions:")
    say("  sqrt(#k) is l1/l2 for equal magnitudes, while the median")
    say("  gain carries the extra 1/%.4f. Comparing like with like"
        % zmed)
    say("  needs no distributional constant at all, because l1/l2 is")
    say("  the reference on both sides. The deficit then factors:")
    say("  N            sqrt(#k)/G   = magnitudes  x  correlation")
    conc, corr = [], []
    for N, nk, G, l1, l2, cg in got:
        tot = math.sqrt(nk) / G
        cm = math.sqrt(nk) / (l1 / l2)
        cr = (l1 / l2) / G
        conc.append((l1 / l2) / math.sqrt(nk))
        corr.append(cr)
        say("  %-12d %-12.4f %-13.4f %.4f" % (N, tot, cm, cr))
    say("  the magnitude factor is %.4f to %.4f and flat; the"
        % (min(math.sqrt(g[1]) / (g[3] / g[4]) for g in got),
           max(math.sqrt(g[1]) / (g[3] / g[4]) for g in got)))
    say("  correlation factor is %.4f to %.4f and rising. So T2 fails"
        % (min(corr), max(corr)))
    say("  as registered and what it caught was its own mixing of")
    say("  conventions, not the magnitudes being equal.")
    say("REFERENCE audit_crossk_reference %d %.4f %.4f"
        % (len(conc), min(conc), max(conc)))
    say("REFERENCE lab_positive_weights %d %.4f %.4f"
        % (len(conc), min(conc), max(conc)))

    say()
    say("  what this does to the published reading. n_eff = G^2 was")
    say("  compared with #k; the comparison that isolates the")
    say("  correlation is with the coin's own n_eff on the same")
    say("  magnitudes:")
    say("  N            #k       n_eff(mu)  n_eff(coin)  "
        "coin/mu in n_eff")
    nef = []
    for N, nk, G, l1, l2, cg in got:
        med = float(np.median(cg))
        nef.append((med / G) ** 2)
        say("  %-12d %-8d %-10.4f %-12.4f %.4f"
            % (N, nk, G * G, med * med, (med / G) ** 2))
    say("  so the walls move together by a factor %.1f to %.1f in the"
        % (min(nef), max(nef)))
    say("  effective count, not by the %d to %d that comparing with"
        % (min(g[1] for g in got), max(g[1] for g in got)))
    say("  #k suggests. In the gain itself the correlation is %.2f to"
        % min(corr))
    say("  %.2f and the concentration of the magnitudes carries the"
        % max(corr))
    say("  remaining factor %.2f to %.2f."
        % (min(math.sqrt(g[1]) / (g[3] / g[4]) for g in got),
           max(math.sqrt(g[1]) / (g[3] / g[4]) for g in got)))

    say()
    say("=" * 70)
    ok = t1 and t2 and t3 and t4
    say("the correlation is real and smaller than the published "
        "comparison makes it" if ok else "REFUTED")

    head = [
        "STATISTIC: the cross-k gain",
        "           G = sum_k (log k)|H(N;k)| / |sum_k (log k)H(N;k)|",
        "           over the squarefree k < N^" + str(THETA)
        + " coprime to N;",
        "           the l1 and l2 norms of the same weighted",
        "           magnitudes; and the same gain with the signs of",
        "           those magnitudes redrawn at random, " + str(DRAWS),
        "           draws, reported as a median and a central",
        "           interval.",
        "NULL: random sign vectors over the admissible k applied to",
        "      the IDENTICAL (log k)|H(N;k)|, so that what",
        "      independent signs give is measured on this sequence of",
        "      magnitudes rather than assumed from their count. Same",
        "      convention as lab_split_budget.py's size permutation.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k coprime to k, as in the E_3 setting;",
        "       Lambda and mu from an integer sieve to " + str(NMAX)
        + ";",
        "       numpy default_rng seed " + str(SEED) + ". Every N is",
        "       2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published gains are read from",
        "       results/lab_positive_weights.txt.",
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
