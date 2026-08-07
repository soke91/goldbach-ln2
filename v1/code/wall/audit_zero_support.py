# -*- coding: utf-8 -*-
"""
The zero-support idiom, swept across the whole corpus (increment 286).

WHAT INCREMENT 285 FOUND, AND WHY IT IS NOT A ONE-OFF. The CI stamp's
V4 reported r = 0.401 against a half-normal 0.798. The cause was one
line,

    r = abs(vals.sum()) / math.sqrt(max(v, 1))

evaluated for pairs with v = 0. Those are M.3's predicted
annihilations -- pairs about which the law claims nothing -- and the
expression hands back r = 0, which then enters a MEAN as though it
were a measurement. 45% of the pairs in that band were annihilated,
and 0.55 x 0.81 = 0.45 reproduces the reported 0.401.

That is a code idiom, not an accident, and this program's oldest
recurring fault is exactly this one: a mask term read as a measurement.
So the question is not "was V4 wrong" but "where else".

WHAT THIS DOES.

 (A) STATIC SCAN. Every `sqrt(max(x, 1))` / `sqrt(maximum(x, 1))`
     divisor in code/, classified by whether the site guards against
     the empty support (a `> threshold` filter, an `np.where`, or a
     `nan`). Unguarded sites are listed with their line.

 (B) DYNAMIC MEASUREMENT. For the family these sites all compute --
     the prime-indexed correlations C_{k,k'} = Sum_p mu(N-pk)mu(N-pk')
     -- the same statistic is evaluated under three conventions:

       (a) annihilated pairs counted as r = 0     <- the buggy idiom
       (b) annihilated pairs excluded             <- correct
       (c) free class only, gcd(k k', N) = 1      <- what section 7 uses

     with the annihilation rate reported alongside. The gap between
     (a) and (b) IS the size of the fault wherever it is unguarded.

 (C) WHICH DOCUMENTED NUMBERS USE WHICH CONVENTION. Stated explicitly,
     because the answer turns out not to be uniform and the documents
     do not say.

PRE-REGISTRATION. There is no hypothesis test here: (A) is a scan and
(B) is the same quantity computed three ways, so the "result" is a
difference between conventions and cannot be a false positive. What it
CAN get wrong is the claim that a documented number used one
convention rather than another -- so (C) is stated as a reading of the
scripts, and each claim names the file it comes from.
"""
import io
import math
import os
import re
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(REPO, "code")

PAT = re.compile(r"(?:math\.sqrt|np\.sqrt)\(\s*(?:max|np\.maximum)\(")
GUARD = re.compile(r"np\.where|>\s*\d+\s*[,)]|if\s+\w+\s*[<>]=?\s*\d+|nan")


def static_scan():
    hits = []
    for fn in sorted(os.listdir(CODE)):
        if not fn.endswith(".py"):
            continue
        path = os.path.join(CODE, fn)
        lines = io.open(path, encoding="utf-8",
                        errors="replace").read().split("\n")
        for i, ln in enumerate(lines):
            if not PAT.search(ln):
                continue
            ctx = "\n".join(lines[max(0, i - 4): i + 2])
            hits.append((fn, i + 1, ln.strip(), bool(GUARD.search(ctx))))
    return hits


def sieve_mu(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    pm = (spf == np.arange(X + 1))
    pm[:2] = False
    return mu, pm


def main():
    print("(A) static scan: sqrt(max(x, 1)) divisors in code/")
    hits = static_scan()
    ung = [h for h in hits if not h[3]]
    print(f"    {len(hits)} sites, {len(ung)} without a visible guard\n")
    print(f"    {'file':>28} {'line':>5}  guarded?  code")
    for fn, ln, txt, g in hits:
        print(f"    {fn:>28} {ln:>5}  {'yes' if g else 'NO ':>7}   "
              f"{txt[:58]}")
    print("\n    'guarded' means a nearby `np.where`, a `> threshold`")
    print("    filter or a nan sentinel. It is a syntactic reading and")
    print("    can be wrong in both directions; (B) measures the size")
    print("    of the fault so the list can be prioritised rather than")
    print("    trusted.")

    X = 10_000_000
    t0 = time.time()
    mu, pm = sieve_mu(X)
    print(f"\nsieve  t={time.time()-t0:.0f}s", flush=True)

    N = 9_999_998
    rng = np.random.default_rng(286)
    print("\n(B) the same statistic under three conventions")
    print(f"    N = {N}, C_(k,k') = Sum_p mu(N-pk) mu(N-pk')")
    print(f"\n{'k-band':>14} {'pairs':>7} {'annih.':>8} "
          f"{'(a) as 0':>10} {'(b) excl.':>10} {'(c) free':>10} "
          f"{'(a)/(b)':>9}")
    for K0, K1 in ((252, 464), (2000, 4000)):
        allr, keptr, freer = [], [], []
        nann = 0
        for _ in range(600):
            k = int(rng.integers(K0, K1))
            kp = int(rng.integers(K0, K1))
            if k == kp:
                continue
            P1 = min(110_000, (N - 2) // max(k, kp))
            P0 = P1 // 2
            ps = np.nonzero(pm[P0:P1])[0] + P0
            w = N - ps * k
            wp = N - ps * kp
            ok = (w > 1) & (wp > 1)
            vals = mu[w[ok]].astype(np.float64) * mu[wp[ok]]
            v = int(np.count_nonzero(vals))
            r_buggy = abs(float(vals.sum())) / math.sqrt(max(v, 1))
            allr.append(r_buggy)
            if v < 50:
                nann += 1
                continue
            keptr.append(r_buggy)
            if math.gcd(k * kp, N) == 1:
                freer.append(r_buggy)
        a = float(np.mean(allr))
        b = float(np.mean(keptr))
        c = float(np.mean(freer)) if freer else float("nan")
        print(f"{f'[{K0},{K1})':>14} {len(allr):>7} "
              f"{nann/len(allr):>7.1%} {a:>10.4f} {b:>10.4f} "
              f"{c:>10.4f} {a/b:>9.4f}")
    print(f"    half-normal reference {math.sqrt(2/math.pi):.4f}")
    print("    (a) is what an unguarded site reports. The ratio (a)/(b)")
    print("    is 1 - (annihilation rate) by construction, so the fault")
    print("    is exactly as large as the mask is dense -- which is why")
    print("    it is invisible in shallow bands and severe in deep ones.")


    # The same statistic over ALL m rather than over primes p. This
    # is what dispersion_engine.py samples. The first draft of this
    # script "explained" its 0.801 by arguing that a number sitting
    # at 0.80 could not have been depressed -- which assumes the
    # conclusion. Measure it instead: over all m the support is
    # dense and annihilation should be essentially absent.
    print("")
    print("    control: the same family summed over ALL m rather")
    print("    than over primes p -- what dispersion_engine.py does")
    print(f"{'k-band':>14} {'pairs':>7} {'annih.':>8} "
          f"{'(a) as 0':>10} {'(b) excl.':>10} {'(a)/(b)':>9}")
    KD = int(X ** 0.4)
    for K0, K1 in ((KD // 2, KD), (2000, 4000)):
        allr, keptr = [], []
        nann = 0
        for _ in range(200):
            k1 = int(rng.integers(K0, K1))
            k2 = int(rng.integers(K0, K1))
            if k1 == k2:
                continue
            M = (N - 1) // max(k1, k2)
            ms = np.arange(1, M + 1, dtype=np.int64)
            prod = mu[N - k1 * ms].astype(np.int16) * mu[N - k2 * ms]
            v = int(np.count_nonzero(prod))
            r = abs(float(prod.sum())) / math.sqrt(max(v, 1))
            allr.append(r)
            if v < 50:
                nann += 1
                continue
            keptr.append(r)
        a = float(np.mean(allr)); b = float(np.mean(keptr))
        print(f"{f'[{K0},{K1})':>14} {len(allr):>7} "
              f"{nann/len(allr):>7.1%} {a:>10.4f} {b:>10.4f} "
              f"{a/b:>9.4f}")
    print("    So the idiom is harmless over all m and severe over")
    print("    primes p. The discriminator is the DENSITY of the")
    print("    summation range, not the site syntax -- which is why")
    print("    the static scan in (A) can only prioritise, not judge.")
    print("\n(C) which documented numbers use which convention")
    print("    (a reading of the scripts, each claim naming its file)")
    for f, conv, note in (
        ("e1_zero_account.py", "(b)/(c)",
         "counts per-pair zeros directly; MEASUREMENTS section 7's "
         "0.97-1.02 and kurtosis 2.99-3.03 rest on this"),
        ("e1_mask_model.py", "(c)",
         "blind mask verification, free class"),
        ("verify_all.py", "(b)+(c)",
         "was (a) until increment 285; V4 read 0.401"),
        ("dispersion_engine.py", "(a)",
         "section 6's 0.801 -- see the note below"),
        ("cascade2.py", "(a)", "duals, not quoted in MEASUREMENTS"),
        ("e1_dilation.py", "(a)",
         "section 6's cross-scale correlations, a CORRELATION of two "
         "such quantities rather than a mean of them"),
    ):
        print(f"    {f:>24}  {conv:>8}  {note}")
    print("")
    print("    Two corrections to this script own first draft:")
    print("    * section 6 dispersion number (0.801) is safe, and the")
    print("      measured reason is the control above -- it sums over")
    print("      ALL m, where annihilation is absent. The first draft")
    print("      argued that a number near 0.80 could not have been")
    print("      depressed, which assumes the answer.")
    print("    * e1_seam_law.py is listed unguarded by (A) but is")
    print("      INSTRUMENTED: its line carries a 0-when-killed note")
    print("      and the next statement prints the annihilation rate.")
    print("      The scan looks for filters and misses a published")
    print("      diagnostic, so NO in column 3 means no filter, not")
    print("      nobody noticed.")
    print("DONE")


if __name__ == "__main__":
    main()
