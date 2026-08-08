# -*- coding: utf-8 -*-
r"""
The wall's bracket at N = 1e8, tested out of sample by measuring the
octave it forecasts.

WHAT IS AT STAKE

audit_margin_bracket.py replaced audit_margin.py's single number with a
bracket for the reciprocal 1/max|C|/N at N = 1e8, and Remark
{#rem:marginbracket} grades that bracket as the safest of the three the
repository carries, on the grounds that it is not made by assuming a
wobble: it IS the measured spread of the decay exponent b across both
octave grid anchors and every leave-one-out subset, so the wobble and
the drift are the same number by construction.

Every one of those refits lives inside the same data. The octaves run
from 3e4 to 1.6e7 and the bracket forecasts a factor 6.25 past the top
of them; nothing in the repository has ever measured on the other side.
A bracket that is honest about its own construction is still only as
good as the assumption that the fitted shape continues, and that
assumption has never been put at risk here.

It can be. The forecast is about a directly measurable quantity -- the
maximum of |C(N)|/N over even N in the octave (5e7, 1e8], read the way
audit_margin.py reads every other octave. Computing C = mu * Lambda to
1.28e8 needs a length-2^28 transform, which is the end of what this
structure reaches, and it settles the question without any fit at all.

THE MECHANISM THIS EXPECTS TO FIND, REGISTERED BEFORE THE RUN

The published table's argmax N are 16170, 30030, 60060, 120120, 300300,
510510, 1021020, 2042040, 4084080, 9699690 -- every one a multiple of a
primorial, because max|C|/N over an octave is taken at the N whose
singular series is largest and that is the N of largest odd radical the
octave contains. Two of the nine octave-to-octave steps do not fall at
all, they RISE, and both are the octave in which a new primorial first
becomes available: 30030 at top 6e4 and 510510 at top 1e6.

So the sequence is a staircase, not a power law, and the fitted b
averages over jumps. The next primorial after 9699690 is 223092870,
which is past 1.28e8. The stretch the bracket extrapolates across
therefore contains no new primorial at all, and should fall faster than
the fit that was calibrated on stretches that did.

BACKS: Remark {#rem:marginoos} in paper/wall_v3.md, which is Remark
{#rem:marginbracket} tested out of sample.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The length-2^28 transform reproduces the published max|C|/N on
      the octave (8e6, 1.6e7], at the published argmax, to 1e-6. A
      transform of a different length is a different computation and
      this is the only place a silent error in it would show.
  W2  Direct summation C(N) = sum_{n<N} Lambda(n) mu(N-n) at the argmax
      of the target octave reproduces the transform's value to 1e-6,
      the control audit_margin_bracket.py ran on every published
      octave and this script owes on every new one.
  W3  The measured reciprocal 1/max|C|/N over the octave (5e7, 1e8]
      lies inside the published bracket, read from
      results/audit_margin_bracket.txt. This is the whole point of the
      script and the first out-of-sample test that bracket has had.
  W4  Every octave-to-octave exponent over the new stretch lies inside
      the published spread of b, also read rather than typed.
  W5  The exponent over the whole new stretch 1.6e7 -> 1.28e8 exceeds
      every fitted b, because no new primorial enters it.

REFUTATION RULE (fixed before the run)

  W1  REFUTED at 1e-6, which would mean the two transforms disagree
      and the published octave table is not reproducible at this
      length.
  W2  REFUTED at 1e-6 at any new octave.
  W3  REFUTED if the measured reciprocal falls outside the published
      bracket. That is the outcome worth having: it would mean the
      safest bracket in the repository does not survive the first
      measurement made on the other side of it, and that the grade
      "safest" was about the bracket's construction and not about its
      reach.
  W4  REFUTED if any new octave step falls outside the published b
      range.
  W5  REFUTED if the new stretch's exponent is not greater than the
      largest fitted b, which would leave the staircase reading
      without support.

  W1 and W2 gate: if either fails nothing else here means anything.
  W3, W4 and W5 are the measurement and do not gate -- a refuted W3 is
  the finding, not an error.

  NO NULL IS RUN and none applies. A published deterministic
  extrapolation is compared against a direct measurement of the very
  quantity it extrapolates; there is no detection against a background.
  The control on the measurement is W2's second route, and the control
  on the transform is W1's reproduction of the published octave.
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
OUT = os.path.join(ROOT, "results", "audit_margin_oos.txt")

XMAX = 128_000_000        # the end of what a length-2^28 transform reaches
TOL = 1e-6                # the tolerance audit_margin_bracket.py used


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
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
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def pow2(n):
    L = 1
    while L < n:
        L <<= 1
    return L


def read_published():
    """the bracket, the b-range, the base and the reach -- read, not typed

    Every number this script judges against comes from
    audit_margin_bracket.py and audit_margin.py. G25 is the reason it is
    parsed: a copied number is a dependency no check can see, and G22
    then keeps the run order honest.
    """
    pb = os.path.join(ROOT, "results", "audit_margin_bracket.txt")
    pm = os.path.join(ROOT, "results", "audit_margin.txt")
    sb = io.open(pb, encoding="utf-8").read()
    sm = io.open(pm, encoding="utf-8").read()

    m = re.search(r"^BRACKET wall_reciprocal_at_1e8 "
                  r"([\d.]+) ([\d.]+) ([\d.]+)\s*$", sb, re.M)
    point, lo, hi = (float(m.group(1)), float(m.group(2)),
                     float(m.group(3)))
    m = re.search(r"b runs ([\d.]+) to ([\d.]+), a relative spread "
                  r"of ([\d.]+)", sb)
    blo, bhi, bdr = (float(m.group(1)), float(m.group(2)),
                     float(m.group(3)))
    base = float(re.search(r"1/top = ([\d.]+)", sm).group(1))

    # the published top octave and its argmax, from the halving grid
    j = sm.index("On that grid:")
    ln = sm[j:].splitlines()[2].split()
    top, val, arg = int(ln[0]), float(ln[2]), int(ln[3])
    return point, lo, hi, blo, bhi, bdr, base, top, val, arg


def octave_max(C, top):
    """max |C(N)|/N over even N in (top/2, top], and where it is taken"""
    lo = top // 2
    Ns = np.arange(lo + 2 - (lo % 2), top + 1, 2, dtype=np.int64)
    r = np.abs(C[Ns]) / Ns
    j = int(np.argmax(r))
    return float(r[j]), int(Ns[j]), int(Ns.size)


def odd_radical(n):
    r, m = 1, n
    while m % 2 == 0:
        m //= 2
    p = 3
    while p * p <= m:
        if m % p == 0:
            r *= p
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        r *= m
    return r


def factorise(n):
    out, m, p = [], n, 2
    while p * p <= m:
        while m % p == 0:
            out.append(p)
            m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        out.append(m)
    return out


def primorials(limit):
    out, v = [], 1
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        v *= p
        if v > limit:
            break
        out.append(v)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    (point, blo_v, bhi_v, blo, bhi, bdr, base,
     ptop, pval, parg) = read_published()
    say("read from results/audit_margin_bracket.txt and "
        "results/audit_margin.txt:")
    say("  the published BRACKET  [%.4f, %.4f]  point %.4f"
        % (blo_v, bhi_v, point))
    say("  the published b range  [%.6f, %.6f]  relative spread %.4f"
        % (blo, bhi, bdr))
    say("  the base 1/max at the top published octave  %.3f" % base)
    say("  the top published octave  (%d, %d]  max %.6f  at N = %d"
        % (ptop // 2, ptop, pval, parg))

    reach = 1e8 / ptop
    say("  so the forecast reaches a factor %.2f past the data, to the"
        % reach)
    say("  octave (%d, %d], which this script measures directly."
        % (int(5e7), int(1e8)))

    say()
    say("sieving to %d ..." % XMAX)
    lam, mu = sieves(XMAX)

    say("convolving C = mu * Lambda by a length-%d transform ..."
        % pow2(2 * (XMAX + 1)))
    n = pow2(2 * (XMAX + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:XMAX + 1] = lam
    F = np.fft.rfft(a)
    a[:] = 0.0
    a[:XMAX + 1] = mu
    F *= np.fft.rfft(a)
    del a
    C = np.fft.irfft(F, n)[:XMAX + 1].copy()
    del F

    # ------------------------------------------------------------- W1
    say()
    say("W1  the new transform against the published octave")
    v16, a16, c16 = octave_max(C, ptop)
    d1 = abs(v16 - pval)
    w1 = d1 < TOL
    say("  octave                    here         published    diff")
    say("  (%d, %d]  %-12.8f %-12.6f %.2e"
        % (ptop // 2, ptop, v16, pval, d1))
    say("  argmax N                  %-12d %-12d %s"
        % (a16, parg, "same" if a16 == parg else "DIFFERENT"))
    say("  W1 %s" % ("holds" if w1 and a16 == parg else "REFUTED"))
    w1 = w1 and a16 == parg

    # ------------------------------- the new octaves, and the target
    grid = []
    t = ptop
    while t * 2 <= XMAX:
        t *= 2
        grid.append(t)
    tgt = int(1e8)

    say()
    say("the new octaves, and the octave the bracket forecasts")
    say("  octave top     count      max |C|/N     argmax N     "
        "odd radical")
    rows = [(ptop, v16, a16, c16)]
    for t in grid + [tgt]:
        v, ar, cnt = octave_max(C, t)
        rows.append((t, v, ar, cnt))
    for t, v, ar, cnt in rows:
        say("  %-14d %-10d %-13.8f %-12d %d"
            % (t, cnt, v, ar, odd_radical(ar)))

    prim = primorials(XMAX)
    fresh = [p for p in prim if ptop < p <= XMAX]
    say("  primorials up to %d: %s" % (XMAX, ", ".join(str(p)
                                                       for p in prim)))
    say("  of these, new inside (%d, %d]: %s"
        % (ptop, XMAX, ", ".join(str(p) for p in fresh) if fresh
           else "none"))
    say("  and what each argmax actually is:")
    for t, v, ar, cnt in rows:
        say("  %-12d = %s" % (ar, " * ".join(str(q)
                                             for q in factorise(ar))))

    # ------------------------------------------------------------- W2
    say()
    say("W2  every new octave maximum, by direct summation")
    say("  C(N) = sum_{n<N} Lambda(n) mu(N-n), one pass per argmax N")
    say("  argmax N       direct         from the transform  diff")
    w2 = True
    seen = {}
    for t, v, ar, _ in rows[1:]:
        if ar in seen:
            direct = seen[ar]
        else:
            direct = abs(float(np.dot(
                lam[1:ar], mu[ar - 1:0:-1].astype(np.float64)))) / ar
            seen[ar] = direct
        d = abs(direct - v)
        if d >= TOL:
            w2 = False
        say("  %-14d %-14.8f %-19.8f %.2e" % (ar, direct, v, d))
    say("  W2 %s" % ("holds" if w2 else "REFUTED"))
    del lam, mu

    if not (w1 and w2):
        say()
        say("the controls do not hold, so nothing below is read")
        _write(lines, ptop, XMAX)
        raise SystemExit(1)

    # ------------------------------------------------------------- W3
    vtgt = rows[-1][1]
    meas = 1.0 / vtgt
    w3 = blo_v <= meas <= bhi_v
    say()
    say("W3  the bracket, against the octave it forecasts")
    say("  measured 1/max|C|/N on (%d, %d]  : %.4f"
        % (int(5e7), tgt, meas))
    say("  the published bracket                   : [%.4f, %.4f]"
        % (blo_v, bhi_v))
    say("  the published point estimate            : %.4f" % point)
    say("  W3 %s" % ("holds" if w3 else "REFUTED"))
    need = math.log(meas / base) / math.log(reach)
    say("  the exponent that would have landed it  : %.6f" % need)
    say("  against the fitted range                : [%.6f, %.6f]"
        % (blo, bhi))
    say("  the measured value is %.4f times the top of the bracket"
        % (meas / bhi_v))

    # ------------------------------------------------------------- W4
    say()
    say("W4  the octave-to-octave exponent over the new stretch")
    say("  step                      b_local     inside [%.6f, %.6f]?"
        % (blo, bhi))
    w4 = True
    locs = []
    for i in range(len(rows) - 2):
        t0, v0 = rows[i][0], rows[i][1]
        t1, v1 = rows[i + 1][0], rows[i + 1][1]
        b = -math.log(v1 / v0) / math.log(t1 / t0)
        locs.append(b)
        ok = blo <= b <= bhi
        w4 = w4 and ok
        say("  %-9d -> %-12d %-11.6f %s"
            % (t0, t1, b, "yes" if ok else "NO"))
    say("  W4 %s" % ("holds" if w4 else "REFUTED"))

    # ------------------------------------------------------------- W5
    v_end = rows[-2][1]
    bstr = -math.log(v_end / v16) / math.log(rows[-2][0] / ptop)
    w5 = bstr > bhi
    say()
    say("W5  the whole new stretch, against every fitted b")
    say("  %d -> %d is a factor %d and no new primorial enters it"
        % (ptop, rows[-2][0], rows[-2][0] // ptop))
    say("  max|C|/N ~ N^{-%.6f} over the new stretch" % bstr)
    say("  the largest b the published grids fit  : %.6f" % bhi)
    say("  W5 %s" % ("holds" if w5 else "REFUTED"))
    sp = max(locs) - min(locs)
    say("  local exponents run %.6f to %.6f" % (min(locs), max(locs)))
    say("SWEPT wall_max_decay_oos N-range %.6f" % sp)

    say()
    say("  DIAGNOSTIC (post hoc). The top new octave does not find its")
    say("  own maximum. Two of the octaves above share an argmax, so")
    say("  the step between them measures where the maximum stopped")
    say("  moving and not how fast the wall falls:")
    shared = 0
    for i in range(len(rows) - 1):
        for j2 in range(i + 1, len(rows)):
            if rows[i][2] == rows[j2][2]:
                shared += 1
                say("  (%d, %d] and (%d, %d] both peak at N = %d"
                    % (rows[i][0] // 2, rows[i][0], rows[j2][0] // 2,
                       rows[j2][0], rows[i][2]))
    say("  octave pairs sharing an argmax: %d" % shared)
    say("  So the last local exponent above is not an independent")
    say("  measurement of the decay, and the whole-stretch exponent")
    say("  W5 reads is a lower bound on it: the octave (%d, %d] would"
        % (XMAX // 2, XMAX))
    say("  need a maximum above the one it inherited to fall further.")

    say()
    say("  DIAGNOSTIC (post hoc on the published grids, pre-registered")
    say("  only for the new one). The published tables already contain")
    say("  steps on which max|C|/N rises rather than falls, so a local")
    say("  exponent there is negative and the fitted b is an average")
    say("  over a staircase. Read off results/audit_margin.txt, on both")
    say("  grids the bracket was fitted on:")
    pm = io.open(os.path.join(ROOT, "results", "audit_margin.txt"),
                 encoding="utf-8").read()
    nneg = nstep = 0
    for anchor, skip in (("octave top     count      max |C|/N     "
                          "argmax N", 1), ("On that grid:", 2)):
        j = pm.index(anchor)
        tab = []
        for ln in pm[j:].splitlines()[skip:]:
            f = ln.split()
            if len(f) != 4 or not f[0].isdigit():
                if tab:
                    break
                continue
            tab.append((int(f[0]), float(f[2]), int(f[3])))
        tab.sort()
        nstep += len(tab) - 1
        for i in range(len(tab) - 1):
            b = -math.log(tab[i + 1][1] / tab[i][1]) / math.log(
                tab[i + 1][0] / tab[i][0])
            if b < 0:
                nneg += 1
                say("  %-9d -> %-10d b_local %-11.6f  argmax %d enters"
                    % (tab[i][0], tab[i + 1][0], b, tab[i + 1][2]))
    say("  negative local steps inside the fitted range: %d of %d"
        % (nneg, nstep))

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s  W5 %s"
        % tuple("holds" if v else "REFUTED"
                for v in (w1, w2, w3, w4, w5)))
    if w3:
        say("the bracket survives its first out-of-sample test")
    else:
        say("the bracket does not contain the octave it forecasts")

    _write(lines, ptop, XMAX)
    return 0


def _write(lines, ptop, xmax):
    head = [
        "STATISTIC: max |C(N)|/N over even N in each octave (T/2, T] with",
        "           C = mu * Lambda by one length-2^28 FFT convolution,",
        "           for the octaves above the published table and for the",
        "           octave (5e7, 1e8] that the published bracket",
        "           forecasts; each new maximum recomputed by direct",
        "           summation at its argmax; the octave-to-octave and",
        "           whole-stretch decay exponents over the new range; and",
        "           the reciprocal 1/max|C|/N on the forecast octave",
        "           against the bracket itself.",
        "NULL: none is run and none applies. A published deterministic",
        "      extrapolation is compared against a direct measurement of",
        "      the quantity it extrapolates, so there is no background to",
        "      detect against. The control on the measurement is the",
        "      direct summation of W2 and the control on the transform is",
        "      W1's reproduction of the published octave at a different",
        "      transform length.",
        "FIELD: even N; Lambda and mu from an integer sieve to 1.28e8;",
        "       octaves with tops doubling from the published top up to",
        "       1.28e8, plus the octave (5e7, 1e8]; the bracket, the",
        "       range of the fitted exponent and the base are read from",
        "       results/audit_margin_bracket.txt and",
        "       results/audit_margin.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    sys.exit(main())
