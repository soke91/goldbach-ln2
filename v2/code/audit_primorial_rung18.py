# -*- coding: utf-8 -*-
r"""
The nineteenth rung: the two shapes disagree about 0.56.

WHAT IS AT STAKE

[rem:rung17] measured rung 17 at 0.555925, short of 0.56 by 0.004075
against a floor of 0.0037, and put the eighteen-rung crossing at
log10 N = 9.6667 with bracket [9.6009, 9.8386].  Rung 18 sits at
9.8961, above the upper end of that bracket.

Fitted on the eighteen uniform rungs the quadratic predicts 0.5668
here and the line 0.5565, against a prediction standard error of
0.0043.  **The quadratic clears 0.56 by 1.60 of that error and the
line does not clear it at all.**  Every earlier rung asked the two
shapes for a number; this one asks them a yes-or-no question about a
level, and they answer differently.

And the answer is what the branch has been walking toward.  Rung 17
is below 0.56 by more than the floor; if rung 18 is above it, then
0.56 is bracketed by the sign of two measurements rather than located
by any fitted curve -- which is what [rem:primorialgap] did for 1/2,
and the only kind of statement about theta' that {#rem:shapepower}
does not forbid.

THE PACKING, AND WHY THE PRIME IS NOT STORED

Rung 18 breaks {#rem:rung17}'s route twice: N = 7872184320 exceeds
2^32, so the prime carrying Lambda no longer fits a uint32, and the
three half-index arrays would be 19.68 GB.  Both go away by not
storing the prime at all.

On the support of Lambda every v is a prime power, and **when v is
prime the prime is v** -- there is nothing to store.  What is left is
the powers p^j with j >= 2, and those have p <= sqrt(N), so there are
fewer than ten thousand odd ones in the whole range: a table, not an
array.  So one uint8 per half-index carries which of the three cases
v is in, the logarithm is taken at the point of use of the same
float64 integer, and the 2^32 ceiling disappears with the prime.

The residue mask sheds a byte for a second reason special to this
ladder.  k is coprime to N = 30030*2^18, so none of 3, 5, 7, 11, 13
can divide k, and their five bits are always required-zero together
-- one bit, "v is coprime to 15015".  Only 17, 19, 23, 29 can divide
k and need bits of their own.  Nine bits become five and the mask
fits a uint8.

Three arrays of one byte per half-index: 11.81 GB, against 19.68 GB
for {#rem:rung17}'s packing at this N and 39.36 GB for
{#rem:rung16}'s.

BACKS: Remark {#rem:rung18} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The sieve control.  At n = 20000000 this route reproduces the
      production route of code/audit_primorial_rung11.py at every
      odd index: Moebius equal, and Lambda -- rebuilt as log v where
      the kind byte says prime and from the power table where it says
      prime power -- equal bit for bit.
  C2  The mask control.  At n = 20000000 the five-bit mask returns
      the same keep decision as the nine-bit mask of
      code/audit_primorial_rung17.py at every odd v, for every
      admissible k below 1000.
  C3  The ladder control.  At the uniform cap 10^6 rungs 16 and 17
      return the K*_R that results/audit_ladder_cap.txt and
      results/audit_primorial_rung17.txt printed, 126079 and 215843,
      exactly, with exponents 0.548808 and 0.555925 to the six
      decimals printed.
  H1  The margin keeps growing: the margin over 1/2 exceeds rung
      17's 0.0559.
  H2  The curvature predicts a sixth time: the departure from the
      quadratic fitted on the eighteen uniform rungs is inside that
      fit's prediction standard error 0.0043.
  H3  And it still beats the line: the quadratic's departure is
      smaller in absolute value than the line's.
  H4  **0.56 is crossed here.**  The exponent at rung 18 is at or
      above 0.56, as the quadratic says by 1.60 of its own prediction
      error and against the line, which says it is not.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a single odd index where either array disagrees.
      No tolerance: one is an integer array and the other is the same
      arithmetic on the same inputs.  THIS ONE GATES.
  C2  REFUTED by a single (v, k) where the two masks decide
      differently.  The five-bit packing rests on k being coprime to
      15015, which is a fact about this ladder and not about the
      statistic, so it is checked rather than asserted.  THIS ONE
      GATES.
  C3  REFUTED by a single K*_R that differs, or an exponent that
      differs in the six decimals printed.  THIS ONE GATES.
  H1  REFUTED if the margin does not grow.  Eight rungs have grown in
      a row; a ninth that does not would end the escalation.
  H2  REFUTED if the departure exceeds the prediction standard error.
      Five out-of-sample hits are what the quadratic's standing
      rests on.
  H3  REFUTED if the line is closer.
  H4  REFUTED if the exponent is below 0.56.  Then the line has won
      the one question the two shapes answer differently, the
      crossing is further out than eighteen rungs place it, and the
      bracketing of 0.56 by two measurements does not happen at this
      rung.  **That is the refutation worth having, and it is the
      reason this rung is worth its two hours**: the level is not
      being confirmed, it is being asked.

  C1, C2 and C3 gate.  H1 to H4 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  A deterministic curve is located
  against a computed threshold; there is no background to detect
  against.  The coin arms for this statistic were run in
  lab_primorial_ladder.py and lab_primorial_share.py, and the scatter
  they left is the floor the margin is judged against.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_primorial_rung18.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
SEED = 20260823
DRAWS = 4000
TARGET = 0.56
BLOCK = 1 << 24                     # the block, in half-indices
SIEVECHECK = 20_000_000             # where C1 and C2 compare
CAP = 1_000_000                     # the uniform cap of {#rem:laddercap}
SMALL = frozenset((3, 5, 7, 11, 13))   # the primes that always divide N
BIG = (17, 19, 23, 29)                 # the ones that may divide k


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
R17 = module("audit_primorial_rung17")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL1 = BASE * (1 << 16)         # 1968046080, the rung 16 point
CONTROL2 = BASE * (1 << 17)         # 3936092160, the rung 17 point
NEW = BASE * (1 << 18)              # 7872184320

KIND_ZERO, KIND_PRIME, KIND_POWER = 0, 1, 2


def power_table(n):
    """the odd prime powers p^j, j >= 2, below n, and their log p

    There are fewer than sqrt(n) of them because p <= sqrt(n), so
    this is what replaces an array of order n: Lambda is log v when v
    is prime, and the only values it cannot read off v itself are
    these.
    """
    vs, lg = [], []
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        if p == 2:
            continue
        lp = math.log(float(p))
        q = p * p
        while q <= n:
            vs.append(q)
            lg.append(lp)
            if q > n // p:
                break
            q *= p
    o = np.argsort(np.array(vs, dtype=np.int64))
    return (np.array(vs, dtype=np.int64)[o],
            np.array(lg, dtype=np.float64)[o])


def kind_and_mu_odd(n, block=BLOCK):
    """the kind byte and Moebius on odd indices, addressed by v >> 1

    kind[v >> 1] is 1 when the odd v is prime, 2 when it is a higher
    power of an odd prime, and 0 otherwise, so that
    Lambda(v) = log v, log p from power_table, or 0 respectively.
    The prime is never stored, which is what lifts the 2^32 ceiling
    of {#rem:rung17}'s route.
    """
    half = n // 2 + 1
    root = int(math.isqrt(n))
    pr = [int(p) for p in primes_upto(root) if p > 2]
    kind = np.zeros(half, dtype=np.uint8)
    mu = np.empty(half, dtype=np.int8)
    for t0 in range(0, half, block):
        t1 = min(t0 + block, half)
        w = t1 - t0
        vals = 2 * np.arange(t0, t1, dtype=np.int64) + 1
        rem = vals.copy()
        m = np.ones(w, dtype=np.int8)
        for p in pr:
            s = ((p - 1) // 2 - t0) % p
            if s < w:
                m[s::p] = -m[s::p]
                rem[s::p] //= p
            q = p * p
            if q <= n:
                s = ((q - 1) // 2 - t0) % q
                if s < w:
                    m[s::q] = 0
            pk = p * p
            while pk <= n:
                s = ((pk - 1) // 2 - t0) % pk
                if s < w:
                    rem[s::pk] //= p
                if pk > n // p:
                    break
                pk *= p
        big = rem > 1
        m[big] = -m[big]
        mu[t0:t1] = m
        idx = np.flatnonzero(big & (rem == vals))
        if idx.size:
            kind[t0 + idx] = KIND_PRIME
        del vals, rem, m, big, idx
    for p in pr:
        kind[p >> 1] = KIND_PRIME
    pv, _ = power_table(n)
    kind[pv >> 1] = KIND_POWER
    return kind, mu


def mask5_odd(n, qs):
    """the residue mask in five bits instead of nine

    Bit 0 says v shares a factor with 15015 = 3*5*7*11*13, which is
    enough because k is coprime to N and those five divide N, so
    their nine-bit counterparts are always required-zero together.
    Bits 1 to 4 are 17, 19, 23, 29, the only ones that can divide k.
    """
    assert set(qs) == SMALL | set(BIG), "the sieve weight moved"
    m = np.zeros(n // 2 + 1, dtype=np.uint8)
    for q in sorted(SMALL):
        m[(q - 1) // 2::q] |= np.uint8(1)
    for i, q in enumerate(BIG):
        m[(q - 1) // 2::q] |= np.uint8(1 << (i + 1))
    return m


def keepbits5(k):
    """the bits of mask5 that must be clear for a given k"""
    b = 1                                   # 3,5,7,11,13 never divide k
    for i, q in enumerate(BIG):
        if k % q:
            b |= 1 << (i + 1)
    return np.uint8(b)


def measure_kind(N, kind, mu, vmask, pv, plg, qs, artin, twin, cap,
                 block=BLOCK):
    """the ladder's K*_R, with Lambda read off the kind byte

    The statistic is that of {#rem:laddercap}; only the storage
    differs.  Where kind says prime the term is log of the index
    itself, and where it says power the logarithm comes from the
    table -- in both cases the logarithm of the same float64 integer
    the production route took, so each term is the identical float.
    """
    PN = R11.factor_set(N)
    assert SMALL <= PN, "this packing needs 15015 | N"
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, cap):
        if any(k % q == 0 for q in PN):
            continue
        if mu[k >> 1] == 0:
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ck = 1.0
        for q in qs:
            if k % q:
                ck *= q / (q - 1.0)
        kb = keepbits5(k)
        drop = [q for q in R11.factor_set(k) if q > 2]
        h, pp, seen = 0.0, 0.0, 0
        for lo in range(1, M + 1, block):
            hi = min(lo + block, M + 1)
            ms = np.arange(lo if lo % 2 else lo + 1, hi, 2,
                           dtype=np.int64)
            if ms.size == 0:
                continue
            ms = ms[mu[ms >> 1] != 0]
            for q in drop:
                ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            seen += ms.size
            vals = N - ms * k
            hv = vals >> 1
            g = mu[ms >> 1].astype(np.float64)
            kd = kind[hv]
            pri = kd == KIND_PRIME
            if pri.any():
                h += float((np.log(vals[pri].astype(np.float64))
                            * g[pri]).sum())
            pwr = kd == KIND_POWER
            if pwr.any():
                w = vals[pwr]
                h += float((plg[np.searchsorted(pv, w)]
                            * g[pwr]).sum())
            keep = (vmask[hv] & kb) == 0
            pp += float(g[keep].sum())
        if seen == 0:
            continue
        ks.append(k)
        Hs.append(h)
        Ps.append(ck * pp)
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())
    cum = np.cumsum(np.log(ks.astype(np.float64))
                    * np.abs(H - beta * P))
    thr = S_ * (1.0 - A_) * N
    j = int(np.searchsorted(cum, thr))
    if j >= ks.size:
        return None
    kk = int(ks[j])
    return kk, math.log(kk) / math.log(N), thr / N, beta, ks.size


def read_ladder():
    """the eighteen uniform rungs at cap 10^6, and the scatter"""
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    ns, ex, star = [], [], {}
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if not m:
            continue
        N = int(m.group(2))
        f = [t for t in m.group(3).split() if t != "cap-invariant"]
        if len(f) < 8 or f[6] == "none":
            continue
        ns.append(N)
        ex.append(float(f[7]))
        star[N] = int(f[6])
    s17 = io.open(os.path.join(RES, "audit_primorial_rung17.txt"),
                  encoding="utf-8").read()
    N17 = BASE * (1 << 17)
    m = re.search(r"^  N = " + str(N17) +
                  r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                  r"K\*_R (\d+)\s+exp ([\d.]+)\s*$", s17, re.M)
    ns.append(N17)
    ex.append(float(m.group(2)))
    star[N17] = int(m.group(1))
    scat = float(re.search(r"^FLOOR primorial_rung17 ([\d.]+)\s*$",
                           s17, re.M).group(1))
    prev = float(re.search(r"^BRACKET ladder_quadratic18_theta_prime "
                           r"([\d.]+) [\d.]+ [\d.]+\s*$",
                           s17, re.M).group(1))
    return ns, ex, star, scat, prev


def quadfit(x, y):
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, s2, float(np.sqrt((r ** 2).mean()))


def cross(c, level):
    a2, b2, c2 = c[0] - level, c[1], c[2]
    if abs(c2) < 1e-18:
        return None if abs(b2) < 1e-18 else -a2 / b2
    disc = b2 * b2 - 4.0 * c2 * a2
    if disc < 0:
        return None
    rs = [r for r in ((-b2 + math.sqrt(disc)) / (2.0 * c2),
                      (-b2 - math.sqrt(disc)) / (2.0 * c2)) if r > 0]
    return min(rs) if rs else None


def main():
    lines = []

    def say(s=""):
        print(s)
        sys.stdout.flush()
        lines.append(s)

    ns, ex, star, scat, prev56 = read_ladder()
    say("read %d uniform rungs at cap %d, their 0.56 crossing %.4f, "
        "and the" % (len(ns), CAP, prev56))
    say("  ladder's scatter %.4f, from results/audit_ladder_cap.txt "
        "and" % scat)
    say("  results/audit_primorial_rung17.txt")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]

    # -------------------------------------------------------------- C1
    say()
    say("C1  the sieve control at n = %d, on odd indices" % SIEVECHECK)
    la, ma = R11.lambda_and_mu(SIEVECHECK)
    kd, mb = kind_and_mu_odd(SIEVECHECK)
    pv, plg = power_table(SIEVECHECK)
    odd = np.arange(1, SIEVECHECK + 1, 2, dtype=np.int64)
    hv = odd >> 1
    dmu = int((ma[odd] != mb[hv]).sum())
    reb = np.zeros(odd.size, dtype=np.float64)
    k1 = kd[hv] == KIND_PRIME
    reb[k1] = np.log(odd[k1].astype(np.float64))
    k2 = kd[hv] == KIND_POWER
    if k2.any():
        reb[k2] = plg[np.searchsorted(pv, odd[k2])]
    dlam = int((la[odd] != reb).sum())
    nz = int((la[odd] != 0).sum())
    c1 = dmu == 0 and dlam == 0
    say("  Moebius disagreements %d, Lambda disagreements %d, over "
        "%d odd indices" % (dmu, dlam, odd.size))
    say("  Lambda support on the odd %d: %d prime, %d prime power, "
        "and the table holds %d"
        % (nz, int(k1.sum()), int(k2.sum()), pv.size))
    say("  the production route is the one code/audit_sieve.py "
        "compared with")
    say("  explicit factorisation (W1, W2), so agreement here is "
        "agreement with that")
    say("  C1 %s   (cap: elementwise equality, no tolerance)"
        % ("hold" if c1 else "REFUTED"))
    del la, ma, reb, k1, k2

    # -------------------------------------------------------------- C2
    say()
    say("C2  the mask control at n = %d" % SIEVECHECK)
    m9 = R17.residue_mask_odd(SIEVECHECK, qs)
    m5 = mask5_odd(SIEVECHECK, qs)
    PNc = R11.factor_set(NEW)
    bad, tried = 0, 0
    for k in range(2, 1000):
        if any(k % q == 0 for q in PNc):
            continue
        if mb[k >> 1] == 0:
            continue
        tried += 1
        kb9 = 0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb9 |= 1 << i
        a = (m9[hv] & np.uint16(~kb9 & 0xFFFF)) == 0
        b = (m5[hv] & keepbits5(k)) == 0
        bad += int((a != b).sum())
    c2 = bad == 0
    say("  %d admissible k below 1000, %d disagreements over %d odd v "
        "each" % (tried, bad, odd.size))
    say("  C2 %s   (cap: identical keep decisions, no tolerance)"
        % ("hold" if c2 else "REFUTED"))
    del m9, m5, kd, mb, odd, hv
    if not (c1 and c2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    say()
    say("the rungs: the controls %d and %d, and the new %d at "
        "log10 N = %.4f" % (CONTROL1, CONTROL2, NEW, math.log10(NEW)))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in R11.factor_set(N) if q > 2))
                  for N in (CONTROL1, CONTROL2, NEW))))
    say()
    say("sieving to %d on half-indices, sieve weight over the odd "
        "primes %s" % (NEW, ", ".join(map(str, qs))))
    kind, mu = kind_and_mu_odd(NEW)
    vmask = mask5_odd(NEW, qs)
    pv, plg = power_table(NEW)
    resident = kind.nbytes + mu.nbytes + vmask.nbytes
    half = NEW // 2 + 1
    say("BYTES resident_arrays %d" % resident)
    say("BYTES rung17_packing %d" % (4 * half + half + 2 * half))
    say("  the three resident arrays are %.2f GB, against %.2f GB for "
        "the packing of" % (resident / 2.0 ** 30,
                            (4 * half + half + 2 * half) / 2.0 ** 30))
    say("  {#rem:rung17} at this N, whose uint32 would not hold the "
        "prime here either")
    say("  the power table holds %d entries" % pv.size)
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    got = {}
    say()
    for N in (CONTROL1, CONTROL2, NEW):
        out = measure_kind(N, kind, mu, vmask, pv, plg, qs, artin,
                           twin, CAP)
        if out is None:
            say("  N = %-12d no crossing below k = %d" % (N, CAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e)
        say("  N = %-12d thr %.6f  #k %-7d beta %.6f  K*_R %-8d "
            "exp %.6f" % (N, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")
    if NEW not in got:
        say()
        say("the statistic has no value at rung 18 below the uniform "
            "cap; H1 to H4 are not evaluable")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)
    e18 = got[NEW][1]

    # -------------------------------------------------------------- C3
    say()
    say("C3  the ladder control at rungs 16 and 17")
    c3 = True
    for j, N in ((16, CONTROL1), (17, CONTROL2)):
        pubk, pube = star[N], ex[ns.index(N)]
        k_, e_ = got[N]
        ok = (k_ == pubk) and abs(e_ - pube) < 5e-7
        c3 = c3 and ok
        say("  rung %d  K*_R here %d against %d, exponent %.6f "
            "against %.6f, %s"
            % (j, k_, pubk, e_, pube, "equal" if ok else "DIFFERENT"))
    say("  C3 %s   (cap: exact on K*_R, six decimals on the exponent)"
        % ("hold" if c3 else "REFUTED"))

    # -------------------------------------------------------------- H1
    say()
    say("H1  does the margin keep growing?")
    marg = e18 - 0.5
    prev = ex[-1] - 0.5
    h1 = marg > prev
    say("  the new exponent is %.6f, margin %.4f, against rung 17's "
        "%.4f and the scatter %.4f" % (e18, marg, prev, scat))
    say("MARGIN audit_primorial_rung18 %.4f %.4f" % (marg, scat))
    if marg <= scat:
        say("INSIDE FLOOR audit_primorial_rung18")
    say("FLOOR primorial_rung18 %.4f" % scat)
    say("  H1 %s   (cap: rung 17's margin)"
        % ("hold" if h1 else "REFUTED"))

    # --------------------------------------------------------- H2, H3
    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    c, cov, s2, rms = quadfit(x, y)
    a, b = np.polyfit(x, y, 1)
    xn = math.log(NEW)
    v = np.array([1.0, xn, xn * xn])
    pq = float(v.dot(c))
    sp = math.sqrt(s2 + float(v.dot(cov).dot(v)))
    pl = a * xn + b
    dq, dl = e18 - pq, e18 - pl
    h2 = abs(dq) <= sp
    h3 = abs(dq) < abs(dl)
    say()
    say("H2/H3  does the curvature predict this one too?")
    say("  fitted on the %d uniform rungs:" % len(ns))
    say("  shape        predicts   measured   departure   pred s.e.  "
        " ratio")
    say("  quadratic    %-10.4f %-10.4f %+-11.4f %-11.4f %.2f"
        % (pq, e18, dq, sp, abs(dq) / sp))
    say("  line         %-10.4f %-10.4f %+-11.4f" % (pl, e18, dl))
    say("  the two shapes are %.4f apart here, %.2f prediction "
        "standard errors" % (abs(pq - pl), abs(pq - pl) / sp))
    say("  H2 %s   (cap: the prediction standard error)"
        % ("hold" if h2 else "REFUTED"))
    say("  H3 %s   (cap: the line's departure)"
        % ("hold" if h3 else "REFUTED"))

    # -------------------------------------------------------------- H4
    say()
    say("H4  is 0.56 crossed here?")
    h4 = e18 >= TARGET
    say("  the exponent is %.6f against the level %.2f, %s by %.6f"
        % (e18, TARGET, "OVER" if h4 else "short", abs(e18 - TARGET)))
    say("MARGIN audit_rung18_target %.6f %.4f"
        % (abs(e18 - TARGET), scat))
    if abs(e18 - TARGET) <= scat:
        say("INSIDE FLOOR audit_rung18_target")
    say("FLOOR rung18_target %.4f" % scat)
    say("  the quadratic predicted %.4f and the line %.4f, so the two "
        "shapes" % (pq, pl))
    say("  answered this question differently and the measurement "
        "settles it")
    say("  H4 %s   (cap: the level 0.56)"
        % ("hold" if h4 else "REFUTED"))
    say("SHAPES 2")

    # ---------------------------------------------- bracketing 0.56
    say()
    say("is 0.56 bracketed by measurement?")
    below = ex[-1]
    lo_n = math.log10(ns[-1])
    hi_n = math.log10(NEW)
    both = below < TARGET <= e18
    say("  rung 17 at log10 N = %.4f reads %.6f, rung 18 at %.4f "
        "reads %.6f" % (lo_n, below, hi_n, e18))
    if both:
        say("  the level is between them: bracketed by the sign of "
            "two measurements,")
        say("  in (10^%.4f, 10^%.4f], a width of %.4f decades"
            % (lo_n, hi_n, hi_n - lo_n))
        say("BRACKET theta_prime_056_measured %.4f %.4f %.4f"
            % ((lo_n + hi_n) / 2.0, lo_n, hi_n))
        say("DRIFT theta_prime_056_measured %.4f"
            % abs((lo_n + hi_n) / 2.0 - prev56))
        say("SCATTER slope_audit_primorial_rung18 %.4f" % rms)
    else:
        say("  the level is not between them; no measured bracket "
            "follows from this pair")
    say("SCALES 2")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v_ else "REFUTED"
                for v_ in (c1, c2, c3, h1, h2, h3, h4)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^18 = 7872184320 and, as controls,",
        "           at N = 30030*2^16 and 30030*2^17; the margin over",
        "           1/2 against rung 17's; the exponent against the",
        "           level 0.56, which the quadratic and the line",
        "           fitted on the eighteen uniform rungs answer",
        "           differently; and whether rungs 17 and 18 put the",
        "           level between them.  C1 compares the kind-byte",
        "           sieve with the production one at every odd index",
        "           below 20000000 and C2 compares the five-bit mask",
        "           with the nine-bit one.  The beta fit and the",
        "           truncation search are both at the uniform cap",
        "           1000000.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py.",
        "FIELD: N = 1968046080, 3936092160 and 7872184320, the odd",
        "       radical 3*5*7*11*13 fixed so the threshold is",
        "       constant; k squarefree and coprime to N with",
        "       2 <= k < 1000000, beta fitted on the same range; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd radical,",
        "       as the RADICALS line declares. The k-cap, the sieve",
        "       weight and the Euler bound are imported from",
        "       code/audit_primorial_rung11.py, whose sieve C1",
        "       compares against, and the nine-bit mask from",
        "       code/audit_primorial_rung17.py; the eighteen uniform",
        "       rungs come from results/audit_ladder_cap.txt and",
        "       results/audit_primorial_rung17.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (c1 and c2 and c3):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
