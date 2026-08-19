# -*- coding: utf-8 -*-
r"""
The eighteenth rung: the one that lands on 0.56.

WHAT IS AT STAKE

[rem:laddercap] left a uniform ladder of seventeen rungs at cap 10^6
whose quadratic reaches 0.56 at log10 N = 9.6200 with bracket
[9.4574, 9.8277], and observed that rung 17 sits at 9.5951 and rung
18 at 9.8961 -- the crossing between them.  This is rung 17.

Fitted on those seventeen rungs the quadratic predicts 0.5592 here
against a prediction standard error of 0.0044.  **That is 0.17
standard errors below 0.56.**  The fit cannot tell whether this rung
clears the level or not; it puts the rung on it.  So the measurement
decides something the fit is not entitled to decide, which is the
only kind of statement about theta' that {#rem:shapepower} does not
forbid -- [rem:primorialgap] did exactly this for 1/2, bracketing the
crossing by the sign of two measurements rather than by a shape.

The line predicts 0.5491 on the same seventeen, so the shapes are
0.0101 apart, 2.29 prediction standard errors.  Rung 16 was the first
rung that could discriminate at 1.97; this one discriminates harder.

THE MEMORY, AND WHY THIS ROUTE IS NEW

Rung 17 does not fit in what rung 16's packing needs.  The three
resident arrays there are 23.62 GB at this N, and the machine's
commit headroom at the time of writing was 19.7 GB -- less than its
free physical memory, which is what the first rung-16 attempt died
on.  The route here halves it again by the same observation that
already halved the residue mask: **N is even, k is coprime to N and
hence odd, and m is odd, so N - mk is odd at every index the
statistic reads.**  The even half of Lambda's support and of Moebius
is never addressed either, so neither is stored.  Everything is kept
on half-indices, v -> v >> 1, and the block sieve runs in that space:
the odd multiples of an odd p are p(2s+1), whose half-indices start
at (p-1)/2 and step by p.

Three arrays of 7.87, 1.97 and 3.94 GB, 13.78 GB in all.  The values
are the values rung 16 used -- C1 is what says so, elementwise, at
every odd index below 2*10^7.

There is one place the addressing bites.  The k-loop tests
squarefreeness through Moebius before it tests coprimality to N, and
with Moebius stored on odd indices only the test must come second;
since 2 divides N, the coprimality test removes every even k first,
so the reordering is exact rather than merely safe.

BACKS: Remark {#rem:rung17} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The sieve control.  At n = 20000000 the half-index route
      reproduces the production route of
      code/audit_primorial_rung11.py at every odd index: Moebius
      equal, and the logarithm of the stored prime equal to the
      stored Lambda, bit for bit.
  C2  The ladder control.  At cap 10^6 rungs 15 and 16 return the
      K*_R that results/audit_ladder_cap.txt printed for them, 75253
      and 126079, exactly, and their exponents to the six decimals
      printed there.
  H1  The margin keeps growing: the margin over 1/2 exceeds rung
      16's 0.0488 on the uniform ladder.
  H2  The curvature predicts a fifth time: the departure from the
      quadratic fitted on the seventeen uniform rungs is inside that
      fit's prediction standard error 0.0044.
  H3  And it still beats the line: the quadratic's departure is
      smaller in absolute value than the line's.
  H4  **0.56 is not yet crossed.** The exponent at rung 17 is below
      0.56, as the seventeen-rung quadratic says by 0.17 of its own
      prediction error.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a single odd index where either array disagrees.
      One is an integer array and the other is the same arithmetic
      on the same inputs; there is no tolerance to argue about.  If
      C1 fails nothing below is a measurement.  THIS ONE GATES.
  C2  REFUTED by a single K*_R that differs, or an exponent that
      differs in the six decimals printed.  The arrays are addressed
      differently from the run that produced those numbers, so this
      is the check that the addressing does not move them.  THIS ONE
      GATES.
  H1  REFUTED if the margin does not grow.  Seven rungs have grown
      in a row; an eighth that does not would end the escalation.
  H2  REFUTED if the departure exceeds the prediction standard
      error.  Four out-of-sample hits are what the quadratic's
      standing rests on, and this is the longest reach yet asked of
      it.
  H3  REFUTED if the line is closer.  At 2.29 prediction standard
      errors apart this is the sharpest discrimination the ladder
      has offered.
  H4  REFUTED if the exponent reaches 0.56.  **This is the
      refutation worth having.**  It would mean the ladder crosses
      0.56 at or below log10 N = 9.5951 -- measured, not
      extrapolated -- one rung earlier than the uniform fit places
      the crossing, and the branch would hold its first observed
      crossing of the level rather than a fitted one.  Either
      outcome is a measurement; neither is a forecast for theta'.

  C1 and C2 gate.  H1 to H4 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  A deterministic curve is located
  against a computed threshold; there is no background to detect
  against.  The coin arms for this statistic were run in
  lab_primorial_ladder.py and lab_primorial_share.py, and the
  scatter they left is the floor the margin is judged against.
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
OUT = os.path.join(RES, "audit_primorial_rung17.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
SEED = 20260823
DRAWS = 4000
TARGET = 0.56
BLOCK = 1 << 24                     # the block, in half-indices
SIEVECHECK = 20_000_000             # where C1 compares the two routes
CAP = 1_000_000                     # the uniform cap of {#rem:laddercap}


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL1 = BASE * (1 << 15)         # 984023040, the rung 15 point
CONTROL2 = BASE * (1 << 16)         # 1968046080, the rung 16 point
NEW = BASE * (1 << 17)              # 3936092160


def prime_and_mu_odd(n, block=BLOCK):
    """Lambda's prime and Moebius on odd indices, addressed by v >> 1

    Returns (lp, mu) of length n//2 + 1 with lp[v >> 1] = p when the
    odd v is a power of the prime p and 0 otherwise, and
    mu[v >> 1] = mu(v).  Even indices are never stored because the
    statistic never reads them: N is even and k and m are odd, so
    N - mk is odd.

    The block sieve runs in half-index space.  The odd multiples of
    an odd p are p(2s+1), whose half-indices are ps + (p-1)/2 -- an
    arithmetic progression of step p starting at (p-1)/2.
    """
    if n >= 2 ** 32:
        raise ValueError("lp is uint32; this route stops below 2^32")
    half = n // 2 + 1
    root = int(math.isqrt(n))
    pr = [int(p) for p in primes_upto(root) if p > 2]
    lp = np.zeros(half, dtype=np.uint32)
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
            lp[t0 + idx] = vals[idx].astype(np.uint32)
        del vals, rem, m, big, idx
    for p in pr:
        lp[p >> 1] = p
        q = p * p
        while q <= n:
            lp[q >> 1] = p
            if q > n // p:
                break
            q *= p
    return lp, mu


def residue_mask_odd(n, qs):
    """bit i of mask[v >> 1] is set exactly when qs[i] divides odd v"""
    m = np.zeros(n // 2 + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[(q - 1) // 2::q] |= np.uint16(1 << i)
    return m


def measure_odd(N, lp, mu, vmask, qs, artin, twin, cap, block=BLOCK):
    """the ladder's K*_R with both the fit and the search at one cap

    Identical to the statistic of {#rem:laddercap}; only the arrays
    are addressed on half-indices.  The coprimality test to N runs
    before the squarefree test, because Moebius is not stored at even
    k -- and since 2 divides N the first test removes them all.
    """
    PN = R11.factor_set(N)
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
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        keepbits = np.uint16(~kb & 0xFFFF)
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
            g = mu[ms >> 1].astype(np.float64)
            pv = lp[vals >> 1]
            nz = pv != 0
            if nz.any():
                h += float((np.log(pv[nz].astype(np.float64))
                            * g[nz]).sum())
            keep = (vmask[vals >> 1] & keepbits) == 0
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


def read_uniform():
    """the seventeen uniform rungs at cap 10^6, the scatter, and the
    crossing they put 0.56 at"""
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
    s15 = io.open(os.path.join(RES, "audit_primorial_rung15.txt"),
                  encoding="utf-8").read()
    scat = float(re.search(r"^FLOOR primorial_rung15 ([\d.]+)\s*$",
                           s15, re.M).group(1))
    prev = float(re.search(r"^BRACKET ladder_cap_theta_prime "
                           r"([\d.]+) [\d.]+ [\d.]+\s*$",
                           src, re.M).group(1))
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

    ns, ex, star, scat, prev56 = read_uniform()
    say("read %d uniform rungs at cap %d from "
        "results/audit_ladder_cap.txt," % (len(ns), CAP))
    say("  their 0.56 crossing %.4f, and the ladder's scatter %.4f "
        "from" % (prev56, scat))
    say("  results/audit_primorial_rung15.txt" % ())
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)

    # -------------------------------------------------------------- C1
    say()
    say("C1  the sieve control at n = %d, on odd indices" % SIEVECHECK)
    la, ma = R11.lambda_and_mu(SIEVECHECK)
    pb, mb = prime_and_mu_odd(SIEVECHECK)
    odd = np.arange(1, SIEVECHECK + 1, 2, dtype=np.int64)
    dmu = int((ma[odd] != mb[odd >> 1]).sum())
    pv = pb[odd >> 1]
    reb = np.zeros(odd.size, dtype=np.float64)
    on = pv != 0
    reb[on] = np.log(pv[on].astype(np.float64))
    dlam = int((la[odd] != reb).sum())
    nz = int((la[odd] != 0).sum())
    non = int(on.sum())
    del la, ma, pb, mb, pv, reb, on
    c1 = dmu == 0 and dlam == 0
    say("  Moebius disagreements %d, Lambda disagreements %d, over "
        "%d odd indices" % (dmu, dlam, odd.size))
    say("  Lambda support on the odd %d against the stored primes %d"
        % (nz, non))
    say("  the production route is the one code/audit_sieve.py "
        "compared with")
    say("  explicit factorisation (W1, W2), so agreement here is "
        "agreement with that")
    say("  C1 %s   (cap: elementwise equality, no tolerance)"
        % ("hold" if c1 else "REFUTED"))
    del odd
    if not c1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    say()
    say("the rungs: the controls %d and %d, and the new %d at "
        "log10 N = %.4f" % (CONTROL1, CONTROL2, NEW, math.log10(NEW)))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in R11.factor_set(N) if q > 2))
                  for N in (CONTROL1, CONTROL2, NEW))))

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d on half-indices, sieve weight over the odd "
        "primes %s" % (NEW, ", ".join(map(str, qs))))
    lp, mu = prime_and_mu_odd(NEW)
    vmask = residue_mask_odd(NEW, qs)
    resident = lp.nbytes + mu.nbytes + vmask.nbytes
    say("BYTES resident_arrays %d" % resident)
    say("BYTES wholeindex_arrays %d"
        % (4 * (NEW + 1) + (NEW + 1) + 2 * (NEW // 2 + 1)))
    say("  the three resident arrays are %.2f GB against the %.2f GB "
        "the whole-index" % (resident / 2.0 ** 30,
                             (4 * (NEW + 1) + (NEW + 1)
                              + 2 * (NEW // 2 + 1)) / 2.0 ** 30))
    say("  packing of {#rem:rung16} would hold at this N")
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
        out = measure_odd(N, lp, mu, vmask, qs, artin, twin, CAP)
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
        say("the statistic has no value at rung 17 below the uniform "
            "cap; H1 to H4 are not evaluable")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)
    e17 = got[NEW][1]

    # -------------------------------------------------------------- C2
    say()
    say("C2  the ladder control at rungs 15 and 16")
    c2 = True
    for j, N in ((15, CONTROL1), (16, CONTROL2)):
        pubk, pube = star[N], ex[ns.index(N)]
        k_, e_ = got[N]
        ok = (k_ == pubk) and abs(e_ - pube) < 5e-7
        c2 = c2 and ok
        say("  rung %d  K*_R here %d against %d, exponent %.6f "
            "against %.6f, %s"
            % (j, k_, pubk, e_, pube, "equal" if ok else "DIFFERENT"))
    say("  C2 %s   (cap: exact on K*_R, six decimals on the exponent)"
        % ("hold" if c2 else "REFUTED"))

    # -------------------------------------------------------------- H1
    say()
    say("H1  does the margin keep growing?")
    marg = e17 - 0.5
    prev = ex[-1] - 0.5
    h1 = marg > prev
    say("  the new exponent is %.6f, margin %.4f, against rung 16's "
        "%.4f and the scatter %.4f" % (e17, marg, prev, scat))
    say("MARGIN audit_primorial_rung17 %.4f %.4f" % (marg, scat))
    if marg <= scat:
        say("INSIDE FLOOR audit_primorial_rung17")
    say("FLOOR primorial_rung17 %.4f" % scat)
    say("  H1 %s   (cap: rung 16's margin)"
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
    dq, dl = e17 - pq, e17 - pl
    h2 = abs(dq) <= sp
    h3 = abs(dq) < abs(dl)
    say()
    say("H2/H3  does the curvature predict this one too?")
    say("  fitted on the %d uniform rungs:" % len(ns))
    say("  shape        predicts   measured   departure   pred s.e.  "
        " ratio")
    say("  quadratic    %-10.4f %-10.4f %+-11.4f %-11.4f %.2f"
        % (pq, e17, dq, sp, abs(dq) / sp))
    say("  line         %-10.4f %-10.4f %+-11.4f" % (pl, e17, dl))
    say("  the two shapes are %.4f apart here, %.2f prediction "
        "standard errors" % (abs(pq - pl), abs(pq - pl) / sp))
    say("  H2 %s   (cap: the prediction standard error)"
        % ("hold" if h2 else "REFUTED"))
    say("  H3 %s   (cap: the line's departure)"
        % ("hold" if h3 else "REFUTED"))

    # -------------------------------------------------------------- H4
    say()
    say("H4  is 0.56 crossed here?")
    h4 = e17 < TARGET
    say("  the exponent is %.6f against the level %.2f, %s by %.6f"
        % (e17, TARGET, "short" if h4 else "OVER", abs(e17 - TARGET)))
    say("MARGIN audit_rung17_target %.6f %.4f"
        % (abs(e17 - TARGET), scat))
    if abs(e17 - TARGET) <= scat:
        say("INSIDE FLOOR audit_rung17_target")
    say("FLOOR rung17_target %.4f" % scat)
    say("  the seventeen-rung quadratic predicted %.4f, which is "
        "%.2f prediction" % (pq, abs(TARGET - pq) / sp))
    say("  standard errors from the level, so the fit could not call "
        "this one")
    say("  H4 %s   (cap: the level 0.56)"
        % ("hold" if h4 else "REFUTED"))

    # ------------------------------------------------- the crossing
    say()
    say("where the crossing stands on eighteen rungs")
    x18 = np.append(x, xn)
    y18 = np.append(y, e17)
    c18, cov18, s218, rms18 = quadfit(x18, y18)
    p18 = cross(c18, TARGET)
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(c18, cov18, size=DRAWS)
    vv = [cross(dd, TARGET) for dd in draws]
    vv = [w / math.log(10.0) for w in vv
          if w is not None and w > x18.max()]
    here = p18 / math.log(10.0) if p18 else None
    say("  the eighteen-rung quadratic is %+.8f in (log N)^2, r.m.s. "
        "%.4f" % (c18[2], rms18))
    if here is None or not vv:
        say("  the refitted quadratic has no crossing above the top "
            "rung; the level is behind us")
        say("SHAPES 2")
        say("SCATTER slope_audit_primorial_rung17 %.4f" % rms18)
    else:
        lo = float(np.percentile(vv, 2.5))
        hi = float(np.percentile(vv, 97.5))
        say("  it reaches 0.56 at log10 N = %.4f, bracket [%.4f, "
            "%.4f] from %d of %d draws"
            % (here, lo, hi, len(vv), DRAWS))
        say("CENSORED ladder_quadratic18_theta_prime %d %d"
            % (DRAWS - len(vv), DRAWS))
        if DRAWS - len(vv):
            say("  %d draws put the crossing at or below the top rung "
                "and are dropped;" % (DRAWS - len(vv)))
            say("  the bracket is conditioned on the rest")
            say("TRUNCATION BIAS ladder_quadratic18_theta_prime")
        say("BRACKET ladder_quadratic18_theta_prime %.4f %.4f %.4f"
            % (here, lo, hi))
        say("DRIFT ladder_quadratic18_theta_prime %.4f"
            % abs(here - prev56))
        say("SHAPES 2")
        say("SCATTER slope_audit_primorial_rung17 %.4f" % rms18)
        say("  the seventeen-rung value was %.4f, from "
            "results/audit_ladder_cap.txt" % prev56)
    say("  rung 18 would sit at log10 N = %.4f, and needs a uint32 "
        "wider than this route has"
        % math.log10(BASE * (1 << 18)))
    say("  no forecast is made from this; {#rem:shapepower} is why.")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v_ else "REFUTED"
                for v_ in (c1, c2, h1, h2, h3, h4)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^17 = 3936092160 and, as controls,",
        "           at N = 30030*2^15 and 30030*2^16; the margin over",
        "           1/2 against rung 16's; the exponent against the",
        "           level 0.56; the quadratic and the line fitted on",
        "           the seventeen uniform rungs of",
        "           results/audit_ladder_cap.txt and asked for this",
        "           one; and the eighteen-rung quadratic's 0.56",
        "           crossing with a bracket from its own parameter",
        "           covariance.  C1 compares the half-index sieve",
        "           with the production one at every odd index below",
        "           n = 20000000.  The beta fit and the truncation",
        "           search are both at the uniform cap 1000000.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py. The bracket is drawn from the",
        "      fit's own parameter covariance with the fixed SEED.",
        "FIELD: N = 984023040, 1968046080 and 3936092160, the odd",
        "       radical 3*5*7*11*13 fixed so the threshold is",
        "       constant; k squarefree and coprime to N with",
        "       2 <= k < 1000000, beta fitted on the same range; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd radical,",
        "       as the RADICALS line declares. The k-cap, the sieve",
        "       weight and the Euler bound are imported from",
        "       code/audit_primorial_rung11.py, whose sieve C1",
        "       compares against; the seventeen uniform rungs come",
        "       from results/audit_ladder_cap.txt and the scatter",
        "       from results/audit_primorial_rung15.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (c1 and c2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
