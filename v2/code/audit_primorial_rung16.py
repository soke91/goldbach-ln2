# -*- coding: utf-8 -*-
r"""
The seventeenth rung: the wall is real, and it is not the one named.

WHAT IS AT STAKE

Remark {#rem:rung15} ends by naming the limit: "the sieve at 2^16
needs sixteen gigabytes for Lambda alone and 2^17 twice that, so the
next doubling is at the memory wall of this machine and the one after
is past it."  Both halves of that sentence are wrong, and they are
wrong in opposite directions.

**The memory is not the wall.**  The routine the ladder has used
since rung 11 holds four arrays of order N where the information fits
in about half that, and three reductions recover it without moving a
number:

  (i)  Lambda is stored as the *prime*, not its logarithm.  On the
       support of Lambda every value is log p for a single p, and
       p < 2^32 for every N this ladder can reach, so a uint32 array
       carries the same information in half the bytes.  The logarithm
       is taken at the point of use, of the same float64 integer, so
       the float that enters the sum is the identical float.
  (ii) Moebius is built blockwise.  The cofactor trick of
       code/audit_primorial_rung11.py -- which code/audit_sieve.py
       compared against explicit factorisation -- needs an int32
       array of order N only because it is run on the whole range at
       once.  Run on a block the remainder is block-local, and the
       full list of primes below N is not needed either: the entries
       of a block that survive division by every prime below sqrt(N)
       are exactly 1 and the primes above sqrt(N), which is where the
       rest of Lambda's support lies.
  (iii) The residue mask is stored on odd indices only.  N is even
       and k is coprime to N, hence odd, and m is odd, so N - mk is
       odd at every index the statistic ever reads.  Half the mask
       was never addressable.

The squarefree flag array is dropped as well; it was a copy of
"Moebius is nonzero".  The per-k sums are streamed over blocks so
that no temporary of order N/k is alive either.  Rung 16's arrays are
half of what the production route would hold, and the BYTES lines
below are the measurement.

**The k-cap is the wall, and it binds here.**  Run at the ladder's
published cap, rung 16 has no crossing below k = 100000: the
truncation K*_R has left the range the ladder searches.  It was
always going to.  K*_R reads 43171 at rung 14 and 72857 at rung 15,
a factor of 1.69 per doubling against a cap that does not move, so
rung 15 was already at 73% of it.  **The ladder stopped one rung
short of its own definition running out, and the remark attributed
the stop to the wrong resource.**

That cap cannot simply be raised.  It is not a search bound: beta is
fitted by least squares over the same k-range, so widening the range
moves beta, moves R = H - beta*P, and moves K*_R at every rung
already published.  Measured in exploration at the two rungs where it
is cheap, doubling the cap moves rung 11's exponent from 0.5099 to
0.5102 and rung 12's from 0.5178 to 0.5180 -- below the ladder's
scatter 0.0037, but above the fourth decimal the ladder prints.

So the definition is extended in the one direction that leaves every
published number untouched: **beta keeps the published estimation
window k < 100000, and only the truncation search is widened.**  On
every published rung the threshold is reached below 100000, so the
fitted beta is the same fit and the cumulative sum below the cap is
the same sum, and the integer is the same integer.  C2 is what says
so, and C3 declares the size of the convention by refitting beta on
the wider window and reporting what that would have done instead.

**And the rung is the one that separates the two shapes.**
30030*2^16 = 1968046080 at log10 N = 9.2940.  Fitted on the sixteen
published rungs the quadratic predicts 0.5499 there and the line
0.5408 -- a gap of 0.0091 against a prediction standard error of
0.0046.  Every earlier rung was predicted by both shapes to within
about one standard error of each other; this is the first one where
they are two apart, so it is the first rung that can choose.

What is downstream of it is worth stating too, because it is close.
The sixteen-rung quadratic puts the 0.56 crossing at log10 N = 9.6358.
Rung 17 sits at 9.5951 and rung 18 at 9.8961.  The crossing is
between them.  If the ladder can be carried two more rungs it stops
being a crossing that is extrapolated to and becomes one that is
bracketed by measurement, which is a different kind of statement and
the only kind {#rem:shapepower} does not forbid.

BACKS: Remark {#rem:rung16} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1 to H4 were registered before any rung-16 exponent existed.  The
  run that established the k-cap binds produced no exponent at rung
  16 -- it reported no crossing and stopped -- so nothing below has
  been seen at the time of writing.

  C1  The sieve control.  At n = 20000000 the blockwise route
      reproduces the production route of
      code/audit_primorial_rung11.py: Moebius elementwise equal, and
      the logarithm of the stored prime equal to the stored Lambda --
      not to a tolerance, bit for bit, because both are the logarithm
      of the same float64 integer.
  C2  The statistic control.  With the arrays built by the new route,
      the sums streamed, and the truncation search widened while
      beta keeps the published window, rungs 14 and 15 reproduce the
      published K*_R of results/audit_primorial_rung14.txt and
      rung15.txt exactly, and their exponents inside the bound their
      printing forces.
  C3  The convention control.  Refitting beta on the widened window
      instead moves the rung 14 and rung 15 exponents by less than
      the ladder's scatter 0.0037, and moves rung 15 by less than
      rung 14 -- the convention costs less the further up it is
      applied.
  H1  The margin keeps growing: the new exponent's margin over 1/2
      exceeds rung 15's 0.0407.
  H2  The curvature predicts a fourth time: the departure from the
      quadratic fitted on the sixteen published rungs is inside that
      fit's prediction standard error at this abscissa.
  H3  And it still beats the line: the quadratic's departure is
      smaller in absolute value than the line's.
  H4  The crossing settles further: refitted on seventeen rungs the
      0.56 crossing moves less than the drift rung 15 declared,
      0.0290.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a single index where either array disagrees.  There
      is no tolerance to argue about: one is an integer array and the
      other is the same arithmetic on the same inputs.  If C1 fails
      nothing below is a measurement of anything.  THIS ONE GATES.
  C2  REFUTED outside the printing bound at either rung, or by a
      single K*_R that differs.  The sums are split differently from
      the published run and the search is wider, so this is the check
      that neither moves the integer the statistic is.  THIS ONE
      GATES.
  C3  REFUTED if either shift reaches the scatter, or if rung 15's
      shift is not smaller than rung 14's.  Then the window is not a
      convention with a decaying cost but a parameter the ladder
      depends on, and the extension below would have to be withdrawn
      in favour of recomputing every published rung at one cap.
      THIS ONE GATES.
  H1  REFUTED if the margin does not grow.  Six rungs have grown in a
      row; a seventh that does not would end the escalation and would
      be the first evidence in this branch that the rise is not
      monotone.
  H2  REFUTED if the departure exceeds the prediction standard error.
      That is the outcome that would matter most: three out-of-sample
      hits are what promoted the quadratic from a fit to a shape that
      predicts, and a miss at the first rung that can discriminate
      would withdraw that reading rather than weaken it.
  H3  REFUTED if the line is closer.  With the shapes 2.0 prediction
      standard errors apart, this is the first rung at which "closer"
      is not a coin flip, and a line that wins here would put the
      curvature back inside its own noise.
  H4  REFUTED if the crossing moves by more than 0.0290.  The
      settling of the extrapolation is one rung old; a second point
      is what says whether it settled or paused.

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
OUT = os.path.join(RES, "audit_primorial_rung16.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
SEED = 20260823
DRAWS = 4000
TARGET = 0.56
BLOCK = 1 << 24                     # the block, in indices
SIEVECHECK = 20_000_000             # where C1 compares the two routes
KSEARCH = 400_000                   # the widened truncation search


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL1 = BASE * (1 << 14)         # 492011520, the rung 14 point
CONTROL2 = BASE * (1 << 15)         # 984023040, the rung 15 point
NEW = BASE * (1 << 16)              # 1968046080


def prime_and_mu_block(n, block=BLOCK):
    """the support of Lambda as its prime, and Moebius, blockwise

    Returns (lp, mu) with lp[v] = p when v is a power of the prime p
    and 0 otherwise, so that Lambda(v) = log(lp[v]) on the support.
    Nothing of order n is alive but those two: the cofactor of the
    production route is block-local here, and the only prime list
    held is the one below sqrt(n).
    """
    if n >= 2 ** 32:
        raise ValueError("lp is uint32; this route stops below 2^32")
    root = int(math.isqrt(n))
    pr = [int(p) for p in primes_upto(root)]
    lp = np.zeros(n + 1, dtype=np.uint32)
    mu = np.empty(n + 1, dtype=np.int8)
    for lo in range(0, n + 1, block):
        hi = min(lo + block, n + 1)
        w = hi - lo
        vals = np.arange(lo, hi, dtype=np.int64)
        rem = vals.copy()
        m = np.ones(w, dtype=np.int8)
        for p in pr:
            s = (-lo) % p
            if s < w:
                m[s::p] = -m[s::p]
                rem[s::p] //= p
            q = p * p
            s = (-lo) % q
            if s < w:
                m[s::q] = 0
            pk = p * p
            while pk <= n:
                s = (-lo) % pk
                if s < w:
                    rem[s::pk] //= p
                if pk > n // p:
                    break
                pk *= p
        big = rem > 1
        m[big] = -m[big]
        mu[lo:hi] = m
        idx = np.flatnonzero(big & (rem == vals))
        if idx.size:
            lp[lo + idx] = vals[idx].astype(np.uint32)
        del vals, rem, m, big, idx
    mu[0] = 0
    for p in pr:
        lp[p] = p
        q = p * p
        while q <= n:
            lp[q] = p
            if q > n // p:
                break
            q *= p
    return lp, mu


def residue_mask_odd(n, qs):
    """bit i of mask[v >> 1] is set exactly when qs[i] divides odd v

    Every index the statistic reads is N - mk with N even and m, k
    odd, so only the odd half of the mask was ever addressable.  The
    odd multiples of an odd q are q(2t+1), whose half-indices start
    at (q-1)/2 and step by q.
    """
    m = np.zeros(n // 2 + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[(q - 1) // 2::q] |= np.uint16(1 << i)
    return m


def measure_block(N, lp, mu, vmask, qs, artin, twin, kfit, ksearch,
                  block=BLOCK):
    """K*_R with beta fitted on k < kfit and the search run to ksearch

    The published statistic takes both bounds equal.  They are split
    here because K*_R has left the published range at rung 16 and the
    cap is not a search bound: beta is a least-squares fit over the
    same k, so widening it moves R = H - beta*P at every rung already
    printed.  Fitting on the published window and searching past it
    is the one extension that leaves those integers alone, and C2 is
    what says it does.

    Both conventions are returned: `frozen` keeps the published
    window, `wide` refits beta on everything searched.  C3 reads the
    distance between them.

    The admissible m for each k are the same numbers in the same
    order and each term of each sum is the identical float; only the
    association is different, because the sum is split into blocks
    and the terms Lambda kills are not carried.
    """
    PN = R11.factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, ksearch):
        if mu[k] == 0:
            continue
        if any(k % q == 0 for q in PN):
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
            ms = ms[mu[ms] != 0]
            for q in drop:
                ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            seen += ms.size
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            pv = lp[vals]
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
    lk = np.log(ks.astype(np.float64))
    thr = S_ * (1.0 - A_) * N
    inw = ks < kfit

    out = {"thr": thr / N, "nfit": int(inw.sum()), "nsearch": ks.size}
    for tag, m in (("frozen", inw), ("wide", np.ones_like(inw))):
        beta = float((H[m] * P[m]).sum() / (P[m] * P[m]).sum())
        cum = np.cumsum(lk * np.abs(H - beta * P))
        j = int(np.searchsorted(cum, thr))
        out[tag] = None if j >= ks.size else (
            int(ks[j]), math.log(float(ks[j])) / math.log(N), beta)
    return out


def read_all():
    """the sixteen rung exponents, the published K*_R, and the crossing"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ns, ex, dec = [], [], 0
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
        dec = max(dec, len(f[2].split(".")[1]))
    star = {}
    for j in (11, 12, 13, 14, 15):
        s = io.open(os.path.join(RES, "audit_primorial_rung%d.txt" % j),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R (\d+)\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(2)))
        star[N] = int(m.group(1))
    s15 = io.open(os.path.join(RES, "audit_primorial_rung15.txt"),
                  encoding="utf-8").read()
    marg15 = float(re.search(r"the new exponent is [\d.]+, margin "
                             r"([\d.]+),", s15).group(1))
    scat = float(re.search(r"^FLOOR primorial_rung15 ([\d.]+)\s*$",
                           s15, re.M).group(1))
    m = re.search(r"^BRACKET ladder_quadratic16_theta_prime ([\d.]+) "
                  r"([\d.]+) ([\d.]+)\s*$", s15, re.M)
    dr = float(re.search(r"^DRIFT ladder_quadratic16_theta_prime "
                         r"([\d.]+)\s*$", s15, re.M).group(1))
    return (ns, ex, dec, star, marg15, scat, float(m.group(1)),
            float(m.group(2)), float(m.group(3)), dr)


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

    (ns, ex, dec, star, marg15, scat, pt56, lo56, hi56,
     drift56) = read_all()
    say("read %d rung exponents, rung 15's margin %.4f, the scatter "
        "%.4f," % (len(ns), marg15, scat))
    say("  and the quadratic's 0.56 crossing %.4f [%.4f, %.4f] with "
        "drift %.4f" % (pt56, lo56, hi56, drift56))
    say("  from results/audit_primorial_rung10 through rung15")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)

    # -------------------------------------------------------------- C1
    say()
    say("C1  the sieve control at n = %d" % SIEVECHECK)
    la, ma = R11.lambda_and_mu(SIEVECHECK)
    pb, mb = prime_and_mu_block(SIEVECHECK)
    dmu = int((ma != mb).sum())
    rebuilt = np.zeros_like(la)
    on = pb != 0
    rebuilt[on] = np.log(pb[on].astype(np.float64))
    dlam = int((la != rebuilt).sum())
    nz = int((la != 0).sum())
    non = int(on.sum())
    del la, ma, pb, mb, rebuilt, on
    c1 = dmu == 0 and dlam == 0
    say("  Moebius disagreements %d, Lambda disagreements %d, over "
        "%d indices" % (dmu, dlam, SIEVECHECK + 1))
    say("  Lambda support %d against the stored primes %d" % (nz, non))
    say("  the production route is the one code/audit_sieve.py "
        "compared with")
    say("  explicit factorisation (W1, W2), so agreement here is "
        "agreement with that")
    say("  C1 %s   (cap: elementwise equality, no tolerance)"
        % ("hold" if c1 else "REFUTED"))
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
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NEW, ", ".join(map(str, qs))))
    say("  beta keeps the published window k < %d; the truncation "
        "search runs to k < %d" % (R11.KCAP, KSEARCH))
    lp, mu = prime_and_mu_block(NEW)
    vmask = residue_mask_odd(NEW, qs)
    resident = lp.nbytes + mu.nbytes + vmask.nbytes
    was = 8 * (NEW + 1) + 2 * (NEW + 1) + 2 * (NEW + 1)
    say("BYTES resident_arrays %d" % resident)
    say("BYTES production_arrays %d" % was)
    say("  the three resident arrays are %.2f GB against the four the "
        "production route" % (resident / 2.0 ** 30))
    say("  would hold, %.2f GB, and the int32 cofactor it builds them "
        "with, %.2f GB more"
        % (was / 2.0 ** 30, 4.0 * (NEW + 1) / 2.0 ** 30))
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
        out = measure_block(N, lp, mu, vmask, qs, artin, twin,
                            R11.KCAP, KSEARCH)
        got[N] = out
        if out["frozen"] is None:
            say("  N = %-12d no crossing below k = %d"
                % (N, KSEARCH))
            continue
        kstar, e, beta = out["frozen"]
        say("  N = %-12d thr %.6f  #k %-7d/%-7d beta %.6f  "
            "K*_R %-8d exp %.4f"
            % (N, out["thr"], out["nfit"], out["nsearch"], beta,
               kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, out["thr"]))
        if out["wide"] is not None:
            kw, ew, bw = out["wide"]
            say("      refitting beta on all %d gives beta %.6f, "
                "K*_R %d, exp %.4f"
                % (out["nsearch"], bw, kw, ew))
    say("RADICALS 1")
    if got[NEW]["frozen"] is None:
        say()
        say("the statistic has no value at rung 16 even at the "
            "widened search; H1 to H4 are not evaluable")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)
    e16 = got[NEW]["frozen"][1]

    # -------------------------------------------------------------- C2
    say()
    say("C2  the statistic control at rungs 14 and 15")
    rnd = 0.5 * 10.0 ** (-dec)
    say("PRINTBOUND audit_primorial_rung16 %d %.8f" % (dec, rnd))
    c2 = True
    shift = {}
    for j, N in ((14, CONTROL1), (15, CONTROL2)):
        pubk, pube = star[N], ex[ns.index(N)]
        k_, e_, _b = got[N]["frozen"]
        d = abs(e_ - pube)
        ok = d <= rnd and k_ == pubk
        c2 = c2 and ok
        say("  rung %d  exponent here %.4f against the published "
            "%.4f, departure %.6f" % (j, e_, pube, d))
        say("          K*_R here %d against the published %d, %s"
            % (k_, pubk, "equal" if k_ == pubk else "DIFFERENT"))
        shift[j] = abs(got[N]["wide"][1] - e_)
    say("  the bound from %d decimals is %.8f" % (dec, rnd))
    say("  C2 %s   (cap: the printing bound, and exact on K*_R)"
        % ("hold" if c2 else "REFUTED"))

    # -------------------------------------------------------------- C3
    say()
    say("C3  what the window convention costs")
    c3 = (shift[14] < scat and shift[15] < scat
          and shift[15] < shift[14])
    for j in (14, 15):
        say("  rung %d  refitting beta on the widened window moves the "
            "exponent by %.6f" % (j, shift[j]))
    say("WINDOW audit_primorial_rung16 %d %d %.6f %.6f"
        % (R11.KCAP, KSEARCH, shift[14], shift[15]))
    say("  against the ladder's scatter %.4f, and %s at rung 15 than "
        "at rung 14" % (scat, "smaller" if shift[15] < shift[14]
                        else "LARGER"))
    say("  C3 %s   (cap: the scatter, and monotone in N)"
        % ("hold" if c3 else "REFUTED"))

    # -------------------------------------------------------------- H1
    say()
    say("H1  does the margin keep growing?")
    marg = e16 - 0.5
    h1 = marg > marg15
    say("  the new exponent is %.4f, margin %.4f, against rung 15's "
        "%.4f and the scatter %.4f" % (e16, marg, marg15, scat))
    say("MARGIN audit_primorial_rung16 %.4f %.4f" % (marg, scat))
    if marg <= scat:
        say("INSIDE FLOOR audit_primorial_rung16")
    say("FLOOR primorial_rung16 %.4f" % scat)
    say("  H1 %s   (cap: rung 15's margin)"
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
    dq, dl = e16 - pq, e16 - pl
    h2 = abs(dq) <= sp
    h3 = abs(dq) < abs(dl)
    say()
    say("H2/H3  does the curvature predict this one too?")
    say("  fitted on the %d published rungs:" % len(ns))
    say("  shape        predicts   measured   departure   pred s.e.  "
        " ratio")
    say("  quadratic    %-10.4f %-10.4f %+-11.4f %-11.4f %.2f"
        % (pq, e16, dq, sp, abs(dq) / sp))
    say("  line         %-10.4f %-10.4f %+-11.4f" % (pl, e16, dl))
    say("  the two shapes are %.4f apart here, %.2f prediction "
        "standard errors" % (abs(pq - pl), abs(pq - pl) / sp))
    say("  H2 %s   (cap: the prediction standard error)"
        % ("hold" if h2 else "REFUTED"))
    say("  H3 %s   (cap: the line's departure)"
        % ("hold" if h3 else "REFUTED"))

    # -------------------------------------------------------------- H4
    say()
    say("H4  does the crossing settle further?")
    x17 = np.append(x, xn)
    y17 = np.append(y, e16)
    c17, cov17, s217, rms17 = quadfit(x17, y17)
    p17 = cross(c17, TARGET)
    moved = abs(p17 / math.log(10.0) - pt56)
    h4 = moved < drift56
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(c17, cov17, size=DRAWS)
    vals = [cross(dd, TARGET) for dd in draws]
    vals = [w / math.log(10.0) for w in vals
            if w is not None and w > x17.max()]
    lo = float(np.percentile(vals, 2.5))
    hi = float(np.percentile(vals, 97.5))
    say("  the seventeen-rung quadratic is %+.8f in (log N)^2, r.m.s. "
        "%.4f" % (c17[2], rms17))
    say("  it reaches 0.56 at log10 N = %.4f, bracket [%.4f, %.4f] "
        "from %d of %d draws"
        % (p17 / math.log(10.0), lo, hi, len(vals), DRAWS))
    say("BRACKET ladder_quadratic17_theta_prime %.4f %.4f %.4f"
        % (p17 / math.log(10.0), lo, hi))
    say("DRIFT ladder_quadratic17_theta_prime %.4f" % moved)
    say("  the sixteen-rung value was %.4f, so it moved %.4f against "
        "the declared drift %.4f" % (pt56, moved, drift56))
    say("SHAPES 2")
    say("SCATTER slope_audit_primorial_rung16 %.4f" % rms17)
    for j in (17, 18):
        say("  rung %d would sit at log10 N = %.4f"
            % (j, math.log10(BASE * (1 << j))))
    say("  H4 %s   (cap: the declared drift)"
        % ("hold" if h4 else "REFUTED"))
    say("  no forecast is made from this; {#rem:shapepower} is why.")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v_ else "REFUTED"
                for v_ in (c1, c2, c3, h1, h2, h3, h4)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^16 = 1968046080 and, as controls,",
        "           at N = 30030*2^14 and 30030*2^15; the margin over",
        "           1/2 against rung 15's; the quadratic and the line",
        "           fitted on the sixteen published rungs and asked",
        "           for this one, with the prediction standard error",
        "           at its abscissa and the gap between the two",
        "           shapes there; and the seventeen-rung quadratic's",
        "           0.56 crossing with a bracket from its own",
        "           parameter covariance, against the sixteen-rung",
        "           value.  C1 compares the blockwise sieve with the",
        "           production one elementwise at n = 20000000.  beta",
        "           is fitted on the published window k < 100000 and",
        "           the truncation search is run to k < 400000; C3",
        "           reports what refitting beta on the wider window",
        "           would move the exponent by.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py. The bracket is drawn from the",
        "      fit's own parameter covariance with the fixed SEED.",
        "FIELD: N = 492011520, 984023040 and 1968046080, the odd",
        "       radical 3*5*7*11*13 fixed so the threshold is",
        "       constant; k squarefree and coprime to N with",
        "       2 <= k < 400000, beta fitted on the 2 <= k < 100000",
        "       of the published ladder; m odd, squarefree and",
        "       coprime to k,",
        "       m < N/k; the sieve weight over the odd primes below",
        "       30; the Euler products at the fixed bound 4000000.",
        "       One odd radical, as the RADICALS line declares. The",
        "       k-cap, the sieve weight and the Euler bound are",
        "       imported from code/audit_primorial_rung11.py, whose",
        "       sieve C1 compares against and whose statistic C2",
        "       reproduces; the sixteen published rungs come from",
        "       results/audit_primorial_rung10.txt and rung11",
        "       through rung15, and the published crossing from",
        "       results/audit_primorial_rung15.txt.",
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
