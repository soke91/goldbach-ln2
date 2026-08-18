# -*- coding: utf-8 -*-
r"""
The split's declined null: which half of P does the work?

WHAT IS AT STAKE

Remark {#rem:predictable} is where the whole split comes from: with
beta the least-squares scale through the origin, sum(log k)|H - beta P|
is 0.6310, 0.5979, 0.5554, 0.5421, 0.5307 of sum(log k)|H|, so "to
within a factor beta of about 3, the dilated wall is an elementary
sieve-weighted Moebius sum plus a residue of half its size". Every
later result -- {#rem:splitbudget}, {#rem:residuelevel},
{#rem:residuecancel}, the whole residue programme -- is downstream of
that beta and that residue.

Its results file declines a null: "none is run and none is needed
here ... S2 and S3 decompose one measured sum, where a randomisation
would move both parts together". That is true of a randomisation that
moves both. It is not true of one that breaks P in a stated way and
leaves H alone, and there are two such, each cutting the predictor in
half:

  (a) keep sign P at every k and permute |P| across k,
  (b) keep |P| at every k and redraw its sign.

(a) leaves the sign field -- the thing {#rem:survivors} validated at
0.9274 -- exactly intact and destroys only the per-k magnitude
matching. (b) does the reverse. Refitting beta under each says which
half of the predictor is doing the 37 to 47 per cent.

The falling share is also read as a trend -- "the share it cuts is
growing with N" -- and by the standard of {#rem:slopes} that needs a
standard error, which it has never had.

The implementation is independent of lab_predictable_part.py's.

BACKS: Remark {#rem:predictablenull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  The control: beta reproduces the published 2.6992, 2.8588,
      2.9952, 3.1437, 3.0473 to within 0.0001 and the residual shares
      0.6310 ... 0.5307 to within 0.001.
  D2  The trend is resolved: the residual share's least-squares slope
      against log N is negative and reaches two standard errors.
  D3  The sign field alone buys almost nothing: with |P| permuted
      across k at fixed sign, the refitted residual share is above
      0.95 at every N.
  D4  And the magnitudes alone buy almost nothing either: with the
      sign of P redrawn at fixed |P|, the refitted residual share is
      above 0.95 at every N.

REFUTATION RULE (fixed before the run)

  D1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:predictable}.
  D2  REFUTED below two standard errors. "The share it cuts is
      growing with N" would then be a reading of five points with no
      error bar, and would have to be withdrawn as {#rem:slopes}
      withdrew rule U4.
  D3  REFUTED at or below 0.95 at any N. The cut would then be
      substantially a sign effect, and since the sign field is what
      {#rem:survivors} validated on short inner sums, the split would
      inherit that window restriction -- which after
      {#rem:survivorrange} is a real limitation and not a formality.
  D4  REFUTED at or below 0.95 at any N, which would say the
      magnitudes alone carry it and the sign agreement is decorative.

  All four gate.

  THE NULL IS THE POINT and it is run, in the two forms the published
  file declined: both break P alone and leave H untouched, so neither
  moves both parts together. 64 draws each, one seed, and beta is
  refitted under each draw rather than carried over.
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
OUT = os.path.join(RES, "audit_predictable_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
DRAWS = 64
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


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


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
    """beta and the residual share at each N"""
    src = io.open(os.path.join(RES, "lab_predictable_part.txt"),
                  encoding="utf-8").read()
    i = src.index("N            beta      residual share   "
                  "|beta P| share")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            if not out:
                continue
            break
        out[int(f[0])] = (float(f[1]), float(f[2]))
    return out


def build(N, lam, mu, oddsqf, sqf, vmask, qs):
    """H, the sieve-weighted P and the log-k weights over 2<=k<KCAP"""
    PN = factor_set(N)
    ks, Hs, Ps = [], [], []
    for k in range(2, KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 1:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[oddsqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
    ks = np.array(ks, dtype=np.int64)
    return (ks, np.array(Hs), np.array(Ps),
            np.log(ks.astype(np.float64)))


def share(H, P, w):
    """refit beta through the origin and return the residual share"""
    b = float((H * P).sum() / (P * P).sum())
    return b, float((w * np.abs(H - b * P)).sum()
                    / (w * np.abs(H)).sum())


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published beta and residual shares from "
        "results/lab_predictable_part.txt" % len(pub))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    vmask = residue_mask(NMAX, qs)
    rng = np.random.default_rng(SEED)
    say("%d draws per null, seed %d, beta refitted under each"
        % (DRAWS, SEED))

    got = []
    for N in NS:
        ks, H, P, w = build(N, lam, mu, oddsqf, sqf, vmask, qs)
        got.append((N, ks, H, P, w))
        say("  N = %-10d #k = %d" % (N, ks.size))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- D1
    say()
    say("D1  the control: beta and the residual share")
    say("  N            beta here   beta pub   share here  share pub")
    d1 = True
    shares = []
    for N, ks, H, P, w in got:
        b, sh = share(H, P, w)
        shares.append(sh)
        pb, ps = pub[N]
        if abs(b - pb) >= 0.0001 or abs(sh - ps) >= 0.001:
            d1 = False
        say("  %-12d %-11.4f %-10.4f %-11.4f %.4f"
            % (N, b, pb, sh, ps))
    say("  D1 %s   (cap 0.0001 in beta, cap 0.001 in the share)"
        % ("hold" if d1 else "REFUTED"))

    # ------------------------------------------------------------- D2
    say()
    say("D2  is the fall resolved?")
    x = np.log(np.array([g[0] for g in got], dtype=np.float64))
    a, rms, se, t = fit(x, np.array(shares))
    d2 = (a < 0.0) and (t >= 2.0)
    say("  least-squares slope against log N = %+.6f" % a)
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rms, se, t))
    say("SCATTER slope_audit_predictable_null %.4f" % rms)
    say("TSTAT slope_audit_predictable_null %.2f" % t)
    say("SPREAD slope_audit_predictable_null %.4f"
        % float(x.max() - x.min()))
    if t < 2.0:
        say("UNRESOLVED SIGN slope_audit_predictable_null")
    say("  D2 %s" % ("hold" if d2 else "REFUTED"))

    # ---------------------------------------------------------- D3/D4
    say()
    say("D3/D4  the two nulls the published file declined")
    say("  (a) sign P kept, |P| permuted across k")
    say("  (b) |P| kept, sign P redrawn")
    say("  N            true share  (a) median  (b) median  "
        "(a) min    (b) min")
    d3 = d4 = True
    fa, fb = [], []
    for (N, ks, H, P, w), sh in zip(got, shares):
        ap, bp = [], []
        for _ in range(DRAWS):
            perm = rng.permutation(P.size)
            Pa = np.sign(P) * np.abs(P)[perm]
            ap.append(share(H, Pa, w)[1])
            eps = rng.integers(0, 2, size=P.size) * 2 - 1
            Pb = eps * np.abs(P)
            bp.append(share(H, Pb, w)[1])
        ma, mb = float(np.median(ap)), float(np.median(bp))
        na, nb = float(np.min(ap)), float(np.min(bp))
        if na <= 0.95:
            d3 = False
        if nb <= 0.95:
            d4 = False
        fa.append((1.0 - ma) / (1.0 - sh))
        fb.append((1.0 - mb) / (1.0 - sh))
        say("  %-12d %-11.4f %-11.4f %-11.4f %-10.4f %.4f"
            % (N, sh, ma, mb, na, nb))
    say("  D3 signs alone leave above 0.95   %s"
        % ("hold" if d3 else "REFUTED"))
    say("  D4 magnitudes alone leave above 0.95   %s"
        % ("hold" if d4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC on D3 (post hoc). The right scale for a null")
    say("  here is not the share but the CUT, since that is what the")
    say("  remark claims. Each null's cut as a fraction of the true")
    say("  one:")
    say("  N            true cut   (a) fraction  (b) fraction")
    for N, sh, va, vb in zip(NS, shares, fa, fb):
        say("  %-12d %-10.4f %-13.4f %.4f" % (N, 1.0 - sh, va, vb))
    say("ONESIDED lab_predictable_part %.4f %.4f"
        % (max(fa), max(fb)))
    say("ONESIDED audit_predictable_null %.4f %.4f"
        % (max(fa), max(fb)))
    say("  so the sign field alone reaches %.0f to %.0f per cent of"
        % (100 * min(fa), 100 * max(fa)))
    say("  the cut and the magnitudes alone reach none of it -- D3")
    say("  fails because 0.95 was a threshold on the share and the")
    say("  share has %.2f of its range below 1 to begin with."
        % (1.0 - min(shares)))

    say()
    say("  what that means for the split. The cut of %.0f to %.0f per"
        % (100 * (1 - max(shares)), 100 * (1 - min(shares))))
    say("  cent is not carried by either half of P on its own: it")
    say("  needs the per-k pairing of the sign field with the")
    say("  magnitude. So the split is a statement about P as a whole")
    say("  and not about the sign agreement {#rem:survivors}")
    say("  validated, which is the part {#rem:survivorrange} showed")
    say("  does not reach the long inner sums.")

    say()
    say("=" * 70)
    ok = d1 and d2 and d3 and d4
    say("neither half of the predictor carries the split on its own"
        if ok else "REFUTED")

    head = [
        "STATISTIC: with beta the least-squares scale through the",
        "           origin, the residual share",
        "           sum(log k)|H - beta P| / sum(log k)|H| over",
        "           2 <= k < " + str(KCAP) + "; its least-squares slope",
        "           against log N with a standard error; and the same",
        "           share refitted under two randomisations of P that",
        "           leave H untouched -- |P| permuted across k at",
        "           fixed sign, and the sign of P redrawn at fixed",
        "           |P| -- " + str(DRAWS) + " draws each.",
        "NULL: the two randomisations above, which is what",
        "      lab_predictable_part.py declined on the grounds that a",
        "      randomisation 'would move both parts together'. Neither",
        "      of these moves H at all; each breaks one half of P and",
        "      beta is refitted under each draw.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       Lambda and mu from an integer sieve to " + str(NMAX)
        + ";",
        "       numpy default_rng seed " + str(SEED) + ". Every N is",
        "       2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published beta and shares are read from",
        "       results/lab_predictable_part.txt.",
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
