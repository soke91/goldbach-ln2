# -*- coding: utf-8 -*-
"""
Conjecture L's ORIGINAL half, at increment 283's precision (inc. 284).

WHAT WAS CLAIMED, AND HOW WELL IT WAS CHECKED. MEASUREMENTS section 7:

    mu-field = (deterministic local mask, computable by finite modular
    enumeration) x (exactly Gaussian fluctuation at half-normal scale)

with the stamp "E|C|^2/support = 0.97-1.02 with kurtosis 2.99-3.03:
exactly Gaussian" on ~12,000 pairs. At n = 12,000 the standard error of
an excess kurtosis is sqrt(24/n) = 0.045, so "3.03" means +0.03 +/-
0.045 -- consistent with Gaussian, and consistent with a good deal
else. Increment 283 reached +/-0.0017 on the wall and used that
precision to show the displayed formula was false at z = 98. The
original half has never been tested at anything like that resolution.

WHY IT MATTERS MORE THAN A PRECISION EXERCISE. Increment 283's finding
was not "the conjecture is wrong" but "the normaliser in the documents
is a FITTED STAND-IN for an exact quantity, and the fit is right in the
mean and wrong in the fluctuation". The same trap is available here.
For D(k) = Sum_m mu(m) mu(N-mk) the exact second moment is

    support(k) = Sum_m mu^2(m) mu^2(N-mk),

which is directly countable, exactly, per k. The alternative is the
MODELLED support of the mask (section 7 builds it from units mod q^2
for q <= 50 plus a tail factor 0.99228). If the campaign normalised by
the model where the exact count was available, this is the wall's story
again on the other half of the conjecture.

PRE-REGISTRATION (fixed before the run).

  A. Z = D(k)/sqrt(support(k)) with support counted EXACTLY, pooled
     over 400 even N and 1000 values of k. Report skewness, excess
     kurtosis and E|Z|/sd(Z), each against its own Gaussian standard
     error sqrt(6/n), sqrt(24/n), sqrt((1-2/pi)/n). Deviations are
     quoted as z, never as "small".
     DECISION RULE: "exactly Gaussian" survives at this precision iff
     |excess kurtosis| < 3 standard errors.

  B. CLASS STRUCTURE. The conjecture says "no class structure". Split
     by gcd(k, N) and report each class separately. A class-dependent
     variance ratio at this n would be the D(k) analogue of what
     increment 283 found for C(N).

  C. THE NORMALISER TEST, the point of the increment. Repeat (A) with
     the band-mean support in place of the per-k exact support -- the
     naive stand-in. If the exact one is Gaussian and the stand-in is
     not, the lesson of 283 generalises; if both are Gaussian, the
     original half is clean where the extension was not, and that is
     worth knowing precisely because it would be the first time this
     audit found nothing.

  D. HALF-NORMAL SCALE. E|Z| against 0.798, the number section 7
     quotes, with its standard error.

  Zero-support k are excluded and COUNTED, not silently dropped: they
  are M.3's predicted annihilations (q^2 | gcd(k,N)) and their number
  is a check on the mask, not a nuisance.
"""
import numpy as np
import math
import time


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
    return mu


def moments(z, label):
    n = len(z)
    s = float(z.std(ddof=1))
    zz = (z - z.mean()) / s
    sk = float((zz ** 3).mean())
    ku = float((zz ** 4).mean()) - 3.0
    ea = float(np.abs(zz).mean())
    se_sk = math.sqrt(6.0 / n)
    se_ku = math.sqrt(24.0 / n)
    se_ea = math.sqrt((1.0 - 2.0 / math.pi) / n)
    d_ea = ea - math.sqrt(2.0 / math.pi)
    return dict(label=label, n=n, sd=s, sk=sk, z_sk=sk / se_sk,
                ku=ku, z_ku=ku / se_ku, ea=ea, d_ea=d_ea,
                z_ea=d_ea / se_ea, se_ku=se_ku)


def show(d):
    print(f"{d['label']:>22} {d['n']:>9} {d['sk']:>8.4f} "
          f"{d['z_sk']:>7.1f} {d['ku']:>9.4f} {d['z_ku']:>8.1f} "
          f"{d['ea']:>9.5f} {d['z_ea']:>8.1f}")


def main():
    X = 4_000_000
    t0 = time.time()
    mu = sieve_mu(X)
    mu2 = (mu != 0).astype(np.int8)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(3_000_000, 4_000_001, 2500)      # 400 even N
    Ns = Ns[Ns % 2 == 0]
    ks = np.arange(500, 1500)                        # 1000 values of k
    print(f"{len(Ns)} values of N, {len(ks)} values of k, "
          f"{len(Ns)*len(ks)} pairs", flush=True)

    D = np.empty(len(Ns) * len(ks), dtype=np.float64)
    SUP = np.empty_like(D)
    GK = np.empty(len(D), dtype=np.int64)
    KK = np.empty(len(D), dtype=np.int64)
    i = 0
    t1 = time.time()
    for N in Ns:
        N = int(N)
        for k in ks:
            k = int(k)
            M = (N - 1) // k
            m = np.arange(1, M + 1, dtype=np.int64)
            w = N - k * m
            a = mu[m].astype(np.int32)
            b = mu[w].astype(np.int32)
            D[i] = float(np.dot(a, b))
            SUP[i] = float(np.dot(mu2[m].astype(np.int32),
                                  mu2[w].astype(np.int32)))
            GK[i] = math.gcd(k, N)
            KK[i] = k
            i += 1
    print(f"field built  t={time.time()-t1:.0f}s", flush=True)

    nz = SUP > 0
    print(f"\nzero-support k: {int((~nz).sum())} of {len(SUP)} "
          f"({100*float((~nz).mean()):.3f}%) -- M.3's predicted")
    print("annihilations, excluded below and counted here.")

    # M.3 says q^2 | gcd(k,N) => D(k) = 0 identically. That is an
    # exact prediction of WHICH pairs vanish, so it can be checked
    # against the observed zeros rather than merely quoted. Both
    # directions matter: a predicted zero that is nonzero refutes M.3,
    # and an unpredicted zero is a second mechanism.
    pred = np.zeros(len(SUP), dtype=bool)
    j = 0
    for N in Ns:
        N = int(N)
        for k in ks:
            g = math.gcd(int(k), N)
            q = 2
            hit = False
            gg = g
            while q * q <= gg:
                if gg % (q * q) == 0:
                    hit = True
                    break
                q += 1
            pred[j] = hit
            j += 1
    obs = ~nz
    fp = int((pred & ~obs).sum())     # predicted zero, observed nonzero
    fn = int((obs & ~pred).sum())     # observed zero, not predicted
    print(f"  M.3 check: predicted {int(pred.sum())}, observed "
          f"{int(obs.sum())}, predicted-but-nonzero {fp}, "
          f"unpredicted zeros {fn}")
    print(f"  {'M.3 EXACT' if fp == 0 and fn == 0 else 'M.3 MISMATCH'}")
    print(f"  NOTE ON THE SAMPLE: every N here is divisible by 4 "
          f"(3e6 and the step 2500 both are), so 4 | gcd(k,N) alone")
    print(f"  annihilates a quarter of the k. That is a property of "
          f"this N-grid, not of the field; it is stated because a")
    print(f"  29% zero rate would otherwise read as a general fact.")
    D, SUP, GK, KK = D[nz], SUP[nz], GK[nz], KK[nz]
    Z = D / np.sqrt(SUP)

    print(f"\n{'set':>22} {'n':>9} {'skew':>8} {'z':>7} "
          f"{'exc.kurt':>9} {'z':>8} {'E|X|/sd':>9} {'z':>8}")
    dA = moments(Z, "exact support")
    show(dA)
    print(f"    Gaussian: 0, 0, {math.sqrt(2/math.pi):.5f}")
    print(f"    variance ratio E[D^2]/support = "
          f"{float((D*D).sum()/SUP.sum()):.5f}")
    verdictA = abs(dA['z_ku']) < 3.0
    print(f"    pre-registered rule |z_kurt| < 3: "
          f"{'SURVIVES' if verdictA else 'FAILS'}")

    print("\n(B) class structure by gcd(k, N) -- the conjecture says none")
    print(f"{'set':>22} {'n':>9} {'skew':>8} {'z':>7} "
          f"{'exc.kurt':>9} {'z':>8} {'E|X|/sd':>9} {'z':>8}")
    for g in (1, 2, 3, 4, 6):
        sel = (GK == g) if g != 6 else (GK >= 5)
        if int(sel.sum()) < 5000:
            continue
        lab = f"gcd = {g}" if g != 6 else "gcd >= 5"
        show(moments(Z[sel], lab))
    vr = []
    for g in sorted(set(GK.tolist()))[:12]:
        sel = GK == g
        if int(sel.sum()) < 5000:
            continue
        vr.append((int(g), float((D[sel]**2).sum()/SUP[sel].sum()),
                   int(sel.sum())))
    print("    variance ratio E[D^2]/support by class:")
    for g, v, n in vr:
        print(f"      gcd = {g:>3}: {v:.5f}   (n = {n})")
    spread = max(v for _, v, _ in vr) - min(v for _, v, _ in vr)
    print(f"    spread across classes: {spread:.5f}")

    print("\n(C) the normaliser test: exact per-k support versus a")
    print("    band-mean stand-in -- increment 283's question, asked")
    print("    of the other half of the conjecture")
    print(f"{'set':>22} {'n':>9} {'skew':>8} {'z':>7} "
          f"{'exc.kurt':>9} {'z':>8} {'E|X|/sd':>9} {'z':>8}")
    show(dA)
    Zb = D / math.sqrt(float(SUP.mean()))
    dC = moments(Zb, "band-mean support")
    show(dC)
    print(f"    support varies: mean {float(SUP.mean()):.1f}, "
          f"sd {float(SUP.std()):.1f} "
          f"({100*float(SUP.std()/SUP.mean()):.1f}%), "
          f"range [{float(SUP.min()):.0f}, {float(SUP.max()):.0f}]")

    print("\n(D) half-normal scale: E|Z| against section 7's 0.798")
    ez = float(np.abs(Z).mean() / Z.std(ddof=1))
    print(f"    E|Z|/sd = {ez:.5f}   half-normal "
          f"{math.sqrt(2/math.pi):.5f}   "
          f"z = {dA['z_ea']:+.1f}")
    print("DONE")


if __name__ == "__main__":
    main()
