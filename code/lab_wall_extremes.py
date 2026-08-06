# -*- coding: utf-8 -*-
"""
The wall's TAIL: max |C(N)| and the margin the chain actually has
(increment 290).

WHY THE TAIL AND NOT THE BULK. The chain needs C(N) = o(N) for EVERY
large N, not for typical N. Increment 283 settled the bulk (Gaussian
under the exact scale V), 288 the variance ratio rho, 289 what rho - 1
is. None of them says anything about the largest C(N), which is the
only quantity the requirement actually constrains. And the extreme is
where non-Gaussianity shows up first: a bulk can look Gaussian to five
decimals while the tail is heavy.

This could not be asked properly before increment 287. The extreme of
C/sqrt(kappa*S*N*log N) is the extreme of a badly normalised variable,
and increment 283 measured that normaliser producing excess kurtosis
+0.17 at z = 98 -- a fake heavy tail. With V(N) exact and the location
mask removed, the question is well posed for the first time.

PRE-REGISTRATION (fixed before the run).

  Z(N) = (C(N) - m(cell)) / sqrt(V(N)), cells keyed by which of
  {3,...,23} divide N, per-cell means removed within each band.

  (A) EXTREME VALUE LAW. For n iid standard Gaussians the maximum
      concentrates at
          a_n = sqrt(2 ln n) - (ln ln n + ln 4pi) / (2 sqrt(2 ln n)),
      with Gumbel fluctuations of scale 1/sqrt(2 ln n) (about 0.18 at
      these n). DECISION RULE: the Gaussian tail law survives if
      max|Z| lies within 3 Gumbel scales of a_n in every band. A
      systematic excess is a heavy tail; a systematic deficit is
      sub-Gaussian and would be just as interesting.

  (B) TAIL COUNTS. #{|Z| > t} against the Gaussian expectation
      2n(1-Phi(t)) for t = 3, 4, 5, reported as a ratio. This is a
      finer instrument than the single maximum and it can disagree
      with (A).

  (C) WHERE THE EXTREME SITS. Is the largest |Z| attained at an N with
      many small prime factors? If so the mask removal is leaking and
      (A) is measuring the mask, not the tail. The radical of the
      argmax N is printed for exactly this reason.

  (D) THE MARGIN, stated in the chain's own terms. The requirement is
      |C(N)| = o(N). If max|Z| ~ sqrt(2 log X), then
          max |C| ~ sqrt(2 log X) * sqrt(A(N) N log N),
      so the margin over the requirement is ~ sqrt(N)/polylog. That
      number is what the program should be quoting, and it has been
      quoting a typical-N figure instead.

WHAT THIS CANNOT DO. Reachable X is 1.6*10^7 and the requirement is
asymptotic; an extreme-value law verified here constrains nothing at
10^480. It is reported as a measurement of the observable range, which
is what increment 279's correction #66 said every statement in this
program should be read as.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]


def sieve(X):
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
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam, primes, spf


def conv(X, a, b):
    n = 1
    while n < 2 * (X + 1):
        n *= 2
    A = np.zeros(n); A[: X + 1] = a
    B = np.zeros(n); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n)[: X + 1]


def phi_tail(t):
    return 0.5 * math.erfc(t / math.sqrt(2.0))

EULER = 0.5772156649015329


def a_n(n):
    L = math.log(n)
    s = math.sqrt(2.0 * L)
    return s - (math.log(L) + math.log(4.0 * math.pi)) / (2.0 * s)


def main():
    X = 16_000_000
    t0 = time.time()
    mu, lam, primes, spf = sieve(X)
    C = conv(X, mu, lam)
    V = conv(X, (mu != 0).astype(np.float64), lam ** 2)
    print(f"sieve + 2 convolutions  t={time.time()-t0:.0f}s", flush=True)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    print("\n(A) extreme value law, per octave band")
    print(f"{'band':>21} {'n':>9} {'max|Z|':>8} {'a_n':>7} "
          f"{'(max-a_n)/scale':>16} {'argmax rad(N)':>22}")
    okA = True
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        n = int(sel.sum())
        if n <= 1000:
            b = hi
            continue
        NN = Ns[sel]
        c = C[NN]
        uniq, inv = np.unique(key[sel], return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
        tot = np.bincount(inv, weights=c, minlength=len(uniq))
        z = (c - (tot / cnt)[inv]) / np.sqrt(V[NN])
        z = z / z.std(ddof=1)          # unit variance, so (A) tests SHAPE
        i = int(np.argmax(np.abs(z)))
        an = a_n(n)
        scale = 1.0 / math.sqrt(2.0 * math.log(n))
        # a_n is the LOCATION of the Gumbel limit, not its mean:
        # E[max] = a_n + gamma*scale with gamma = 0.5772. Centering
        # on a_n alone makes a normal result look like a systematic
        # +0.58 excess. Center on E[max].
        dev = (float(np.abs(z).max()) - (an + EULER * scale)) / scale
        okA &= abs(dev) <= 3.0
        w = int(NN[i])
        rad = []
        t = w
        while t > 1:
            q = int(spf[t]); rad.append(q)
            while t % q == 0:
                t //= q
        rows.append((n, float(np.abs(z).max()), an, dev, z))
        print(f"{b:>9}-{hi:>11} {n:>9} {rows[-1][1]:>8.3f} {an:>7.3f} "
              f"{dev:>+16.2f} {'.'.join(map(str, rad[:6])):>22}")
        b = hi
    devs = [r[3] for r in rows]
    mdev = sum(devs) / len(devs)
    se = 1.2825 / math.sqrt(len(devs))   # Gumbel sd = pi/sqrt(6)
    print(f"    Gumbel scale ~ {scale:.3f}; deviations are from E[max] = ")
    print(f"    a_n + {EULER:.4f}*scale, so a centred result has mean 0")
    print(f"    mean deviation {mdev:+.2f} +/- {se:.2f} (Gumbel sd/sqrt(8))")
    print(f"    pre-registered |dev| <= 3 in every band  ->  {'PASS' if okA else 'FAIL'}")

    print("\n(B) tail counts against the Gaussian expectation")
    print(f"{'band':>21} " + " ".join(f"{'t='+str(t):>14}" for t in (3, 4, 5)))
    for (n, mx, an, dev, z), b0 in zip(rows, range(len(rows))):
        cells = []
        for t in (3, 4, 5):
            obs = int((np.abs(z) > t).sum())
            exp = 2.0 * n * phi_tail(t)
            cells.append(f"{obs:>6}/{exp:>7.1f}")
        print(f"{'band ' + str(b0 + 1):>21} " + " ".join(f"{c:>14}" for c in cells))
    print("    observed / expected. Per band these are noisy; the")
    print("    aggregate is the statistic that discriminates:")
    for tt in (3, 4, 5):
        o = sum(int((np.abs(r[4]) > tt).sum()) for r in rows)
        e = sum(2.0 * r[0] * phi_tail(tt) for r in rows)
        print(f"      t = {tt}: observed {o:>6}, expected {e:>8.1f}, ratio {o/e:>6.3f}")

    print("\n(D) the margin, in the chain's own terms")
    for k in (7, 12, 20, 50):
        Xk = 10.0 ** k
        mx = a_n(Xk / 2.0)
        print(f"    X = 1e{k:<3}: max|Z| ~ {mx:5.2f}, so "
              f"max|C| ~ {mx:.2f}*sqrt(A*N*log N), and the margin over")
        print(f"{'':>15} the requirement |C| = o(N) is "
              f"~ sqrt(N)/({mx:.1f}*sqrt(log N)) = "
              f"1e{0.5*k - math.log10(mx*math.sqrt(k*math.log(10))):.1f}")
    print("    The program has been quoting a TYPICAL-N margin. This is")
    print("    the margin at the extreme, which is the one the")
    print("    requirement actually constrains.")
    print("\n    Scope: X <= 1.6e7 here. An extreme-value law verified in")
    print("    that range constrains nothing at the sizes where the")
    print("    no-go theorems bite (correction #66).")
    print("DONE")


if __name__ == "__main__":
    main()
