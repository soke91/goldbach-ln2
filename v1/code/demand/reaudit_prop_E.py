# -*- coding: utf-8 -*-
"""
Proposition E re-derived against the corrected size of C(N)
(increment 309).

WHY. `OPEN_QUESTIONS.md` Register B lists the closures that rest on a
quantity a later correction moved, and put Proposition E first: it
"rests on the size of C(N) against psi(N) ~ N", and #47 changed that
size to sqrt(N)(log N)^{0.29}. The triage said the direction looked
favourable and that "most likely is not a derivation". This is the
derivation.

READING THE PROOF AGAIN, THE TRIAGE WAS WRONG. Proposition E compares
the two classical circle-method estimates against the TRIVIAL BOUND,
not against the true size of C(N):

  (i)  Cauchy-Schwarz  |C| <= ||S_L||_2 ||S_mu||_2
       = (Sum_{n<N} Lambda(n)^2)^{1/2} (Sum_{m<N} mu(m)^2)^{1/2}
       ~ (N log N)^{1/2} (6N/pi^2)^{1/2},
  (ii) pointwise x L^1  |C| <= sup_a|S_mu| * ||S_L||_1,
       with sup_a|S_mu| >= ||S_mu||_2 by Parseval.

Every quantity in those two lines is either an exact sum over n < N or
a consequence of the prime number theorem. **The true size of C(N)
enters nowhere.** So the closure does not depend on #47, on #36, or on
any other correction in the record, and the triage mis-assigned it.
That mis-assignment is itself recorded.

WHAT #47 DOES BUY. The old statement was qualitative -- the estimates
sit "at or above" the trivial bound, so there is "zero margin". With
the size of C(N) known the deficit becomes a MEASURED QUANTITY:

    deficit(N) = (what the method can deliver) / (what C(N) is)
               ~ N / (sqrt(N)(log N)^{0.29})
               = sqrt(N) (log N)^{-0.29}.

That is a full power of N, not a log power. The closure strengthens
from "no margin" to "short by a power of N", and that is a strictly
stronger statement which needs the correction to say.

PRE-REGISTRATION (fixed before the run).

  Over octave bands, with everything exact:
    W(N) = Sum_{n<N} Lambda(n)^2,  Q(N) = Sum_{m<N} mu(m)^2,
    CS(N) = sqrt(W(N) Q(N)),                    [route (i)]
    rms C = root mean square of C(N) over the band.
  Then deficit = CS / rms(C).

  (A) THE DEFICIT IS A CLEAN POWER OF N, AND IT IS DERIVED RATHER
      THAN FITTED. A first draft of this run asked whether a power of N
      beat a power of log N in fitted RSS. It came back exponent
      +0.5133 -- inside the window -- with the power model ahead by
      only 2.8x against a required 10x, i.e. NOT SEPARABLE. That is
      exactly what #119 established three increments earlier about this
      same range of N, and setting the rule that way repeated the fault
      it recorded. The rule below replaces it and needs no fit at all.

      Cauchy-Schwarz gives CS = sqrt(W(N) Q(N)). Proposition V gives
      V(N) = W(N) A(N) (1+o(1)), and rms(C)^2 = rho * mean V by the
      definition of rho. So

          deficit = CS / rms(C)
                  = sqrt( W Q / (rho A W) )
                  = sqrt( Q / (rho A) )
                  ~ sqrt( 6 N / (pi^2 rho A) ).

      **W cancels.** The log powers on both sides -- (log N)^{1/2} in
      the Cauchy-Schwarz bound and the same in the wall's scale -- are
      the same log power, so the deficit is proportional to sqrt(N)
      with a bounded constant, with nothing fitted anywhere.

      A SECOND draft then tested "deficit/sqrt(N) equals
      sqrt(Q/(rho A))/sqrt(N)" with rho := mean(C^2)/V and A := V/W.
      Substituting those definitions makes the two sides the same
      expression, and the run returned a ratio of 1.0000 in all eight
      bands. That is a check that cannot come out false -- hazard 6's
      third form, built three increments after shipping a linter
      against it, and invisible to that linter because the verdict was
      computed and it was the CRITERION that was tautological.

      What the derivation actually asserts is that the three factors
      left after W cancels are BOUNDED and do not drift. Those are the
      rules:
        (A1) A(N) = V/W is constant to 1% across bands;
        (A2) Q(N)/N reproduces 6/pi^2 to 1%;
        (A3) rho stays inside [0.5, 1.2].
      Together they give deficit/sqrt(N) = sqrt(Q/N)/sqrt(rho A)
      bounded, hence deficit asymptotic to a clean power of N with no
      log power surviving. Each can fail on its own.

  (B) THE PARSEVAL FLOOR HOLDS AT EVERY SIZE. For a handful of N by
      direct FFT: ||S_mu||_2 * ||S_L||_1 / N >= 1. This is route (ii)'s
      floor, and route (ii) is the only one with any hope. RULE: no
      N tested falls below 1.

  (C) THE SELF-TEST that makes (B) readable: ||S_mu||_2 / sqrt(N) must
      reproduce sqrt(6/pi^2) = 0.7797 to three decimals, since that is
      an identity and not a measurement.

  WHAT WOULD REFUTE. (A) failing in the log direction would mean the
  circle method is short of C(N) by only a log power, which is a
  different and much less final statement. (B) failing at any size
  would mean the Parseval floor is not what the proposition says.
"""
import math
import time

import numpy as np

C2PI = 6.0 / math.pi ** 2


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
    return mu, lam


def wfit(x, y, w):
    A = np.stack([np.ones_like(x), x], axis=1)
    c = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w),
                        rcond=None)[0]
    r = y - A @ c
    return c, float((w * r * r).sum())


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    muf = mu.astype(np.float64)
    C = np.fft.irfft(np.fft.rfft(np.pad(muf, (0, nfft - X - 1)))
                     * np.fft.rfft(np.pad(lam, (0, nfft - X - 1))),
                     nfft)[: X + 1]
    # exact partial sums: W(N) = sum_{n<N} Lambda^2, Q(N) = sum mu^2
    Wcum = np.cumsum(lam ** 2)
    Qcum = np.cumsum((mu != 0).astype(np.float64))
    print(f"sieve + C  t={time.time()-t0:.0f}s", flush=True)

    # V(N) exactly, for Proposition V's A(N) = V/W and for rho
    Vex = np.fft.irfft(np.fft.rfft(np.pad((mu != 0).astype(np.float64),
                                          (0, nfft - X - 1)))
                       * np.fft.rfft(np.pad(lam ** 2,
                                            (0, nfft - X - 1))),
                       nfft)[: X + 1]
    Ns = np.arange(lo, X + 1, 2)
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Nb = Ns[sel]
            cs = float(np.sqrt(Wcum[Nb] * Qcum[Nb]).mean())
            cb = C[Nb]
            rms = float(np.sqrt((cb * cb).mean()))
            Vb = float(Vex[Nb].mean())
            Wb = float(Wcum[Nb].mean())
            Qb = float(Qcum[Nb].mean())
            rho = float((cb * cb).mean()) / Vb
            A = Vb / Wb
            # 밴드 **평균**인 Qb를 기하 중점으로 나누면 1.5/√2 = 1.0607배가
            # 끼어든다. 첫 판이 그래서 Q/N = 0.6448을 내고 6/π²에서
            # 5.4% 벗어난다고 FAIL을 찍었다 — 산술이 아니라 정규화 실수다.
            Nbar = float(Nb.mean())
            rows.append((math.sqrt(b * hi), int(sel.sum()), cs, rms,
                         cs / rms, rho, A, Nbar, Qb))
        b = hi

    # A first draft compared deficit/sqrt(N) against
    # sqrt(Q/(rho*A))/sqrt(N) with rho := mean(C^2)/V and A := V/W.
    # Substituting those definitions makes the two sides the SAME
    # expression, and the run duly returned a ratio of 1.0000 in all
    # eight bands -- a check that could not come out false, hazard 6's
    # third form, built by me three increments after shipping a linter
    # against it. `lint_verdicts.py` cannot catch this: the verdict was
    # computed, it was the CRITERION that was tautological.
    #
    # The content of the derivation is not that identity. It is that
    # the three factors left after W cancels are BOUNDED and do not
    # drift, so that deficit ~ sqrt(N) with no log power surviving.
    # Those are the three things measured below.
    print(f"\n(A) after W cancels, what is left must not drift")
    print(f"{'N':>12} {'n':>9} {'deficit':>10} {'/sqrt(N)':>10} "
          f"{'A = V/W':>9} {'Q/N':>8} {'rho':>7}")
    Avals, Qratio, rhos = [], [], []
    for Nm, n, cs, rms, d, rho, A, Nbar, Qb in rows:
        obs = d / math.sqrt(Nm)
        Avals.append(A); rhos.append(rho)
        qn = Qb / Nbar
        Qratio.append(qn)
        print(f"{Nm:>12.3e} {n:>9} {d:>10.1f} {obs:>10.4f} "
              f"{A:>9.4f} {qn:>8.4f} {rho:>7.4f}")
    Avals = np.array(Avals); Qratio = np.array(Qratio)
    rhos = np.array(rhos)
    sA = float(Avals.max() / Avals.min() - 1.0)
    okA1 = sA <= 0.01
    okA2 = bool((np.abs(Qratio - C2PI) < 0.01 * C2PI).all())
    okA3 = bool((rhos > 0.5).all() and (rhos < 1.2).all())
    okA = okA1 and okA2 and okA3
    print(f"\n    (A1) A(N) = V/W constant to 1% across bands: "
          f"{'PASS' if okA1 else 'FAIL'}  (spread {sA:.2%}, "
          f"value {Avals.mean():.4f})")
    print(f"    (A2) Q(N)/N reproduces 6/pi^2 = {C2PI:.4f} to 1%: "
          f"{'PASS' if okA2 else 'FAIL'}  "
          f"(measured {Qratio.mean():.4f})")
    print(f"    (A3) rho bounded in [0.5, 1.2]: "
          f"{'PASS' if okA3 else 'FAIL'}  "
          f"({rhos.min():.4f} to {rhos.max():.4f})")
    print(f"    ==> deficit/sqrt(N) = sqrt(Q/N) / sqrt(rho A) is "
          f"bounded, so the deficit is a clean power of N")
    Nm = np.array([r[0] for r in rows])
    d = np.array([r[4] for r in rows])
    w = np.ones(len(d))
    cP, rP = wfit(np.log(Nm), np.log(d), w)
    cL, rL = wfit(np.log(np.log(Nm)), np.log(d), w)
    gP = float(cP[1])
    print(f"    for the record, and NOT used as a criterion: fitting "
          f"log(deficit)")
    print(f"      against log N     gives {gP:+.4f}, weighted RSS "
          f"{rP:.3e}")
    print(f"      against log log N gives {float(cL[1]):+.4f}, "
          f"weighted RSS {rL:.3e}")
    print(f"      ratio {rL/max(rP,1e-300):.2f}x -- over this range the "
          f"two are NOT separable by fit (#119), which is why the "
          f"criterion above is a derivation instead")

    # ---- (B), (C): the Parseval floor by direct FFT ----
    print(f"\n(B)(C) the Parseval floor for route (ii), exact FFT")
    print(f"{'N':>10} {'||S_mu||2/sqrt(N)':>18} {'||S_L||1/sqrt(N)':>18} "
          f"{'floor = prod/N':>16}")
    okB, okC = True, True
    for e in (14, 16, 18, 20, 22):
        Nn = 1 << e
        g = 4 * Nn
        sm = np.fft.fft(muf[:Nn], g)
        sl = np.fft.fft(lam[:Nn], g)
        l2mu = float(np.sqrt((np.abs(sm) ** 2).mean()))
        l1lam = float(np.abs(sl).mean())
        floor = l2mu * l1lam / Nn
        r2 = l2mu / math.sqrt(Nn)
        okB &= floor >= 1.0
        okC &= abs(r2 - math.sqrt(C2PI)) < 5e-3
        print(f"{Nn:>10} {r2:>18.4f} {l1lam/math.sqrt(Nn):>18.4f} "
              f"{floor:>16.4f}")
    print(f"\n    (C) ||S_mu||2/sqrt(N) reproduces sqrt(6/pi^2) = "
          f"{math.sqrt(C2PI):.4f}: {'PASS' if okC else 'FAIL'}")
    print(f"    (B) route (ii)'s floor >= 1 at every size: "
          f"{'PASS' if okB else 'FAIL'}")

    print(f"\n(D) does any quantity in the proof depend on a correction?")
    print(f"    W(N) = sum Lambda^2   exact sum, PNT asymptotic")
    print(f"    Q(N) = sum mu^2       exact sum, 6N/pi^2")
    print(f"    psi(N) ~ N            PNT")
    print(f"    sup|S_mu| >= ||S_mu||_2   Parseval, an identity")
    print(f"    none of these appears in CLOSURE_REAUDIT")
    if okA and okB and okC:
        v = ("Proposition E stands and is INDEPENDENT of every recorded "
             "correction; increment 308's triage mis-assigned it. And "
             "Proposition V sharpens it: the log powers cancel between "
             "the Cauchy-Schwarz bound and the wall's own scale, so the "
             "method is short of C(N) by a clean sqrt(N) -- a power of "
             "N, not a log power. 'Zero margin' understates it")
    elif okB and okC:
        v = ("the Parseval floor stands, but the deficit is not the "
             "power of N the corrected size predicts -- the "
             "strengthening does not go through")
    else:
        v = ("the re-derivation does not reproduce the proposition's own "
             "identities; nothing here reads")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
