# -*- coding: utf-8 -*-
"""
The wall in the units a proof is written in (increment 320)

WHY. Increment 319 reported "the margin is a power of N and the circle
method loses the same power", from Proposition E's loss
sqrt(Q/(rho*A)) ~ N^0.498 set beside a margin trivial/truth ~ N^0.454.
Both numbers are right and the comparison is wrong: **a proof needs
bound <= target, not bound ~ truth.** Comparing an upper bound against
the measured truth answers "how lossy is the method", which is not the
question; the question is "does the bound reach the target".

Done in the right units the picture changes and gets sharper. Write
the target as |C(N)| << N (log N)^{-A}.

    target                N (log N)^{-A}
    trivial   psi(N)      N                      -> (log N)^{A}   above
    Cauchy-Schwarz        sqrt(W(N) Q(N))
                          = N sqrt(6 (log N - 1) / pi^2)
                                                 -> (log N)^{A+1/2} above
    measured truth        sqrt(rho * A(N) * W)
                          ~ sqrt(rho A N log N)  -> BELOW target by ~sqrt(N)

So the truth sits a power of N *below* the target while every available
bound sits a log power *above* it. **The entire difficulty is a log
power**, and that is the parity problem stated in its own units: no
method gives ANY saving over trivial for a Mobius-weighted sum.

AND IT UNIFIES THE TWO ROUTES, which 319 called "two different
reasons". Both must supply a saving of (log N)^A over trivial:

    circle method     supplies -(1/2) of a log power (it LOSES
                      sqrt(log N) against trivial)
    divisor switch    loses exp(c sqrt(log N)), which exceeds EVERY
                      fixed power of log (Theorem D)

Neither supplies any. Same currency, same failure, one wall.

Note what does NOT enter: rho and A(N) cancel out of the
Cauchy-Schwarz-versus-trivial comparison entirely -- the arithmetic of
the wall's own scale is irrelevant to whether the method reaches the
target. That is why sharpening the Mobius input cannot help, and it is
Proposition E's content restated so the units are visible.

PRE-REGISTRATION (fixed before the run).

  (V1) THE CAUCHY-SCHWARZ BOUND IS N sqrt(6 (log N - 1) / pi^2), with
       nothing else in it. Measured: sqrt(W(N) Q(N)) / N divided by
       sqrt(6 (log N - 1)/pi^2) must be 1 within 3% at every band, and
       must not drift. If a rho or an A(N) appeared there, the
       cancellation claimed above would be false.

  (V2) THE TRUTH IS BELOW THE TARGET BY A POWER OF N, the target's
       log power notwithstanding: rms(C) / (N (log N)^{-A}) must fall
       like N^{-1/2} for A = 1 and A = 2. RULE: the fitted exponent of
       rms(C)/N against N is -1/2 within 0.05.

  (V3) NO RULE: the table of ratios, printed in the units of the
       target, which is the deliverable.

  WHAT WOULD REFUTE. (V1) failing means the Cauchy-Schwarz bound
  carries arithmetic that could in principle be exploited, and
  Proposition E's "no improvement in Mobius technology can help" would
  need qualifying.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

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


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam, (0, nf - X - 1))),
                     nf)[: X + 1]
    W = np.cumsum(lam ** 2)
    Q = np.cumsum((mu != 0).astype(np.float64))
    Ns = np.arange(lo, X + 1, 2)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(V1)(V3) every quantity in units of N, and of the target")
    print(f"{'N':>11} {'CS/N':>9} {'sqrt(6(lnN-1)/pi2)':>19} {'ratio':>7} "
          f"{'rms C/N':>10} {'trivial/CS':>11}")
    okV1 = True
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Nb = Ns[sel]
            Nbar = float(Nb.mean())
            cs = float(np.sqrt(W[Nb] * Q[Nb]).mean())
            rms = float(np.sqrt((C[Nb] ** 2).mean()))
            # W(N) = sum Lambda(n)^2 ~ N(log N - 1), not N log N:
            # sum_{p<=N}(log p)^2 ~ int log t dt = N log N - N. A first
            # draft used the leading term alone and missed by 4%, with
            # the ratio climbing monotonically toward 1 -- the
            # signature of a missing secondary term, not of arithmetic
            # hiding inside the bound.
            pred = math.sqrt(C2PI * (math.log(Nbar) - 1.0))
            r = (cs / Nbar) / pred
            okV1 &= abs(r - 1.0) <= 0.03
            rows.append((Nbar, cs / Nbar, pred, r, rms / Nbar))
            print(f"{Nbar:>11.3e} {cs/Nbar:>9.4f} {pred:>19.4f} "
                  f"{r:>7.4f} {rms/Nbar:>10.3e} {Nbar/cs:>11.4f}")
        b = hi
    print(f"\n    (V1) CS/N equals sqrt(6(ln N-1)/pi^2) to 3%, with no "
          f"rho and no A(N) in it: {'PASS' if okV1 else 'FAIL'}")

    Nm = np.array([r[0] for r in rows])
    rn = np.array([r[4] for r in rows])
    g = float(np.polyfit(np.log(Nm), np.log(rn), 1)[0])
    okV2 = abs(g + 0.5) <= 0.05
    print(f"    (V2) rms(C)/N falls like N^(-1/2): "
          f"{'PASS' if okV2 else 'FAIL'}  (fitted {g:+.4f})")

    print(f"\n(V3) the wall in the units a proof is written in, "
          f"at N = 1e8")
    Nq = 1e8
    lg = math.log(Nq)
    for A in (1, 2):
        tgt = Nq * lg ** (-A)
        triv = Nq
        cs = Nq * math.sqrt(C2PI * (lg - 1.0))
        truth = math.sqrt(0.81 * 0.8106 * Nq * lg)
        print(f"    A = {A}:")
        print(f"      target   N(log N)^-{A}      {tgt:>12.4e}   "
              f"ratio 1")
        print(f"      trivial  psi(N)            {triv:>12.4e}   "
              f"ratio {triv/tgt:>10.3e}  = (log N)^{math.log(triv/tgt)/math.log(lg):.2f}")
        print(f"      Cauchy-Schwarz             {cs:>12.4e}   "
              f"ratio {cs/tgt:>10.3e}  = (log N)^{math.log(cs/tgt)/math.log(lg):.2f}")
        print(f"      measured truth             {truth:>12.4e}   "
              f"ratio {truth/tgt:>10.3e}  = N^{math.log(truth/tgt)/math.log(Nq):.3f}")
    print(f"\n    The truth is a POWER OF N below the target and every")
    print(f"    bound is a LOG POWER above it. The whole difficulty is")
    print(f"    a log power -- which is the parity problem in its own")
    print(f"    units: no method gives ANY saving over trivial for a")
    print(f"    Mobius-weighted sum.")
    print(f"    Both routes are then the same failure in one currency:")
    print(f"    the circle method must supply (log N)^A and instead")
    print(f"    LOSES sqrt(log N); the divisor switch loses")
    print(f"    exp(c sqrt(log N)), beyond every fixed log power.")
    if okV1 and okV2:
        v = ("increment 319's comparison was against the wrong "
             "reference -- bound versus truth, where a proof needs "
             "bound versus target. In the right units the margin to be "
             "closed is a LOG POWER, not a power of N, and the two "
             "routes fail in the same currency rather than for "
             "different reasons. rho and A(N) cancel out of it "
             "entirely, which is why sharpening the Mobius input "
             "cannot help")
    elif okV1:
        v = ("the Cauchy-Schwarz bound is arithmetic-free as claimed, "
             "but rms(C)/N does not fall like N^(-1/2), so the "
             "truth-to-target statement needs its own re-derivation")
    else:
        v = ("the Cauchy-Schwarz bound carries arithmetic after all; "
             "Proposition E's reach would need qualifying")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
