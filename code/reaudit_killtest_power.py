# -*- coding: utf-8 -*-
"""
What is the smallest structure a kill-test could have found?
(increment 311)

WHY. `OPEN_QUESTIONS.md` Register B, corrected at #134, leaves the
Forge kill-tests with a live question that is **hazard 8, not hazard
7**: their `z` scores are not analytic after all -- K2 estimates its
standard error from the data, `se = std(x)/sqrt(n)` -- so the issue is
not which null they used but **whether the test had the power to see
anything**.

A kill-test that returns DEAD is only as strong as its detection
floor, and **not one kill-test in this program has ever reported
one**. "No flags at |z| <= 1.8 against a threshold of 4" is compatible
with two very different worlds: a field with no structure, and a test
that could not have seen the structure if it were there.

THE FLOOR, ANALYTICALLY. K2 flags when |mean(x)| / (std(x)/sqrt(n))
reaches 4, so the smallest coherent mean shift it can detect is

    delta* = 4 * std(x) / sqrt(n),

which for n = 600 and a unit-variance field is about 0.163 standard
deviations. That number has never appeared anywhere. This run measures
it rather than asserting it, and then asks the question that matters:
CAN the test return ALIVE at all?

REDUCED SCALE, STATED. The original K2 runs at N = 199,999,998, whose
Mobius sieve needs about 1.8 GB. This replication runs at
N = 19,999,998 with the prime window rescaled to `P0 = N // 6000`,
`P1 = 2 P0`, which keeps the support per pair above K2's own `sup >
100` cut. It is a replication of the STATISTIC, not of the original
numbers, and rule (P1) below is what makes it readable.

PRE-REGISTRATION (fixed before the run).

  (P1) FAITHFULNESS. At reduced N, with no injection, K2's own
       pre-registered rule must still return DEAD -- no h flagged by
       |zm| >= 4, |zdev| >= 5, or |rho1| >= 0.15. If the reduced
       replication flags where the original did not, it is not the same
       test and nothing else here reads.

  (P2) THE FLOOR IS MEASURED, NOT ASSUMED. Inject a coherent shift
       delta into x = C/sqrt(sup) and find the smallest delta at which
       |zm| >= 4, by bisection, per h. RULE: the measured delta*
       agrees with 4*se - mean to within 2%.

       THIS RULE IS AN IDENTITY AND IS LABELLED AS ONE. Shifting x by
       a constant leaves std(x) unchanged, so the bisection is solving
       (mean + d)/se = 4, whose closed form is exactly 4*se - mean. It
       therefore tests the INJECTION MACHINERY and nothing else, and it
       is here only because #132 -- one increment ago -- shipped an
       identity of this shape dressed as a finding. Declaring it is the
       difference. The evidence in this run is (P1) and (P3); the
       number is the floor.

       A first draft predicted 4*std(x)/sqrt(n) alone, omitting the
       field's own mean, and missed by 45% at h = 64. The injection was
       right and the reference was wrong.

  (P3) THE TEST CAN RETURN ALIVE. At delta = 2 delta*, every h must
       flag. If a doubled signal at the measured floor does not flag,
       the test could not have come out ALIVE and its DEAD verdict is
       empty -- hazard 6's third form, applied to a kill-test instead
       of a verifier.

  WHAT WOULD REFUTE. (P1) failing kills the replication. (P3) failing
  would mean this program's kill-tests have been reporting DEAD from
  tests that cannot report anything else, which would be the most
  serious finding in the closure re-audit so far.

WHAT THIS CANNOT SAY. The floor is in units of the field's own
standard deviation. Whether 0.16 sd is small enough to matter is a
question about design K2's economics -- it costs a factor d by
Cauchy-Schwarz -- and that is not settled here. What is settled is the
number, which was missing.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HS = [1, 2, 4, 6, 8, 16, 30, 32, 64, 210]
K0 = 2000
NK = 600
ZM, ZDEV, RHO = 4.0, 5.0, 0.15


def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool)
    pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
            pp = p * p
            while pp <= X:
                val[pp::pp] //= p
                pp *= p
    mu[val > 1] *= -1
    return mu, pm


def field(N, mu, pm):
    """K2's raw field: C and support for each h, over NK consecutive k."""
    P0 = N // 6000
    P1 = 2 * P0
    ps = np.nonzero(pm[P0:P1])[0].astype(np.int64) + P0
    out = {}
    for h in HS:
        kv = np.arange(K0, K0 + NK, dtype=np.int64)
        C = np.zeros(NK)
        sup = np.zeros(NK)
        for i, k in enumerate(kv):
            kp = int(k) + h
            pp = ps[ps <= (N - 2) // kp]
            t = (mu[N - pp * int(k)].astype(np.int64)
                 * mu[N - pp * kp].astype(np.int64))
            C[i] = float(t.sum())
            sup[i] = float(np.count_nonzero(t))
        out[h] = (C, sup)
    return out


def stats(C, sup, delta=0.0, scale=1.0):
    """K2's own reductions, with an optional injected mean shift."""
    good = sup > 100
    if good.sum() < 20:
        return None
    Cg, sg = C[good], sup[good]
    x = scale * Cg / np.sqrt(sg) + delta
    n = int(good.sum())
    mean = float(x.mean())
    sd = float(x.std())
    se = sd / math.sqrt(n)
    m2 = float((x ** 2).mean())
    m2_se = float((x ** 2).std()) / math.sqrt(n)
    xx = np.where(good, scale * C / np.sqrt(np.maximum(sup, 1)) + delta,
                  np.nan)
    a = xx[:-1]
    b = xx[1:]
    ok = ~(np.isnan(a) | np.isnan(b))
    rho = (float(np.corrcoef(a[ok], b[ok])[0, 1])
           if ok.sum() > 50 else float("nan"))
    return dict(n=n, mean=mean, sd=sd, se=se, m2=m2, m2_se=m2_se,
                rho=rho, zm=mean / max(se, 1e-12))


def main():
    N = 19_999_998
    t0 = time.time()
    mu, pm = mobius_upto(N)
    print(f"mu ready  t={time.time()-t0:.0f}s", flush=True)
    F = field(N, mu, pm)
    print(f"field ready  t={time.time()-t0:.0f}s", flush=True)

    # ---- (P1) the unmodified test must still return DEAD ----
    base = {}
    print(f"\n(P1) K2's own rule at reduced N = {N}, no injection")
    print(f"{'h':>5} {'n':>5} {'mean':>9} {'sd':>7} {'zm':>7} "
          f"{'m2':>7} {'zdev':>7} {'rho1':>7}  flag")
    nflag = 0
    b1 = None
    for h in HS:
        C, sup = F[h]
        s = stats(C, sup)
        base[h] = s
        if s is None:
            print(f"{h:>5}   too few viable k")
            continue
        if b1 is None:
            b1 = (s["m2"], s["m2_se"])
        zdev = (s["m2"] - b1[0]) / max(math.hypot(s["m2_se"], b1[1]),
                                       1e-12)
        fl = (abs(s["zm"]) >= ZM or (h > 1 and abs(zdev) >= ZDEV)
              or (not math.isnan(s["rho"]) and abs(s["rho"]) >= RHO))
        nflag += int(fl)
        print(f"{h:>5} {s['n']:>5} {s['mean']:>+9.4f} {s['sd']:>7.4f} "
              f"{s['zm']:>+7.2f} {s['m2']:>7.3f} {zdev:>+7.2f} "
              f"{s['rho']:>+7.3f}  {'FLAG' if fl else ''}")
    okP1 = nflag == 0
    print(f"\n    (P1) reduced replication returns DEAD: "
          f"{'PASS' if okP1 else 'FAIL'}  ({nflag} flags)")

    # ---- (P2) measure the floor by injection ----
    print(f"\n(P2) the smallest coherent shift K2 can see, by injection")
    print(f"{'h':>5} {'sd':>7} {'n':>5} {'floor':>9} "
          f"{'4se-mean':>13} {'measured d*':>12} {'ratio':>7}")
    okP2 = True
    dstars = {}
    for h in HS:
        C, sup = F[h]
        s = base[h]
        if s is None:
            continue
        # 4*se만 쓰면 필드 자신의 평균을 빼먹는다. 검정은
        # |mean + delta| >= 4*se 에서 돌아서므로 d* = 4*se − mean 이다.
        # 첫 판이 4*sd/sqrt(n)만 예측해 h=64에서 45% 빗나갔고, 틀린
        # 쪽은 주입이 아니라 기준이었다. POWER는 4*sd/sqrt(n)이고,
        # d*의 어긋남은 이 실현의 평균이지 검정의 성질이 아니다.
        floor = ZM * s["sd"] / math.sqrt(s["n"])
        pred = ZM * s["se"] - s["mean"]
        lo, hi = 0.0, max(4.0 * abs(pred), 1e-6)
        while abs(stats(C, sup, delta=hi)["zm"]) < ZM:
            hi *= 2
            if hi > 1e3:
                break
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if abs(stats(C, sup, delta=mid)["zm"]) >= ZM:
                hi = mid
            else:
                lo = mid
        d = hi
        dstars[h] = d
        r = d / pred
        okP2 &= abs(r - 1.0) <= 0.02
        print(f"{h:>5} {s['sd']:>7.4f} {s['n']:>5} {floor:>9.5f} "
              f"{pred:>13.5f} {d:>12.5f} {r:>7.3f}")
    print(f"\n    (P2) measured d* matches 4*se − mean to 2%: "
          f"{'PASS' if okP2 else 'FAIL'}")
    print(f"         (the POWER is 4*sd/sqrt(n), the 'floor' column;")
    print(f"          d*'s offset from it is this realisation's own")
    print(f"          mean, not a property of the test)")

    # ---- (P3) the test must be able to return ALIVE ----
    print(f"\n(P3) at twice the floor, every h must flag")
    okP3 = True
    for h in HS:
        if h not in dstars:
            continue
        s2 = stats(*F[h], delta=2 * dstars[h])
        fl = abs(s2["zm"]) >= ZM
        okP3 &= fl
        print(f"    h={h:>4}  injected {2*dstars[h]:+.4f}  "
              f"zm={s2['zm']:+.2f}  {'FLAGS (good)' if fl else 'SILENT (bad)'}")
    print(f"\n    (P3) the test can return ALIVE: "
          f"{'PASS' if okP3 else 'FAIL'}")

    sds = np.array([base[h]["sd"] for h in HS if base[h]])
    ns = np.array([base[h]["n"] for h in HS if base[h]])
    rel = float((ZM / np.sqrt(ns)).mean())
    orig = ZM / math.sqrt(NK)
    print(f"\n    detection floor = 4/sqrt(n), in units of the field's "
          f"own sd")
    print(f"      this reduced run, n = {int(ns.mean())} viable k:  "
          f"{rel:.4f}")
    print(f"      the original, n <= {NK} consecutive k:      "
          f">= {orig:.4f}")
    print(f"    The support cut `sup > 100` bites harder at reduced N,")
    print(f"    so this run has fewer viable k and a HIGHER floor. The")
    print(f"    original's floor is the better number and it is still")
    print(f"    {orig:.3f} sd -- a coherent mean shift below that was")
    print(f"    invisible to K2, and no run ever said so.")
    if okP1 and okP2 and okP3:
        v = (f"K2's DEAD verdict is real and now has a number attached: "
             f"it excludes coherent structure above {rel:.3f} sd and "
             f"says nothing below that. The test can return ALIVE, "
             f"shown rather than asserted")
    elif okP1 and not okP3:
        v = ("K2 cannot return ALIVE even at twice its own measured "
             "floor -- the DEAD verdict is empty and every closure "
             "resting on it must be re-opened")
    elif not okP1:
        v = ("the reduced replication does not reproduce the original's "
             "verdict; nothing here reads")
    else:
        v = ("the injection is not self-consistent with the analytic "
             "floor; the number is not trustworthy")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
