# -*- coding: utf-8 -*-
"""
Fast exact checks of every proposition proved in this campaign
(increment 275).

verify_all.py is the heavy CI stamp for the measurement corpus and runs
at X = 10^8. The propositions added since -- M.1 to M.3 for the
location mask and P.1 to P.7 for transform P -- are exact statements,
so they can be checked in seconds and should be, because the proved
part of a program is the part that has no excuse for being unverifiable.

Every check below is an identity or a counting statement, and every one
is written so that it CAN come out false: nothing is rearranged from
the quantity it is meant to confirm (the fault of increment 272), each
compares an independently computed left and right side, and the stated
exceptions of M.1 and P.5 are classified rather than ignored (the fault
of increment 271).

Usage:  python code/verify_propositions.py [N]
"""
import numpy as np
import math
import sys


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    om = np.zeros(X + 1, dtype=np.int8)
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
        om[i] = om[j] + (0 if j % p == 0 else 1)
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p); lg = math.log(int(p))
        while q <= X:
            lam[q] = lg; q *= int(p)
    return mu, om, lam, primes, spf


def radset(N, spf):
    s = set()
    n = N
    while n > 1:
        p = int(spf[n]); s.add(p)
        while n % p == 0:
            n //= p
    return s


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 210210
    if N % 2:
        N -= 1
    X = N + 10
    mu, om, lam, primes, spf = sieve(X)
    R = radset(N, spf)
    ps = primes[primes < N]
    lp = np.log(ps.astype(np.float64))
    v = np.arange(1, N, dtype=np.int64)
    mv = mu[1:N].astype(np.float64)
    lamr = lam[N - v]
    logv = np.log(np.maximum(v, 2).astype(np.float64))
    t = mv * lamr / logv; t[0] = 0.0          # summand of P.1
    a = np.abs(mv) * lamr / logv; a[0] = 0.0
    print(f"N = {N},  rad(N) = {sorted(R)}\n")
    res = []

    # --- M.1: Lambda(N-v) != 0 and q | N force q not| v, except when
    #          N - v is a power of a prime dividing N.
    live = np.nonzero((mv != 0) & (lamr != 0.0))[0] + 1
    bad = exc = 0
    for q in R:
        hit = live[live % q == 0]
        for vv in hit:
            w = N - int(vv)
            qq = int(spf[w]) if w > 1 else 0
            if qq in R and lam[w] != 0.0:
                exc += 1
            else:
                bad += 1
    res.append(("M.1 forcing: v coprime to rad(N)", bad == 0,
                f"{len(live)} live v, {bad} violations, "
                f"{exc} stated exceptions"))

    # --- M.2: C(N) = (Sum_j T_j) * R_A(N), by grouping on omega
    JM = 12
    jj = om[1:N]
    aw = np.abs(mv) * lamr; aw[0] = 0.0
    tw = mv * lamr; tw[0] = 0.0
    T = np.array([float(aw[jj == j].sum()) for j in range(JM + 1)])
    alt = float(sum((-1) ** j * T[j] for j in range(1, JM + 1)))
    Cdirect = float(tw.sum())
    res.append(("M.2 identity C = (Sum T_j) R_A",
                abs(alt - Cdirect) < 1e-6 * max(1.0, abs(Cdirect)),
                f"grouped {alt:+.6f} vs direct {Cdirect:+.6f}"))

    # --- M.3: q^2 | gcd(k,N)  =>  D(k) = 0 identically
    SQ = int(N ** 0.5)
    m3bad = 0; m3n = 0
    for k in range(2, min(4000, N // (SQ + 2))):
        g = 0
        for q in R:
            if k % (q * q) == 0 and N % (q * q) == 0:
                g = q; break
        if not g:
            continue
        m3n += 1
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        if len(ms) and float(np.dot(mu[ms].astype(np.float64),
                                    mu[N - k * ms].astype(np.float64))):
            m3bad += 1
    res.append(("M.3 support mask: q^2|gcd(k,N) => D(k)=0",
                m3bad == 0, f"{m3n} such k, {m3bad} nonzero"))

    # --- P.1: C(N) - Lambda(N-1) = Sum_p log p D_p
    Dp = np.array([float(t[p - 1::p].sum()) for p in ps])
    lhs = float(np.dot(lp, Dp))
    rhs = float(tw[1:].sum())                 # Sum_{v>=2} mu Lam
    res.append(("P.1 split identity",
                abs(lhs - rhs) < 1e-6 * max(1.0, abs(rhs)),
                f"{lhs:+.6f} vs {rhs:+.6f}"))

    # --- P.2: Sum_p log p M_p = Sum_{v>=2} mu^2 Lam
    Mp = np.array([float(a[p - 1::p].sum()) for p in ps])
    l2 = float(np.dot(lp, Mp)); r2 = float(aw[1:].sum())
    res.append(("P.2 losslessness", abs(l2 - r2) < 1e-6 * r2,
                f"{l2:.6f} vs {r2:.6f}"))

    # --- P.3: Sum_w mu^2(w) G_w = the trivial bound, and
    #          -Sum_w mu(w) G_w = C(N) - Lambda(N-1).
    # The first version of this check returned a hardcoded True, i.e.
    # it could not come out false -- the fault named at increment 272,
    # shipped inside a verifier. Both identities are now summed over
    # every w and compared against independently computed right-hand
    # sides.
    Gabs = 0.0
    Gsig = 0.0
    for w in range(1, N // 2 + 1):
        if mu[w] == 0:
            continue
        hi = (N - 1) // w
        j = int(np.searchsorted(ps, hi, side='right'))
        if j == 0:
            continue
        pw = ps[:j] * w
        keep = (w % ps[:j] != 0)
        if not keep.any():
            continue
        terms = (lp[:j][keep] * lamr[pw[keep] - 1]
                 / np.log(pw[keep].astype(np.float64)))
        s = float(terms.sum())
        Gabs += s
        Gsig += float(mu[w]) * s
    res.append(("P.3a  Sum_w mu^2(w) G_w = trivial bound",
                abs(Gabs - r2) < 1e-6 * r2,
                f"{Gabs:.6f} vs {r2:.6f}"))
    res.append(("P.3b  -Sum_w mu(w) G_w = C(N) - Lambda(N-1)",
                abs(-Gsig - rhs) < 1e-6 * max(1.0, abs(rhs)),
                f"{-Gsig:+.6f} vs {rhs:+.6f}"))

    # --- P.5 + P.7: terms of D_p are v = mp with (m, rad N) = 1, and
    #     for (N-1)/p < r^2 their signs are -1 at m=1 and +1 at m prime
    r = min(q for q in primes if N % int(q) != 0)
    p5bad = p7bad = p5exc = nchk = 0
    for p in ps:
        p = int(p)
        M = (N - 1) // p
        if M < 1 or M >= r * r:
            continue
        for m in range(1, M + 1):
            vv = m * p
            if mu[vv] == 0 or lamr[vv - 1] == 0.0:
                continue
            nchk += 1
            w = N - vv
            qq = int(spf[w]) if w > 1 else 0
            isexc = qq in R and lam[w] != 0.0
            if math.gcd(m, N) != 1:
                if isexc:
                    p5exc += 1
                else:
                    p5bad += 1
                continue
            term = float(mu[vv]) * lamr[vv - 1] / math.log(vv)
            want = -1.0 if m == 1 else 1.0
            if term * want < 0 and not isexc:
                p7bad += 1
    res.append(("P.5 effective terms: (m, rad N) = 1",
                p5bad == 0,
                f"{nchk} terms, {p5bad} violations, {p5exc} exceptions"))
    res.append(("P.7 forced signs: -1 at m=1, +1 at m prime",
                p7bad == 0, f"{p7bad} violations"))

    print(f"{'proposition':>44} {'verdict':>8}   evidence")
    allok = True
    for name, ok, ev in res:
        allok &= bool(ok)
        print(f"{name:>44} {'PASS' if ok else 'FAIL':>8}   {ev}")
    print(f"\n{'ALL PASS' if allok else 'SOMETHING FAILED'}")

    # --- sensitivity: a check that cannot come out false is not a
    # check. Four times in this run a verification was shipped or
    # nearly shipped that could not fail (increments 259, 272, 274,
    # 275). Asserting that these can fail is worth nothing; showing it
    # is worth something. Each identity is re-evaluated with one side
    # perturbed by 1 part in 1000, and must flip to FAIL.
    print("\nsensitivity check: perturb one side by 1e-3, expect FAIL")
    pairs = [("M.2", alt, Cdirect), ("P.1", lhs, rhs),
             ("P.2", l2, r2), ("P.3a", Gabs, r2), ("P.3b", -Gsig, rhs)]
    sens_ok = True
    for nm, L, Rr in pairs:
        tol = 1e-6 * max(1.0, abs(Rr))
        flips = abs(L * (1 + 1e-3) - Rr) >= tol
        sens_ok &= flips
        print(f"  {nm:>5}  perturbed |lhs-rhs| = "
              f"{abs(L*(1+1e-3)-Rr):.6g}  tol = {tol:.3g}   "
              f"{'flips to FAIL (good)' if flips else 'STILL PASSES (BAD)'}")
    print(f"  counting checks (M.1, M.3, P.5, P.7) test "
          f"'violations == 0', so a single injected violation flips "
          f"them by construction")
    print(f"\n{'SENSITIVITY OK' if sens_ok else 'A CHECK CANNOT FAIL'}")

    # 증분 307. 이 파일은 CLOSURE_REAUDIT #61이 "실패할 수 없는 검사"의
    # **답**으로 내세우는 파일이고 STATUS가 재현 명령으로 인용하는데,
    # 정작 자기는 'SOMETHING FAILED'와 'A CHECK CANNOT FAIL'을 찍고도
    # **종료코드 0으로 끝났다**. 민감도 블록은 있고 실패 경로가 없었다 —
    # 위험 6번 셋째 형태가 그 형태의 수리로 지목된 파일 안에 있었다.
    # `code/lint_gates.py`가 기계로 잡는다.
    if not (allok and sens_ok):
        print("DONE (failed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
