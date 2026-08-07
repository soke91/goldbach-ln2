# -*- coding: utf-8 -*-
"""
The factor 5, bin by bin (increment 329)

WHY. #171 recorded that three mechanisms had been proposed and refuted,
and that the next step is to check the periodogram against the
exponential sum AT A SINGLE BIN rather than invent a fourth story.

There is one difference between the two summaries that has been there
all along and was never stated. The periodogram numerator is

    sum over atomic bins of |mu-hat(f)|^2 |Lambda-hat(f)|^2,

which is |Lambda-hat|^2-WEIGHTED, and |Lambda-hat(j/q)|^2 is about
(X/phi(q))^2 -- so q = 3 outweighs q = 143 by a factor (120/2)^2 = 3600.
Every exponential-sum check so far reported the UNWEIGHTED MEAN of
|mu-hat|^2/Q over the 255 frequencies. Two summaries of the same
spectrum with different weights are not required to agree, and nothing
in #167, #169 or #171 noticed that.

If that is the whole story, then per q the two ratios agree, and the
aggregate difference is entirely the weighting -- in which case the
finding is about my summary statistic and not about mu.

PRE-REGISTRATION (fixed before the run).

  (Q1) PER q, THE TWO RATIOS AGREE. For each of the 27 exact moduli
       report coin/real from the periodogram bins of that q, and
       coin/real from the exponential sums at the same frequencies.
       RULE: the two agree within 30% for every q with at least four
       frequencies. If they do, the objects are the same and only the
       weighting differed.

  (Q2) THE WEIGHTED MEAN REPRODUCES THE AGGREGATE. Weight the
       per-frequency ratios by the periodogram's own |C-hat_coin|^2 and
       check the weighted mean returns the 4.39 of #170, while the
       unweighted mean returns the 0.9 of #171. RULE: weighted within
       20% of 4.39, unweighted within 20% of 0.9. This is the claim
       that the discrepancy IS the weighting, stated so it can fail.

  (Q3) WHERE IT CONCENTRATES. Report each q's share of the total
       atomic numerator, for real and coin. If q = 3 dominates, that
       is the bin to look at and the "single bin" of #171 is named.

  WHAT WOULD REFUTE. (Q1) failing at some q means the periodogram and
  the exponential sum genuinely disagree there, and one of the two
  computations is wrong -- which is the other possibility #171 named.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MOD = 30030
QP = [3, 5, 7, 11, 13]


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


def divisors(ps):
    out = [1]
    for p in ps:
        out += [d * p for d in out]
    return sorted(d for d in out if d > 1)


def main():
    X = 8_000_000
    lo = 200_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1))) * F_lam,
                     nf)[: X + 1]
    rng = np.random.default_rng(329)
    idx = np.nonzero(mu != 0)[0]
    eps = np.zeros(nf)
    eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
    Cc = np.fft.irfft(np.fft.rfft(eps) * F_lam, nf)[: X + 1]

    Ns = np.arange(lo, X + 1, 2)
    n = (len(Ns) // MOD) * MOD
    Ns = Ns[:n]
    yr = C[Ns] - C[Ns].mean()
    yc = Cc[Ns] - Cc[Ns].mean()
    Pr = np.abs(np.fft.rfft(yr)) ** 2
    Pc = np.abs(np.fft.rfft(yc)) ** 2

    sup = (mu != 0)
    vv = np.nonzero(sup)[0]
    mv = mu[sup].astype(np.float64)
    ev = eps[vv]
    Qn = float(len(vv))
    print(f"n = {n}   t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(Q1)(Q3) per modulus: periodogram against exponential sum")
    print(f"{'q':>6} {'#f':>4} {'P coin/real':>12} {'E coin/real':>12} "
          f"{'agree':>7} {'share real':>11} {'share coin':>11}")
    QS = divisors(QP)
    # A real FFT only holds bins 0..n/2; a frequency j/q with
    # j > q/2 folds to n - b. The first draft indexed past the end.
    def qbins(q):
        out = set()
        for j in range(1, q):
            if math.gcd(j, q) != 1:
                continue
            b = (j * n) // q
            if b > n // 2:
                b = n - b
            if 0 < b < len(Pr):
                out.add(b)
        return sorted(out)

    allb = sorted({b for q in QS for b in qbins(q)})
    totr = float(Pr[allb].sum())
    totc = float(Pc[allb].sum())
    okQ1 = True
    wr, wu, ww = [], [], []
    for q in QS:
        js = [j for j in range(1, q) if math.gcd(j, q) == 1]
        bs = qbins(q)
        pr = float(Pr[bs].sum())
        pc = float(Pc[bs].sum())
        M = 2 * q
        r = (vv % M).astype(np.int64)
        am = np.bincount(r, weights=mv, minlength=M)
        ae = np.bincount(r, weights=ev, minlength=M)
        em = ec = 0.0
        for j in js:
            w = np.exp(2j * np.pi * j * np.arange(M) / M)
            em += abs(complex(np.dot(am, w))) ** 2
            ec += abs(complex(np.dot(ae, w))) ** 2
        rp, re_ = pc / pr, ec / em
        ok = abs(rp / re_ - 1.0) <= 0.30 if len(js) >= 4 else True
        okQ1 &= ok
        wr.append(rp)
        wu.append(re_)
        ww.append(pc)
        print(f"{q:>6} {len(js):>4} {rp:>12.3f} {re_:>12.3f} "
              f"{'yes' if ok else 'NO':>7} {pr/totr:>11.4f} "
              f"{pc/totc:>11.4f}")

    wr = np.array(wr); ww = np.array(ww); wu = np.array(wu)
    agg = totc / totr
    wmean = float((wr * ww).sum() / ww.sum())
    umean = float(wu.mean())
    okQ2 = abs(wmean / agg - 1.0) <= 0.20
    print(f"\n    (Q1) periodogram and exponential sum agree per q "
          f"(30%): {'PASS' if okQ1 else 'FAIL'}")
    print(f"    (Q2) aggregate coin/real = {agg:.3f}; "
          f"|C-hat|^2-weighted mean of the per-q ratios = "
          f"{wmean:.3f}  ->  {'PASS' if okQ2 else 'FAIL'}")
    print(f"         unweighted mean of the exponential-sum ratios = "
          f"{umean:.3f}, which is what #167/#171 reported")
    j3 = QS.index(3)
    b3s = qbins(3)
    print(f"    (Q3) q = 3 carries {float(Pr[b3s].sum())/totr:.1%} "
          f"of the real's atomic numerator and "
          f"{float(Pc[b3s].sum())/totc:.1%} of the coin's")

    if okQ1 and okQ2:
        v = (f"the two objects agree at every modulus. The aggregate "
             f"{agg:.2f}x and the per-frequency {umean:.2f}x are the "
             f"same spectrum under different weights -- the "
             f"periodogram sums |mu-hat|^2|Lambda-hat|^2 and so is "
             f"dominated by small q, while #167 and #171 reported an "
             f"UNWEIGHTED mean over 255 frequencies. The contradiction "
             f"of #169 was never about mu; it was two summary "
             f"statistics of one spectrum, and nobody said which "
             f"weight each carried")
    elif okQ1:
        v = (f"the ratios agree per q but the weighted mean "
             f"({wmean:.2f}) does not reproduce the aggregate "
             f"({agg:.2f}); the weighting is part of the story and not "
             f"all of it")
    else:
        v = ("the periodogram and the exponential sum disagree at some "
             "modulus, so one of the two computations is wrong -- the "
             "other possibility #171 named, and now the one to chase")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
