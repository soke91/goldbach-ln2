# -*- coding: utf-8 -*-
"""
Is #166's suppression real, or is it the 1/sqrt(V) rescaling?
(increment 327)

WHY. Increment 326 (#166) found that a coin carries 2.1x the real's
atomic energy share in Z = C/sqrt(V) -- mu apparently SUPPRESSING the
periodic covariance that Lambda supplies. That was the first place in
this program where mu looked unlike a coin in a direction that is not
an estimator defect.

Its own mechanism check (#167) then measured, on the UNWINDOWED and
UNRESCALED object,

    |sum_v mu(v) e(v j/q)|^2 / Q  =  0.9492
    |sum_v eps(v) e(v j/q)|^2 / Q =  0.9975

-- a factor 1.05. But the covariance of the raw C is exactly
|mu-hat(f)|^2 |Lambda-hat(f)|^2, and the coin's is
|eps-hat(f)|^2 |Lambda-hat(f)|^2, so at every frequency the ratio of
the two spectra IS that 1.05. **On the raw C there is no suppression at
all.**

Which leaves two candidates for where the 2.1x came from, and both are
mine rather than mu's:

    the 1/sqrt(V(N)) rescaling, which is an N-dependent multiplication
        and therefore a CONVOLUTION in frequency; and
    the window, N in [2e5, 8e6] rather than the whole range.

This program's most frequent fault by a wide margin is
scale-normalisation drift -- #30, #36, RV #3, #133 -- so the
normalisation is the first thing to suspect, not the last.

PRE-REGISTRATION (fixed before the run).

  Compute the atomic energy share four ways on the same N grid, for
  both the real and a coin:

      raw       C(N)                      no rescaling
      scaled    C(N)/sqrt(V(N))           the Z of #166

  (N1) THE RAW OBJECTS AGREE. On the raw C the real and coin atomic
       shares must agree within 20%, since their spectra differ by the
       measured factor 1.05 at every frequency. RULE: |ratio - 1| <=
       0.20. If they do not agree, the 1.05 measurement or this
       reasoning is wrong and nothing else reads.

  (N2) THE GAP IS THE RESCALING. If (N1) holds and the scaled objects
       still differ by about 2x, then dividing by sqrt(V) is what
       creates the gap, and #166's "mu suppresses" is an artefact of
       the normalisation. RULE: report the scaled ratio beside the raw
       one; the finding is which of the two carries the gap.

  (N3) AND IF IT IS THE RESCALING, WHY. sqrt(V(N)) is itself a
       function with strong periodic structure -- V = W * A(N) and
       A(N) = prod_{q not| N}(1 - 1/(q(q-1))) depends only on which
       small primes divide N. Multiplying by 1/sqrt(V) therefore mixes
       the atomic frequencies with each other. Report the atomic energy
       share of 1/sqrt(V(N)) itself: if it is large, the mechanism is
       named.

  WHAT WOULD REFUTE. (N1) failing means the raw objects already differ
  and #166 survives at the level of C itself. That would make it a
  finding about mu rather than about my normalisation.
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
    return sorted(out)


def main():
    X = 8_000_000
    lo = 200_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    V = np.fft.irfft(np.fft.rfft(np.pad((mu != 0).astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1))),
                     nf)[: X + 1]
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1))) * F_lam,
                     nf)[: X + 1]
    rng = np.random.default_rng(327)
    idx = np.nonzero(mu != 0)[0]
    eps = np.zeros(nf)
    eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
    Cc = np.fft.irfft(np.fft.rfft(eps) * F_lam, nf)[: X + 1]

    Ns = np.arange(lo, X + 1, 2)
    n = (len(Ns) // MOD) * MOD
    Ns = Ns[:n]
    atom = np.zeros(n // 2 + 1, dtype=bool)
    for q in divisors(QP):
        if q == 1:
            continue
        for j in range(1, q):
            if math.gcd(j, q) == 1:
                b = (j * n) // q
                if 0 < b < len(atom):
                    atom[b] = True
    print(f"n = {n}, {int(atom.sum())} atomic bins   "
          f"t={time.time()-t0:.0f}s", flush=True)

    def share(y):
        y = y - y.mean()
        P = np.abs(np.fft.rfft(y)) ** 2
        return float(P[atom].sum() / P[1:].sum())

    sr_raw = share(C[Ns])
    sc_raw = share(Cc[Ns])
    iv = 1.0 / np.sqrt(V[Ns])
    sr_sc = share(C[Ns] * iv)
    sc_sc = share(Cc[Ns] * iv)
    s_iv = share(iv)

    print(f"\n(N1)(N2) atomic energy share, four ways")
    print(f"{'object':<26} {'real':>10} {'coin':>10} {'coin/real':>11}")
    print(f"{'raw  C(N)':<26} {sr_raw:>10.5f} {sc_raw:>10.5f} "
          f"{sc_raw/sr_raw:>11.3f}")
    print(f"{'scaled  C(N)/sqrt(V(N))':<26} {sr_sc:>10.5f} "
          f"{sc_sc:>10.5f} {sc_sc/sr_sc:>11.3f}")
    okN1 = abs(sc_raw / sr_raw - 1.0) <= 0.20
    print(f"\n    (N1) the RAW objects agree within 20%: "
          f"{'PASS' if okN1 else 'FAIL'}  "
          f"(coin/real = {sc_raw/sr_raw:.3f})")
    print(f"    (N2) the scaled ratio is {sc_sc/sr_sc:.3f} against a "
          f"raw {sc_raw/sr_raw:.3f}")

    print(f"\n(N3) the multiplier itself")
    print(f"    atomic energy share of 1/sqrt(V(N)) alone: "
          f"{s_iv:.5f}")
    print(f"    (V = W*A(N), and A depends only on which small primes")
    print(f"     divide N, so 1/sqrt(V) is nearly periodic mod 30030")
    print(f"     and multiplying by it CONVOLVES the spectrum with its")
    print(f"     own atoms.)")

    # (N4) THE CONTRADICTION WITH #167, resolved. That run evaluated
    # the exponential sums at j/q and found the real and coin within
    # 1.05x. But the periodogram runs over EVEN N with step 2, so a
    # bin b/n is the N-frequency b/(2n): the atoms sit at j/(2q), not
    # j/q. #167 measured the right quantity at the wrong frequencies.
    print("\n(N4) the exponential sums at the frequencies the atoms "
          "actually sit at")
    supm = (mu != 0)
    vv = np.nonzero(supm)[0]
    mv = mu[supm].astype(np.float64)
    ev = eps[vv]
    Qn = float(len(vv))
    print(f"{'evaluated at':<16} {'mu':>12} {'coin':>12} {'coin/mu':>10}")
    for tag, den in (("j/q", 1), ("j/(2q)", 2)):
        rm, re_ = [], []
        for q in divisors(QP):
            if q == 1:
                continue
            M = den * q
            r = (vv % M).astype(np.int64)
            am = np.bincount(r, weights=mv, minlength=M)
            ae = np.bincount(r, weights=ev, minlength=M)
            for jj in range(1, q):
                if math.gcd(jj, q) != 1:
                    continue
                w = np.exp(2j * np.pi * jj * np.arange(M) / M)
                rm.append(abs(complex(np.dot(am, w))) ** 2 / Qn)
                re_.append(abs(complex(np.dot(ae, w))) ** 2 / Qn)
        rm, re_ = np.array(rm), np.array(re_)
        print(f"{tag:<16} {rm.mean():>12.4f} {re_.mean():>12.4f} "
              f"{re_.mean()/rm.mean():>10.2f}")
        if den == 2:
            gap2 = float(re_.mean() / rm.mean())
    okN4 = gap2 > 2.0
    print(f"    (N4) at j/(2q) the coin's sum exceeds mu's by more "
          f"than 2x: {'PASS' if okN4 else 'FAIL'}  ({gap2:.2f}x)")

    if okN1 and abs(sc_sc / sr_sc - 1.0) > 0.5:
        v = (f"#166 is an artefact of the normalisation. On the raw "
             f"C(N) the real and the coin carry the same atomic share "
             f"to {abs(sc_raw/sr_raw-1):.1%}, exactly as the 1.05 "
             f"exponential-sum ratio of #167 requires. The 2.1x gap "
             f"appears only after dividing by sqrt(V(N)) -- an "
             f"N-dependent multiplication, hence a convolution in "
             f"frequency, by a function whose own atomic share is "
             f"{s_iv:.3f}. mu does NOT suppress the periodic "
             f"covariance; my normalisation redistributed it")
    elif okN1:
        v = ("the raw objects agree and the scaled ones do too, so "
             "#166's gap does not reproduce here and its own "
             "measurement needs re-checking before anything is said")
    else:
        v = (f"the RAW objects already differ by "
             f"{sc_raw/sr_raw:.2f}x, so #166 survives at the level of "
             f"C(N) itself: it is a statement about mu and NOT about "
             f"the normalisation, which contributes only "
             f"{abs(sc_sc/sr_sc - sc_raw/sr_raw)/(sc_raw/sr_raw):.1%}. "
             f"And #167's contradiction is resolved: it evaluated the "
             f"exponential sums at j/q, but the periodogram runs over "
             f"EVEN N, so its atoms sit at j/(2q). At those "
             f"frequencies the coin's sum exceeds mu's by "
             f"{gap2:.1f}x -- which is the mechanism #167 proposed and "
             f"mis-measured"
             if okN4 else
             f"the RAW objects already differ by "
             f"{sc_raw/sr_raw:.2f}x, so #166 survives at the level of "
             f"C(N) itself, and #167's 1.05 was measured at j/q while "
             f"the atoms sit at j/(2q) -- but even there the sums "
             f"differ by only {gap2:.2f}x, so the mechanism is still "
             f"not identified")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
