"""H1'' 파괴시험: 분절 체(segmented sieve) r(n) 챔피언 사냥 — 10^10급.

메모리에 전체 체를 못 올리므로 구간 [lo, hi)마다:
  1) 기저 소수(≤ √LIMIT)로 구간 체질
  2) 구간 내 짝수들을 p 오름차순 벡터 마스크로 해소 (r(n) 극값만 수집)
H1'' 포락선 r < 0.28·ln³n 검정 + H2(3∤n) 검정.
"""

import sys
import time
from math import log

import numpy as np

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000_000
SEG = 40_000_000          # 구간 크기 (짝수 2천만 개)
PMARGIN = 8000            # r(n) 후보 소수 상한 (0.28·ln³(10^10) ≈ 3420 의 여유)

t0 = time.time()
base_lim = int(LIMIT ** 0.5) + 1
base = np.ones(base_lim + 1, dtype=bool)
base[:2] = False
for p in range(2, int(base_lim ** 0.5) + 1):
    if base[p]:
        base[p * p :: p] = False
base_primes = np.flatnonzero(base)
small = base_primes[base_primes < PMARGIN].astype(np.int64)
print(f"기저 체 완료 (소수 {len(base_primes):,}개, {time.time()-t0:.0f}s)", flush=True)

records, best = [], 0
seg_starts = list(range(0, LIMIT, SEG))
for si, s0 in enumerate(seg_starts):
    lo = max(4, s0)
    hi = min(s0 + SEG, LIMIT)
    left = max(2, lo - PMARGIN)          # n-p 조회용 여유
    size = hi - left
    seg = np.ones(size, dtype=bool)      # seg[i] = (left+i) 소수 여부
    if left <= 2:
        seg[: 3 - left] = False          # 0,1,2 처리 (2는 소수)
        if left <= 2 <= hi - 1:
            seg[2 - left] = True
    for p in base_primes:
        p = int(p)
        st = max(p * p, ((left + p - 1) // p) * p)
        if st >= hi:
            break
        seg[st - left :: p] = False
    nvals = np.arange(lo + (lo % 2), hi, 2, dtype=np.int64)
    work = nvals
    hi_n, hi_r = [], []
    best0 = best
    for p in small:
        if len(work) == 0:
            break
        hit = seg[work - p - left]
        if p > best0 and hit.any():
            hi_n.append(work[hit]); hi_r.append(np.full(int(hit.sum()), p))
        work = work[~hit]
    assert len(work) == 0, f"미해결 {len(work)}개 (구간 {lo:,}) — PMARGIN 확대 필요"
    if hi_n:
        cn, cr = np.concatenate(hi_n), np.concatenate(hi_r)
        o = np.argsort(cn)
        for n, r in zip(cn[o], cr[o]):
            if r > best:
                best = int(r)
                records.append((int(n), best))
                v = best / log(n) ** 3
                print(f"  ★ 챔피언 n={n:,} r={best} r/ln³n={v:.4f} mod3={n%3}", flush=True)
    if si % 25 == 0:
        print(f"..{hi:,} ({time.time()-t0:.0f}s) best={best}", flush=True)

print(f"\n총 소요 {time.time()-t0:.0f}s | 범위 ≤ {LIMIT:,}")
worst = max(r / log(n) ** 3 for n, r in records if n > 100)
mods = [n % 3 for n, _ in records if n > 30]
print(f"H1'' max r/ln³n = {worst:.4f} → {'생존 (0.28 이내)' if worst < 0.28 else '사망! (0.28 돌파)'}")
print(f"H2   3의 배수 챔피언 {mods.count(0)}/{len(mods)} → {'생존' if mods.count(0)==0 else '사망!'}")
