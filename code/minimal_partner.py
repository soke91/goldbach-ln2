"""최소 골드바흐 파트너 r(n) 기록 사냥.

r(n) = min{p 소수 : n - p 소수}.  골드바흐 ⟺ 모든 짝수 n ≥ 4에서 r(n) 존재.
증명 급소: r(n) ≤ f(n) 꼴의 '오차 없는 상계' 하나면 충분하다.
후보 포락선: r(n) < C·ln²n (크라메르식 사고에서 나오는 자연스러운 스케일).

기록 갱신 지점(champion)만 추적 — 포락선을 뚫는 놈이 있는지 본다.
"""

import sys
from math import log

import numpy as np

from goldbach.sieve import primes_upto, sieve

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 20_000_000

is_p = sieve(LIMIT)
small = [int(p) for p in primes_upto(100_000)]

records = []          # (n, r(n))
best = 0
freq = {}             # r 값 분포
for n in range(4, LIMIT + 1, 2):
    for p in small:
        if is_p[n - p]:
            if p > best:
                best = p
                records.append((n, p))
            freq[p] = freq.get(p, 0) + 1
            break
    else:
        print(f"!!! r({n}) > {small[-1]} — 탐색 한계 초과 (반례 후보 아님, 한계 늘릴 것)")
        break

print(f"범위: 짝수 4..{LIMIT:,} ({(LIMIT-2)//2:,}개)")
print(f"\n[r(n) 분포 상위 — 대부분의 짝수는 아주 작은 파트너로 해결]")
total = sum(freq.values())
cum = 0
for p in sorted(freq)[:8]:
    cum += freq[p]
    print(f"  r(n)={p:>3}: {freq[p]/total:>7.2%} (누적 {cum/total:.2%})")

print(f"\n[기록 갱신 지점 — r(n)의 최악 사례들]")
print(f"{'n':>14} {'r(n)':>6} {'ln²n':>8} {'r/ln²n':>7}")
for n, p in records:
    l2 = log(n) ** 2
    print(f"{n:>14,} {p:>6} {l2:>8.1f} {p/l2:>7.3f}")

worst = max(p / log(n) ** 2 for n, p in records)
print(f"\n포락선 검정: max r(n)/ln²n = {worst:.3f}")
print(f"→ 후보 법칙 'r(n) < {worst + 0.05:.2f}·ln²n' — 이 범위에서 위반 0건")
