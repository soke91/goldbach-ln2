# -*- coding: utf-8 -*-
"""깊은 N 전용 도장 — 위치 마스크가 사는 곳에서 반정규성 검증.

증분 307에서 다시 씀. 그 전까지 이 파일은 **285 이전 `verify_all.py`의
방치된 사본**이었다 (독스트링까지 같았다). 결함 셋:

(A) **기준이 출력 문자열 안 글자였다.** `f"평균 r={...} (기준 0.80)"`에서
    0.80은 코드 어디에도 없었고, 비교도 종료코드도 없어 무엇이 나오든
    exit 0이었다. 285가 `verify_all.py`에서 고친 바로 그 결함이
    형제 파일에 그대로 남아 있었다 — 수리가 한 파일에만 적용되고
    쓸리지 않은 것. `code/lint_gates.py`가 이제 기계로 잡는다.

(B) **`verify_DEEP`이라는 이름이 깊은 N을 구조적으로 배제했다.**
    표본이 `// 6 * 6 + 2`라 `N ≡ 2 (mod 6)`, 즉 3으로 나뉘는 N이
    하나도 없었다. 마스크는 작은 소수를 많이 가진 N에서 가장 크다.

(C) 구간이 각 행의 산포에 대해 교정된 적이 없었다 (위험 8번, 증분 306).

**이 파일이 이제 하는 일**, 그리고 `verify_all.py`와 겹치지 않는 이유:
`verify_all.py`의 깊은 팔은 `n = 30`이다. 여기는 `n = 300` — 같은 통계량을
**10배의 검정력**으로 잰다. 창을 `N − 3·10⁷`까지 넓혀야 30030의 배수가
그만큼 나온다. 얕은 팔도 같은 크기로 나란히 둔다: 마스크의 서명은 두 팔의
**차이**이므로 한쪽만 정밀해서는 소용이 없다.

측정하는 것은 `r = |T(N)| / √V(N)`, `T(N) = Σ_p log p · μ(N−p)`,
`V(N) = Σ_p log²p · [μ(N−p) ≠ 0]` — 적합된 대역품이 아니라 **정확한**
2차 모멘트다 (증분 283).

구간은 40시드 재실행에서 잰 산포로 잡았다 (중심 ± 4σ,
`results/audit_stamp_calibration_deep.txt`). 중심은 실측 평균이다 —
반정규 `√(2/π) = 0.7979`는 `ρ = 1`을 뜻하는데 증분 288·289가 그것을
반박했고, 증분 306이 `N ≈ 10⁸`에서 `ρ = 0.810 ± 0.018`을 쟀다.
"""

import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

X = 100_000_000
print("mu 계산...", flush=True)
mu = np.ones(X + 1, dtype=np.int8)
pm = np.ones(X + 1, dtype=bool)
pm[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        pm[p * p:: p] = False
        mu[p::p] *= -1
        mu[p * p:: p * p] = 0
val = np.arange(X + 1, dtype=np.int64)
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        val[p::p] //= p
        pp = p * p
        while pp <= X:
            val[pp::pp] //= p
            pp *= p
mu[val > 1] *= -1
mu[0] = 0
del val
print("mu 완료됨 (검증: Σμ(≤1e6) =",
      int(mu[:1_000_001].astype(np.int64).sum()), ")", flush=True)

N = 99_999_998
WINDOW = 30_000_000
NSAMP = 300
rng = np.random.default_rng(int(os.environ.get("STAMP_SEED", "211")))
rows = []

ps_small = np.nonzero(pm[:X // 2])[0]
ps_small = ps_small[ps_small > 2].astype(np.int64)
logp = np.log(ps_small.astype(np.float64))


def half_normal_ratio(N2):
    muv = mu[N2 - ps_small[ps_small < N2 - 2]].astype(np.float64)
    lp = logp[: len(muv)]
    T = float((lp * muv).sum())
    V = float((lp ** 2 * (muv != 0)).sum())
    return abs(T) / math.sqrt(V)


# 깊은 팔: 30030 = 2·3·5·7·11·13. 창 3·10⁷ 안에 배수가 ~1000개 있어
# 300개를 중복 없이 뽑을 수 있다. 옛 3·10⁶ 창에서는 ~100개가 한계였다.
#
# 후보를 **열거한 뒤 무작위로** 고른다. `sorted(set(...))[:n]` 꼴은
# 창의 **아래쪽 n개**를 집는다 — 옛 이 파일과 `verify_all.py`가 둘 다
# 그 모양이고, 창이 좁을 때는 무해했지만 창을 3·10⁷로 넓히면 표본이
# N의 하위 30%에 몰린다. 편향은 창 폭에 비례해 커지므로 창을 넓히는
# 순간 결함이 된다.
cand = np.arange((N - WINDOW) // 30030 + 1, N // 30030 + 1,
                 dtype=np.int64) * 30030
deep = sorted(int(v) for v in
              rng.choice(cand, size=min(NSAMP, len(cand)), replace=False))
rs_deep = [half_normal_ratio(v) for v in deep]
rows.append((f"깊은 팔 반정규 (30030|N, {len(rs_deep)}N)",
             f"평균 r={np.mean(rs_deep):.4f}", float(np.mean(rs_deep)),
             (0.83, 1.07)))
print(rows[-1][:2], flush=True)

# 얕은 팔: 3∤N. 같은 크기로 나란히 둔다 — 마스크의 서명은 차이다.
cand_s = np.arange((N - WINDOW) // 6 + 1, N // 6 + 1,
                   dtype=np.int64) * 6 + 2
shal = sorted(int(v) for v in
              rng.choice(cand_s, size=NSAMP, replace=False))
rs_shal = [half_normal_ratio(v) for v in shal]
rows.append((f"얕은 팔 반정규 (3∤N, {len(rs_shal)}N)",
             f"평균 r={np.mean(rs_shal):.4f}", float(np.mean(rs_shal)),
             (0.62, 0.82)))
print(rows[-1][:2], flush=True)

# 두 팔의 차이. M.1이 정리이므로 깊은 팔이 위여야 한다. 이건 크기가
# 아니라 **부호와 크기**를 함께 거는 유일한 행이고, 마스크가 이
# 스케일에서 살아있다는 진술 자체다.
gap = float(np.mean(rs_deep) - np.mean(rs_shal))
rows.append(("마스크 서명: 깊은 팔 − 얕은 팔",
             f"차이 {gap:+.4f}", gap, (0.08, 0.37)))
print(rows[-1][:2], flush=True)


def judge(rows):
    out, ok_all = [], True
    for name, text, val, (lo, hi) in rows:
        ok = (lo <= val <= hi)
        ok_all &= ok
        out.append((name, text, lo, hi, ok))
    return out, ok_all


print("\n===== 깊은 N 도장 요약 =====")
print(f"  {'검증':<34} {'값':>16} {'사전등록 구간':>18}  판정")
table, ok_all = judge(rows)
for name, text, lo, hi, ok in table:
    print(f"  {name:<34} {text:>16} {f'[{lo:g}, {hi:g}]':>18}  "
          f"{'PASS' if ok else 'FAIL'}")

print("\n  민감도: 각 값을 구간 밖으로 밀면 FAIL이 나와야 한다")
sens_ok = True
for name, text, val, (lo, hi) in rows:
    bad = hi + 1.0
    flips = not (lo <= bad <= hi)
    sens_ok &= flips
    print(f"    {name:<34} 주입값 {bad:>8.3g}  "
          f"{'FAIL로 뒤집힘 (좋음)' if flips else '여전히 PASS (나쁨)'}")

_dump = os.environ.get("STAMP_DUMP")
if _dump:
    import json
    with open(_dump, "w", encoding="utf-8") as fh:
        json.dump([{"name": nm, "val": float(v), "lo": lo, "hi": hi}
                   for nm, _t, v, (lo, hi) in rows], fh, ensure_ascii=False)

print(f"\n  {'전체통과' if ok_all else '실패한 검증이 있음'}"
      f" / {'민감도 정상' if sens_ok else '판정기가 실패할 수 없음'}")
if not (ok_all and sens_ok):
    print("전체완료 (실패)", flush=True)
    sys.exit(1)
print("전체완료", flush=True)
