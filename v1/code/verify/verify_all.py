"""재현성 종합검증 — 핵심 엔진 축약판 일괄 실행 (공개 리포용 CI-도장).

각 검증의 축약판 (표본 축소)을 순차 실행, 요약표 출력.

증분 285 감사에서 두 결함이 나와 고쳤다.

(A) **이 도장은 거짓이 될 수 없었다.** `assert`도 `sys.exit`도 없었고,
    "기준 0.80" 같은 문구는 **코드가 아니라 출력 문자열 안의 글자**여서
    어떤 결과가 나오든 exit 0이었다. 증분 272에서 명명하고
    `verify_propositions.py`에서 고친 바로 그 결함(위험 6번 셋째 형태)이
    STATUS가 "the heavy CI stamp"라 부르는 이 파일에 그대로 있었다.
    이제 각 검증이 사전 등록된 구간을 갖고, 벗어나면 FAIL을 찍고
    종료코드 1로 죽는다. 그리고 끝에서 **일부러 어긋난 값**을 같은
    판정기에 넣어 FAIL이 실제로 나오는지 보인다 — 검사가 실패할 수
    있다고 주장하는 것과 보이는 것은 다르다(증분 276).

(B) **V1이 얕은 N만 봤다.** `// 6 * 6 + 2`는 표본을 전부
    `N ≡ 2 (mod 6)`로 만들어 **3으로 나뉘는 N을 통째로 배제**했다.
    위치 마스크는 작은 소수를 많이 가진 N에서 가장 크고(증분 240:
    깊은 셀에서 −5~−7 sd), 그 N이 구조적으로 표본에서 빠져 있었다.
    이제 얕은 팔과 **깊은 팔(N ≡ 0 mod 30030)**을 함께 돌린다.

정규화에 대해서는 감사 결과가 깨끗했다: V1의 `Σ log²p·[μ≠0]`,
V5·V6의 직접 센 받침은 모두 **정확한** 2차 모멘트다. 증분 283이
문서에서 찾아낸 `𝔖N` 대역품은 이 도장에 들어온 적이 없다.

증분 306에서 **위험 8번**을 이 도장 자신에게 겨눴고, 세 번째 결함이
나왔다 (`code/audit_stamp_calibration.py`).

(C) **구간 폭이 각 행의 산포에 대해 교정된 적이 없었다.** 285가 고친
    것은 "판정기가 실패할 수 있는가"였지 "구간이 옳은 폭인가"가
    아니다. 시드를 40개로 바꿔 돌리자 도장 전체가 **25%의 시드에서
    exit 1**을 냈다 — 시드 211의 초록불은 운이었다. 얕은 팔은
    자기 구간 안쪽으로 **1.14σ**밖에 안 들어와 있었고(혼자 ~13%
    실패), E1은 1.65σ였다.

    고친 방식: 폭은 **실측 산포에서** 잡는다(중심 ± 4σ, 40시드).
    중심은 이론값이 살아있는 행에서는 이론값 그대로 두고(사다리 0,
    인수분해 법칙 1, E1 1, 반정규 0.7979), **이론값이 이미 반박된
    두 행에서만** 실측 평균으로 옮긴다 — V1의 두 팔이다. 통과하도록
    넓힌 것이 아니라 폭의 근거를 바꾼 것이고, 결과 여유는 3.4~4.1σ로
    사전등록 창 [2, 20] 안에 든다.

(D) **V1의 옛 중심은 이 프로그램이 이미 반박한 예측이었다.**
    구간 [0.60, 1.00]은 반정규 `√(2/π) = 0.7979`, 즉 `ρ = 1`을
    중심으로 잡은 것인데 증분 288·289가 `ρ = Var T / V < 1`을
    확정했다. 40시드 평균은 얕은 팔에서 `r̄ = 0.7181 ± 0.0163`,
    이론값에서 **4.9σ 아래**다. `r̄ = √(2ρ/π)`로 풀면
    **`N ≈ 10⁸`에서 `ρ = 0.810 ± 0.018`** — FFT 작업이 닿는
    `1.6·10⁷`보다 한 자리 위에서의 독립 측정이다. 깊은 팔이
    `0.8018`로 더 높은 것은 마스크가 결정론적 성분을 얹기 때문이고,
    그래서 두 팔의 차이 자체가 마스크의 서명이다.

    한편 **유도해 적어둔 표준오차 둘은 옳았다**: 이음새의
    `SE = 0.0953`(반정규, n=40) 대 실측 `0.0904`(0.95배), E1의
    "실SE ~0.15" 대 실측 `0.1320`(0.88배). 오차 분석은 멀쩡했고,
    그것을 구간 폭에 대조하는 단계가 없었을 뿐이다.
"""

import math
import os
import sys

import numpy as np
from math import gcd as gcd0

# 이 파일은 한국어와 수학 기호를 찍는다. 콘솔이 cp949면 U+2014 같은
# 문자에서 UnicodeEncodeError로 죽는데, 그러면 도장이 "실패"가 아니라
# "예외"로 끝나 판정 자체가 사라진다. 출력 인코딩을 고정한다.
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
SQ = int(N ** 0.5)
# 시드는 기본 211로 고정된다 — 기본 실행은 종전과 바이트 동일하다.
# 환경변수 STAMP_SEED는 증분 306의 교정 감사(위험 8번)가 이 도장을
# 여러 표본에서 다시 돌려 **각 행의 산포**를 재기 위한 것이다.
# 사전등록 구간은 그 산포에 대해 교정된 적이 없었다.
rng = np.random.default_rng(int(os.environ.get("STAMP_SEED", "211")))
rows = []

# V1: 핵 총합 T(N) 반정규 (30 N)
# 정규화는 V = Σ log²p·[μ≠0], 즉 **정확한** 2차 모멘트다. 적합된
# 𝔖N 대역품이 아니다 (증분 283 참조).
ps_small = np.nonzero(pm[:X // 2])[0]
ps_small = ps_small[ps_small > 2].astype(np.int64)
logp = np.log(ps_small.astype(np.float64))


def half_normal_ratio(N2):
    muv = mu[N2 - ps_small[ps_small < N2 - 2]].astype(np.float64)
    lp = logp[: len(muv)]
    T = float((lp * muv).sum())
    V = float((lp ** 2 * (muv != 0)).sum())
    return abs(T) / math.sqrt(V)


# 얕은 팔: 옛 표본. `// 6 * 6 + 2`가 N ≡ 2 (mod 6)을 강제해 3|N을
# 전부 배제하므로, 이건 **얕은 팔이라고 이름 붙여** 남긴다.
rs_shallow = [half_normal_ratio(int(v)) for v in
              sorted(rng.integers(N - 3_000_000, N, 30) // 6 * 6 + 2)]
# 깊은 팔: 마스크가 사는 곳. N = 30030·j 는 3·5·7·11·13 으로 나뉜다.
deep = sorted({int(v) // 30030 * 30030 for v in
               rng.integers(N - 3_000_000, N, 60)})
deep = [v for v in deep if v > 2][:30]
rs_deep = [half_normal_ratio(v) for v in deep]
rows.append(("핵 T(N) 반정규 — 얕은 팔 (3∤N, 30N)",
             f"평균 r={np.mean(rs_shallow):.3f}", np.mean(rs_shallow),
             (0.30, 1.14)))
print(rows[-1][:2], flush=True)
rows.append((f"핵 T(N) 반정규 — 깊은 팔 (30030|N, {len(rs_deep)}N)",
             f"평균 r={np.mean(rs_deep):.3f}", np.mean(rs_deep),
             (0.40, 1.20)))
print(rows[-1][:2], flush=True)

# V2: 사다리 항등식 (10쌍, p=3)
errs = []
for _ in range(10):
    k = int(rng.integers(500, 1200))
    ms = np.arange(SQ + 1, (N - 1) // k + 1, dtype=np.int64)
    sel = ms[ms % 3 == 0]
    A = int((mu[sel].astype(np.int16) * mu[N - k * sel]).sum(dtype=np.int64))
    m2 = np.arange(SQ // 3 + 1, (N - 1) // (3 * k) + 1, dtype=np.int64)
    m2 = m2[m2 % 3 != 0]
    B = int((mu[m2].astype(np.int16) * mu[N - 3 * k * m2]).sum(dtype=np.int64))
    errs.append(abs(A + B))
rows.append(("사다리 항등식 (10쌍)", f"최대 오차 {max(errs)}",
             float(max(errs)), (0.0, 0.0)))
print(rows[-1][:2], flush=True)

# V3: 분산 엔진 (50쌍)
K0 = int(X ** 0.4)
rr = []
for _ in range(50):
    k1, k2 = int(rng.integers(K0 // 2, K0)), int(rng.integers(K0 // 2, K0))
    if k1 == k2:
        continue
    M = (N - 1) // max(k1, k2)
    msv = np.arange(1, M + 1, dtype=np.int64)
    prod = (mu[N - k1 * msv].astype(np.int16) * mu[N - k2 * msv]).astype(float)
    v = int(np.count_nonzero(prod))
    rr.append(abs(prod.sum()) / math.sqrt(max(v, 1)))
rows.append(("분산 비대각 (50쌍)", f"평균 r={np.mean(rr):.3f}",
             float(np.mean(rr)), (0.51, 1.08)))
print(rows[-1][:2], flush=True)

# V4: 이음새 대역
#
# 증분 285에서 이 검증이 처음 FAIL(r=0.401)을 냈고, 진단을 두 번
# 틀렸다. 둘 다 그럴듯했고 둘 다 데이터가 반박했다.
#
#   틀린 진단 1: "구간을 유비로 잡았다" — 사실이지만(위험 4번)
#     값 0.401의 원인은 아니었다.
#   틀린 진단 2: "V5와 달리 공유 클래스를 섞어 재서 낮다" —
#     분리해 보니 공유 클래스도 r=0.808로 반정규였다.
#
# 진짜 원인: 옷 코드가 `abs(sum)/sqrt(max(v,1))`로 **받침이 0인
# 쌍을 r = 0으로 평균에 넣고** 있었다. 그건 M.3가 예측하는
# 소멸이고, 법칙이 아무 주장도 하지 않는 쌍이다. 소멸된 쌍을
# "측정된 0"으로 세면 평균이 내려간다 — 마스크 항목을
# 측정치로 오독한 것이고, 이 캠페인이 반복해서 잡아온 종류다.
#
# 그래서: 받침 50 미만을 제외하고 그 개수를 찍는다. 자유/공유
# 분리는 유지한다 — 원인은 아니었지만 V5와 규율을 맞추는 것이
# 옳고, 둘을 따로 보는 덗분에 진단 2가 반박됐다.
# 표본 크기: N이 짝수라 짝수 k는 전부 공유 클래스로 빠진다.
# 40쌍을 모을 때까지 돌린다 — 6쌍짜리 게이트는 거의 검정력이 없다.
N_FREE = 40
n_low4 = 0
rr2 = []
rr2_shared = []
tries4 = 0
while len(rr2) < N_FREE and tries4 < 4000:
    tries4 += 1
    k = int(rng.integers(252, 464))
    kp = int(rng.integers(252, 464))
    if k == kp:
        continue
    free = gcd0(k * kp, N) == 1
    P1 = min(110_000, (N - 2) // max(k, kp))
    P0 = P1 // 2
    ps2 = np.arange(P0, P1, dtype=np.int64)
    ps2 = ps2[pm[P0:P1]]
    w = N - ps2 * k
    wp = N - ps2 * kp
    ok = (w > 1) & (wp > 1)
    vals = mu[w[ok]].astype(np.float64) * mu[wp[ok]]
    v = int(np.count_nonzero(vals))
    if v < 50:
        n_low4 += 1
        continue
    r = abs(vals.sum()) / math.sqrt(v)
    (rr2 if free else rr2_shared).append(r)
# 구간은 유비가 아니라 **표본 크기에서 유도**한다.
# 반정규 표준편차 sqrt(1-2/pi)=0.6028, n=40 ⇒ SE=0.0953,
# 0.798 ± 2.5·SE = [0.56, 1.04]. n은 보기 전에 고정됐다.
rows.append((f"이음새 대역 — 자유 클래스 ({len(rr2)}쌍)",
             f"평균 r={np.mean(rr2):.3f}", float(np.mean(rr2)),
             (0.43, 1.16)))
print(rows[-1][:2], flush=True)
print(f"  [진단] 받침 50 미만이라 제외된 쌍: {n_low4} "
      f"({100*n_low4/max(n_low4+len(rr2)+len(rr2_shared),1):.0f}%) — "
      f"옛 코드는 이걸 r=0으로 평균에 넣었다. 산수가 맞는다: "
      f"살아남은 비율 × 반정규 ≈ 0.55×0.81 ≈ 0.45, 옛 값 0.401.",
      flush=True)
print(f"  [참고, 게이트 아님] 이음새 공유 클래스 "
      f"({len(rr2_shared)}쌍): 평균 r="
      f"{np.mean(rr2_shared) if rr2_shared else float('nan'):.3f}"
      f"  — §6·§7이 기록한 분산 억제 클래스. "
      f"깨끗한 예측이 없으므로 판정하지 않는다.", flush=True)

# V5: 인수분해 법칙 (144~154차분) — 자유-클래스 가우시안 + 마스크 블라인드
from math import gcd
K0f, K1f = 2000, 4000
P0f = N // (2 * K1f)
ps_f = np.nonzero(pm[P0f:2 * P0f])[0] + P0f
m2s = []
n_free = 0
tries = 0
while n_free < 400 and tries < 20000:
    tries += 1
    k1 = int(rng.integers(K0f, K1f)); k2 = int(rng.integers(K0f, K1f))
    if k1 == k2 or gcd(k1 * k2, N) != 1:
        continue
    ppf = ps_f[ps_f <= (N - 2) // max(k1, k2)]
    if len(ppf) < 200:
        continue
    t = mu[N - ppf * k1].astype(np.int64) * mu[N - ppf * k2]
    nz = int(np.count_nonzero(t))
    if nz > 50:
        m2s.append(float(t.sum()) ** 2 / nz)
        n_free += 1
m2f = float(np.mean(m2s))
rows.append(("인수분해 법칙: 자유-클래스 분산비",
             f"m2_eff={m2f:.3f} ({n_free}쌍)", m2f, (0.72, 1.28)))
print(rows[-1][:2], flush=True)

# V6: E1 비율 미니 도장 (100 k, 대역 [3000,4000))
ks6 = rng.choice(np.arange(3000, 4000), size=100, replace=False)
S2 = 0.0; SS = 0.0
for k6 in ks6:
    k6 = int(k6)
    ms6 = np.arange(SQ + 1, N // k6 + 1, dtype=np.int64)
    t6 = mu[ms6].astype(np.int64) * mu[N - k6 * ms6]
    S2 += float(t6.sum()) ** 2
    SS += int(np.count_nonzero(t6))
rows.append(("E1 비율 (100k, 실SE ~0.15)",
             f"{S2/SS:.3f}", float(S2/SS), (0.47, 1.53)))
print(rows[-1][:2], flush=True)

def judge(rows):
    """사전 등록 구간에 대해 각 행을 판정. 반환: (표, 전체통과?)"""
    out = []
    ok_all = True
    for name, text, val, (lo, hi) in rows:
        ok = (lo <= val <= hi)
        ok_all &= ok
        out.append((name, text, lo, hi, ok))
    return out, ok_all


print("\n===== 종합검증 요약 =====")
print(f"  {'검증':<38} {'값':>18} {'사전등록 구간':>18}  판정")
table, ok_all = judge(rows)
for name, text, lo, hi, ok in table:
    print(f"  {name:<38} {text:>18} {f'[{lo:g}, {hi:g}]':>18}  "
          f"{'PASS' if ok else 'FAIL'}")

# 이 판정기가 실제로 FAIL을 낼 수 있는지 보인다. 주장하는 것과
# 보이는 것은 다르다 — 증분 272·276, 위험 6번 셋째 형태.
print("\n  민감도: 각 값을 구간 밖으로 밀면 FAIL이 나와야 한다")
sens_ok = True
for name, text, val, (lo, hi) in rows:
    bad = hi + 1.0 if hi > lo else 1.0
    flips = not (lo <= bad <= hi)
    sens_ok &= flips
    print(f"    {name:<38} 주입값 {bad:>8g}  "
          f"{'FAIL로 뒤집힘 (좋음)' if flips else '여전히 PASS (나쁨)'}")

# STAMP_DUMP가 있으면 판정 결과를 JSON으로 남긴다. 감사용이며
# 기본 실행에서는 아무 일도 하지 않는다.
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
