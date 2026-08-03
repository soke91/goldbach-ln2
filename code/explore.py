"""골드바흐 탐구 플랫폼 — 통합 데모 실행.

사용법:
    python explore.py            # 전체 데모
    python explore.py verify 10000
    python explore.py comet 20000
"""

import sys

from goldbach import bases, methods, primes, stats, verify, viz


def sec(title):
    print(f"\n{'=' * 62}\n■ {title}\n{'=' * 62}")


def demo_verify(limit=100_000):
    sec(f"1. 검증 엔진 — 4..{limit:,} 전수 검증")
    checked, fails = verify.verify_range(4, limit)
    print(f"검증한 짝수: {checked:,}개 / 반례: {fails if fails else '없음'}")
    p = verify.find_partition(limit)
    print(f"예시: {p.n:,} = {p.p:,} + {p.q:,}")


def demo_bases():
    sec("2. 진법 탐구 — 같은 명제, 다른 옷")
    print(bases.show_partition(100, 3, 97))
    print()
    for row in bases.base_pattern_report([2, 8, 10, 12, 16]):
        top = list(row["실제 분포"].items())[:6]
        print(f"{row['진법']:>2}진법 | 짝수판별: {row['짝수판별']:<4} | "
              f"소수 끝자리 {row['허용 끝자리 φ(b)']} ({row['압축률']}) | 분포 {top}")
    print("\n[홀수 진법 — 패턴이 자릿수합으로 숨는다]")
    for b in (3, 9):
        r = bases.odd_base_patterns(b)
        print(f"{b}진법: 끝자리 {r['허용 끝자리']}, "
              f"홀수소수 자릿수합 홀수 비율 = {r['홀수 소수 중 자릿수합이 홀수인 비율']}")


def demo_primes():
    sec("3. 소수의 증가·간격 패턴")
    print("[소수 정리 — pi(x)의 성장]")
    for r in primes.pnt_table([10_000, 100_000, 1_000_000]):
        print(f"  x={r['x']:>9,} pi(x)={r['pi(x)']:>7,} "
              f"x/lnx 오차 {r['오차%']:+.2f}% | Li(x) 오차 {r['Li오차%']:+.3f}%")
    print("[n번째 소수 ~ n·ln n]")
    for r in primes.nth_prime_growth([100, 10_000, 100_000]):
        print(f"  p_{r['n']} = {r['p_n']:,} vs n·lnn = {r['n·ln n']:,} (비율 {r['비율']})")
    g = primes.gap_stats(1_000_000)
    print(f"[간격] 평균 {g['평균 간격']} (이론 ln x = {g['이론 평균(ln x)']}), "
          f"최빈 {g['최빈 간격']}, 쌍둥이 {g['쌍둥이(간격2)']:,}쌍, 최대 {g['최대 간격']}")
    print(f"       분포 상위: {g['간격 분포 상위']}")


def demo_randomness():
    sec("4. 소수는 랜덤인가? — 세 가지 실험")
    c = primes.chebyshev_bias(1_000_000)
    print(f"[체비쇼프 편향] 4k+1: {c['4k+1']:,} vs 4k+3: {c['4k+3']:,} "
          f"→ {c['우세']} 우세 (+{abs(c['차이'])}) — 랜덤이면 없어야 할 치우침")
    r = primes.last_digit_repulsion(1_000_000)
    print(f"[끝자리 반발] 연속 소수 같은 끝자리 비율 {r['같은 끝자리 연속 비율']} "
          f"(랜덤 기대 {r['랜덤 기대치']}) → {r['결론']}")
    e = primes.residue_equidistribution(1_000_000, 12)
    live = {k: v for k, v in e["분포"].items() if v > 10}
    print(f"[디리클레 균등분포] mod 12 분포: {live} — 허용 잉여류 안에선 거의 균등")


def demo_stats(limit=50_000):
    sec(f"5. g(n) 통계 — 왜 참일 수밖에 없어 보이는가 (≤{limit:,})")
    for r in stats.minimum_g_growth(limit, buckets=8):
        print(f"  {r['구간']:>18} min g(n) = {r['min g(n)']:>4} (n={r['달성 n']:,}) "
              f"| HL 예측 {r['HL 예측(구간초)']}")
    print("[Hardy-Littlewood 예측 정확도]")
    for r in stats.hl_accuracy([1_000, 10_000, 50_000]):
        print(f"  n={r['n']:>7,} 실제 {r['실제 g(n)']:>4} vs 예측 {r['HL 예측']:>7} "
              f"(비율 {r['비율']})")


def demo_methods():
    sec("6. 증명 방법론 시뮬레이션")
    c = methods.circle_method_demo(10_000)
    print(f"[원 방법/FFT] g(10000) = {c['지수합 기반 g(n)']} — {c['설명']}")
    s = methods.schnirelmann_demo(50_000)
    print(f"[슈니렐만] 짝수 커버율 {s['짝수 커버율']} → {s['결론']}")
    ch = methods.chen_demo(10_000)
    for row in ch["표본"][-3:]:
        print(f"[천징룬] n={row['n']:,}: 소수+소수 {row['소수+소수']}개 "
              f"vs 소수+(소수|P2) {row['소수+(소수|P2)']}개")
    print("[Cramér 확률 모델] 반례가 나올 확률:")
    for r in methods.cramer_heuristic([100, 1_000, 10_000, 100_000]):
        print(f"  n={r['n']:>7,} 기대 g(n)={r['기대 g(n)']:>7} → P(분할없음) ≈ {r['P(분할 없음)']}")
    w = methods.weak_goldbach_demo(7, 5_000)
    print(f"[약한 골드바흐(2013 증명됨)] {w['범위']}: 반례 {w['반례'] or '없음'} — {w['결론']}")


def demo_viz(limit=20_000):
    sec("7. 시각화 생성")
    print(viz.goldbach_comet(limit, "comet.png"))
    print(viz.gap_histogram(1_000_000, "gaps.png"))
    print(viz.min_g_trend(limit, "min_g.png"))


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        demo_verify(); demo_bases(); demo_primes(); demo_randomness()
        demo_stats(); demo_methods(); demo_viz()
    elif args[0] == "verify":
        demo_verify(int(args[1]) if len(args) > 1 else 100_000)
    elif args[0] == "comet":
        demo_viz(int(args[1]) if len(args) > 1 else 20_000)
