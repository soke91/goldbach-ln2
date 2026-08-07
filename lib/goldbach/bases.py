"""진법 관점 탐구 — 골드바흐 추측의 진법 불변성.

핵심 사실: 짝수성·소수성은 '수 자체'의 성질이므로 추측의 참/거짓은
진법과 무관하다. 다만 진법마다 "어떻게 보이는지"가 다르다:

- 짝수 밑(2,8,10,12,16...): 끝자리만으로 짝수 판별 가능
- 홀수 밑(3,5,7...): 끝자리로 판별 불가 — 자릿수 합의 홀짝으로 판별
- 2진수: 2를 제외한 모든 소수는 끝이 1 → 분할 p+q는 항상 끝자리 1+1
- 각 진법 b에서 소수의 끝자리는 gcd(끝자리, b)=1 인 것만 가능 (유한 예외 제외)
"""

from math import gcd

import numpy as np

from .sieve import primes_upto

DIGITS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def to_base(n: int, b: int) -> str:
    if n == 0:
        return "0"
    out = []
    while n:
        n, r = divmod(n, b)
        out.append(DIGITS[r])
    return "".join(reversed(out))


def parity_rule(b: int) -> str:
    """진법 b에서 짝수를 표기만 보고 판별하는 규칙 설명."""
    if b % 2 == 0:
        evens = ",".join(DIGITS[d] for d in range(0, b, 2))
        return f"{b}진법(짝수 밑): 끝자리가 [{evens}] 중 하나면 짝수"
    return f"{b}진법(홀수 밑): 끝자리로 판별 불가 — 자릿수 합이 짝수면 짝수"

def show_partition(n: int, p: int, q: int, bases: list[int] = [2, 8, 10, 12, 16]) -> str:
    """같은 분할 n = p + q 를 여러 진법으로 표시 — 명제는 하나, 옷만 다름."""
    lines = [f"n = p + q  (10진수: {n} = {p} + {q})"]
    for b in bases:
        lines.append(
            f"  {b:>2}진법: {to_base(n, b)} = {to_base(p, b)} + {to_base(q, b)}"
        )
    return "\n".join(lines)


def allowed_last_digits(b: int, limit: int = 100_000) -> dict[str, int]:
    """진법 b에서 소수의 끝자리 분포. gcd(d,b)=1 인 끝자리에만 몰림을 확인."""
    counts: dict[str, int] = {}
    for p in primes_upto(limit):
        d = DIGITS[int(p) % b]
        counts[d] = counts.get(d, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def coprime_digits(b: int) -> list[str]:
    """진법 b에서 (큰) 소수가 가질 수 있는 끝자리 = b와 서로소인 자릿수."""
    return [DIGITS[d] for d in range(b) if gcd(d, b) == 1]


def odd_base_patterns(b: int, limit: int = 100_000) -> dict:
    """홀수 진법 b에서만 보이는 패턴 탐구.

    홀수 진법의 고유 성질:
    - 짝수성이 끝자리가 아닌 '자릿수 합의 홀짝'에 나타남
      (b ≡ 1 (mod 2) 이므로 n ≡ 자릿수합 (mod 2))
    - 홀수 소수는 자릿수 합이 항상 홀수
    - 골드바흐 분할 p+q: 홀수합 + 홀수합 = 짝수합 구조
    - 소수 끝자리는 여전히 gcd(d,b)=1 자리에만: 3진법이면 {1,2} 둘뿐
    """
    assert b % 2 == 1, "홀수 진법 전용"
    ps = [int(p) for p in primes_upto(limit) if p > b]
    digit_sum = lambda n: sum(int(DIGITS.index(c)) for c in to_base(n, b))
    odd_sum = sum(1 for p in ps if digit_sum(p) % 2 == 1)
    dist = allowed_last_digits(b, limit)
    return {
        "진법": b,
        "홀수 소수 중 자릿수합이 홀수인 비율": f"{odd_sum/len(ps):.1%} (이론값 100%)",
        "허용 끝자리": coprime_digits(b),
        "끝자리 분포": dist,
        "비고": "짝수 판별이 끝자리가 아닌 자릿수합으로 이동 — 패턴이 '숨는' 진법",
    }


# ── 왜 12진법인가? ────────────────────────────────────────────

def compression_scan(max_base: int = 36) -> list[dict]:
    """모든 진법의 소수 끝자리 압축률 φ(b)/b 스캔.

    핵심 정리: φ(b)/b = ∏_{p|b} (1 - 1/p) — b의 '서로 다른 소인수'로만
    결정된다. 지수는 무관: 6, 12, 24는 모두 1/3로 동률.
    최소화하려면 작은 소수를 많이 나눠야 함 → 소수계승(primorial) 밑
    (2, 6, 30, 210...)이 국소 최솟값.
    """
    rows = []
    for b in range(2, max_base + 1):
        phi = len([d for d in range(b) if gcd(d, b) == 1])
        rows.append({"b": b, "φ(b)": phi, "압축률": phi / b,
                     "소인수": _prime_factors(b)})
    return sorted(rows, key=lambda r: (r["압축률"], r["b"]))


def _prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def digit_information_gain(b: int, limit: int = 1_000_000) -> dict:
    """진법 b에서 '끝자리 하나를 아는 것'이 소수 판별에 주는 정보량.

    I(소수; 끝자리) = H(소수) - H(소수|끝자리)  [비트]
    압축률이 낮을수록 끝자리가 소수 후보를 강하게 걸러줌 → 정보량 큼.
    """
    from math import log2
    ps = primes_upto(limit)
    n_total = limit - 1  # 2..limit
    n_prime = len(ps)
    p_prime = n_prime / n_total

    def h(p):
        return 0.0 if p in (0.0, 1.0) else -p * log2(p) - (1 - p) * log2(1 - p)

    h_prior = h(p_prime)
    # 끝자리별 조건부 엔트로피
    h_cond = 0.0
    for d in range(b):
        cnt_all = len(range(2 + ((d - 2) % b), limit + 1, b))
        if cnt_all == 0:
            continue
        cnt_p = int(np.sum(ps % b == d)) if cnt_all else 0
        h_cond += (cnt_all / n_total) * h(cnt_p / cnt_all)
    return {"진법": b, "사전 엔트로피": round(h_prior, 4),
            "조건부 엔트로피": round(h_cond, 4),
            "정보 이득(비트)": round(h_prior - h_cond, 4),
            "후보 압축": f"{len(coprime_digits(b))}/{b}"}


def goldbach_digit_pairs(n_mod: int, b: int = 12) -> list[tuple[str, str]]:
    """짝수 n ≡ n_mod (mod b)일 때 분할 p+q의 가능한 끝자리 쌍.

    골드바흐 혜성의 '띠'가 갈라지는 이유: n mod b에 따라 허용되는
    끝자리 쌍 개수가 달라진다 (많을수록 g(n)이 큼).
    """
    allowed = [d for d in range(b) if gcd(d, b) == 1]
    pairs = []
    for dp in allowed:
        dq = (n_mod - dp) % b
        if dq in allowed and dp <= dq:
            pairs.append((DIGITS[dp], DIGITS[dq]))
    return pairs


def base_pattern_report(bases: list[int] = [2, 8, 10, 12, 16], limit: int = 100_000) -> list[dict]:
    """진법별 소수 패턴 비교 — 진법마다 '보이는 패턴'이 어떻게 다른가.

    핵심 지표: 허용 끝자리 개수 = φ(b) (오일러 피 함수).
    φ(b)/b 가 작을수록 소수가 소수의 끝자리에 압축되어 패턴이 또렷하다.
    (12진법: 4/12 = 33%로 10진법 4/10 = 40%보다 압축적 — 12진법이
    소수 관찰에 유리하다고 말하는 이유)
    """
    rows = []
    for b in bases:
        allowed = coprime_digits(b)
        dist = allowed_last_digits(b, limit)
        rows.append({
            "진법": b,
            "짝수판별": "끝자리" if b % 2 == 0 else "자릿수합",
            "허용 끝자리 φ(b)": f"{len(allowed)}/{b} ({','.join(allowed)})",
            "압축률": f"{len(allowed)/b:.0%}",
            "실제 분포": dist,
        })
    return rows
