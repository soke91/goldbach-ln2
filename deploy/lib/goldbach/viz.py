"""시각화 — 골드바흐 혜성, 소수 간격, 진법별 끝자리 분포."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .primes import gap_stats
from .sieve import primes_upto
from .stats import comet_data, hardy_littlewood_estimate

plt.rcParams["font.family"] = ["Malgun Gothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False


def goldbach_comet(limit: int, path: str) -> str:
    """골드바흐 혜성: n vs g(n) 산점도 + Hardy-Littlewood 하한 곡선."""
    ns, gs = comet_data(limit)
    fig, ax = plt.subplots(figsize=(12, 7), dpi=110)
    # n mod 6 에 따라 색 구분 — 혜성의 '띠' 구조가 드러남
    for r, c, lbl in [(0, "#d62728", "n≡0 (mod 6)"), (2, "#1f77b4", "n≡2 (mod 6)"),
                      (4, "#2ca02c", "n≡4 (mod 6)")]:
        m = ns % 6 == r
        ax.scatter(ns[m], gs[m], s=2, c=c, label=lbl, alpha=0.5)
    xs = np.linspace(6, limit, 300)
    hl = [hardy_littlewood_estimate(int(x) + int(x) % 2) for x in xs]
    ax.plot(xs, hl, "k--", lw=1.2, label="Hardy–Littlewood 예측(6∤n 기준선)")
    ax.set_xlabel("짝수 n")
    ax.set_ylabel("골드바흐 분할 개수 g(n)")
    ax.set_title(f"골드바흐 혜성 (n ≤ {limit:,}) — g(n)은 0에서 멀어지기만 한다")
    ax.legend(markerscale=4)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def gap_histogram(limit: int, path: str) -> str:
    """소수 간격 히스토그램 — 6의 배수 간격이 우세한 패턴."""
    ps = primes_upto(limit)
    gaps = np.diff(ps)
    fig, ax = plt.subplots(figsize=(10, 5), dpi=110)
    vals, cnts = np.unique(gaps, return_counts=True)
    colors = ["#d62728" if v % 6 == 0 else "#1f77b4" for v in vals]
    ax.bar(vals, cnts, color=colors)
    ax.set_xlabel("연속 소수 간격")
    ax.set_ylabel("빈도")
    ax.set_title(f"소수 간격 분포 (≤{limit:,}) — 빨강 = 6의 배수 간격 (우세)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def min_g_trend(limit: int, path: str, window: int = 500) -> str:
    """이동 최솟값으로 본 g(n)의 하한 추세 — '반례 가능성'의 시각화."""
    ns, gs = comet_data(limit)
    k = max(1, window // 2)
    mins = np.array([gs[max(0, i - k): i + k].min() for i in range(len(gs))])
    fig, ax = plt.subplots(figsize=(11, 5), dpi=110)
    ax.plot(ns, mins, lw=1, color="#d62728", label=f"이동 최솟값 (창 {window})")
    ax.plot(ns, gs, ",", color="#bbbbbb", alpha=0.3)
    ax.set_xlabel("짝수 n")
    ax.set_ylabel("g(n)")
    ax.set_title("g(n)의 하한 추세 — 반례(=0)에서 계속 멀어진다")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
