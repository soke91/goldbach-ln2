# SPEC — 배포용 논문 작성 규약

이 폴더(`deploy/`)에 들어가는 문서의 **양식·규칙·수록 기준**을 못박는다.
연구 저장소(`v1/`, `v2/`)는 실험실 노트다. 이 폴더는 **독자에게 내보내는
것**이고, 두 문서는 목적이 다르므로 규칙도 다르다.

원칙 한 줄: **저장소는 어떻게 알게 됐는지를 남기고, 배포본은 무엇이
참인지만 남긴다.**

---

## 0. 저장소 문서와 배포본의 차이

| | `v2/paper/*.md` (저장소) | `deploy/papers/*.md` (배포본) |
|---|---|---|
| 독자 | 저자 자신, 재검증자 | 논문 심사자, 외부 수학자 |
| 목적 | 무엇을 어떻게 알게 됐는지의 완전한 기록 | 무엇이 참이고 무엇이 증명됐는지 |
| 정정 이력 | **남긴다** (v3가 뭘 틀렸는지 포함) | **지운다** (고쳐진 결과만 싣는다) |
| 실패한 사전등록 | 규칙 번호까지 남긴다 | 싣지 않는다 |
| 철회된 수치 | 철회 사실과 함께 남긴다 | 싣지 않는다 |
| Remark 개수 | 218 | 0 (필요한 것은 본문·주석으로 흡수) |
| 형식 | 마크다운 | **LaTeX (`.tex`)** — 투고 형식 |
| 게이트 | `gate.py` 77개 검사 통과 필수 | 게이트 마커(`<!-- evidence: -->`) 제거, 대신 `MANIFEST.md`로 대응 |

배포본은 저장소의 **요약이 아니라 사영(projection)** 이다. 내용을 줄이는
것이 아니라 한 축(역사)을 없앤다.

---

## 1. 파일 양식

### 1.1 형식 — LaTeX

- **LaTeX 1파일 = 논문 1편. 확장자 `.tex`.** 마크다운은 저장소의 형식이고
  배포본의 형식이 아니다.
- 문서클래스 `\documentclass[11pt,a4paper]{article}`, `amsmath`·`amssymb`·
  `amsthm`·`hyperref` 만 쓴다. 외부 스타일 파일에 의존하지 않는다
  (arXiv 업로드에서 그대로 컴파일돼야 한다).
- 정리 환경은 `amsthm` 의 `\newtheorem` 으로 선언하고, **번호는
  하나의 연속 카운터**를 공유한다:

  ```latex
  \theoremstyle{plain}
  \newtheorem{theorem}{Theorem}
  \newtheorem{corollary}[theorem]{Corollary}
  \newtheorem{proposition}[theorem]{Proposition}
  \newtheorem{lemma}[theorem]{Lemma}
  \newtheorem{conjecture}[theorem]{Conjecture}
  \theoremstyle{definition}
  \newtheorem{observation}[theorem]{Observation}
  \newtheorem{measurement}[theorem]{Measurement}
  ```

- 상호참조는 **전부 `\label`/`\ref`**. 번호를 손으로 적지 않는다.
  저장소 라벨(`thm:A`, `prop:V` …)을 `\label{thm:A}` 로 **그대로 승계**해
  추적성을 유지한다.
- 매크로는 **프리앰블에서 한 번만** 정의하고 그 뒤에는 정의를 쓰지 않는다:
  `\newcommand{\SS}{\mathfrak{S}}`, `\newcommand{\AAA}{\mathfrak{A}}`,
  `\newcommand{\Emu}{E_\mu}`, `\DeclareMathOperator{\rad}{rad}`.
- 참고문헌은 `thebibliography` 환경(외부 `.bib` 없음), 키는 저장소와
  동일(`HL`, `GY`, `Bom76`, `Li23` …).
- 한 줄 72–76자로 접는다. 표는 접지 않는다.
- 주석(`%`)에 저장소 메타(evidence 마커, 게이트 번호)를 남기지 않는다.

### 1.2 문서 골격 (고정)

```latex
\title{...}  \author{...}  \date{}
\begin{abstract} ... \end{abstract}
\section{Introduction}
  \subsection{The Huang--Li reduction}
  \subsection{Results}
  \subsection{What is not claimed}      % 필수
\section{Notation and the mechanism}
\section{...}                            % 본문: 진술 — 증명 — 필요한 수치
\section{Relation to the literature}
\section{Summary}
\begin{thebibliography}{99} ... \end{thebibliography}
\appendix
\section{Reproducibility}                % 스크립트 ↔ 결과 대응표
```

`\subsection{What is not claimed}` 는 **필수 절**이다. 이 프로그램의 결과
대부분이 음성이거나 골드바흐-중립이라, 범위를 본문 앞에서 못박지 않으면
과대주장으로 읽힌다.

### 1.3 진술 번호

한 논문 안에서 **종류를 가리지 않고 하나의 연속 번호**를 쓴다
(위 프리앰블의 `[theorem]` 공유 카운터). 결과:

```
Theorem 1, Corollary 2, Theorem 3, Proposition 4, Lemma 5, ...
```

독자가 "Lemma 5"를 찾을 때 Lemma만 세지 않아도 된다.

---

## 2. 진술 등급 — 이 규약의 핵심

배포본은 **증명된 것과 측정된 것을 이름으로 구분한다.** 저장소는 둘 다
`Proposition` 으로 부르는 자리가 있고, 그것이 배포본에서 가장 위험한
결함이다.

| 등급 | 조건 | 예 |
|---|---|---|
| **Theorem / Corollary / Lemma / Proposition** | 완전한 증명이 본문에 있다. 계산은 확인용이지 근거가 아니다. | Theorem 1 (thm:A), Lemma (lem:MP), Proposition (prop:V) |
| **Observation** | 유도는 있으나 상수·오차항이 통제되지 않았거나, 결론이 점근이 아니라 측정 범위 위의 진술이다. | prop:coh (오차막대가 $(\log N)^{-1/2}$), prop:W 의 예산 초과 부분 |
| **Measurement** | 순수 측정. 사전등록된 판정규칙과 널이 있고, 재현 스크립트가 있다. | prop:placebo(마스크 플라시보), prop:scaleinv |
| **Conjecture** | 증거는 있으나 증명이 없다. 증거의 등급을 함께 적는다. | Conjecture (conj:L) |

**규칙 T1.** `Theorem`/`Lemma`/`Proposition` 으로 부르는 진술에 증명이
없으면 등급을 내린다. 예외 없음.

**규칙 T2.** `Observation`·`Measurement` 는 **측정 범위를 진술 안에**
적는다. "over even $N \in [10^5, 1.6\cdot10^7]$" 이 진술의 일부다.

**규칙 T3.** 측정으로 뒷받침되는 모든 진술은 통계량(`STATISTIC`)과
필드(`FIELD`)를 문장으로 고정한다. "셀 평균"은 진술이 아니다. "무엇으로
색인된 셀의, 어느 범위의, 평균만 뺀 것"이 진술이다.

---

## 3. 넣을 것

### 3.1 반드시 넣는다

1. **모든 정리·보조정리·명제와 그 증명 전문.** 축약하지 않는다.
2. **범위 진술** (§1.3). "골드바흐를 향한 순진전은 0" 같은 문장은
   겸양이 아니라 정확한 진술이므로 남긴다.
3. **고쳐진 뒤의 수치.** 감사에서 재계산된 값만.
4. **재현 정보.** 스크립트 이름, 실행 대상 $N$, 결과 파일명.
5. **관례 고정.** "dyadic" 이 $K<k\le 2K$ 인지 $K\le k<2K$ 인지처럼,
   값이 달라지는 관례는 본문에 적는다. (저장소 `rem:supp`의 교훈.)
6. **덫(trap) 세 개.** $(q,N)>1$ 주항 덫, 적은 뽑기로 만든 널,
   틀린 모듈러스에서 잰 국소인자 — 이것들은 후속 연구자가 다시 밟을
   함정이므로 방법론 절에 남긴다. **단, "우리가 밟았다"가 아니라
   "여기에 덫이 있다"로 쓴다.**
7. **널과 대조군의 존재.** 어떤 널을 돌렸는지는 결과의 일부다.

### 3.2 넣되 형태를 바꾼다

| 저장소의 형태 | 배포본의 형태 |
|---|---|
| "버전 3은 X라고 적었다. 그건 규칙 X5이고 실패한다. 옳은 값은 Y다." | "Y." (각주도 없이) |
| "우리 예측이 반증됐다: 실제로는 두 조각이 같은 크기다." | "두 조각은 같은 크기다." |
| "이 널은 사양됐다가 나중에 돌렸다." | "이 주장의 대조군은 …이며, 결과는 …" |
| "K1은 재개방한다" | 열린 항목 목록에 그대로 열린 것으로 적는다 |
| Remark 105개 | 본문 문단 또는 진술 뒤 한 문단 주석 |

### 3.3 넣지 않는다 (오류·정정 계층)

배포본에서 **삭제**하는 것 — 이것이 "정정 없음"의 정의다.

1. **버전 간 정정 서술.** `"Version 3 said …"`, `"An earlier version
   claimed …"`, `"An earlier draft stated …"` 로 시작하는 모든 문단.
   → 고쳐진 진술만 남긴다.
2. **사전등록 규칙의 실패 보고.** `"the audit's rule X5 and it fails"`,
   `"rules M2 and M4, and both fail"` 등 규칙 이름(X5, C2, C3, C4, I1–I4,
   L3, M1, N4, Z1–Z6, H4, U4, B5, Y1, Y2 …)이 등장하는 모든 문장.
   → 판정이 바뀐 결과만 남긴다.
3. **철회된 수치.** `"withdrawn"`, `"not supportable"`,
   `"restated rather than corrected"` 가 붙은 값. 예: `prop:V` 의
   `0.000582`, R1·R2의 정밀도, `rem:e1row`의 "exact square-root
   cancellation", `rem:levelmeas`의 레벨 측정값.
   → 아예 인쇄하지 않는다. 대신 살아남은 약한 진술을 싣는다.
4. **자기 예측의 반증 서사.** `rem:toprdom`, `rem:cap`, `rem:band`,
   `rem:artifact` 류.
5. **게이트·저장소 운영 메타.** `<!-- evidence: -->` 마커, `gate.py`,
   G1–G77, M1–M9, `OPEN.md`·`DECISIONS.md` 참조, 접두사 규약
   (`lab_`/`audit_`/`verify_`).
   → 재현 정보는 부록의 표 하나로 대체.
6. **개인적 서술.** `"두 번 물렸다"`, `"내 출력이 걸렸다"`,
   `"the cost objection was mine"` 등 1인칭 실패 고백.
7. **한국어.** 배포본 본문은 영어. (이 SPEC과 INVENTORY만 한국어.)

### 3.4 삭제하면 **안 되는** 것 (정정처럼 보이지만 아닌 것)

이것을 지우면 논문이 틀린다. §3.3과 혼동 금지.

- **`rem:bound`의 결론.** Theorem A 의 오차항은
  $N e^{-c\sqrt{\log N}}$ 이 **아니라** $N(\log N)^{-A}$ 다.
  정정 서사는 지우고 **고쳐진 진술을 싣는다.**
- **`rem:sign`의 부호.** $A(N)\widetilde G(1) = -\mathfrak{S}(N)$.
  부호가 반대면 Theorem C 가 $2\mathfrak{S}(N)N$ 만큼 틀린다.
- **`rem:threeway`의 세 갈래.** $E_3 \ll_A N(\log N)^{-A}$ 는 골드바흐를
  주지만 점근 $\tilde r \sim \mathfrak{S}N$ 은 주지 않는다. 둘을 하나의
  동치로 적으면 틀린다.
- **`rem:trap`의 $(q,N)=1$ 제한.** 이것은 증명의 일부다.
- **`prop:onesided`가 `thm:D`를 재개방하지 않는다는 문장.**
- **`lem:coin`의 단방향성.** "코인에 재현되면 μ를 재는 게 아니다"의 역은
  성립하지 않는다. 이 단서를 빼면 이후의 모든 널 판정이 과대해석된다.
- **`rem:secondorder`의 2차 형태.** $V \sim \mathfrak{A}(N)(N\log N - N)$.
  1차 형태만 싣고 수치를 함께 실으면 5% 어긋난다.

---

## 4. 수치 인용 규칙

- **N1.** 인쇄하는 모든 소수(소수점 아래 3자리 이상)는 `deploy/results/`의
  결과 파일에서 온다. 손으로 옮겨 적지 않는다.
- **N2.** 모든 수치에는 **측정 범위**가 붙는다. 범위 없는 수는 인쇄하지
  않는다.
- **N3.** 개수는 오차 막대가 아니다. 상관된 항의 평균은 $1/\sqrt{n}$ 로
  떨어지지 않는다 (P4의 주제). 표본 크기를 정밀도처럼 쓰지 않는다.
- **N4.** 비를 인쇄할 때는 분자·분모의 **가중치**를 먼저 말한다.
  평균의 비 $\ne$ 비의 평균.
- **N5.** 두 관례가 다른 값을 주면 둘 다 인쇄하거나, 하나를 고정하고
  고정했다고 적는다.
- **N6.** 외삽은 구간(bracket)과 함께 적고, 계산 범위 밖임을 명시한다.

---

## 5. 증명 서술 규칙

- **P1.** 증명은 `\begin{proof} ... \end{proof}` 환경으로 쓴다.
  증명이 없는 진술은 `Observation`/`Measurement`/`Conjecture` 로 부른다.
- **P2.** "표준적이다", "잘 알려져 있다" 로 끝내지 않는다. 어느 문헌의
  무엇인지 적는다. (Theorem A의 재료는 전부 고전이며, 그렇게 적는다.)
- **P3.** 상수 의존성을 명시한다: $\ll_{A,\theta'}$ 처럼.
- **P4.** 절단(truncation)은 $A$ 에 의존해도 된다는 것을 적는다.
  고정 절단은 유계 $A$ 만 덮는다.
- **P5.** 증명 안에서 계산 결과를 근거로 쓰지 않는다. 계산은 확인이다.
  예외: 유한 검증(예 `lem:density`의 유리수 산술 검증)은 그렇게 적는다.

---

## 6. 인용·귀속 규칙

- **A1.** Huang–Li [HL]의 결함(식 (18)의 $n$-의존 제약 누락)은 사실로
  적되, 그들의 Theorem 1·Corollary 1이 **그대로 성립함**을 같은 문단에서
  적는다. 결함 보고가 결과 취소로 읽히면 안 된다.
- **A2.** Theorem A의 메커니즘은 고전적이다. "제공하는 것은 적용이다"라고
  적는다. 새로움을 부풀리지 않는다.
- **A3.** Lichtman [Li23]의 $66/107$ 은 **소수**에 대한 것이고 μ에 대한
  것이 아니다. $3/5$ 는 [Li20] 계열이며 원문 확인 전에는 그렇게 적는다.
- **A4.** Bombieri의 점근 체 [Bom76]와의 관계를 Theorem D 에 반드시 붙인다.
  장르가 고전임을 밝히고, 새로운 것은 **정량적 형태**($N^{1/2}$ 분리)임을
  적는다.

---

## 7. 검수 체크리스트 (배포 전 필수)

- [ ] `Version`, `earlier`, `withdrawn`, `fails`, `rule [A-Z]\d` 문자열이
      본문에 0회
- [ ] 모든 `theorem`/`lemma`/`proposition`/`corollary` 환경에
      `\begin{proof}` 존재 (규칙 T1)
- [ ] 모든 `observation`/`measurement` 에 측정 범위 존재 (규칙 T2)
- [ ] 인쇄된 모든 소수가 `deploy/results/` 의 파일에 존재 (규칙 N1)
- [ ] `\subsection{What is not claimed}` 존재
- [ ] §3.4 목록의 7개 항목이 전부 살아 있음
- [ ] `\ref` 가 전부 실재하는 `\label` 을 가리킴 (미해결 참조 0)
- [ ] `pdflatex` 2회 통과, 경고 외 오류 0
- [ ] `MANIFEST.md` 의 스크립트·결과가 `deploy/code`·`deploy/results` 에 실재
- [ ] 한국어 0자 (`papers/` 안)
