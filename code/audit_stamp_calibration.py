# -*- coding: utf-8 -*-
"""
Is the CI stamp's own gate calibrated? (increment 306)

HAZARD 8, named at increment 304 and not yet swept. A pre-registered
tolerance means nothing until the target's own spread is measured. It
was named because increment 304c's RULE B "passed 7 of 8 bands within
1.5x" on a statistic whose per-draw spread turned out to be 86.9%: at
R = 8 that is a standard error near 31%, so the rule could not have
discriminated anything, and the pass was noise.

THE OBVIOUS PLACE TO POINT IT. `code/verify_all.py` is this program's
own gate -- STATUS calls it the CI stamp, it exits 1 on failure, and
every claim downstream is "verified" by it. Increment 285 repaired it
after finding it had no assertions and no failure path, and added a
sensitivity block showing each verdict CAN flip. But a sensitivity
block shows a test is capable of failing; it says nothing about
whether the interval is the right width. Every row draws its sample
from one seed, `default_rng(211)`, and the spread of any row across
seeds has never been measured.

Two ways that goes wrong, and they are opposite:

  TOO TIGHT   the row's own sampling spread is comparable to the
              interval, so the stamp fails at random. A green stamp is
              then luck, and a red one means nothing.
  TOO LOOSE   the spread is minute against the interval, so the row
              would pass through any plausible change in the
              underlying mathematics. The row is decoration.

Both are invisible from a single run, which is the whole content of
hazard 8.

PRE-REGISTRATION (fixed before the run).

  Re-run the stamp under R = 40 seeds, unchanged in every other
  respect, and collect each row's value. Per row:

      margin = min(mean - lo, hi - mean) / sd     [in sigmas]

  (S) SELF-TEST. The ladder identity is an exact algebraic identity
      with interval [0, 0]. Its sd across seeds must be exactly 0. If
      an identity row varies with the seed it is not an identity, and
      that finding outranks everything else here.

  (C) CALIBRATION. A statistical row is well-calibrated iff
      2 <= margin <= 20. Below 2 it fails by chance more than 5% of
      the time; above 20 it cannot detect anything smaller than a
      20-sigma move. RULE: the stamp is calibrated iff every
      statistical row lands in that window.

  (F) FALSE-FAILURE RATE. The fraction of the 40 seeds at which the
      stamp as a whole would exit 1. RULE: below 10%.

  (D) TWO STATED STANDARD ERRORS get checked against the measured
      spread, since both were derived rather than observed:
        the seam row's interval is written as 0.798 +/- 2.5 SE with
        SE = 0.0953 from a half-normal at n = 40;
        the E1 row's header says "real SE ~ 0.15".
      No pass/fail is attached -- they are reported beside the
      measured sd, which is the first time either has been checked.

  WHAT WOULD REFUTE. (S) failing means an identity row is not exact.
  (C) failing in the tight direction means recorded PASSes were luck;
  in the loose direction it means those rows verify nothing. Either is
  a real finding about the gate every other claim leans on.
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

# 이 스크립트는 도장의 한국어 행 이름을 그대로 찍는다. 콘솔이 cp949면
# U+2014에서 UnicodeEncodeError로 죽는다 — verify_all.py가 자기 머리에
# 적어둔 바로 그 결함이고, 첫 실행에서 그대로 밟았다. 40회 재실행을
# 전부 끝내고 **인쇄 단계에서** 죽어 결과가 사라졌다.
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

R = 40
HERE = os.path.dirname(os.path.abspath(__file__))
# 감사 대상은 인자로 받는다. 기본은 verify_all.py — 증분 307에서
# verify_deep.py가 게이트를 갖게 되면서 같은 교정이 필요해졌다.
TARGET = sys.argv[1] if len(sys.argv) > 1 else "verify_all.py"
STAMP = os.path.join(HERE, TARGET)
CACHE = os.path.join(os.path.dirname(HERE), "results",
                     "audit_stamp_calibration_raw_"
                     + TARGET.replace(".py", "") + ".json")


def collect():
    """40회 재실행. 인쇄 전에 원자료를 먼저 디스크에 남긴다 —
    표시 단계의 예외가 계산을 날리지 않도록."""
    t0 = time.time()
    tmpdir = tempfile.mkdtemp(prefix="stampcal_")
    runs = []
    for r in range(R):
        seed = 1000 + r
        dump = os.path.join(tmpdir, f"s{seed}.json")
        env = dict(os.environ)
        env["STAMP_SEED"] = str(seed)
        env["STAMP_DUMP"] = dump
        p = subprocess.run([sys.executable, STAMP], capture_output=True,
                           env=env)
        if not os.path.exists(dump):
            print(f"  seed {seed}: no dump, exit {p.returncode}", flush=True)
            continue
        with open(dump, encoding="utf-8") as fh:
            runs.append((seed, p.returncode, json.load(fh)))
        if (r + 1) % 10 == 0:
            print(f"  {r+1}/{R} seeds  t={time.time()-t0:.0f}s", flush=True)
    with open(CACHE + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(runs, fh, ensure_ascii=False)
    os.replace(CACHE + ".tmp", CACHE)
    return runs


def main():
    if os.environ.get("STAMPCAL_REUSE") and os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as fh:
            runs = json.load(fh)
        print(f"reusing {len(runs)} cached runs from {CACHE}")
    else:
        runs = collect()

    nrows = min(len(x[2]) for x in runs)
    names = [runs[0][2][i]["name"] for i in range(nrows)]
    los = np.array([runs[0][2][i]["lo"] for i in range(nrows)])
    his = np.array([runs[0][2][i]["hi"] for i in range(nrows)])
    if os.environ.get("STAMPCAL_REUSE"):
        # 캐시된 값은 그대로 쓰되 **구간은 지금 파일에서 다시 읽는다**.
        # 구간을 고친 뒤 재감사할 때 옛 구간으로 판정하면 아무 의미가 없다.
        dump = os.path.join(tempfile.mkdtemp(prefix="stampcur_"), "cur.json")
        env = dict(os.environ)
        env.pop("STAMPCAL_REUSE", None)
        env["STAMP_SEED"] = "211"
        env["STAMP_DUMP"] = dump
        subprocess.run([sys.executable, STAMP], capture_output=True, env=env)
        with open(dump, encoding="utf-8") as fh:
            cur = json.load(fh)
        los = np.array([cur[i]["lo"] for i in range(nrows)])
        his = np.array([cur[i]["hi"] for i in range(nrows)])
        print(f"intervals re-read from {TARGET}: "
              + ", ".join(f"[{a:g}, {b:g}]" for a, b in zip(los, his)))
    vals = np.array([[x[2][i]["val"] for i in range(nrows)] for x in runs])
    exits = np.array([x[1] for x in runs])

    mean = vals.mean(axis=0)
    sd = vals.std(axis=0, ddof=1)
    ident = (his - los) == 0.0

    print(f"\nthe CI stamp under {len(runs)} seeds "
          f"(default is one seed, 211)")
    print(f"{'row':<40} {'mean':>9} {'sd':>9} {'min':>9} {'max':>9} "
          f"{'interval':>16} {'margin':>8}")
    margin = np.full(nrows, np.nan)
    for i in range(nrows):
        if sd[i] > 0:
            margin[i] = min(mean[i] - los[i], his[i] - mean[i]) / sd[i]
        nm = names[i][:39]
        mtxt = "identity" if ident[i] else (
            f"{margin[i]:.1f}" if np.isfinite(margin[i]) else "sd=0")
        print(f"{nm:<40} {mean[i]:>9.4f} {sd[i]:>9.4f} "
              f"{vals[:, i].min():>9.4f} {vals[:, i].max():>9.4f} "
              f"{f'[{los[i]:g}, {his[i]:g}]':>16} {mtxt:>8}")

    # 항등식 행이 **없는** 대상이 있다 (verify_deep.py). 첫 판은 그때
    # okS = False를 내고 요약문이 "항등식 행이 시드에 따라 변한다"는
    # 거짓말을 찍었다 — 변한 행은 없고 그런 행이 0개였을 뿐이다.
    # 없는 것을 실패로 세는 것은 판정이 아니라 결함이다.
    okS = bool((sd[ident] == 0).all()) if ident.any() else None
    stat = ~ident
    tight = stat & (margin < 2.0)
    loose = stat & (margin > 20.0)
    okC = not (tight.any() or loose.any())
    # 실패율은 기록된 종료코드가 아니라 **값과 현재 구간**에서 다시
    # 계산한다. 값은 구간과 무관하게 시드만으로 정해지므로, 이러면
    # 구간을 고친 뒤 40회를 다시 돌리지 않아도 된다 (STAMPCAL_REUSE).
    inside = (vals >= los[None, :]) & (vals <= his[None, :])
    fail_rate = float((~inside.all(axis=1)).mean())
    okF = fail_rate < 0.10

    print(f"\n    (S) identity rows have sd exactly 0: "
          f"{'n/a' if okS is None else ('PASS' if okS else 'FAIL')}  "
          f"({int(ident.sum())} identity row(s))")
    print(f"    (C) every statistical row has 2 <= margin <= 20: "
          f"{'PASS' if okC else 'FAIL'}")
    if tight.any():
        for i in np.nonzero(tight)[0]:
            print(f"        TOO TIGHT  {names[i][:44]}  "
                  f"margin {margin[i]:.2f} sigma")
    if loose.any():
        for i in np.nonzero(loose)[0]:
            print(f"        TOO LOOSE  {names[i][:44]}  "
                  f"margin {margin[i]:.1f} sigma")
    print(f"    (F) stamp-wide failure rate across seeds: "
          f"{fail_rate:.1%}  ->  {'PASS' if okF else 'FAIL'}")

    print(f"\n    (D) two stated standard errors, measured for the "
          f"first time")
    for i in range(nrows):
        if "이음새" in names[i]:
            print(f"        seam:  stated SE 0.0953 (half-normal, n=40) "
                  f"vs measured sd {sd[i]:.4f}  "
                  f"ratio {sd[i]/0.0953:.2f}")
        if "E1" in names[i]:
            print(f"        E1:    stated SE ~0.15 "
                  f"vs measured sd {sd[i]:.4f}  "
                  f"ratio {sd[i]/0.15:.2f}")

    nt, nl = int(tight.sum()), int(loose.sum())
    okS_eff = (okS is None) or okS
    if okS_eff and okC and okF:
        v = ("the gate is calibrated: every statistical row sits between "
             "2 and 20 sigma inside its interval, and no seed fails")
    elif okS_eff and okF and nl and not nt:
        v = (f"{nl} of {int(stat.sum())} statistical rows are DECORATION "
             f"-- their interval is more than 20 sigma wide, so they "
             f"would pass through any plausible change in the "
             f"mathematics they claim to check")
    elif okS_eff and nt:
        v = (f"{nt} of {int(stat.sum())} statistical rows are TOO TIGHT "
             f"-- their own sampling spread reaches their interval, so "
             f"a green stamp there is luck")
    elif okS is False:
        v = ("an identity row varies with the seed; it is not the exact "
             "identity the stamp says it is, and that outranks the "
             "calibration question")
    else:
        v = "the gate fails its calibration in more than one way"
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
