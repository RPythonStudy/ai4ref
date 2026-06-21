#!/usr/bin/env python3
# src/alert/run.py
"""ai4ref alert — 랜드마크 감시 오케스트레이터.

흐름:  수집(PubMed 증분) → 기계필터 → efetch → LLM게이트 판정 → dedup(seen-set) → sink fan-out

  python src/alert/run.py --demo                  # 가짜 랜드마크로 배관 검증
  python src/alert/run.py --collect-only          # 실 PubMed 후보만 확인(전송 안 함)
  python src/alert/run.py --collect-only --reldate 60 --limit 20
  python src/alert/run.py                         # 전체 실행 (claude_cli 백엔드 권장)

증분 검색: (검색식) AND pdat≥T(지침 이후) · datetype=edat reldate=N(최근 색인).
dedup 은 seen-set 이 권위 — edat 창은 재조회량을 줄이는 효율 장치.
"""
import os
import sys
import json
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))   # src/ 를 import 경로에

from common.features import load_features, is_enabled, log, project_root   # noqa: E402
from common.state import SeenStore                                  # noqa: E402
from common.pubmed import esearch, efetch_meta                      # noqa: E402
from notify.base import LandmarkItem                                # noqa: E402
from notify.registry import build_sinks, fan_out                    # noqa: E402
from alert.llm_gate import judge                                    # noqa: E402


def _demo_candidates() -> list:
    return [
        LandmarkItem(kq="ERAS 수액전략", pmid="29742967",
                     title="Restrictive vs Liberal Fluid Therapy in Major Abdominal Surgery (RELIEF)",
                     journal="NEJM", year="2018"),
        LandmarkItem(kq="ERAS 수액전략", pmid="24842135",
                     title="Effect of a Perioperative Goal-Directed Therapy Protocol (OPTIMISE)",
                     journal="JAMA", year="2014"),
    ]


def _load_kqs() -> list:
    """search_collections.json 에서 감시 대상 KQ(enabled + term + guideline)만."""
    path = project_root() / "config" / "search_collections.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for r in data.get("search_collections", []):
        if r.get("enabled") and r.get("term") and r.get("guideline"):
            out.append(r)
    return out


def _collect_candidates(features, reldate: int = 30, limit: int = None) -> list:
    """PubMed 증분 검색 → efetch(제목·초록) → 후보 LandmarkItem 목록."""
    kqs = _load_kqs()
    if not kqs:
        log.warning("[collect] 감시 KQ 없음 (search_collections.json: enabled+term+guideline 필요)")
        return []
    items = []
    for kq in kqs:
        term, T = kq["term"], kq["guideline"]["date"]
        Ty = int(T[:4])
        term_dated = f"({term}) AND ({Ty + 1}:3000[dp])"            # 지침 이후(pdat)
        pmids = esearch(term_dated, datetype="edat", reldate=reldate)   # 최근 색인(edat)
        log.info(f"[collect] {kq['kq']}: 최근 {reldate}일 신규 색인 {len(pmids)}건")
        if limit:
            pmids = pmids[:limit]
        metas = efetch_meta(pmids)
        for pmid in pmids:
            m = metas.get(pmid, {})
            items.append(LandmarkItem(
                kq=kq["kq"], pmid=pmid, title=m.get("title", ""),
                journal=m.get("journal", ""), year=m.get("year", ""),
                abstract=m.get("abstract", ""),
            ))
    return items


def run(demo: bool = False, collect_only: bool = False,
        reldate: int = 30, limit: int = None, features_path: str = None) -> None:
    features = load_features(features_path)
    if not is_enabled(features, "alert"):
        log.warning("[run] alert 기능 비활성(features.yml). 종료.")
        return

    # --collect-only: 실수집 배관만 확인 (판정·전송 없음)
    if collect_only:
        cands = _collect_candidates(features, reldate=reldate, limit=limit)
        log.info(f"[collect-only] 후보 {len(cands)}건 (전송 안 함)")
        for c in cands[:20]:
            print(f"- {c.pmid}  {c.title[:88]}  ({c.journal} {c.year})")
        if len(cands) > 20:
            print(f"... 외 {len(cands) - 20}건")
        return

    sinks = build_sinks(features)
    backend = (features.get("alert", {}) or {}).get("llm_backend", "stub")

    # 안전장치: stub 백엔드로 실수집 전체 전송은 모두 '랜드마크'가 되어 스팸 → telegram 차단
    if not demo and backend == "stub":
        sinks = [s for s in sinks if s.name == "stdout"]
        log.warning("[run] stub 백엔드 실수집 → 안전상 stdout 만 출력(telegram 차단). "
                    "실판정은 features.yml 의 llm_backend: claude_cli 설정 후.")

    if not sinks:
        log.warning("[run] 활성 sink 없음 — features.yml 의 notify.stdout 를 켜세요.")

    state_path = (features.get("state", {}) or {}).get("path", "state/alerted.jsonl")
    seen = SeenStore(state_path)

    # KQ → 현행 지침 라벨 (게이트 프롬프트에 '무엇과 다른가' 기준 제공)
    gl_map = {}
    for kq in _load_kqs():
        g = kq.get("guideline", {})
        gl_map[kq["kq"]] = f"{g.get('name','')} ({g.get('date','')})".strip()

    candidates = _demo_candidates() if demo else _collect_candidates(features, reldate=reldate, limit=limit)

    sent = skipped = 0
    for c in candidates:
        if seen.is_seen(c.key):
            skipped += 1
            continue
        verdict = judge(
            {"title": c.title, "abstract": c.abstract},
            {"kq": c.kq, "guideline": gl_map.get(c.kq, "")},
            backend=backend,
        )
        if not verdict.get("is_landmark"):
            continue
        c.reason = verdict.get("reason", "")
        results = fan_out(c, sinks)
        seen.mark(c.key, {"pmid": c.pmid, "title": c.title, "sinks": results})
        sent += 1

    log.info(f"[run] 완료 — 전송 {sent} · 중복 skip {skipped} · seen 총 {len(seen)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="ai4ref 랜드마크 감시")
    ap.add_argument("--demo", action="store_true", help="가짜 랜드마크로 배관 검증")
    ap.add_argument("--collect-only", action="store_true", help="실 PubMed 후보만 확인(전송 안 함)")
    ap.add_argument("--reldate", type=int, default=30, help="edat 최근 N일 창 (기본 30)")
    ap.add_argument("--limit", type=int, default=None, help="후보 수 상한(테스트용)")
    ap.add_argument("--features", default=None, help="features.yml 경로 override")
    args = ap.parse_args()
    run(demo=args.demo, collect_only=args.collect_only,
        reldate=args.reldate, limit=args.limit, features_path=args.features)
