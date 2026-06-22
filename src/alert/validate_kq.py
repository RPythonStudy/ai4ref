#!/usr/bin/env python3
# src/alert/validate_kq.py
"""ai4ref — key_questions.yml KQ 레코드 검증 (수정 개념 MVP)

  [0단계] 사용자가 입력한 question_type ↔ classify_qtype 보조 제안 대조 (결정형 검증)
  [검증 A] term(검색식)이 guideline_refs(지침 근거)를 재현하나 (relative recall)
  [검증 B] alert 가 지침 이후(T 이후) 랜드마크를 포착하나
"""
from __future__ import annotations

import os
import sys

# 경로 세팅
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from classify_qtype import classify_question_type
from common.query import build_query
from common.kq import load_kqs
from common.pubmed import esearch


def validate_single_kq(r: dict) -> dict:
    """KQ 레코드 1건에 대해 유형 검증 및 Recall A/B 채점을 수행합니다."""
    # [0단계] 유형 분류 대조
    classify_res = classify_question_type(r["kq"])
    sug = classify_res.get("question_type") or ""
    user = r.get("question_type") or ""
    type_match = (sug == user)

    # 시간축 (T = 지침 제정 시점)
    T = (r.get("guideline") or {}).get("date")          # "YYYY-MM" or None
    Ty = int(T[:4]) if T else None
    query = build_query(r["pico"])                       # (I) AND (P), OR 블록

    # [검증 A] 지침근거 재현율 (T 이전)
    pmids = r.get("guideline_refs", [])
    refs_results = []
    refs_hit = 0
    for pmid in pmids:
        term = f"({query}) AND {pmid}[uid]"
        if Ty:
            res_pmids = esearch(term, retmax=1, datetype="pdat", mindate="1900/01/01", maxdate=f"{Ty}/12/31")
        else:
            res_pmids = esearch(term, retmax=1)
        ok = (pmid in res_pmids)
        refs_hit += ok
        refs_results.append((pmid, ok))

    # [검증 B] 랜드마크 포착 (T 이후)
    lms = r.get("post_guideline_landmarks", [])
    landmarks_results = []
    landmarks_hit = 0
    if Ty and lms:
        for pmid in lms:
            term = f"({query}) AND {pmid}[uid]"
            res_pmids = esearch(term, retmax=1, datetype="pdat", mindate=f"{Ty + 1}/01/01", maxdate="2030/12/31")
            ok = (pmid in res_pmids)
            landmarks_hit += ok
            landmarks_results.append((pmid, ok))

    return {
        "kq": r["kq"],
        "type_match": type_match,
        "user_type": user,
        "suggested_type": sug,
        "refs_results": refs_results,
        "refs_hit": refs_hit,
        "refs_total": len(pmids),
        "landmarks_results": landmarks_results,
        "landmarks_hit": landmarks_hit,
        "landmarks_total": len(lms) if Ty else 0
    }


def main():
    recs = load_kqs(only_enabled=True)
    checked = 0
    for r in recs:
        if not r.get("guideline_refs"):
            continue                                   # 검증셋 없는 KQ(general 등) skip
        checked += 1
        res = validate_single_kq(r)
        print(f"\n=== KQ: {res['kq']} ===")
        
        # 0단계 출력
        flag = "✅ 일치" if res["type_match"] else "⚠️ 불일치"
        print(f"  [0단계] 유형: 사용자(결정형)={res['user_type']} · classify제안={res['suggested_type']} → {flag}")
        
        # 검증 A 출력
        for pmid, ok in res["refs_results"]:
            print(f"    {'✅' if ok else '❌'} {pmid}")
        print(f"  [검증 A] 지침근거 재현 = {res['refs_hit']}/{res['refs_total']}")

        # 검증 B 출력
        if res["landmarks_total"] > 0:
            for pmid, ok in res["landmarks_results"]:
                print(f"    {'✅' if ok else '❌'} {pmid}  (landmark)")
            print(f"  [검증 B] 랜드마크 포착 = {res['landmarks_hit']}/{res['landmarks_total']}")
            
    if checked == 0:
        print("(검증셋(guideline_refs) 있는 KQ 없음)")


if __name__ == "__main__":
    main()
