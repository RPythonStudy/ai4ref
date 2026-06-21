#!/usr/bin/env python3
"""
ai4ref — kq_pico_collection.yml KQ 레코드 검증 (수정 개념 MVP)
================================================================
  [0단계] 사용자가 입력한 question_type ↔ classify_qtype 보조 제안 대조 (결정형 검증)
  [검증 A] term(검색식)이 guideline_refs(지침 근거)를 재현하나 (relative recall)
  .venv/bin/python src/alert/validate_kq.py
"""
from __future__ import annotations
import os, sys, json, time, urllib.parse, urllib.request
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from classify_qtype import classify_question_type

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch_count(term: str, mindate: str | None = None, maxdate: str | None = None) -> int:
    q = {"db": "pubmed", "term": term, "retmode": "json",
         "tool": "ai4ref-kqval", "email": "r.python.ai@gmail.com"}
    if mindate and maxdate:
        q.update({"datetype": "pdat", "mindate": mindate, "maxdate": maxdate})
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(q)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        d = json.load(r)
    time.sleep(0.34)
    return int(d["esearchresult"]["count"])

def main():
    cfg = os.path.join(os.path.dirname(__file__), "..", "..", "config", "kq_pico_collection.yml")
    recs = yaml.safe_load(open(cfg, encoding="utf-8"))["kqs"]
    checked = 0
    for r in recs:
        if not r.get("guideline_refs"):
            continue                                   # 검증셋 없는 KQ(general 등) skip
        checked += 1
        print(f"\n=== KQ: {r['kq']} ===")
        # [0단계] 결정형 유형 ↔ classify 보조 제안
        sug = classify_question_type(r["kq"])["question_type"]
        user = r.get("question_type")
        flag = "✅ 일치" if sug == user else f"⚠️ 불일치"
        print(f"  [0단계] 유형: 사용자(결정형)={user} · classify제안={sug} → {flag}")
        # 시간축 (T = 지침 제정 시점) — 있으면 A=T이전 / B=T이후 윈도우
        T = (r.get("guideline") or {}).get("date")          # "YYYY-MM" or None
        Ty = int(T[:4]) if T else None

        # [검증 A] 검색식이 지침 근거(T 이전)를 재현하나
        pmids = r["guideline_refs"]
        hit = 0
        for pmid in pmids:
            term = f"({r['term']}) AND {pmid}[uid]"
            ok = (esearch_count(term, "1900/01/01", f"{Ty}/12/31") if Ty else esearch_count(term)) >= 1
            hit += ok
            print(f"    {'✅' if ok else '❌'} {pmid}")
        print(f"  [검증 A] 지침근거 재현 = {hit}/{len(pmids)}")

        # [검증 B] alert 가 지침 이후(T 이후) 랜드마크를 포착하나
        lms = r.get("post_guideline_landmarks", [])
        if lms and Ty:
            bhit = 0
            for pmid in lms:
                ok = esearch_count(f"({r['term']}) AND {pmid}[uid]", f"{Ty + 1}/01/01", "2030/12/31") >= 1
                bhit += ok
                print(f"    {'✅' if ok else '❌'} {pmid}  (landmark)")
            print(f"  [검증 B] 랜드마크 포착 = {bhit}/{len(lms)}")
    if checked == 0:
        print("(검증셋(guideline_refs) 있는 KQ 없음)")

if __name__ == "__main__":
    main()
