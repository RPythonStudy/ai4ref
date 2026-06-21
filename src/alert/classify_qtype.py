#!/usr/bin/env python3
"""
ai4ref 0단계 — KQ 질문유형 판별 (heuristic STUB)
=================================================
주제(자연어) → 유형: intervention / diagnostic / prognostic / predictive.
유형이 PICO 구조·이상적 연구설계·검색 설계필터를 결정 (Vella&Yao 2026, J Clin Epidemiol).
  .venv/bin/python src/alert/classify_qtype.py "복부수술 수액요법 효과"
  .venv/bin/python src/alert/classify_qtype.py            # 4유형 데모

설계: 검색·선별 5단계 *앞*의 0단계. 유형이 정해져야 question_type·DESIGN_HEDGE·PICO 틀이 결정됨.

역할 = *보조*. 권위는 **사용자가 설정파일(search_collections.json)에 직접 입력한
question_type (결정형)**. 0단계 정확이 전체 성능을 좌우하고 전문가에겐 유형 판별이 자명하므로,
사람이 결정형으로 박는 게 최선(재현가능·garbage-in 방지). 이 모듈은 사용자가 *빠뜨렸거나 모호할 때*
제안·검증(입력값과 대조해 경고)만 한다. LLM 은 "의미 판단이 꼭 필요한 4단계 선별"에만.
"""
from __future__ import annotations
import sys, json

# 유형별 신호어 (heuristic). 검사형(진단/예후/예측)이 더 구체적 → 신호 있으면 우선.
SIGNALS = {
    "diagnostic": ["진단", "정확도", "민감도", "특이도", "검출", "감별", "선별검사",
                   "diagnos", "sensitivity", "specificity", "accuracy", "index test",
                   "reference standard", "roc", "auc"],
    "prognostic": ["예후", "생존", "재발", "진행", "경과", "사망위험", "위험인자",
                   "prognos", "survival", "recurrence", "mortality", "cohort", "follow-up"],
    "predictive": ["치료반응", "반응예측", "반응 예측", "예측마커", "예측 마커", "예측인자",
                   "바이오마커", "동반진단", "서브그룹",
                   "predictive marker", "predictive biomarker", "treatment response",
                   "effect modif", "companion diagnostic", "subgroup"],
    "intervention": ["치료", "중재", "효과", "효능", "수술", "약물", "요법", "투여", "비교",
                     "vs", "versus", "therapy", "treatment", "intervention", "efficacy",
                     "rct", "randomized", "무작위"],
}

# 유형 → PubMed 설계 필터 (검색식의 설계 게이트)
DESIGN_HEDGE = {
    "intervention": "randomized controlled trial[pt]",
    "diagnostic":   '(sensitivity[tiab] OR specificity[tiab] OR "diagnostic accuracy"[tiab] OR "predictive value"[tiab])',
    "prognostic":   '(prognos*[tiab] OR cohort[tiab] OR "follow-up studies"[Mesh] OR survival[tiab])',
    "predictive":   '(subgroup[tiab] OR "effect modif*"[tiab] OR "predictive marker"[tiab] OR "treatment interaction"[tiab])',
}

# 유형 → PICO 의 I/C/O 자리 의미 (Vella&Yao)
PICO_FRAME = {
    "intervention": {"I": "중재(치료)",      "C": "대조 치료",        "O": "임상 결과"},
    "diagnostic":   {"I": "Index test",     "C": "Reference standard", "O": "정확도(민감도·특이도)"},
    "prognostic":   {"I": "예후 인자",       "C": "(인자 유무)",       "O": "경과(생존·재발)"},
    "predictive":   {"I": "예측 마커",       "C": "(마커 유무+치료)",  "O": "차별적 치료반응"},
}

def classify_question_type(topic: str) -> dict:
    t = topic.lower()
    scores = {k: sum(1 for kw in kws if kw.lower() in t) for k, kws in SIGNALS.items()}
    test_types = ["diagnostic", "prognostic", "predictive"]
    test_hit = {k: scores[k] for k in test_types if scores[k] > 0}
    if test_hit:                                   # 검사형 신호 우선 (더 구체적)
        best = max(test_hit, key=test_hit.get)
    elif scores["intervention"] > 0:
        best = "intervention"
    else:
        best = "intervention"                      # 기본값 (가장 흔함)
    return {
        "topic": topic,
        "question_type": best,
        "pico_frame": PICO_FRAME[best],
        "design_filter": DESIGN_HEDGE[best],
        "scores": scores,
        "method": "heuristic-stub",                # TODO: Claude 분류로 교체
    }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("=== 0단계 유형 판별 데모 (heuristic STUB) ===")
        for d in [
            "복부수술 환자의 제한적 vs 자유로운 수액요법 효과",
            "갑상선결절 PET 검사의 진단 정확도(민감도·특이도)",
            "대장암 환자에서 종양표지자의 예후 예측(생존)",
            "유방암 HER2 검사로 trastuzumab 치료반응 예측",
        ]:
            r = classify_question_type(d)
            print(f"\n  주제: {d}")
            print(f"   → 유형: {r['question_type']}  (I={r['pico_frame']['I']}, O={r['pico_frame']['O']})")
            print(f"   → 설계필터: {r['design_filter']}")
    else:
        print(json.dumps(classify_question_type(sys.argv[1]), ensure_ascii=False, indent=2))
