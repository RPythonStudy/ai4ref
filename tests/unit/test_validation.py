# tests/unit/test_validation.py
"""검증기 핵심 로직(validate_single_kq) 단위테스트.

classify_question_type 및 esearch 모킹을 적용하여 네트워크 호출 없이 채점 로직을 검증합니다.
"""
from unittest.mock import patch
import pytest

from alert.validate_kq import validate_single_kq


def test_validate_single_kq_success():
    """모든 매칭이 성공하는 경우 검증 A/B 및 유형 교차 검사 판정 단언."""
    mock_kq = {
        "kq": "주요 복부수술 환자 수액",
        "question_type": "intervention",
        "pico": {"P": ["p"], "I": ["i"]},
        "guideline": {"name": "g", "date": "2018-05"},
        "guideline_refs": ["111", "222"],
        "post_guideline_landmarks": ["333"]
    }

    # classify_question_type이 사용자 지정과 일치하도록 모킹
    mock_classify = {"question_type": "intervention"}
    # 모든 esearch 질의가 매칭 성공(PMID가 반환값에 포함되도록)하도록 모킹
    # esearch(term, ...) 호출 시 term에 포함된 PMID를 리턴하여 매칭 성공 처리
    def mock_esearch(term, **kwargs):
        # term: "((i) AND (p)) AND 111[uid]" 형태이므로 PMID를 추출하여 리스트로 반환
        for pmid in ["111", "222", "333"]:
            if f"{pmid}[uid]" in term:
                return [pmid]
        return []

    with patch("alert.validate_kq.classify_question_type", return_value=mock_classify) as mock_cl, \
         patch("alert.validate_kq.esearch", side_effect=mock_esearch) as mock_es:

        res = validate_single_kq(mock_kq)
        
        # 0단계 유형 일치 확인
        assert res["type_match"] is True
        assert res["user_type"] == "intervention"
        assert res["suggested_type"] == "intervention"
        
        # 검증 A 채점 확인 (2/2)
        assert res["refs_hit"] == 2
        assert res["refs_total"] == 2
        assert res["refs_results"] == [("111", True), ("222", True)]
        
        # 검증 B 채점 확인 (1/1)
        assert res["landmarks_hit"] == 1
        assert res["landmarks_total"] == 1
        assert res["landmarks_results"] == [("333", True)]


def test_validate_single_kq_partial_match():
    """일부 매칭 실패 및 유형 불일치 시의 리포팅 검증."""
    mock_kq = {
        "kq": "주요 복부수술 환자 수액",
        "question_type": "intervention",
        "pico": {"P": ["p"], "I": ["i"]},
        "guideline": {"name": "g", "date": "2018-05"},
        "guideline_refs": ["111", "222"],
        "post_guideline_landmarks": ["333"]
    }

    mock_classify = {"question_type": "diagnostic"}  # 불일치 유도
    
    # 111은 성공, 222와 333은 실패 처리
    def mock_esearch(term, **kwargs):
        if "111[uid]" in term:
            return ["111"]
        return []

    with patch("alert.validate_kq.classify_question_type", return_value=mock_classify), \
         patch("alert.validate_kq.esearch", side_effect=mock_esearch) as mock_es:

        res = validate_single_kq(mock_kq)
        
        # 0단계 유형 불일치 확인
        assert res["type_match"] is False
        assert res["user_type"] == "intervention"
        assert res["suggested_type"] == "diagnostic"
        
        # 검증 A 채점 확인 (1/2)
        assert res["refs_hit"] == 1
        assert res["refs_total"] == 2
        assert res["refs_results"] == [("111", True), ("222", False)]
        
        # 검증 B 채점 확인 (0/1)
        assert res["landmarks_hit"] == 0
        assert res["landmarks_total"] == 1
        assert res["landmarks_results"] == [("333", False)]
