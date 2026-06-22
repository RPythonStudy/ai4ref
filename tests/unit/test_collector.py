# tests/unit/test_collector.py
"""후보 수집기(_collect_candidates) 단위테스트.

FR-001~FR-011 요구사항을 만족하는지 검증한다.
NCBI E-utilities API 네트워크 호출을 차단하기 위해 esearch와 efetch_meta를 모킹한다.
"""
from unittest.mock import patch, MagicMock
import pytest

from alert.run import _collect_candidates
from notify.base import LandmarkItem


def get_mock_features():
    return {"alert": {"enabled": True}}


def test_collect_candidates_no_kqs():
    """활성 KQ가 없을 때 빈 목록 반환 (FR-009)."""
    with patch("alert.run.load_kqs", return_value=[]):
        cands = _collect_candidates(get_mock_features(), reldate=30)
        assert cands == []


def test_collect_candidates_success():
    """정상적인 수집 흐름 검증 (FR-001, FR-002, FR-003, FR-005, FR-006, FR-007)."""
    # 1. 테스트용 KQ 설정 (T = 2018-05 이므로 pdat 필터는 2019:3000[dp])
    mock_kq = {
        "kq": "테스트 수액전략",
        "question_type": "intervention",
        "pico": {
            "P": ["abdominal surgery[tiab]"],
            "I": ["fluid therapy[tiab]"]
        },
        "guideline": {
            "name": "테스트 가이드라인",
            "pmid": "111111",
            "date": "2018-05"
        },
        "enabled": True
    }

    # 2. 모킹값 준비
    # build_query(pico) = "(fluid therapy[tiab]) AND (abdominal surgery[tiab])"
    # term_dated = "((fluid therapy[tiab]) AND (abdominal surgery[tiab])) AND (2019:3000[dp])"
    expected_pmids = ["12345", "67890"]
    mock_metas = {
        "12345": {
            "title": "Title 1",
            "abstract": "Abstract 1",
            "journal": "NEJM",
            "year": "2020"
        },
        "67890": {
            "title": "Title 2",
            "abstract": "Abstract 2",
            "journal": "JAMA",
            "year": "2021"
        }
    }

    with patch("alert.run.load_kqs", return_value=[mock_kq]) as mock_load, \
         patch("alert.run.esearch", return_value=expected_pmids) as mock_esearch, \
         patch("alert.run.efetch_meta", return_value=mock_metas) as mock_efetch:

        cands = _collect_candidates(get_mock_features(), reldate=15, limit=None)

        # load_kqs 호출 확인
        mock_load.assert_called_once_with(only_enabled=True)

        # esearch 쿼리 확인 (지침년도 T=2018 이므로 T+1=2019 적용, reldate=15 전달)
        expected_query = "((fluid therapy[tiab]) AND (abdominal surgery[tiab])) AND (2019:3000[dp])"
        mock_esearch.assert_called_once_with(expected_query, datetype="edat", reldate=15)

        # efetch 호출 확인
        mock_efetch.assert_called_once_with(expected_pmids)

        # 수집된 후보 검증
        assert len(cands) == 2
        assert all(isinstance(c, LandmarkItem) for c in cands)
        
        # 첫 번째 후보 데이터 검증
        c1 = cands[0]
        assert c1.kq == "테스트 수액전략"
        assert c1.pmid == "12345"
        assert c1.title == "Title 1"
        assert c1.journal == "NEJM"
        assert c1.year == "2020"
        assert c1.abstract == "Abstract 1"


def test_collect_candidates_limit():
    """limit 설정 시 반환되는 후보 개수 제한 및 efetch 호출 개수 제한 검증 (FR-011)."""
    mock_kq = {
        "kq": "테스트 수액전략",
        "pico": {"P": ["p"], "I": ["i"]},
        "guideline": {"name": "g", "date": "2020-01"},
        "enabled": True
    }
    
    mock_pmids = ["1", "2", "3", "4", "5"]
    mock_metas = {pmid: {"title": f"T{pmid}", "journal": "J", "year": "2021"} for pmid in mock_pmids}

    with patch("alert.run.load_kqs", return_value=[mock_kq]), \
         patch("alert.run.esearch", return_value=mock_pmids), \
         patch("alert.run.efetch_meta", return_value=mock_metas) as mock_efetch:

        # limit=2로 제한 호출
        cands = _collect_candidates(get_mock_features(), reldate=30, limit=2)

        # efetch_meta에는 limit 처리된 2개만 전달되어야 함
        mock_efetch.assert_called_once_with(["1", "2"])
        assert len(cands) == 2
        assert [c.pmid for c in cands] == ["1", "2"]
