# tests/unit/test_query.py
"""build_query 등가 단위테스트 — contracts/build_query.md C1~C6 + 헌법 XIII.

ERAS 2012 KQ 의 pico 로 build_query 를 호출한 결과가 acceptance 참값
(recall_baseline.yml.derived_query)과 문자 단위로 일치함을 단언한다.
참값이 곧 등가 가드: build_query 재작성이 기존 검색식을 그대로 재현해야 통과.
"""
from pathlib import Path

import yaml

from common.query import build_query

_ROOT = Path(__file__).resolve().parents[2]
_KQ = _ROOT / "config" / "key_questions.yml"
_BASELINE = _ROOT / "tests" / "acceptance" / "recall_baseline.yml"


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_eras_derived_query_matches_baseline():
    """ERAS KQ pico → build_query 출력 == acceptance derived_query (C3·C5·C6)."""
    pico = _load(_KQ)["kqs"][0]["pico"]
    expected = _load(_BASELINE)["key_questions"][0]["derived_query"]
    assert build_query(pico) == expected


def test_block_order_and_exclusion():
    """C1·C2·C4 — I 가 P 앞, C·O 토큰은 검색식에 등장하지 않음."""
    pico = {"P": ["a[tiab]", "b[tiab]"], "I": ["x[tiab]", "y[tiab]"],
            "C": "control text", "O": "outcome text"}
    q = build_query(pico)
    assert q == "(x[tiab] OR y[tiab]) AND (a[tiab] OR b[tiab])"
    assert "control" not in q and "outcome" not in q


def test_deterministic():
    """C6 — 동일 입력 → 동일 출력(항목 순서 보존)."""
    pico = {"P": ["p1", "p2", "p3"], "I": ["i1", "i2"]}
    assert build_query(pico) == build_query(pico)


def test_items_verbatim():
    """C5 — Mesh/따옴표 항목이 변형 없이 그대로 사용됨."""
    pico = {"P": ['"Surgical Procedures, Operative"[Mesh]'], "I": ['"Fluid Therapy"[Mesh]']}
    assert build_query(pico) == '("Fluid Therapy"[Mesh]) AND ("Surgical Procedures, Operative"[Mesh])'
