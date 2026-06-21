# tests/unit/test_kq.py
"""KQ 로더·검증 단위테스트 — US2(두-검증 앵커) 범위.

검증 앵커(guideline·guideline_refs·post_guideline_landmarks) 필드 보존과
검증 규칙 5(검증 A∩B 중복 PMID 경고, data-model.md)를 단언한다.
US1 기본 검증(규칙 1~4)도 회귀 가드로 함께 둔다.
"""
import warnings
from pathlib import Path

import pytest

from common.kq import (
    load_kqs,
    validate_kq_record,
    effective_design_strictness,
    is_kq_enabled,
)

_ROOT = Path(__file__).resolve().parents[2]
_KQ = _ROOT / "config" / "key_questions.yml"


def _valid_rec() -> dict:
    """검증 통과하는 최소 KQ 레코드(앵커 포함)."""
    return {
        "kq": "테스트 KQ",
        "question_type": "intervention",
        "pico": {"P": ["p[tiab]"], "I": ["i[tiab]"]},
        "guideline": {"name": "테스트 지침", "pmid": "123", "date": "2012-12"},
        "guideline_refs": ["111", "222"],
        "post_guideline_landmarks": ["333", "444"],
    }


def test_load_real_config_preserves_anchors():
    """실제 config 로드 → 앵커 필드(guideline·refs·landmarks) 보존."""
    recs = load_kqs(_KQ)
    rec = recs[0]
    assert rec["guideline"]["date"] == "2012-12"
    assert len(rec["guideline_refs"]) == 9
    assert len(rec["post_guideline_landmarks"]) == 2


def test_duplicate_pmid_across_a_and_b_warns():
    """규칙 5 — guideline_refs ∩ post_guideline_landmarks 중복 PMID 시 경고."""
    rec = _valid_rec()
    rec["post_guideline_landmarks"] = ["222", "555"]   # 222 = refs 와 중복
    with pytest.warns(UserWarning, match="222"):
        validate_kq_record(rec)


def test_no_duplicate_pmid_no_warning():
    """A∩B 중복 없으면 경고 없음(앵커 정상)."""
    rec = _valid_rec()
    with warnings.catch_warnings():
        warnings.simplefilter("error")     # 경고가 뜨면 테스트 실패
        validate_kq_record(rec)


def test_anchor_fields_optional():
    """guideline_refs·post_guideline_landmarks·guideline.date 모두 선택."""
    rec = _valid_rec()
    del rec["guideline_refs"]
    del rec["post_guideline_landmarks"]
    del rec["guideline"]["date"]
    validate_kq_record(rec)   # 예외·경고 없이 통과


def test_anchor_refs_must_be_list():
    """guideline_refs 가 리스트가 아니면 차단."""
    rec = _valid_rec()
    rec["guideline_refs"] = "111"     # 문자열 — 잘못된 형태
    with pytest.raises(ValueError, match="guideline_refs"):
        validate_kq_record(rec)


def test_guideline_date_format_validated():
    """guideline.date 가 YYYY-MM 형식이 아니면 차단."""
    rec = _valid_rec()
    rec["guideline"]["date"] = "2012/12"
    with pytest.raises(ValueError, match="date"):
        validate_kq_record(rec)


# ── US3 감시 제어 (T010, FR-007) ────────────────────────────────────────────

def test_design_strictness_defaults_to_loose():
    """design_strictness 미지정 → 기본 loose(recall 우선)."""
    rec = _valid_rec()
    assert "design_strictness" not in rec
    assert effective_design_strictness(rec) == "loose"


def test_design_strictness_explicit_preserved():
    """design_strictness 명시값은 그대로 사용."""
    rec = _valid_rec()
    rec["design_strictness"] = "strict"
    assert effective_design_strictness(rec) == "strict"


def test_enabled_defaults_to_true():
    """enabled 미지정 → 기본 활성(True)."""
    rec = _valid_rec()
    assert "enabled" not in rec
    assert is_kq_enabled(rec) is True


def test_enabled_false_disables():
    """enabled: false 만 비활성."""
    rec = _valid_rec()
    rec["enabled"] = False
    assert is_kq_enabled(rec) is False


def test_load_only_enabled_excludes_disabled(tmp_path):
    """load_kqs(only_enabled=True) 가 enabled:false KQ 를 제외."""
    import yaml
    on, off = _valid_rec(), _valid_rec()
    on["kq"], off["kq"], off["enabled"] = "on", "off", False
    p = tmp_path / "kq.yml"
    p.write_text(yaml.safe_dump({"kqs": [on, off]}), encoding="utf-8")
    assert [r["kq"] for r in load_kqs(p)] == ["on", "off"]                 # 기본: 전부
    assert [r["kq"] for r in load_kqs(p, only_enabled=True)] == ["on"]     # 비활성 제외


def test_enabled_must_be_bool():
    """enabled 가 bool 이 아니면 차단."""
    rec = _valid_rec()
    rec["enabled"] = "yes"
    with pytest.raises(ValueError, match="enabled"):
        validate_kq_record(rec)


def test_collection_must_be_str():
    """collection 이 문자열이 아니면 차단(있을 때만)."""
    rec = _valid_rec()
    rec["collection"] = ["a", "b"]
    with pytest.raises(ValueError, match="collection"):
        validate_kq_record(rec)


def test_collection_passes_through():
    """collection 문자열은 통과·보존."""
    rec = _valid_rec()
    rec["collection"] = "eras-fluid-abdominal"
    validate_kq_record(rec)   # 예외 없음
