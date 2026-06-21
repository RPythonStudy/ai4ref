# src/common/kq.py
"""KQ 레코드 적재·유효성 검증 (결정론, stdlib + PyYAML).

권위 원본 = `config/key_questions.yml` 의 `kqs:` 리스트. 각 원소 = 1 KQ 레코드.
로더는 검증 규칙(data-model.md 1~4)을 적용해, 잘못된 KQ 가 파이프라인에
진입하기 전에 조기 차단한다(명확한 오류). 검색식은 저장하지 않고 build_query 로
파생하므로 `query` 필드는 금지(헌법 III).

계약: specs/001-kq-pico-authoring/contracts/kq_schema.md.
검색식 파생은 [[query]] (`common.query.build_query`).
"""
from __future__ import annotations

from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[2]          # ai4ref/
DEFAULT_KQ_PATH = _ROOT / "config" / "key_questions.yml"

# 허용 enum (data-model.md 검증 규칙 4)
QUESTION_TYPES = ("intervention", "diagnostic", "prognostic", "predictive")
DESIGN_STRICTNESS = ("strict", "loose")

# 저장 금지 필드 — 검색식은 PICO 에서 파생(헌법 III)
FORBIDDEN_FIELDS = ("query",)


def validate_kq_record(rec: dict) -> None:
    """KQ 레코드 1건을 검증(규칙 1~4). 위반 시 ValueError(명확한 메시지).

    1. `kq`·`question_type`·`pico`·`guideline.name` 필수.
    2. `pico.P`·`pico.I` 가 비어있지 않은 리스트.
    3. 검색식 필드(`query` 등) 부재(레코드·pico 양쪽).
    4. `question_type` ∈ enum, `design_strictness` ∈ {strict, loose} 또는 미지정.
    """
    if not isinstance(rec, dict):
        raise ValueError(f"KQ 레코드가 dict 가 아님: {type(rec).__name__}")
    label = rec.get("kq", "<이름 없음>")

    # 규칙 1 — 필수 필드
    if not rec.get("kq"):
        raise ValueError("KQ 검증 실패: 'kq'(주제) 필수")
    if not rec.get("question_type"):
        raise ValueError(f"KQ '{label}' 검증 실패: 'question_type' 필수")
    pico = rec.get("pico")
    if not isinstance(pico, dict):
        raise ValueError(f"KQ '{label}' 검증 실패: 'pico' 필수(dict)")
    guideline = rec.get("guideline")
    if not isinstance(guideline, dict) or not guideline.get("name"):
        raise ValueError(f"KQ '{label}' 검증 실패: 'guideline.name' 필수")

    # 규칙 2 — P·I 비어있지 않은 리스트
    for blk in ("P", "I"):
        items = pico.get(blk)
        if not isinstance(items, (list, tuple)) or not [t for t in items if str(t).strip()]:
            raise ValueError(f"KQ '{label}' 검증 실패: pico.{blk} 가 비어있음(OR 리스트 필요)")

    # 규칙 3 — 검색식 필드 금지
    for fld in FORBIDDEN_FIELDS:
        if fld in rec or fld in pico:
            raise ValueError(
                f"KQ '{label}' 검증 실패: '{fld}' 필드 금지 — 검색식은 build_query 로 파생(헌법 III)")

    # 규칙 4 — enum
    qt = rec.get("question_type")
    if qt not in QUESTION_TYPES:
        raise ValueError(
            f"KQ '{label}' 검증 실패: question_type='{qt}' 미허용 (허용: {', '.join(QUESTION_TYPES)})")
    ds = rec.get("design_strictness")
    if ds is not None and ds not in DESIGN_STRICTNESS:
        raise ValueError(
            f"KQ '{label}' 검증 실패: design_strictness='{ds}' 미허용 (허용: {', '.join(DESIGN_STRICTNESS)} 또는 미지정)")


def load_kqs(path: str | Path | None = None) -> list[dict]:
    """key_questions.yml 의 `kqs` 리스트를 적재·검증해 반환.

    각 레코드에 validate_kq_record 를 적용 — 하나라도 위반하면 ValueError(차단).
    검증 통과 레코드만 반환(결정론, 파일 순서 보존).
    """
    p = Path(path) if path else DEFAULT_KQ_PATH
    data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    records = data.get("kqs", []) or []
    for i, rec in enumerate(records):
        try:
            validate_kq_record(rec)
        except ValueError as e:
            raise ValueError(f"[{p.name} #{i}] {e}") from e
    return list(records)
