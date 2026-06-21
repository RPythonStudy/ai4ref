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

import re
import warnings
from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[2]          # ai4ref/
DEFAULT_KQ_PATH = _ROOT / "config" / "key_questions.yml"

# 허용 enum (data-model.md 검증 규칙 4)
QUESTION_TYPES = ("intervention", "diagnostic", "prognostic", "predictive")
DESIGN_STRICTNESS = ("strict", "loose")

# 감시 제어 기본값 (US3, FR-007·헌법 VII — recall 우선)
DEFAULT_DESIGN_STRICTNESS = "loose"   # design_strictness 미지정 → loose
DEFAULT_ENABLED = True                # enabled 미지정 → 활성(명시 false 만 비활성)

# 저장 금지 필드 — 검색식은 PICO 에서 파생(헌법 III)
FORBIDDEN_FIELDS = ("query",)

# 검증 앵커 — PMID 리스트 필드(검증 A·B). 둘 다 선택(없으면 채점/포착 제외)
ANCHOR_LIST_FIELDS = ("guideline_refs", "post_guideline_landmarks")

# guideline.date(=T, 제정 시점) 형식 — YYYY-MM
_DATE_RE = re.compile(r"^\d{4}-\d{2}$")


def validate_kq_record(rec: dict) -> None:
    """KQ 레코드 1건을 검증(규칙 1~5). 위반 시 ValueError(명확한 메시지).

    1. `kq`·`question_type`·`pico`·`guideline.name` 필수.
    2. `pico.P`·`pico.I` 가 비어있지 않은 리스트.
    3. 검색식 필드(`query` 등) 부재(레코드·pico 양쪽).
    4. `question_type` ∈ enum, `design_strictness` ∈ {strict, loose} 또는 미지정.
       감시 제어(US3): `enabled` = bool, `collection` = str (있을 때만, 모두 선택).
    5. (경고) `guideline_refs` ∩ `post_guideline_landmarks` 중복 PMID = 정의 오류 표면화.

    검증 앵커(`guideline.date`·`guideline_refs`·`post_guideline_landmarks`)는 모두
    선택 — 있으면 형태(YYYY-MM·PMID 리스트)를 확인하고, 검증 A/B 시간축 분리에 쓰인다.
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

    # 감시 제어 속성(US3, FR-007) — 모두 선택. enabled 는 bool, collection 은 str(있을 때).
    enabled = rec.get("enabled")
    if enabled is not None and not isinstance(enabled, bool):
        raise ValueError(
            f"KQ '{label}' 검증 실패: enabled={enabled!r} 는 bool 이어야 함(true|false 또는 미지정)")
    collection = rec.get("collection")
    if collection is not None and not isinstance(collection, str):
        raise ValueError(
            f"KQ '{label}' 검증 실패: collection={collection!r} 는 문자열이어야 함(컬렉션명)")

    # 검증 앵커 — guideline.date(T) 형식(선택)
    gdate = guideline.get("date")
    if gdate is not None and not _DATE_RE.match(str(gdate)):
        raise ValueError(
            f"KQ '{label}' 검증 실패: guideline.date='{gdate}' 형식 오류 (YYYY-MM 필요)")

    # 검증 앵커 — guideline_refs·post_guideline_landmarks 형태(선택, PMID 리스트)
    for fld in ANCHOR_LIST_FIELDS:
        val = rec.get(fld)
        if val is None:
            continue
        if not isinstance(val, (list, tuple)):
            raise ValueError(
                f"KQ '{label}' 검증 실패: '{fld}' 는 PMID 리스트여야 함(현재 {type(val).__name__})")

    # 규칙 5 — 검증 A∩B 중복 PMID 경고(차단 아님, 정의 오류 표면화)
    refs = {str(x) for x in (rec.get("guideline_refs") or [])}
    landmarks = {str(x) for x in (rec.get("post_guideline_landmarks") or [])}
    dup = sorted(refs & landmarks)
    if dup:
        warnings.warn(
            f"KQ '{label}': PMID {', '.join(dup)} 가 guideline_refs(A)·"
            f"post_guideline_landmarks(B) 양쪽에 존재 — T 기준 시간축 모순(정의 확인 필요)",
            UserWarning, stacklevel=2)


def effective_design_strictness(rec: dict) -> str:
    """KQ 의 유효 `design_strictness` — 미지정이면 기본 `loose`(US3, FR-007·헌법 VII).

    게이트 2단계의 랜드마크 판정 엄격도. 미지정/빈값 → recall 우선의 loose.
    """
    ds = rec.get("design_strictness")
    return ds if ds in DESIGN_STRICTNESS else DEFAULT_DESIGN_STRICTNESS


def is_kq_enabled(rec: dict) -> bool:
    """감시 토글 — `enabled` 미지정이면 기본 활성(True). 명시 `false` 만 비활성(US3).

    감시 대상에서 KQ 를 끄는 단일 기준. `load_kqs(only_enabled=True)` 가 이 함수로 거른다.
    """
    return bool(rec.get("enabled", DEFAULT_ENABLED))


def load_kqs(path: str | Path | None = None, *, only_enabled: bool = False) -> list[dict]:
    """key_questions.yml 의 `kqs` 리스트를 적재·검증해 반환.

    각 레코드에 validate_kq_record 를 적용 — 하나라도 위반하면 ValueError(차단).
    검증 통과 레코드만 반환(결정론, 파일 순서 보존).

    `only_enabled=True` 면 감시 비활성(`enabled: false`) KQ 를 제외한다(US3). 검증은
    전 레코드에 먼저 적용하므로 비활성 KQ 의 스키마 오류도 똑같이 조기 차단된다.
    """
    p = Path(path) if path else DEFAULT_KQ_PATH
    data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    records = data.get("kqs", []) or []
    for i, rec in enumerate(records):
        try:
            validate_kq_record(rec)
        except ValueError as e:
            raise ValueError(f"[{p.name} #{i}] {e}") from e
    if only_enabled:
        records = [r for r in records if is_kq_enabled(r)]
    return list(records)
