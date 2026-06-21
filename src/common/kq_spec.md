# kq_spec.md — `common/kq.py` (L3 명세)

> L3 = 모듈 단위 명세. API 시그니처·구현 상세는 코드(docstring·type hint)로 일원화한다
> (헌법 XI·XVI). 이 문서는 **무엇을·왜**(검증 규칙·단일 출처)와 참조만 둔다.

## 책임

`config/key_questions.yml` 의 KQ 레코드를 **적재·검증**한다. 잘못된 KQ 가
파이프라인에 진입하기 전 조기 차단(결정론, 명확한 오류).

- 모듈: `src/common/kq.py`
- 공개 API (시그니처 = 코드 docstring):
  - `load_kqs(path?, only_enabled=False) -> list[dict]` — `kqs` 리스트 적재·검증·반환
  - `validate_kq_record(rec) -> None | raises` — 레코드 1건 검증
  - `effective_design_strictness(rec) -> str` · `is_kq_enabled(rec) -> bool` — 제어 기본값(US3)
- 권위 원본: `config/key_questions.yml` (단일 출처). 검색식 파생 = [[query]].

## 검증 규칙 (data-model.md 1~4, US1 범위)

1. **필수 필드**: `kq`·`question_type`·`pico`·`guideline.name`.
2. **P·I 비어있지 않음**: `pico.P`·`pico.I` 가 비어있지 않은 OR 리스트.
3. **검색식 필드 금지**: `query` 등 파생 검색식을 레코드/pico 에 저장 금지(헌법 III).
4. **enum**: `question_type` ∈ {intervention, diagnostic, prognostic, predictive},
   `design_strictness` ∈ {strict, loose} 또는 미지정.

검증 실패 = **차단**: `load_kqs` 가 위반 레코드 인덱스·사유를 담아 `ValueError`.

## 검증 앵커 규칙 (data-model.md 5, US2 범위, FR-006)

두-검증(검색 recall 참값)의 시간축 분리 앵커를 보유·검증한다 — 모두 **선택** 필드.

5. **A∩B 중복 PMID 경고**: 동일 PMID 가 `guideline_refs`(검증 A, T 이전 근거)와
   `post_guideline_landmarks`(검증 B, T 이후 랜드마크) 양쪽에 있으면 T 기준 시간축이
   모순 — `UserWarning` 으로 표면화(**차단 아님**, 정의 오류 점검 유도).
- **형태 검증(있을 때만)**: `guideline.date` = `YYYY-MM`(=T), `guideline_refs`·
  `post_guideline_landmarks` = PMID 리스트. 형태 위반은 차단(`ValueError`).
- 앵커는 검색식·게이트에 직접 쓰이지 않음 — 검증 A/B(`validate_kq`)의 참값일 뿐.

## 감시 제어 속성 (US3, FR-007)

KQ 별 감시 제어 속성 — 모두 **선택**, 기본값은 헌법 VII(recall 우선)을 따른다.

- **`design_strictness`** ∈ {strict, loose}, 미지정 → 기본 `loose`. 유효값은
  `effective_design_strictness(rec)` 로 일원화(게이트 2단계 랜드마크 판정 엄격도).
- **`enabled`** = bool, 미지정 → 기본 `True`(활성). 명시 `false` 만 감시 비활성.
  `is_kq_enabled(rec)` 로 판정하고 `load_kqs(only_enabled=True)` 가 비활성 KQ 를 제외.
- **`collection`** = str(Zotero 보관 컬렉션명), 통과·보존(스키마는 형태만 확인).
- **`question_type`** ∈ enum — 검증은 규칙 4 와 동일(필수 필드, US1).

형태 위반(`enabled` 비-bool, `collection` 비-str)은 차단(`ValueError`). 기본값은
**저장하지 않고** 조회 시점에 적용(단일 출처 — YAML 에 명시된 값만 권위).

## 단일 출처 (헌법 XVI)

- 검색식은 저장하지 않고 `build_query(pico)` 로 파생. 스키마에 `query` 필드 부재.
- KQ 레코드의 권위는 YAML 1곳. 코드는 그 스키마를 검증할 뿐 복제하지 않는다.

## 계약·검증

- 계약: `specs/001-kq-pico-authoring/contracts/kq_schema.md`
- 데이터 모델: `specs/001-kq-pico-authoring/data-model.md`

## 관련 FR

- FR-001 KQ 정의 · FR-002·FR-003 PICO/파생 입력 · FR-006 두-검증 앵커 보유·정합 ·
  FR-007 감시 제어 속성(검증) · FR-009 config 확장만으로 KQ 추가(코드 변경 0)
