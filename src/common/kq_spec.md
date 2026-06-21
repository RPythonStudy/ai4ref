# kq_spec.md — `common/kq.py` (L3 명세)

> L3 = 모듈 단위 명세. API 시그니처·구현 상세는 코드(docstring·type hint)로 일원화한다
> (헌법 XI·XVI). 이 문서는 **무엇을·왜**(검증 규칙·단일 출처)와 참조만 둔다.

## 책임

`config/key_questions.yml` 의 KQ 레코드를 **적재·검증**한다. 잘못된 KQ 가
파이프라인에 진입하기 전 조기 차단(결정론, 명확한 오류).

- 모듈: `src/common/kq.py`
- 공개 API (시그니처 = 코드 docstring):
  - `load_kqs(path?) -> list[dict]` — `kqs` 리스트 적재·검증·반환
  - `validate_kq_record(rec) -> None | raises` — 레코드 1건 검증
- 권위 원본: `config/key_questions.yml` (단일 출처). 검색식 파생 = [[query]].

## 검증 규칙 (data-model.md 1~4, US1 범위)

1. **필수 필드**: `kq`·`question_type`·`pico`·`guideline.name`.
2. **P·I 비어있지 않음**: `pico.P`·`pico.I` 가 비어있지 않은 OR 리스트.
3. **검색식 필드 금지**: `query` 등 파생 검색식을 레코드/pico 에 저장 금지(헌법 III).
4. **enum**: `question_type` ∈ {intervention, diagnostic, prognostic, predictive},
   `design_strictness` ∈ {strict, loose} 또는 미지정.

검증 실패 = **차단**: `load_kqs` 가 위반 레코드 인덱스·사유를 담아 `ValueError`.

> 증분 규칙(US2·US3): 검증 앵커 중복 경고(규칙 5)·`design_strictness` 기본 `loose`·
> `enabled` 기본 처리는 후속 태스크(T008·T010)에서 본 명세에 추가한다.

## 단일 출처 (헌법 XVI)

- 검색식은 저장하지 않고 `build_query(pico)` 로 파생. 스키마에 `query` 필드 부재.
- KQ 레코드의 권위는 YAML 1곳. 코드는 그 스키마를 검증할 뿐 복제하지 않는다.

## 계약·검증

- 계약: `specs/001-kq-pico-authoring/contracts/kq_schema.md`
- 데이터 모델: `specs/001-kq-pico-authoring/data-model.md`

## 관련 FR

- FR-001 KQ 정의 · FR-002·FR-003 PICO/파생 입력 · FR-007 감시 제어 속성(검증) ·
  FR-009 config 확장만으로 KQ 추가(코드 변경 0)
