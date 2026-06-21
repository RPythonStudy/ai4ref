# query_spec.md — `common/query.py` (L3 명세)

> L3 = 모듈 단위 명세. API 시그니처·구현 상세는 코드(docstring·type hint)로 일원화한다
> (헌법 XI·XVI). 이 문서는 **무엇을·왜**(비즈니스 규칙)와 참조만 둔다.

## 책임

PICO(OR-블록)에서 PubMed 검색식을 **파생**한다. 검색식 = `(I 블록) AND (P 블록)`.
별도 query 필드를 저장하지 않고 매 호출 파생(단일 출처, 헌법 III).

- 모듈: `src/common/query.py`
- 공개 API: `build_query(pico: dict) -> str` (시그니처·인자 = 코드 docstring)
- 호출처: `alert/run.py`(증분 검색식), `alert/validate_kq.py`(recall 검증)

## 비즈니스 규칙

- **I·P 만 검색식**: I·P 항목을 ` OR ` 로 묶어 블록 → `(I) AND (P)`. 넓은 그물 = recall.
- **C·O 제외**: 대조·결과는 검색식에 넣지 않는다. 게이트가 의미 판정(헌법 II).
  검색식에 넣으면 검증 A 가 붕괴(3/9) — 재논의 금지 결정.
- **정련 = 리스트 성장**: P·I OR 리스트에 synonym 항목을 **추가**하는 것이 정련.
  점수 상승(6/9→8/9)은 블록 확장의 결과. 별도 정련 검색식 필드 없음(파생).
- **항목 원형 보존**: `"Fluid Therapy"[Mesh]` 등 토큰을 변형 없이 사용(C5).
  이스케이프·문법 검증은 검색 단계(003) 책임.
- **결정론**: 동일 PICO → 동일 검색식(항목 순서 보존, C6).

## 계약·검증

- 계약: `specs/001-kq-pico-authoring/contracts/build_query.md` (C1~C6)
- 단위테스트: `tests/unit/test_query.py` — derived_query 등가 + 블록 순서·C/O 제외·결정론
- 등가 참값: `tests/acceptance/recall_baseline.yml.derived_query` (헌법 XIII)

## 관련 FR

- FR-002 검색식 파생 · FR-003 C·O 제외 · FR-004 결정론·재현
