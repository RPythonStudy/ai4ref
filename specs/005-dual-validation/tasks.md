# Tasks: 두-검증 (Dual Validation — Recall)

**Feature**: `specs/005-dual-validation/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: `validate_kq.py` 실행 시 검증 A 8/9 및 검증 B 2/2 재현 결과 콘솔 출력이 정확하게 동일해야 함.

## Phase 1: Setup

- [x] T001 [P] `src/alert/validate_spec.md` 생성 (L3 <= 200줄) - `validate_kq.py` 모듈 및 API 규격 문서화.

## Phase 2: Foundational

- [x] T002 [P] `tests/unit/test_validation.py` 단위 테스트 파일 작성 시작 (esearch 및 classify_question_type 모킹 인프라 구축).

## Phase 3: User Story 1 & 2 & 3 — 검증 엔진 리팩터링 및 API 통합

**Goal**: 중복 API 코드를 폐기하고, 검증 핵심 로직을 테스트 가능한 함수로 추출 및 이식.
**Independent Test**: `test_validation.py` 단위 테스트 및 실제 실행 리포트 결과 일치.

- [x] T003 `src/alert/validate_kq.py` 리팩터링:
  * `EUTILS` 상수 및 `esearch_count` 함수 폐기.
  * `common.pubmed.esearch` 및 `common.kq.load_kqs` 임포트 적용.
  * 개별 KQ 검증 핵심 비즈니스 로직을 `validate_single_kq(kq: dict) -> dict` 함수로 추출하고 포맷팅 및 출력 기능과 역할 분리.
- [x] T004 [P] `tests/unit/test_validation.py` 작성 완료 (모의 응답을 통한 `validate_single_kq` 채점 정밀도 검증).
- [x] T005 `uv run pytest`로 신규 테스트를 포함한 전체 테스트 통과 확인.

## Phase 4: Polish & 최종 검증

- [x] T006 `uv run python -m alert.validate_kq`를 호출하여 실제 PubMed API를 탄 리포트 출력이 기존 포맷과 완벽히 일치하는지 최종 확인.
- [x] T007 [P] 헌법 XI 줄수 점검 — `validate_kq.py` <= 300줄, `validate_spec.md` <= 200줄 확인.

## Dependencies

- Phase 1 → Phase 2 → Phase 3 → Phase 4.
