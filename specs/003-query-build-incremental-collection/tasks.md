# Tasks: 검색식 빌드 & 증분 수집 (Query Build & Incremental Collection)

**Feature**: `specs/003-query-build-incremental-collection/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: `validate_kq` 수행 시 검증 A 8/9, B 2/2가 정확히 달성되며, `--collect-only` 실행 결과가 동일해야 함.

## Phase 1: Setup

- [x] T001 [P] `src/common/pubmed_spec.md` 생성 (L3 <= 200줄) - `pubmed.py` 모듈 및 API 규격 문서화.

## Phase 2: Foundational

- [x] T002 [P] `tests/unit/test_collector.py` 단위 테스트 파일 작성 시작 (NCBI E-utilities API 모킹 인프라 구축).

## Phase 3: User Story 1 — KQ 로더 단일화 & context map 정비

**Goal**: `src/alert/run.py` 내의 중복 로더를 거두고 `common.kq.load_kqs`로 단일화하여 검증 적용 및 조기 오류 차단.
**Independent Test**: `validate_kq_record`에 실패하는 가짜 KQ가 yml에 있을 때 `run.py` 가 초기에 에러를 내며 죽는지 검증.

- [x] T003 `src/alert/run.py`에서 로컬 `_load_kqs()` 제거 및 `common.kq.load_kqs(only_enabled=True)` 임포트 및 활용.
- [x] T004 `src/alert/run.py`의 `run()` 내 `ctx_map` 빌드 시 `effective_design_strictness(kq)` 적용하도록 수정.

## Phase 4: User Story 2 — 후보 수집기 테스트 및 등가성 검증

**Goal**: 모킹된 NCBI API를 이용하여 `_collect_candidates`가 요구사항 사양서(FR-001~011)에 명시된 대로 쿼리를 결합하고 결과를 반환하는지 단언.
**Independent Test**: `test_collector.py` 실행 통과.

- [x] T005 [P] `tests/unit/test_collector.py` 작성 완료:
  - 활성 KQ가 없을 때 빈 목록 반환 여부.
  - 검색식 빌드 및 날짜 필터(`pdat` 조건) 결합이 사양서 규칙대로 올바르게 일치하는지 검증.
  - `limit` 인자 제공 시 후보 개수가 올바르게 제한되는지 검증.
  - `reldate` 인자가 `esearch`에 정상 전달되는지 검증.

## Phase 5: Polish & 등가 검증

- [x] T006 `uv run pytest`로 기존 테스트 및 신규 단위 테스트 통과 확인.
- [x] T007 `uv run python -m alert.validate_kq`를 실행하여 검증 A 8/9, B 2/2가 완벽히 재현되는지 최종 확인.
- [x] T008 [P] 헌법 XI 줄수 점검 — `pubmed.py` <= 300줄, `pubmed_spec.md` <= 200줄 확인.

## Dependencies

- Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5.
- T002(test setup)가 T005(test completion) 이전에 수행되어야 함.

## MVP Scope

**US1 (Phase 3) & US2 (Phase 4)**: KQ 로더의 올바른 통합과 `_collect_candidates` 기능에 대한 모의(mock) 단위 테스트 확보가 완성도의 핵심입니다.
