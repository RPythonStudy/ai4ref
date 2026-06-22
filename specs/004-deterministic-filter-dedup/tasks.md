# Tasks: 기계 필터 & 중복 차단 (Deterministic Filter & Dedup)

**Feature**: `specs/004-deterministic-filter-dedup/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: `SeenStore`가 정상 파일 및 손상된 파일 모두에 대해 완벽하게 기동하며, 중복 Skip 기능 작동이 보장되어야 함.

## Phase 1: Setup

- [x] T001 [P] `src/common/state_spec.md` 생성 (L3 <= 200줄) - `state.py` 모듈 및 SeenStore API 문서화.

## Phase 2: Foundational

- [x] T002 [P] `tests/unit/test_state.py` 단위 테스트 파일 작성 시작.

## Phase 3: User Story 1 & 2 — 상태 저장소 검증 및 예외 복구

**Goal**: `SeenStore`가 비어있거나 손상된 파일에서도 복구되어 기동하고, 마킹 작업이 멱등적이고 정상 작동하는지 확인.
**Independent Test**: `test_state.py` 내 모든 단언문 실행 통과.

- [x] T003 [P] `tests/unit/test_state.py` 내 테스트 케이스 구현 완료:
  - 지정 경로에 파일이 없을 때 빈 set으로 로드.
  - 깨진 줄(파싱 실패) 및 빈 행 무시 및 정상 행만 로드 (fail-soft).
  - `mark` 호출 시 파일 끝 기록 및 `_seen` 집합 반영 검증.
  - 상위 디렉토리 미존재 시 자동 생성 기능 검증.
- [x] T004 `uv run pytest`로 신규 테스트 통과 확인.

## Phase 4: Polish & 정합성 점검

- [x] T005 `git status`를 실행하여 런타임 상태 디렉토리 `state/`가 정상적으로 `.gitignore` 정책 하에 추적 제외(untracked에 보이지 않음)되는지 최종 점검.
- [x] T006 [P] 헌법 XI 줄수 점검 — `state.py` <= 300줄, `state_spec.md` <= 200줄 확인.

## Dependencies

- Phase 1 → Phase 2 → Phase 3 → Phase 4.
