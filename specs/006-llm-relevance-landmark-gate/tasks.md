# Tasks: LLM 게이트 (Relevance & Landmark Gate — Precision)

**Feature**: `specs/006-llm-relevance-landmark-gate/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: `tests/acceptance/gate_baseline.yml` 재현 — `landmark_true` 2/2(is_relevant&
is_landmark=true) · `reject_false_positive` 2/2(is_relevant=false). 기존 `alert.run` 게이트
동작(judge 호출·verdict 소비) 보존. 코드 표면 = `src/alert/llm_gate.py` + `run.py` 호출부.

## Phase 1: Setup

- [x] T001 [P] `tests/unit/test_gate.py` 단위테스트 파일 생성 — 백엔드 모의(monkeypatch
  `subprocess.run`/stub) 인프라 골격 구축. `from alert.llm_gate import judge, _extract_json`.

## Phase 2: Foundational — 출력 스키마 정합 (blocking)

**Goal**: 게이트 출력 키를 acceptance 스키마(`is_relevant`)로 일원화. 모든 US 의 전제.
**근거**: 현행 `judge()`는 `{"relevant",...}` 반환, baseline·§006 은 `is_relevant`. `run.py`는
`verdict.get("is_landmark")`만 소비(grep 확인됨) → 변경은 `llm_gate.py` 내부로 국소.

- [x] T002 `src/alert/llm_gate.py` 출력 키 `relevant` → `is_relevant` 일괄 정합 — `_stub`·
  `_claude_cli`(정상/실패 경로 모두)·docstring·PROMPT JSON 예시·`data.get(...)` 추출부.
- [x] T003 `src/alert/run.py` 게이트 소비부 점검 — `verdict.get("is_landmark")` 유지 확인, 향후
  `is_relevant` 소비 시 정합. 동작 변화 없음 단언.

## Phase 3: User Story 1 — 관련성 게이트 (P1, MVP)

**Goal**: ① 대상(P) AND 중재(I) 동시매칭으로 토픽-인접 FP 차단. ② 는 ①=true 일 때만 평가.
**Independent Test**: `gate_baseline.yml` 의 `reject_false_positive` 2건(신장이식·로봇수술)이
is_relevant=false 로 거부되고, ① reject 시 ②(랜드마크)가 평가되지 않는다.

- [x] T004 [US1] `src/alert/llm_gate.py` 2단계 순서 보존 확인 — ① 관련성 판정, ①=false 면
  is_landmark 평가 없이 즉시 탈락(FR-001·002·SC-005). PROMPT 의 P∧I 동시매칭 규칙 유지.
- [x] T005 [P] [US1] `tests/unit/test_gate.py` — 모의 백엔드로 ① reject(is_relevant=false) 시
  ②(is_landmark) 미평가/false 단언 + 대상불일치·중재불일치 케이스(FR-002).

## Phase 4: User Story 2 — 랜드마크 + 설계 엄격도 (P2)

**Goal**: relevant 일 때 practice-changing 판정. strict=기대설계 요구 / loose=경계선 포함.
**Independent Test**: `landmark_true` 2건(RELIEF·OPTIMISE) is_landmark=true 재현(acceptance).
strict/loose 컨텍스트가 PROMPT 에 정확히 주입되는지 단위 검증.
**Depends on**: US1(2단계 순서).

- [x] T006 [US2] `src/alert/llm_gate.py` strict/loose·`DESIGN_LABEL`(question_type→기대설계)
  컨텍스트 주입 보존 확인 + PROMPT 포맷 인자(C·O·guideline·design·strictness) 정합(FR-003·004·005).
- [x] T007 [P] [US2] `tests/unit/test_gate.py` — `_claude_cli` 의 PROMPT 포맷 인자 구성·초록
  3000자 절단·question_type 미지정 시 기본 설계라벨 대체를 모의로 단언(FR-005·011).

## Phase 5: User Story 3 — 교체형 백엔드 & fail-soft (P3)

**Goal**: 백엔드 교체(stub/claude_cli/api)·무비용 로컬·실패 시 보수적 보류·래퍼 견고 파싱.
**Independent Test**: 엔진 호출 실패/JSON 깨짐 시 보류(is_relevant=false·is_landmark=false)·
로그·계속. 래퍼(`{"result":"...json..."}`) 파싱 성공. api 백엔드 선택 시 stub 안전강등.
**Depends on**: 없음(독립 — 백엔드 디스패치 레이어).

- [x] T008 [US3] `src/alert/llm_gate.py` 백엔드 디스패치(`judge` 분기)·`_stub`·api 미구현 경고
  강등·fail-soft(except→보수적 보류) 보존 확인(FR-007·008·009).
- [x] T009 [P] [US3] `tests/unit/test_gate.py` — `subprocess.run` 모의로 (a) is_error 래퍼 →
  RuntimeError→보류, (b) 비JSON 출력→보류, (c) `_extract_json` 래퍼/중괄호추출 성공 단언(FR-009·010).

## Phase 6: Polish & 등가 검증

- [x] T010 `src/alert/gate_spec.md` 작성 (L3 ≤200줄) — judge 계약(`{is_relevant,is_landmark,
  reason}`)·2단계·백엔드·fail-soft·FR006~012 참조 + classify_qtype 보조(권위=사용자 입력) 명시.
- [x] T011 `uv run python -m pytest -q` 전체 통과 확인(기존 28 + 신규 test_gate).
- [x] T012 [P] 헌법 XI 줄수 점검 — `llm_gate.py` ≤300줄, `gate_spec.md` ≤200줄.
- [x] T013 등가 확인 — `gate_baseline.yml` 의 4 케이스(landmark 2·reject 2)를 실 `claude_cli`
  또는 문서화된 수동 점검으로 대조(LLM 비결정성은 acceptance 참값 책임, 단위테스트와 분리).

## Dependencies

- Phase 1 → Phase 2(정합, blocking) → Phase 3(US1) → Phase 4(US2, US1 의존) · Phase 5(US3 독립) → Phase 6.
- T002(키 정합)는 US 구현·테스트 전 선행(전 출력 경로 영향).
- US2 는 US1(2단계 순서) 위에 선다. US3 는 독립(백엔드 레이어).

## Parallel Opportunities

- T001 [P] (setup, 독립).
- T005·T007·T009 [P] (각 US 테스트, 서로 다른 단언 — 단 같은 파일이면 순차 편집).
- T012 [P] (줄수 점검, 독립).

## MVP Scope

**US1(Phase 3) + Phase 2 정합** = MVP — 관련성 게이트로 토픽-인접 FP 차단(precision 1차 관문).
이것만으로 §007 알림이 잡음 없이 동작한다. US2(랜드마크 정밀)·US3(백엔드·견고성)은 증분.
