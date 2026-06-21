# Tasks: KQ·PICO 정의 (Key Question & PICO Authoring)

**Feature**: `specs/001-kq-pico-authoring/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: build_query 출력 == `tests/acceptance/recall_baseline.yml` 의 `derived_query`
(단위테스트) + recall A 8/9·B 2/2(`validate_kq`). 코드 표면 = `common/query.py` + `common/kq.py`.

## Phase 1: Setup

- [x] T001 [P] `tests/unit/` 디렉토리 생성 + `tests/unit/__init__.py`(빈 파일) 추가

## Phase 2: Foundational

(추가 인프라 없음 — uv·src 레이아웃·acceptance 참값 기존 완비. 다음 단계로 진행)

## Phase 3: User Story 1 — KQ·PICO 정의 + 검색식 파생 (P1, MVP)

**Goal**: KQ 를 PICO(OR-블록)로 정의하고 검색식을 파생한다.
**Independent Test**: ERAS KQ 1건 로드 → `build_query(pico)` 출력이 acceptance `derived_query`
와 문자 일치(`quickstart.md` §1·§4).

- [x] T002 [P] [US1] build_query 단위테스트 작성 `tests/unit/test_query.py` — `key_questions.yml`
  의 ERAS KQ `pico` 로 `build_query` 호출 → `recall_baseline.yml.derived_query` 와 일치 단언
  (contracts/build_query.md C1~C6)
- [x] T003 [US1] `src/common/query.py` 의 `build_query(pico)` 재작성 — `(I OR …) AND (P OR …)`,
  C·O 제외, 항목 원형 보존, 결정론 (contracts/build_query.md). API 는 docstring·type hint (헌법 XVI)
- [x] T004 [US1] `src/common/kq.py` 신규 — `load_kqs(path?)`·`validate_kq_record(rec)`:
  YAML `kqs` 로드 + data-model.md 검증규칙 1~4(필수필드·P/I 비어있지 않음·query 필드 부재·enum)
- [x] T005 [US1] `src/common/query_spec.md` 작성 (L3 ≤200줄) — 모듈·테스트·FR(002·003·004) 참조 +
  비즈니스 규칙(파생·C/O 제외·정련=리스트성장). API 시그니처는 코드로 일원화(헌법 XI)
- [x] T006 [US1] `src/common/kq_spec.md` 작성 (L3 ≤200줄) — 모듈·테스트·FR(001·002·003·007·009) 참조 +
  검증규칙·단일출처
- [x] T007 [US1] `config/key_questions.yml` 스키마 적합 점검 — 기존 ERAS 레코드가 스키마 준수,
  검색식(`query`) 필드 없음 확인 (kq_schema.md)

## Phase 4: User Story 2 — 두-검증 앵커 (P2)

**Goal**: 검증 A/B 앵커(근거·랜드마크 PMID)와 guideline 앵커를 보유·검증.
**Independent Test**: 앵커 달린 KQ 로드 시 필드 보존 + 중복 PMID(A∩B) 경고.
**Depends on**: US1 (`kq.py` 존재) — 동일 파일 확장이라 순차.

- [x] T008 [US2] `src/common/kq.py` 에 guideline 앵커(name/pmid/date=T)·`guideline_refs`·
  `post_guideline_landmarks` 필드 처리 + 검증규칙 5(A∩B 중복 PMID 경고) 추가 (data-model.md)
- [x] T009 [US2] `src/common/kq_spec.md` 에 검증 앵커 규칙·FR-006 매핑 반영

## Phase 5: User Story 3 — 감시 제어 (P3)

**Goal**: enabled·design_strictness·collection·question_type 기본·검증.
**Independent Test**: `enabled:false` KQ 제외, `design_strictness` 미지정 → `loose`.
**Depends on**: US1 (`kq.py`) — 순차.

- [ ] T010 [US3] `src/common/kq.py` 에 `design_strictness` 기본 `loose`·`enabled` 기본 처리·
  `question_type` enum 검증·`collection` 통과 추가 (FR-007, data-model.md)
- [ ] T011 [US3] `src/common/kq_spec.md` 에 제어 속성 기본·FR-007 반영

## Phase 6: Polish & 등가 검증

- [ ] T012 [P] `quickstart.md` §1(검색식 파생)·§2(로드)·§4(단위테스트) 실행 통과 확인
- [ ] T013 recall 참값 재현 — `uv run python -m alert.validate_kq` → 검증 A 8/9·B 2/2 일치 확인
  (`recall_baseline.yml`)
- [ ] T014 [P] 헌법 XI 줄수 점검 — `query.py`·`kq.py` ≤300줄, `*_spec.md` ≤200줄

## Dependencies

- Phase 1 → Phase 3(US1) → US2·US3(US1 의 `kq.py` 확장, 순차) → Phase 6.
- T002(test) 는 T003(구현) 전 작성(TDD, 등가 가드).
- US2·US3 는 동일 `kq.py`·`kq_spec.md` 수정이라 상호 병렬 불가(순차).

## Parallel Opportunities

- T001 [P] (setup, 독립).
- T002 [P] (test 파일, 구현과 다른 파일 — 단 T003 전).
- T012·T014 [P] (검증·점검, 서로 독립).

## MVP Scope

**US1(Phase 3)** = MVP — KQ·PICO 정의 + 검색식 파생 + 등가 검증. 이것만으로 003(검색) 의
입력이 완성된다. US2·US3 는 증분.
