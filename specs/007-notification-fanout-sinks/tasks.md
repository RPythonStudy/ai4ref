# Tasks: Sink 멀티탭 (Notification Fan-out)

**Feature**: `specs/007-notification-fanout-sinks/` | retro-spec (기존 동작 보존 + 헌법 정합 재작성)
**등가 기준**: 기존 `alert.run` 의 sink fan-out 동작 보존 — 토글로 켜진 통로만 활성, prereq
미충족 skip, 전송 실패 격리, stdout 항상 동작. 코드 표면 = `src/notify/*` + `run.py` 호출부.
기존 notify 단위테스트 없음 → 신규 추가.

## Phase 1: Setup

- [ ] T001 [P] `tests/unit/test_notify.py` 단위테스트 파일 생성 — `urllib.request.urlopen`/
  `pyzotero` 모의 인프라 골격. `from notify.registry import build_sinks, fan_out`·
  `from notify.base import LandmarkItem`·`from notify.stdout_sink import render`.

## Phase 2: Foundational

(추가 인프라 없음 — notify 패키지·features.yml 토글·LandmarkItem 기존 완비. 다음 단계로 진행)

## Phase 3: User Story 1 — 랜드마크 fan-out + 화면 렌더 (P1, MVP)

**Goal**: 켜진 모든 통로로 랜드마크 1건 분배. stdout 의존성 0·항상 동작·사람 읽을 형식.
**Independent Test**: 통로 2개 활성 + 랜드마크 1건 → 통로별 결과 반환. render 가 제목·메타·
근거·링크 포함하고 내부 abstract 는 제외.

- [ ] T002 [US1] `src/notify/registry.py`·`base.py` fan-out 경로 동작 보존 확인 — `fan_out(item,
  sinks)` 가 통로별 `{name: ref|None}` 반환, `Sink` 인터페이스(available·send) 일관(FR-001·004).
- [ ] T003 [US1] `src/notify/stdout_sink.py` `render` 보존 확인 — 제목·메타(저널·연도)·한 줄 근거·
  PubMed 링크 포함, 내부 전용 `abstract` 미노출(FR-006). `LandmarkItem.url` 자동 파생(FR-011).
- [ ] T004 [P] [US1] `tests/unit/test_notify.py` — render 출력에 title·journal·reason·url 포함·
  abstract 제외 단언 + `LandmarkItem(pmid=…)` url 자동 파생·`key="kq:pmid"` 단언.

## Phase 4: User Story 2 — 토글 + prereq skip (fail-soft) (P2)

**Goal**: features.yml 토글로 켜진 통로만 활성, 키·의존성 부재 시 자동 skip(오류 아님).
**Independent Test**: 텔레그램 켜되 토큰 비움 → skip·로그, stdout 정상. 꺼진 통로 미인스턴스화.
**Depends on**: US1(build_sinks/fan_out 경로).

- [ ] T005 [US2] `src/notify/registry.py` `build_sinks` 보존 확인 — `enabled` 토글 필터 +
  `available()` prereq 점검 → 미충족 skip·`log.warning`(FR-002·003).
- [ ] T006 [US2] `src/notify/telegram_sink.py`·`zotero_sink.py` `available()` 보존 확인 — 토큰/
  키 부재·`pyzotero` 미설치 시 `(False, 사유)` 반환(fail-soft, FR-003·010).
- [ ] T007 [P] [US2] `tests/unit/test_notify.py` — 토큰 빈 telegram = skip·stdout 활성 단언 +
  `enabled:false` 통로 미인스턴스화 단언(모의 features dict).

## Phase 5: User Story 3 — 전송 격리 + Zotero 보관 + 확장 (P3)

**Goal**: 한 통로 전송 실패 격리(다른 통로·다음 항목 계속). Zotero 컬렉션 보관·근거 메모. 등록 1줄 확장.
**Independent Test**: 한 sink 가 예외 → 결과 None·로그, 다른 sink 정상. 새 sink 등록만으로 fan-out 포함.
**Depends on**: 없음(격리는 fan_out 레이어 독립).

- [ ] T008 [US3] `src/notify/registry.py` `fan_out` 전송 격리 보존 확인 — `try/except` 로 통로
  실패 시 `out[name]=None`·`log.error`, 다른 통로 계속(FR-005).
- [ ] T009 [US3] `src/notify/zotero_sink.py` 보존 확인 — 지정 컬렉션(기본 ai4ref-landmarks) 생성/
  탐색·`reason`→Extra 메모 보강·`collector.zotero_add` 재사용(FR-008). 새 통로 = `SINK_CLASSES`
  한 줄 확장 지점 명시(FR-009).
- [ ] T010 [P] [US3] `tests/unit/test_notify.py` — sink A 예외·sink B 성공 시 결과
  `{A:None, B:ref}`·B 호출됨 단언(전송 격리, FR-005).

## Phase 6: Polish & 등가 검증

- [ ] T011 `src/notify/notify_spec.md` 작성 (L3 ≤200줄) — Sink 계약(available·send)·build_sinks·
  fan_out·fail-soft 2종·플러그 확장(SINK_CLASSES)·LandmarkItem·FR001~011 참조.
- [ ] T012 `src/alert/run.py` 호출부 점검 — `build_sinks`/`fan_out` 사용·실수집+stub 백엔드
  외부통로 차단 안전장치 보존(§006/§007 경계). 동작 변화 없음 단언.
- [ ] T013 `uv run python -m pytest -q` 전체 통과 확인(기존 44 + 신규 test_notify).
- [ ] T014 [P] 헌법 XI 줄수 점검 — `notify/*.py` 각 ≤300줄, `notify_spec.md` ≤200줄.

## Dependencies

- Phase 1 → Phase 3(US1) → Phase 4(US2, US1 의존) · Phase 5(US3 독립) → Phase 6.
- US2 는 US1(build_sinks/fan_out) 위에 선다. US3(격리)는 fan_out 레이어 독립.

## Parallel Opportunities

- T001 [P] (setup, 독립).
- T004·T007·T010 [P] (각 US 테스트 — 단 같은 파일이면 순차 편집).
- T014 [P] (줄수 점검, 독립).

## MVP Scope

**US1(Phase 3)** = MVP — 랜드마크를 켜진 통로로 fan-out + stdout 항상 동작. 이것만으로 감시
파이프라인이 사용자에게 가치를 전달한다(003→004→006→007 완결). US2(토글·skip)·US3(격리·보관·
확장)은 증분 견고성.
