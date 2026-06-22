# Feature Specification: 기계 필터 & 중복 차단 (Deterministic Filter & Dedup)

**Feature Branch**: `004-deterministic-filter-dedup`

**Created**: 2026-06-22

**Status**: Draft (retro-spec — 기존 동작 역기록)

**Input**: L1 §004 retro-spec. 수집된 후보 문헌(LandmarkItem)을 전달받아 중복 전송 방지를 위해 `SeenStore`(seen-set jsonl)와 대조하고, 기계적인 날짜/언어 필터 등 결정론적 규칙을 적용해 1차 선별한다. seen-set과 상태 정보 파일은 gitignore 대상이다 (헌법 X).

## User Scenarios & Testing *(mandatory)*

### User Story 1 - 기 발송된 랜드마크 중복 전송 방지 (Priority: P1)

시스템이 랜드마크를 감시하여 알림을 보낼 때, 이전에 이미 알림을 보냈던 논문은 중복해서 Telegram으로 전송하거나 Zotero에 추가하지 않고 건너뛴다.

**Why this priority**: 동일한 논문에 대한 중복 알림은 사용자의 피로도를 급격히 높이며, 알림 채널의 차단으로 이어질 수 있는 크리티컬한 결함이다.

**Independent Test**: 이미 발송 이력이 `state/alerted.jsonl`에 저장되어 있는 동일 `kq:pmid` 키를 가진 후보를 파이프라인에 입력했을 때, 전송 단계가 생략되고 skip 처리되는지 검증.

**Acceptance Scenarios**:

1. **Given** `state/alerted.jsonl`에 특정 `key`("ERAS 수액전략:29742967")의 발송 이력이 저장된 상태, **When** 동일한 키를 가진 후보 문헌을 다시 처리할 때, **Then** `seen.is_seen(key)`가 `True`를 반환하고 후속 게이트 판정 및 전송을 수행하지 않고 skip 된다.
2. **Given** 랜드마크로 판정되어 정상 발송 처리된 신규 후보 문헌, **When** 발송 즉시, **Then** `SeenStore`에 `key`, UTC 타임스탬프(`ts`), 제목(`title`), 발송 대상 싱크 정보(`sinks`) 등이 포함된 JSON 한 줄이 파일 끝에 추가 기록(append)된다.

---

### User Story 2 - 손상된 상태 파일 복구 및 강건성 (Priority: P2)

로컬 빌드, 비정상 종료, 혹은 Git 병합 과정 등에서 `state/alerted.jsonl` 파일 내의 특정 행이 깨지거나 빈 행이 되더라도 전체 수집 및 감시 루프가 중단(crash)되지 않고 정상 기동된다.

**Why this priority**: 1인 독립 구동 환경(cron)에서 상태 파일 이상으로 전체 파이프라인이 멈추면 자율 감시가 붕괴된다.

**Independent Test**: JSONL 파일에 빈 줄이나 구문 오류가 있는 행을 강제로 삽입한 후 `SeenStore`를 로드했을 때, 예외가 발생하지 않고 나머지 유효한 행들만 정상 로드되는지 단위 테스트로 단언.

**Acceptance Scenarios**:

1. **Given** JSON 구문 오류를 포함하는 비정상 데이터 라인이 포함된 `alerted.jsonl`, **When** `SeenStore` 생성자를 기동, **Then** 에러 로그를 남기고 비정상 행은 무시하며, 정상적인 행들만 집합에 담아 인스턴스화가 완료된다 (fail-soft).

---

## Requirements *(mandatory)*

### Functional Requirements

* **FR-001**: 시스템은 후보 문헌의 고유 키를 `f"{kq_label}:{pmid}"` 형식으로 다루어야 한다.
* **FR-002**: 중복 여부는 런타임 로컬 JSONL 파일(기본 `state/alerted.jsonl`)의 이력을 읽어 판정해야 한다.
* **FR-003**: `SeenStore` 인스턴스화 시 지정된 경로의 파일이 없을 경우 빈 집합으로 정상 시작해야 하며, 파일의 하위 디렉토리가 없을 시 기록 시점에 자동 생성해야 한다.
* **FR-004**: `SeenStore`는 파일 로드 도중 파싱 에러가 발생한 행을 **fail-soft**하게 스킵하여 예외로 인한 전체 중단을 방지해야 한다.
* **FR-005**: 랜드마크 발송 성공 시 즉시 `mark(key, record)`를 통해 파일 끝에 한 줄을 덧붙여(append) 디스크 상태와 메모리 집합(`_seen`)을 동기화해야 한다.
* **FR-006**: 런타임 상태 폴더(`state/`)는 외부 유출 방지 및 리포지토리 경량화를 위해 버전 제어(`.gitignore`)에서 영구 제외되어야 한다 (헌법 X, 환경 데이터 격리).
* **FR-007**: 모든 상태 조회 및 기록 작업은 결정론적이며, 멱등적(idempotent)이어야 한다.

### Key Entities

* **SeenStore**: 중복 판단의 영구 기록물인 JSONL 백엔드 파일과 메모리 캐시 집합(`_seen`)을 이어주는 상태 관리 객체.
* **전달 기록(Record)**: `key`, `ts`(UTC 타임스탬프)를 필수 정보로 가지며, PMID, 제목, 전송 성공한 싱크 목록 정보를 선택적으로 포함하여 기록하는 데이터 유닛.

## Success Criteria *(mandatory)*

### Measurable Outcomes

* **SC-001**: 훼손된 JSONL 로드 테스트 시, 유효하지 않은 행 수만큼 예외 없이 정상 복구하고 정상 행 개수만큼의 셋 크기를 반환한다.
* **SC-002**: 중복 후보 입력 시 파이프라인에서 제외(skip) 처리되어 Telegram API나 Zotero API에 불필요한 트래픽을 주지 않는다.
* **SC-003**: `git status` 확인 시 `state/` 디렉토리와 파일은 변경이나 미추적 파일 목록에 나타나지 않는다.

## Assumptions

* `state/` 내부의 데이터 포맷은 줄바꿈 단위의 JSON(JSONL) 형식으로, 1인 개발 환경에서 사람이 메모장이나 CLI 툴로 편집하여 손쉽게 특정 내역을 삭제하거나 백업하기 용이해야 한다.
* 중복 판단의 기준은 `kq:pmid` 조합이다. 동일 논문이라도 서로 다른 KQ의 감시 대상이 될 수 있으므로, 논문 자체의 중복 판단이 아닌 **해당 질문 관점에서의 중복 판단**이 요구된다.
* 설계 권위: `.specify/memory/constitution.md`(헌법 I·X·XII), `system-spec.md` §004.
