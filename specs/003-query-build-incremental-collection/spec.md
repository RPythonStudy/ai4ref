# Feature Specification: 검색식 빌드 & 증분 수집 (Query Build & Incremental Collection)

**Feature Branch**: `003-query-build-incremental-collection`

**Created**: 2026-06-22

**Status**: Draft (retro-spec — 기존 동작 역기록)

**Input**: L1 §003 retro-spec. KQ 의 PICO 에서 파생한 검색식으로 PubMed 신규분만 증분
수집한다. `build_query`(I AND P, 001 산출)를 입력으로, 지침-이후(pdat) 게이트 + 최근색인
(edat reldate) 창으로 esearch → efetch 로 후보 문헌(PMID·제목·초록·저널·연도)을 확보한다.
설계 하드필터(`[pt]` 등)는 검색식에 넣지 않는다(헌법 IV).

## User Scenarios & Testing *(mandatory)*

### User Story 1 - PICO 검색식으로 PubMed 후보 메타 확보 (Priority: P1)

감시 시스템이 정의된 KQ 의 PICO(OR-블록)에서 파생한 넓은 검색식으로 PubMed 를 질의하고,
일치하는 문헌의 식별자(PMID)와 선별에 필요한 메타데이터(제목·초록·저널·연도)를 가져온다.
이 후보 목록이 이후 기계필터·게이트 선별의 입력이 된다.

**Why this priority**: 후보를 실제로 가져오지 못하면 선별·검증·알림 어떤 후속 단계도
동작할 수 없다. 감시 파이프라인의 데이터 공급원이다.

**Independent Test**: 활성 KQ 1건에 대해 검색식을 파생해 질의하면, PMID 목록과 각 PMID 의
제목·초록·저널·연도가 채워진 후보 항목이 산출되는지로 단독 검증(`--collect-only`).

**Acceptance Scenarios**:

1. **Given** PICO 가 정의된 활성 KQ, **When** 후보 수집을 실행, **Then** PICO 에서 파생된
   검색식으로 PubMed 가 질의되고 매칭 PMID 별 제목·초록·저널·연도가 채워진 후보 목록이 반환된다.
2. **Given** 다수 PMID(>100), **When** 메타 확보, **Then** 100건 단위로 나눠 가져오며 일부
   배치 파싱 실패 시 해당 배치만 건너뛰고 나머지는 계속 수집된다(fail-soft).
3. **Given** 검색식, **When** 질의 문자열 구성, **Then** 설계 하드필터(`[pt]` 등)는 검색식에
   포함되지 않는다(의미판단은 후속 게이트 책임).

---

### User Story 2 - 신규 색인분만 증분 수집 (Priority: P2)

감시 시스템이 매 실행 시 PubMed 전체가 아니라 "최근 N일 사이 새로 색인된 것"만 가져와,
매일 소량(증분)으로 효율적으로 감시한다.

**Why this priority**: 일일 자율 감시(§008)의 비용·시간을 좌우한다. 증분이 없으면 매번
전체를 재조회해 비효율적이다. 단, 최종 중복 차단 권위는 seen-set(§004)이고 이 창은 재조회량을
줄이는 효율 장치다.

**Independent Test**: `reldate` 창을 좁히면(예: 1일) 후보 수가 줄고, 넓히면(예: 60일) 늘어나는지로
검증. 창 길이는 실행 인자로 조절 가능.

**Acceptance Scenarios**:

1. **Given** `reldate=N`, **When** 수집, **Then** 색인일(edat) 기준 최근 N일 내 신규분만 질의된다.
2. **Given** 동일 KQ 를 좁은 창과 넓은 창으로 각각 수집, **When** 비교, **Then** 좁은 창의
   후보가 넓은 창의 부분집합이다(창은 결정론적 시간 필터).

---

### User Story 3 - 지침-이후 한정 & 무비용 운영 (Priority: P3)

감시 시스템이 지침 제정시점(T) 이후 출판물만 후보로 삼고(이전 근거는 이미 지침에 반영됨),
유료 API 없이 rate-limit 을 지켜 안정적으로 동작한다.

**Why this priority**: 랜드마크 감시의 의미상 "지침 이후 새 근거"가 대상이며, 무비용·로컬
운영(헌법 III)을 유지해야 한다. 기본값으로도 동작하나 키가 있으면 더 빠르다.

**Independent Test**: `guideline.date`(T) 가 2018 인 KQ 의 후보가 모두 출판연도 2019 이상인지,
그리고 API 키 없이도(rate-limit sleep 적용) 수집이 완료되는지로 검증.

**Acceptance Scenarios**:

1. **Given** `guideline.date = T` 인 KQ, **When** 검색식 구성, **Then** 출판일(pdat) ≥ T+1년
   조건이 검색식에 결합되어 지침-이전 문헌은 제외된다.
2. **Given** PUBMED_API_KEY 부재, **When** 수집, **Then** 더 긴 rate-limit 간격으로(키 있으면 더
   짧게) 질의가 수행되어 차단 없이 완료된다.
3. **Given** `--collect-only` 실행, **When** 수집, **Then** 후보만 출력되고 선별·전송은 일어나지
   않는다(수집 배관 단독 점검).

---

### Edge Cases

- 활성 KQ 가 0건이면? → 경고 로그 후 빈 후보 목록 반환(실패 아님).
- 검색식이 매칭 0건이면? → 빈 PMID 목록 → 빈 후보(정상 종료).
- efetch 응답 XML 파싱 실패? → 해당 100건 배치만 skip, 로그 남기고 나머지 계속.
- 어떤 PMID 의 초록이 비어 있으면? → 초록 빈 문자열로 후보 포함(누락 아님, 게이트가 판단).
- `guideline.date` 가 없거나 4자리 연도로 파싱 불가? → 지침-이후 게이트 적용 불가(정의 오류로 표면화).
- 테스트 시 후보 폭증? → 실행 인자(`limit`)로 후보 수 상한 지정 가능.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: 시스템은 활성 KQ 의 PICO 에서 검색식을 **파생**해 사용해야 하며(001 의
  `build_query` 재사용, 재구현 금지), 검색식을 별도 저장하지 않는다(단일 출처 = PICO).
- **FR-002**: 시스템은 검색식으로 PubMed 를 질의해 매칭 문헌의 식별자(PMID) 목록을 얻어야 한다.
- **FR-003**: 시스템은 각 PMID 의 선별용 메타데이터 — 제목·초록·저널·출판연도 — 를 확보해야 한다.
- **FR-004**: 메타 확보는 다수 PMID 를 **100건 단위 배치**로 처리하고, 배치 단위 파싱 실패는
  **fail-soft**(해당 배치 skip·로그, 전체 중단 금지)로 다뤄야 한다.
- **FR-005**: 증분 수집은 색인일(edat) 기준 **최근 N일 창**(`reldate`)으로 신규분만 질의해야
  하며, N 은 실행 인자로 조절 가능해야 한다(기본 = 일일 운영 친화값).
- **FR-006**: 검색식에는 지침 제정시점(T) **이후 출판**(pdat ≥ T+1년) 조건을 결합해 지침-이전
  문헌을 제외해야 한다.
- **FR-007**: 검색식에는 연구설계 하드필터(`[pt]` 등)를 **포함해서는 안 된다** — 설계·관련성
  판단은 후속 게이트(§006) 책임이다(헌법 IV, 검증 A recall 보호).
- **FR-008**: 수집은 **무비용·로컬**이어야 한다 — PUBMED_API_KEY 는 선택이며, 키 유무에 따라
  rate-limit 간격만 달라지고(키 있으면 더 빠름) 키 없이도 동작해야 한다(헌법 III).
- **FR-009**: 수집 대상은 `enabled` 이며 `pico`·`guideline` 을 가진 KQ 로 한정되어야 한다.
- **FR-010**: 수집은 **결정론적**이어야 한다 — 동일 입력(검색식·창·시점)에 동일 후보 집합.
  (중복 차단·상태는 본 피처 밖 = §004 seen-set.)
- **FR-011**: 선별·전송 없이 수집 배관만 점검하는 **collect-only 모드**와 후보 수 **상한
  (limit)** 을 제공해야 한다(테스트·운영 점검용).

### Key Entities *(include if feature involves data)*

- **후보 문헌(Candidate)**: 하나의 PubMed 매칭 결과. 식별자(PMID) + 선별용 메타(제목·초록·
  저널·연도) + 소속 KQ 라벨로 구성. 이후 필터·게이트의 입력 단위.
- **증분 창(Window)**: 색인일(edat) 기준 최근 N일. 재조회량을 줄이는 결정론적 시간 필터.
- **지침-이후 게이트(Date Gate)**: 출판일(pdat) ≥ T+1년 조건. 검색식에 결합되는 결정론 필터.
- **파생 검색식(Derived Query)**: `(I 블록) AND (P 블록)` (001 산출) + 지침-이후 게이트. 별도
  저장 없이 매 실행 파생.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 활성 KQ 에 대해 `--collect-only` 수집을 실행하면 후보 목록(PMID + 제목·초록·
  저널·연도)이 산출되며, 기존 `alert.run --collect-only` 동작과 **등가**다(retro-spec 보존 기준).
- **SC-002**: 파생 검색식에 설계 하드필터(`[pt]` 등)가 **0개** 포함된다 → 검증 A recall **8/9
  유지**(헌법 IV 가드, `tests/acceptance/recall_baseline.yml`).
- **SC-003**: PUBMED_API_KEY **없이도** 수집이 차단 없이 완료된다(rate-limit 간격 자동 적용).
- **SC-004**: 100건을 초과하는 후보도 **배치 처리**로 메타가 채워지며, 1개 배치 파싱 실패가
  전체 수집을 **중단시키지 않는다**.
- **SC-005**: `reldate` 창을 좁히면 후보 수가 단조 감소(좁은 창 ⊆ 넓은 창)한다.

## Assumptions

- 후보 수집의 등가 기준은 기존 구현 동작이다(`src/common/pubmed.py` 의 `esearch`/`efetch_meta`,
  `src/alert/run.py` 의 `_collect_candidates`; 재작성 시 동작 보존).
- 검색식 파생은 001 의 `build_query`(`src/common/query.py`)를 그대로 입력으로 사용한다(본 피처에서
  재구현하지 않음).
- 중복 차단(seen-set)·기계 날짜/언어 필터는 본 피처 범위 밖이며 L1 §004 에서 다룬다. 본 피처는
  "검색식 dating + 증분 질의 + 메타 확보"까지다.
- LLM 게이트 판정(§006)·두-검증 채점(§005)·알림 fan-out(§007)은 본 피처 밖이다.
- PubMed E-utilities 응답 형식(esearch JSON, efetch XML)과 가용성은 외부 제약으로 전제한다.
- 기본 `reldate` 는 일일 운영 친화값(현행 CLI 기본 30, 자율 cron §008 은 1)을 가정한다.
- 설계 권위: `.specify/memory/constitution.md`(헌법 III·IV·VI), `system-spec.md` §003.
