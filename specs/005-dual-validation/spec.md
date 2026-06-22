# Feature Specification: 두-검증 (Dual Validation — Recall)

**Feature Branch**: `005-dual-validation`

**Created**: 2026-06-22

**Status**: Draft (retro-spec — 기존 동작 역기록)

**Input**: L1 §005 retro-spec. 지침 앵커와 검증 앵커(PMID 목록)가 설정된 KQ 목록을 받아, 검색식의 상대 재현율(Recall)을 자동으로 측정 및 채점한다.
- 검증 A: 지침 제정 시점(T) 이전 근거문헌(guideline_refs) 재현율 채점 (현행 8/9 목표).
- 검증 B: 지침 제정 시점(T) 이후 랜드마크 문헌(post_guideline_landmarks) 포착 여부 채점 (현행 2/2 목표).
- 0단계 검증: 사용자가 명시한 `question_type`과 로컬 LLM/분류기(`classify_question_type`)가 제안한 유형 간의 일치 여부 비교 출력.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - 지침 근거 재현율(Recall A) 자동 채점 (Priority: P1)

의학 연구자가 작성한 PICO 검색식이 기존에 지침을 제정할 때 근거가 되었던 알려진 문헌들(guideline_refs)을 실제로 검색해낼 수 있는지 컴퓨터가 자동으로 채점하여, 검색 그물의 촘촘함(Recall)을 객관적으로 입증한다.

**Why this priority**: 검색 그물이 성겨서 기존 지침의 중요 논문을 놓친다면, 새로운 랜드마크 또한 모니터링망을 빠져나갈 확률이 높으므로 검색식 검증의 가장 핵심이다.

**Independent Test**: `validate_kq.py`를 실행했을 때 각 KQ별로 "검증 A 지침근거 재현 = N/M"과 개별 PMID별 매칭 결과(✅/❌)가 정확히 콘솔에 출력되는지 확인.

**Acceptance Scenarios**:

1. **Given** 9건의 검증 A PMID를 가진 ERAS 수액전략 KQ, **When** 검증 스크립트 실행, **Then** 각 PMID에 대해 검색식 매칭 여부를 판단하고 최종 재현 점수 `8/9`를 콘솔에 정확히 출력한다.

---

### User Story 2 - 지침 이후 랜드마크 포착율(Recall B) 자동 채점 (Priority: P2)

지침 제정 시점(T) 이후에 발표된 것으로 알려진 주요 랜드마크 논문들(post_guideline_landmarks)을 현재의 검색식이 성공적으로 잡아낼 수 있는지 채점한다.

**Why this priority**: 지침 이후 신규 문헌을 실시간 감시(Alert)하는 시스템의 유효성을 과거 데이터를 기반으로 소급하여 검증한다.

**Independent Test**: `validate_kq.py` 실행 결과로 각 KQ별로 "검증 B 랜드마크 포착 = X/Y" 및 개별 PMID 매칭 결과가 콘솔에 출력되는지 확인.

**Acceptance Scenarios**:

1. **Given** 2건의 랜드마크 PMID를 가진 ERAS 수액전략 KQ, **When** 검증 스크립트 실행, **Then** 각 PMID에 대해 검색식 매칭 여부를 판단하고 최종 포착 점수 `2/2`를 정확히 출력한다.

---

### User Story 3 - 질문 유형 분류 교차 검증 (0단계) (Priority: P3)

사용자가 설정한 질문 유형(`question_type`)과 자연어 분석 분류기가 분류한 유형을 대조하여, 의도된 질문 속성이 올바르게 설정되었는지 교차 검증한다.

**Why this priority**: 질문 유형 설정에 따라 후속 LLM 게이트의 판정 템플릿과 기준이 세부 조정되므로, 설정 오기를 미연에 방지한다.

**Independent Test**: 스크립트 실행 결과에 0단계 유형 일치 여부가 `✅ 일치` 또는 `⚠️ 불일치` 형태로 나타나는지 확인.

---

## Requirements *(mandatory)*

### Functional Requirements

* **FR-001**: 시스템은 활성 KQ 목록을 불러올 때 공인 로더인 `common.kq.load_kqs`를 사용해야 한다.
* **FR-002**: 각 KQ에 대해 `classify_question_type` 분류기를 호출하여 교차 검증 결과를 출력해야 한다 (0단계).
* **FR-003**: 검색식 파생 시 001 피처의 `common.query.build_query` 함수를 재사용해야 한다.
* **FR-004**: 검증 A(지침 이전 문헌) 판단 시 지침 연도(T) 이전 출판물로 제한하기 위해 `pdat` 기간 조건(`1900/01/01` ~ `T년/12/31`)을 쿼리와 결합해야 한다 (지침 제정 시점이 있는 경우).
* **FR-005**: 검증 B(지침 이후 문헌) 판단 시 지침 연도 이후 출판물로 제한하기 위해 `pdat` 기간 조건(`T+1년/01/01` ~ `2030/12/31`)을 쿼리와 결합해야 한다.
* **FR-006**: 개별 PMID의 매칭 판단은 PubMed E-utilities `esearch` API를 통해 검색 조건과 PMID(`AND PMID[uid]`)를 함께 쿼리하여 결과 건수가 1건 이상인지 판단하는 방식을 사용해야 한다.
* **FR-007**: 모든 PubMed API 호출 시 속도 지연(rate limit) 조율 로직을 반영해야 하며, 가능한 한 `common.pubmed` 모듈에 준하는 안전한 호출 방식을 지향하거나 이를 포용해야 한다.

### Key Entities

* **두-검증기 (validate_kq)**: PICO 검색식의 상대 재현율을 검증 A/B 기준으로 자동 채점하고 콘솔에 리포팅하는 실행기 모듈.

## Success Criteria *(mandatory)*

### Measurable Outcomes

* **SC-001**: `validate_kq` 전체 파이프라인 구동 시, ERAS 수액전략 KQ 기준 검증 A 재현 결과가 정확히 `8/9` (88.9%), 검증 B 포착 결과가 정확히 `2/2` (100.0%)를 달성해야 한다 (Acceptance 참값 재현성 보장).
* **SC-002**: `validate_kq` 소스 코드는 중복 코드를 제거하고, 라인 수 제한(<= 300줄) 규격을 준수해야 한다.

## Assumptions

* `validate_kq` 검증 스크립트는 실시간 PubMed API(E-utilities)에 네트워크 요청을 보낸다. 따라서 인터넷 연결과 유효한 `ENTREZ_EMAIL`/`PUBMED_API_KEY` 설정이 필요하다.
* 설계 권위: `.specify/memory/constitution.md`(헌법 VI·XIII), `system-spec.md` §005.
