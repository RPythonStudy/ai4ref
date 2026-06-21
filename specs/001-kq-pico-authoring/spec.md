# Feature Specification: KQ·PICO 정의 (Key Question & PICO Authoring)

**Feature Branch**: `001-kq-pico-authoring`

**Created**: 2026-06-22

**Status**: Draft (retro-spec — 기존 동작 역기록)

**Input**: L1 §001 retro-spec. 진료지침 감시 대상 KQ 를 PICO(OR-블록)로 정의해 config 에
적재한다. 검색식은 PICO 에서 파생, C·O 는 게이트 필터. guideline·두-검증 앵커 보유.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - 진료지침 KQ 를 감시 대상으로 정의 (Priority: P1)

임상 연구자가 감시하고 싶은 진료지침의 핵심 질문(KQ)을 대상(P)·중재(I)·대조(C)·결과(O)로
구조화하고, 대상·중재를 각각 동의어 묶음(OR 리스트)으로 적어 넓은 검색 그물을 구성한다.
이렇게 정의한 KQ 하나가 곧 하나의 감시 단위가 된다.

**Why this priority**: KQ·PICO 정의가 없으면 검색·선별·검증 어떤 후속 단계도 동작할 수 없다.
감시 파이프라인 전체의 입력이자 단일 출처다.

**Independent Test**: KQ 레코드 1건을 작성한 뒤, 그 PICO 로부터 검색식이 파생되어 PubMed
질의가 만들어지는지로 단독 검증 가능(코드 변경 없이 config 추가만으로).

**Acceptance Scenarios**:

1. **Given** 비어있는 감시 설정, **When** 대상·중재를 OR 리스트로, 대조·결과를 필터 기준으로
   하는 KQ 1건을 추가, **Then** 그 KQ 가 감시 대상 목록에 나타나고 검색식이 `(I 항목 OR …)
   AND (P 항목 OR …)` 형태로 파생된다.
2. **Given** 정의된 KQ, **When** 대상(P)·중재(I) 블록에 동의어 항목을 추가(정련), **Then**
   별도 검색식 필드 수정 없이 검색 범위가 넓어진다(검색식은 리스트에서 파생).

---

### User Story 2 - 두-검증 앵커로 검색 품질을 검증 가능하게 함 (Priority: P2)

연구자가 KQ 에 "지침 제정(T) 이전의 근거 문헌"과 "T 이후의 practice-changing 랜드마크"를
식별자(PMID)로 함께 적어, 이 KQ 의 검색식이 *알려진 문헌을 실제로 잡는지* 나중에 자동
채점할 수 있게 한다.

**Why this priority**: 검색식의 recall 을 객관적으로 측정·정련하는 기준점(golden). 없으면
검색식이 좋은지 나쁜지 판단할 수 없다.

**Independent Test**: 검증 앵커가 달린 KQ 에 대해 채점을 돌려 "근거 재현율 / 랜드마크 포착"
점수가 산출되는지로 검증.

**Acceptance Scenarios**:

1. **Given** 근거 PMID 목록(검증 A)과 랜드마크 PMID 목록(검증 B)이 달린 KQ, **When** 채점,
   **Then** "재현 N/M" 형태의 recall 점수가 KQ 별로 산출된다.

---

### User Story 3 - 감시 동작 제어 (Priority: P3)

연구자가 KQ 별로 감시 on/off(`enabled`), 설계 엄격도(`design_strictness`), 보관처(Zotero
`collection`), 질문 유형(`question_type`)을 지정해 감시 동작을 조절한다.

**Why this priority**: 다수 KQ 운영 시 선택적 활성화·튜닝에 필요하나, 단일 KQ MVP 에는
기본값으로도 충분하다.

**Independent Test**: `enabled:false` 인 KQ 가 감시에서 제외되는지로 검증.

**Acceptance Scenarios**:

1. **Given** `enabled:false` KQ, **When** 감시 대상 로드, **Then** 해당 KQ 는 제외된다.
2. **Given** `design_strictness` 미지정 KQ, **When** 로드, **Then** 기본값 `loose` 로 간주된다.

---

### Edge Cases

- 대상(P) 또는 중재(I) 블록이 비어 있으면? → 검색식 파생 불가, 유효하지 않은 KQ 로 거부.
- 검증 앵커(guideline_refs)가 없는 KQ? → 채점 대상에서 제외(검증 없이 감시만).
- 같은 PMID 가 근거(A)와 랜드마크(B)에 중복? → 정의 오류로 표면화.
- `guideline.date`(T)가 없으면? → 시간축 분리(A=T이전/B=T이후) 없이 전체 기간 채점.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: 시스템은 KQ 를 PICO(대상 P·중재 I·대조 C·결과 O)로 정의할 수 있어야 한다.
- **FR-002**: 대상(P)·중재(I)는 각각 **OR 리스트**(복수 검색어 항목)로 표현되어야 한다.
- **FR-003**: 검색식은 P·I 리스트에서 `(I 항목 OR …) AND (P 항목 OR …)` 로 **파생**되어야
  하며, 정련된 검색식을 **별도 필드로 저장해서는 안 된다**(단일 출처 = PICO 리스트).
- **FR-004**: 대조(C)·결과(O)는 검색식에서 **제외**하고 후속 선별(게이트)의 필터 기준으로만
  쓸 수 있도록 보관되어야 한다.
- **FR-005**: 각 KQ 는 감시 대상 지침 앵커(`guideline`: 이름·식별자 PMID·제정시점 T)를
  보유해야 한다.
- **FR-006**: 각 KQ 는 두-검증 앵커 — 검증 A(T 이전 지침 근거 PMID 목록)와 검증 B(T 이후
  랜드마크 PMID 목록) — 를 보유할 수 있어야 한다.
- **FR-007**: 각 KQ 는 질문 유형(intervention|diagnostic|prognostic|predictive),
  설계 엄격도(strict|loose, 기본 loose), 보관처 컬렉션명, 감시 토글(enabled)을 가질 수 있어야
  한다.
- **FR-008**: 정련(검색 품질 개선)은 P·I OR 리스트에 항목을 **추가**하는 방식으로만 이뤄지며,
  리스트 성장이 곧 recall 향상이어야 한다.
- **FR-009**: 다수 KQ 를 추가·운영해도 KQ 간·코드 변경 없이 **설정만으로** 확장 가능해야 한다.

### Key Entities *(include if feature involves data)*

- **KQ 레코드**: 하나의 감시 단위. PICO + 지침 앵커 + 검증 앵커 + 제어 속성(유형·엄격도·
  토글·컬렉션)으로 구성.
- **PICO(OR-블록)**: P·I = 검색어 OR 리스트(검색식 파생원), C·O = 게이트 필터 기준.
- **지침 앵커(guideline)**: 감시 대상 지침의 이름·식별자(PMID)·제정시점(T).
- **검증셋**: 검증 A(근거 PMID, recall 참값) / 검증 B(랜드마크 PMID, 포착 대상).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 정의된 KQ 의 파생 검색식이 검증 A(지침 근거)의 **8/9 이상**을 재현한다
  (현행 acceptance 참값 = `tests/acceptance/recall_baseline.yml`).
- **SC-002**: 검증 B(T 이후 랜드마크)를 **2/2** 포착한다.
- **SC-003**: 새 KQ 1건 추가가 기존 KQ·코드 수정 0회로, 설정 파일 편집만으로 완료된다.
- **SC-004**: PICO OR 리스트에 항목을 추가하는 정련 후 검증 A recall 이 **감소하지 않으며**,
  그 과정에서 추가되는 별도 검색식 필드는 **0개**다.

## Assumptions

- KQ 레코드는 사람이 읽고 편집 가능한 단일 설정 파일(YAML)에 보관된다(기존 구현 참조:
  `config/key_questions.yml`).
- PICO → 검색식 파생은 기존 `build_query` 동작과 등가여야 한다(기존 구현 참조:
  `src/common/query.py`; 재작성 시 동작 보존).
- 검증 채점의 acceptance 참값은 이미 박제됨: `tests/acceptance/recall_baseline.yml`
  (검증 A 8/9 — 누락 19602972 = 의도된 천장 · B 2/2).
- 수집 트랙(L1 §009) 주제를 동일 설정 파일에 `type` 필드로 통합할지 여부는 본 피처 범위
  밖이며 L2(009)에서 확정한다(헌법 X).
- 설계 권위: `.specify/memory/constitution.md`(헌법 II·III·VI·XIII), `system-spec.md` §001.
