# Specification Quality Checklist: 검색식 빌드 & 증분 수집

**Purpose**: Validate specification completeness and quality before proceeding to planning
**Created**: 2026-06-22
**Feature**: [spec.md](../spec.md)

## Content Quality

- [x] No implementation details (languages, frameworks, APIs)
- [x] Focused on user value and business needs
- [x] Written for non-technical stakeholders
- [x] All mandatory sections completed

## Requirement Completeness

- [x] No [NEEDS CLARIFICATION] markers remain
- [x] Requirements are testable and unambiguous
- [x] Success criteria are measurable
- [x] Success criteria are technology-agnostic (no implementation details)
- [x] All acceptance scenarios are defined
- [x] Edge cases are identified
- [x] Scope is clearly bounded
- [x] Dependencies and assumptions identified

## Feature Readiness

- [x] All functional requirements have clear acceptance criteria
- [x] User scenarios cover primary flows
- [x] Feature meets measurable outcomes defined in Success Criteria
- [x] No implementation details leak into specification

## Notes

- retro-spec: 기존 동작(`pubmed.esearch`/`efetch_meta`, `run._collect_candidates`)을 보존 기준으로
  역기록. 등가 기준 = `alert.run --collect-only` 동작(SC-001).
- 범위 경계 명시: 중복 차단(seen-set)·기계필터 = §004, 게이트 = §006, 채점 = §005, 알림 = §007 —
  모두 본 피처 밖. 본 피처 = "검색식 dating + 증분 질의 + 메타 확보".
- SC-001 의 "기존 동작 등가"는 PubMed 라이브 응답에 의존하므로 결정론 단언이 아닌 **배관/형식**
  등가로 검증한다(PMID 목록·메타 필드 채움·배치·fail-soft). golden recall 단언은 §005 책임.
- FR 표면 용어(SC-002 검증 A 8/9, `[pt]` 금지)는 헌법 IV·`recall_baseline.yml` 과 일치.
- 헌법 IV·III 정합: 하드필터 금지(SC-002), 무비용·키선택(SC-003) — 모두 측정 가능 기준으로 반영.
