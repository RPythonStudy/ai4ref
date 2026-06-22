# Specification Quality Checklist: LLM 게이트 (Relevance & Landmark Gate)

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

- retro-spec: 기존 동작(`llm_gate.judge`/`_claude_cli`/`_stub`/`_extract_json`, `run.py` 게이트
  호출·컨텍스트 구성)을 보존 기준으로 역기록. 등가 기준 = `gate_baseline.yml`(landmark 2/2 ·
  reject FP 2/2 · 백필 96건 FP 0).
- 범위 경계 명시: 두-검증(§005, recall)·수집(§003)·기계필터·dedup(§004)·알림(§007) 모두 본 피처 밖.
  본 피처 = "후보 초록 2단계 판정(precision)".
- recall·precision 2층 분리(헌법) 반영: 게이트=precision, 두-검증=recall 별개 레이어로 명시.
- 헌법 IV·V·IX 정합: `[pt]`=의미판단 컨텍스트이지 검색식 하드필터 아님(FR-012, edge case).
- 헌법 III 정합: 무비용·로컬 백엔드(FR-007·008, SC-002).
- 알려진 정합 과제(구현 단계): 출력 키 명명 — 현행 코드 `relevant` vs baseline/§006 `is_relevant`.
  spec 은 acceptance 스키마(`is_relevant`)를 권위로 두고 plan 에서 코드를 맞추도록 Assumptions 에 명시.
- LLM 비결정성은 precision 레이어 특성으로 수용 — acceptance 단언은 baseline 의 명확 사례(극단 FP/
  명백 랜드마크)로 한정(SC-001), 회색지대 정량단언은 두지 않음.
