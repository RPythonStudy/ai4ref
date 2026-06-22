# Implementation Plan: LLM 게이트 (Relevance & Landmark Gate — Precision)

**Branch**: `006-llm-relevance-landmark-gate` | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/006-llm-relevance-landmark-gate/spec.md`

## Summary

후보 논문의 제목·초록을 LLM 으로 읽어 2단계(① 관련성 = 대상 P AND 중재 I 동시매칭 / ② 랜드마크
= practice-changing)로 판정하는 정밀(precision) 레이어를 헌법 정합으로 정비한다. 교체형 백엔드
(stub/claude_cli/api)·무비용 로컬(OAuth Haiku)·fail-soft·견고 JSON 추출을 보존하며, 기존
`llm_gate.judge` 동작과 등가(acceptance `gate_baseline.yml` 재현)임을 단위테스트로 가드한다.
핵심 정비 = 출력 키 명명을 acceptance 스키마(`is_relevant`)에 맞추고 `run.py` 호출부와 정합한다.

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: stdlib (`subprocess`, `json`), `common.features`(log). 운영 백엔드 =
로컬 `claude` CLI(OAuth·Haiku, 추가비용 0). 단위테스트는 CLI 미호출(stub/모의).

**Storage**: None (후보 단위 무상태 판정)

**Testing**: `pytest` 단위테스트(`tests/unit/test_gate.py`) — 2단계 순서·fail-soft·`_extract_json`
래퍼 파싱을 백엔드 모의로 결정론 검증. 실 `claude_cli` 판정 정확도는 acceptance 참값
(`gate_baseline.yml`)으로 분리(LLM 비결정성은 단위테스트 밖).

**Constraints**:
* 게이트 출력 스키마 = `{is_relevant, is_landmark, reason}` (acceptance·§006 권위). 현행 코드
  `relevant` 키를 `is_relevant` 로 맞추고 `run.py` 게이트 호출부와 정합 (FR-006).
* `[pt]` 등 설계 필터는 검색식이 아닌 게이트의 *의미 판단 컨텍스트*로만 사용 (헌법 IV·V·IX, FR-012).
* 판정 실패·파싱 실패 시 보수적 보류(`is_relevant=false`·`is_landmark=false`)·로그·계속 (FR-009).
* 운영 무비용 — 유료 API 미연결, 로컬 OAuth 백엔드로 동작 (헌법 III, FR-007·008).
* ②(랜드마크)는 ①(관련성)이 true 일 때만 평가 (2단계 순서, FR-001·SC-005).

## Constitution Check

| 원칙 | 적용 | 상태 |
|---|---|---|
| II recall·precision 2층 | 게이트=precision, 두-검증(§005)=recall 별개 레이어로 운영·평가. | ✅ |
| III 추가비용 0 / 로컬 | 운영 백엔드 = 로컬 claude CLI(OAuth·Haiku). 유료 API 미도입. | ✅ |
| IV·V pt=게이트 의미판단 | 설계(`[pt]`)는 strict 판정 컨텍스트일 뿐 검색식 하드필터 아님. | ✅ |
| IX LLM 경계 | 의미판단(관련성·랜드마크)만 LLM, 결정론 backbone과 분리. | ✅ |
| XI 3계층·줄수 | L3 명세(`gate_spec.md`) 작성, 코드 ≤300·`*_spec.md` ≤200줄. | ✅ |
| XII 절차적 명확성 | 명시적 제어흐름(백엔드 분기·2단계 순서), 숨은 마법 없음. | ✅ |
| XIII acceptance 참값 | `gate_baseline.yml`(landmark 2/2·reject FP 2/2) = 등가 기준. | ✅ |

## Project Structure

### Documentation (this feature)

```text
specs/006-llm-relevance-landmark-gate/
├── spec.md                 # 피처 요구사항 사양서
├── plan.md                 # 이 파일
├── checklists/requirements.md
└── tasks.md                # [NEXT] /speckit-tasks 산출
```

### Source Code

```text
src/alert/
├── llm_gate.py             # judge·_claude_cli·_stub·_extract_json — 출력키 is_relevant 정합·동작 보존
├── gate_spec.md            # [NEW] llm_gate.py L3 규격서 (judge 계약·2단계·백엔드·fail-soft, ≤200줄)
├── run.py                  # 게이트 호출부(ctx_map·verdict 소비) is_relevant 정합 확인
└── classify_qtype.py       # 0단계 유형판별 보조 — 권위=사용자입력, 본 피처는 컨텍스트 입력으로만 참조

tests/unit/
└── test_gate.py            # [NEW] 2단계 순서·fail-soft·_extract_json 래퍼 파싱 (백엔드 모의)
```

## Phase 0 / 1 Notes (retro-spec — 경량)

기존 동작 역기록이라 신규 research·data-model·contracts 산출물은 만들지 않는다(003/004/005 동일).
- **계약**: 게이트 인터페이스 계약 = `judge(paper, kq, backend) -> {is_relevant, is_landmark, reason}`.
  본문은 `gate_spec.md`(L3)로 일원화(헌법 XI, API 시그니처는 코드·docstring 권위).
- **등가 검증**: `gate_baseline.yml` 재현(landmark_true 2/2·reject_false_positive 2/2). 단위테스트는
  배관(2단계 순서·fail-soft·JSON 추출)만 결정론 가드.
- **정합 리스크**: `relevant`→`is_relevant` 키 변경이 `run.py`·기타 소비처를 깨지 않는지 grep 확인 후 일괄.
