# Implementation Plan: 기계 필터 & 중복 차단 (Deterministic Filter & Dedup)

**Branch**: `004-deterministic-filter-dedup` | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/004-deterministic-filter-dedup/spec.md`

## Summary

감시 알림의 안정성과 단일성을 유지하기 위한 중복 제거용 `SeenStore`(seen-set jsonl)와 결정론적 1차 기계 필터 규격을 정비합니다. 기존 작동 방식(`common.state.SeenStore`)의 안정성 및 예외 처리를 단위 테스트로 밀도 있게 점검하며 설계 헌법에 합치되도록 정렬합니다.

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: stdlib (`json`, `pathlib`, `datetime`)

**Storage**: Append-Only JSONL 텍스트 파일 (기본 `state/alerted.jsonl`, gitignore 적용)

**Testing**: `pytest`를 통한 단위 테스트(`tests/unit/test_state.py`) 및 전체 회귀 테스트.

**Constraints**:
* `state/` 디렉토리는 런타임 데이터 격리 정책에 따라 `.gitignore`로 버전 관리에서 영구 제외되어야 합니다 (헌법 X).
* 파일 손상이나 병합 시 충돌 흔적 등으로 파일 구문 오류가 발생하더라도, 이를 예외로 방치하지 않고 fail-soft하게 회복 기동해야 합니다 (FR-004).
* 모든 중복 제거 기준은 질문 관점인 `f"{kq}:{pmid}"` 형식(key)으로 통일되어야 합니다 (FR-001).
* 코드 모듈 및 L3 명세의 라인 수 규칙을 준수합니다 (헌법 XI).

## Constitution Check

| 원칙 | 적용 | 상태 |
|---|---|---|
| I 결정론 backbone | SeenStore의 중복 차단 여부 판단은 100% 결정론적입니다. | ✅ |
| X Config 단순화 | DB 사용을 완벽히 배제하고 파일 기반(JSONL) 백엔드만을 유지합니다. | ✅ |
| XI 3계층·줄수 | L3 명세(state_spec.md)를 작성하며 각 파일 라인 수 기준을 충족합니다. | ✅ |
| XII 절차적 명확성 | SeenStore 객체는 클래스 상속 없이 명시적이고 읽기 쉬운 구조로 유지됩니다. | ✅ |
| XVI No Duplication | 상태 경로 및 데이터 정의를 단일 진실 원천(`base.LandmarkItem.key`)으로 연동합니다. | ✅ |

## Project Structure

### Documentation (this feature)

```text
specs/004-deterministic-filter-dedup/
├── spec.md              # 요구사항 명세서
├── plan.md              # 이 파일
└── tasks.md             # 작업 목록
```

### Source Code

```text
src/common/
├── state.py             # SeenStore 구현체 (기존 코드 보존)
└── state_spec.md        # [NEW] state.py의 L3 규격서 (줄수 <= 200줄)

tests/unit/
└── test_state.py        # [NEW] SeenStore 예외 회복 및 마킹 단위 테스트
```
