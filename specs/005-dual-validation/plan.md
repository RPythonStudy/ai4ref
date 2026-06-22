# Implementation Plan: 두-검증 (Dual Validation — Recall)

**Branch**: `005-dual-validation` | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/005-dual-validation/spec.md`

## Summary

PICO 검색식의 상대 재현율(Recall)을 검증 A(지침 이전 근거)와 검증 B(지침 이후 랜드마크) 기준으로 계산하고 채점하는 `validate_kq.py`의 구조를 개선합니다. `common.kq.load_kqs` 및 `common.pubmed.esearch` 통합을 통해 외부 호출을 단일 채널로 결합하고, 핵심 판정 로직을 분리해 테스트 코드를 확보합니다.

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: `common.kq`, `common.pubmed`, `common.query`, `alert.classify_qtype`

**Storage**: None (로컬 메모리 채점 엔진)

**Testing**: `pytest`를 통한 단위 테스트(`tests/unit/test_validation.py`) + 실제 리포트 최종 대조.

**Constraints**:
* API 중복 호출을 방지하기 위해 PubMed E-utilities HTTP 요청은 오직 `common.pubmed.esearch`를 매개로만 수행되어야 합니다 (헌법 XVI).
* 사용자 설정이 잘못된 경우의 조기 실패(fail-fast)와 미지정 기본값 처리가 KQ 로더 레벨에서 일괄 적용되어야 합니다.
* 출력 포맷(0단계 일치 여부, ✅/❌ 아이콘 등)은 기존 사용자의 CLI 보고 환경을 정확히 보존해야 합니다 (헌법 XIII).
* 코드 모듈 및 L3 명세의 줄 수 제약을 충족합니다 (헌법 XI).

## Constitution Check

| 원칙 | 적용 | 상태 |
|---|---|---|
| I 결정론 backbone | 입력 조건에 따른 재현 점수 산정은 결정론적입니다. | ✅ |
| VI 두-검증 책임분리 | 본 피처는 LLM 게이트의 선별성과 독립된 순수 검색 재현율 평가 레이어입니다. | ✅ |
| XI 3계층·줄수 | L3 명세(validate_spec.md)를 작성하며 각 모듈의 길이 제약을 준수합니다. | ✅ |
| XVI API 중복 금지 | `esearch_count`를 폐기하고 `common.pubmed.esearch`로 단일 출처를 적용합니다. | ✅ |

## Project Structure

### Documentation (this feature)

```text
specs/005-dual-validation/
├── spec.md              # 요구사항 명세서
├── plan.md              # 이 파일
└── tasks.md             # 작업 목록
```

### Source Code

```text
src/alert/
├── validate_kq.py       # API 교체, 로더 교체 및 validate_single_kq 함수 추출
└── validate_spec.md     # [NEW] validate_kq.py의 L3 규격서 (줄수 <= 200줄)

tests/unit/
└── test_validation.py   # [NEW] validate_single_kq 모킹 단위 테스트
```
