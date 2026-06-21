# Implementation Plan: KQ·PICO 정의 (Key Question & PICO Authoring)

**Branch**: `001-kq-pico-authoring` (main; git extension 미설치) | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/001-kq-pico-authoring/spec.md`

## Summary

진료지침 감시 KQ 를 PICO(OR-블록)로 정의·검증하는 데이터 계약과 검색식 파생을 제공한다.
핵심 산출물 = ① `key_questions.yml` 레코드 스키마(단일 출처) ② `build_query(pico)` 검색식
파생(=(I) AND (P)) ③ KQ 레코드 적재·유효성 검증. 기존 동작(ERAS 2012 1건, 검증 A 8/9·B 2/2)
을 보존하며 헌법 정합으로 재작성한다(retro-spec). DB 없음(YAML only).

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: PyYAML(설정 로드), python-dotenv. `build_query` 는 순수 stdlib.

**Storage**: YAML 파일 `config/key_questions.yml` (DB 없음 — 헌법 X)

**Testing**: pytest(단위) + acceptance YAML(`tests/acceptance/recall_baseline.yml`).
검색 recall 실측은 `alert.validate_kq`(PubMed, 005 영역; 본 피처는 derived_query 등가까지).

**Target Platform**: Linux(WSL2), 로컬 CLI

**Project Type**: 단일 프로젝트 (library/CLI 모듈, src 레이아웃)

**Performance Goals**: 검색식 파생·스키마 검증은 결정론·즉시(<10ms). 성능 목표 N/A.

**Constraints**: 결정론·재현(헌법 I). 검색식은 PICO 에서 파생하고 별도 query 필드 저장 금지
(II·III). C·O 는 검색식 제외(II). 파생 검색식이 acceptance 의 `derived_query` 와 일치하고
recall 참값(A 8/9·B 2/2)을 재현해야 등가(XIII). 코드 ≤300줄/L3 명세 ≤200줄(XI).

**Scale/Scope**: 현행 KQ 1건 → 수십 KQ. config 추가만으로 확장(코드 변경 0, FR-009).

## Constitution Check

*GATE: Phase 0 전 통과, Phase 1 후 재점검.*

| 원칙 | 적용 | 상태 |
|---|---|---|
| I 결정론 backbone | 스키마 로드·build_query = 결정론, LLM 미사용 | ✅ |
| II PICO=OR블록 | 스키마가 P·I 를 OR 리스트로 강제, C·O 는 검색식 제외 | ✅ |
| III 정련=리스트성장 | 스키마에 query 필드 없음(파생 전용) | ✅ |
| IX 추가비용 0 | 외부 유료 API 미사용 | ✅ |
| XII 절차적 명확성 | build_query·loader = 순수 함수, 상속·데코레이터 없음 | ✅ |
| XIII Validation Triad | acceptance recall_baseline(8/9·B 2/2)·derived_query 재현 | ✅ |
| XVI No Duplication | key_questions.yml 단일 출처, API 는 docstring/type hint | ✅ |
| XI 3계층·줄수 | L3 `<모듈>_spec.md` 동거, 코드 ≤300·명세 ≤200 | ✅ (Phase 1 산출) |

**위반 없음** → Complexity Tracking 불요.

## Project Structure

### Documentation (this feature)

```text
specs/001-kq-pico-authoring/
├── plan.md              # 이 파일
├── research.md          # Phase 0
├── data-model.md        # Phase 1 (KQ 레코드 엔티티·검증 규칙)
├── quickstart.md        # Phase 1 (등가 검증 가이드)
├── contracts/           # Phase 1 (kq 스키마 + build_query 계약)
└── tasks.md             # Phase 2 (/speckit-tasks)
```

### Source Code (repository root)

```text
config/
└── key_questions.yml          # KQ 레코드(스키마 권위·단일 출처)

src/common/
├── kq.py                      # KQ 적재·유효성 검증 (load_kqs, validate)
├── kq_spec.md                 # L3 명세 (kq.py)
├── query.py                   # build_query(pico) = (I) AND (P)
└── query_spec.md              # L3 명세 (query.py)

tests/
├── acceptance/recall_baseline.yml   # 참값(기존)
└── unit/test_query.py               # build_query 단위 (derived_query 등가)
```

**Structure Decision**: 단일 프로젝트·src 레이아웃(top-level 패키지). 감시·수집 공유
로직은 `common/`. 001 의 코드 표면 = `common/query.py`(검색식 파생) + `common/kq.py`(적재·검증),
각각 L3 `<모듈>_spec.md` 동거. 데이터 권위는 `config/key_questions.yml`. DB·클래스 계층 없음.
