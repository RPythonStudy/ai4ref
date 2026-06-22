# Implementation Plan: 검색식 빌드 & 증분 수집 (Query Build & Incremental Collection)

**Branch**: `003-query-build-incremental-collection` | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/003-query-build-incremental-collection/spec.md`

## Summary

KQ의 PICO 정보로부터 파생된 검색식에 지침 제정 이후 출판 조건(pdat >= T+1년)을 결합하고, 최근 N일 이내 색인(edat reldate)된 신규 문헌 정보만을 PubMed에서 가져와 후보군을 생성하는 수집 배관을 정비합니다. 기존 작동 방식의 등가성(단위 테스트 통과 및 `validate_kq` 수행 결과 동일성)을 완벽히 보장하면서 설계 헌법에 합치되도록 리팩터링합니다.

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: stdlib (`urllib.request`, `xml.etree.ElementTree`), `common.features`, `common.kq`, `common.query`

**Storage**: None (로컬 메모리 데이터 수집 파이프라인)

**Testing**: `pytest`를 통한 단위 테스트(`tests/unit/test_collector.py`) + `validate_kq.py`를 통한 등가성 검증.

**Constraints**:
* 검색식은 PICO에서 실시간 파생(`build_query`)되어야 하며 따로 물리 파일에 저장하지 않습니다 (헌법 III).
* 연구설계 하드필터(`[pt]` 등)는 검색식에 포함되지 않아야 합니다 (헌법 IV).
* API Key 유무에 따라 속도(sleep)를 자동 조절하여 무비용 작동을 보장합니다 (헌법 III).
* 100건 단위의 efetch 배치 처리를 지원하며, 일부 파싱 오류 발생 시 전체 실패가 아닌 fail-soft로 해당 배치만 넘어가도록 보장합니다 (FR-004).

## Constitution Check

| 원칙 | 적용 | 상태 |
|---|---|---|
| I 결정론 backbone | 입력(PICO, 창, 기준년도)에 따른 수집 쿼리 생성 및 매칭 로직은 결정론적입니다. | ✅ |
| IV pt는 게이트 | 연구설계 필터는 쿼리에 결합하지 않고 후속 게이트에서만 다룹니다. | ✅ |
| IX 추가비용 0 | 외부 상용 API나 유료 서비스를 일절 사용하지 않으며, API Key 유무와 관계없이 지연율만 변경되어 작동합니다. | ✅ |
| XI 3계층·줄수 | L3 명세(pubmed_spec.md)를 작성하며 각 모듈의 길이 제약을 준수합니다. | ✅ |
| XII 절차적 명확성 | 복잡한 상속이나 데코레이터가 없는 명시적 제어 흐름을 유지합니다. | ✅ |
| XVI No Duplication | 중복 정보 제거를 위해 `_load_kqs`를 `common.kq.load_kqs`로 대체합니다. | ✅ |

## Project Structure

### Documentation (this feature)

```text
specs/003-query-build-incremental-collection/
├── spec.md              # 피처 요구사항 사양서
├── plan.md              # 이 파일
└── tasks.md             # 할 일 목록
```

### Source Code

```text
src/alert/
└── run.py               # `_load_kqs` 제거 및 `common.kq.load_kqs` 교체

src/common/
├── pubmed.py            # PubMed E-utilities 클라이언트 (기존 코드 보존)
└── pubmed_spec.md       # [NEW] pubmed.py의 L3 규격서 (줄수 <= 200줄)

tests/unit/
└── test_collector.py    # [NEW] 후보 수집기 단위 테스트 (esearch/efetch 모킹)
```
