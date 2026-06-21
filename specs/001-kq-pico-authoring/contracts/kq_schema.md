# Contract: key_questions.yml 스키마 + KQ 로더

`config/key_questions.yml` 의 레코드 계약 + `src/common/kq.py` 적재·검증 계약.
(필드 상세 = data-model.md, 구현 = docstring/type hint — 헌법 XVI)

## YAML 구조

```yaml
kqs:
  - kq: <str>                       # 필수
    question_type: <enum>           # 필수: intervention|diagnostic|prognostic|predictive
    design_strictness: <enum>       # 선택: strict|loose (기본 loose)
    pico:
      P: [<str>, ...]               # 필수, 비어있지 않음 (OR 리스트)
      I: [<str>, ...]               # 필수, 비어있지 않음 (OR 리스트)
      C: <str>                      # 선택 (검색식 제외)
      O: <str>                      # 선택 (검색식 제외)
    guideline:
      name: <str>                   # 필수
      pmid: <str>                   # 선택
      date: <YYYY-MM>               # 선택 (T = 제정시점)
    guideline_refs: [<PMID>, ...]   # 선택 (검증 A)
    post_guideline_landmarks: [<PMID>, ...]  # 선택 (검증 B)
    collection: <str>               # 선택
    enabled: <bool>                 # 선택
```

## 로더 계약 (`common/kq.py`)

| 함수(개념) | 계약 |
|---|---|
| `load_kqs(path?) -> list[dict]` | YAML 의 `kqs` 리스트 반환. 각 레코드 검증 통과분만 |
| `validate_kq_record(rec) -> None\|raises` | data-model.md 검증 규칙 1~4 적용, 위반 시 명확한 오류 |

## 검증 실패 = 차단(결정론)

- 필수 필드 누락 / `pico.P`·`pico.I` 빈 리스트 / 검색식 필드(`query`) 존재 / enum 위반
  → 로드 거부(명확한 메시지). 잘못된 KQ 가 파이프라인에 진입하지 못하게 조기 차단.

## 단일 출처

- 검색식은 저장하지 않고 `build_query(pico)` 로 파생(헌법 III). 스키마에 `query` 필드 금지.
