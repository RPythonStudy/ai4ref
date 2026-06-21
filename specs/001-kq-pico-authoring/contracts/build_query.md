# Contract: build_query(pico)

`src/common/query.py` 의 검색식 파생 계약. (구현 상세는 docstring/type hint 일원화 — 헌법 XVI)

## 시그니처 (개념)

```
build_query(pico: dict) -> str
```

- **입력** `pico`: `{P: list[str], I: list[str], C?: str, O?: str}`. (C·O 는 무시)
- **출력**: PubMed 검색식 문자열 `(I 항목 OR …) AND (P 항목 OR …)`.

## 규칙

| # | 규칙 |
|---|---|
| C1 | I 블록 = `pico.I` 항목을 ` OR ` 로 결합, 괄호로 감쌈 |
| C2 | P 블록 = `pico.P` 항목을 ` OR ` 로 결합, 괄호로 감쌈 |
| C3 | 최종 = `(<I 블록>) AND (<P 블록>)` |
| C4 | `pico.C`·`pico.O` 는 검색식에 **포함하지 않음**(헌법 II) |
| C5 | 항목 문자열은 변형 없이 그대로 사용(예: `"Fluid Therapy"[Mesh]`) |
| C6 | 결정론 — 동일 입력 → 동일 출력(항목 순서 보존) |

## 등가 기준 (acceptance)

`tests/acceptance/recall_baseline.yml` 의 ERAS 2012 KQ 에 대해:

```
build_query(pico) == recall_baseline.derived_query
```

= `("Fluid Therapy"[Mesh] OR fluid therapy[tiab] OR … OR "central venous"[tiab])
   AND (surgery[tiab] OR surgical[tiab] OR perioperative[tiab] OR abdominal[tiab]
   OR "Surgical Procedures, Operative"[Mesh])`

## 예외

- `pico.P` 또는 `pico.I` 가 비어있으면 → 무효(검색식 파생 불가). 호출 측(kq 검증)이 차단.
