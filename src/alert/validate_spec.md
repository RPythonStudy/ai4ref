# Module Specification: alert.validate_kq (L3)

**File**: [validate_kq.py](file:///home/ben/projects/ai4ref/src/alert/validate_kq.py) | **Spec Version**: 1.0 | **Code Lines Constraint**: <= 300줄 | **Spec Lines Constraint**: <= 200줄

## Summary

PICO 검색식의 상대 재현율(Recall)을 검증 A(지침 이전 근거)와 검증 B(지침 이후 랜드마크) 기준으로 계산하고 채점합니다. 또한, 질문의 의학적 유형과 분류기 판단 결과를 대조하는 교차 검증을 제공합니다.

## API Specifications

### `validate_single_kq`
```python
def validate_single_kq(kq: dict) -> dict:
```
* **역할**: 개별 KQ 레코드 1건을 전달받아 결정론적 유형 대조 및 Recall A/B 검색 여부를 채점합니다.
* **매개변수**:
  * `kq`: `common.kq.load_kqs`를 통해 로드 및 검증된 KQ 레코드 딕셔너리.
* **출력**:
  ```python
  {
      "kq": str,
      "type_match": bool,
      "user_type": str,
      "suggested_type": str,
      "refs_results": list[tuple[str, bool]],  # (pmid, ok_flag) 목록
      "refs_hit": int,
      "refs_total": int,
      "landmarks_results": list[tuple[str, bool]],  # (pmid, ok_flag) 목록
      "landmarks_hit": int,
      "landmarks_total": int,
  }
  ```
* **PubMed 매칭 판단 (FR-006)**:
  * `common.pubmed.esearch`를 재사용하여 검색식과 대상 PMID 조건을 결합한 쿼리를 던져 매칭 여부를 판단합니다.
  * `ok = len(esearch(term, retmax=1, ...)) >= 1`
  * 검증 A(지침 이전): `datetype="pdat", mindate="1900/01/01", maxdate="T년/12/31"`로 기간을 제한합니다 (지침 연도가 있는 경우).
  * 검증 B(지침 이후): `datetype="pdat", mindate="T+1년/01/01", maxdate="2030/12/31"`로 기간을 제한합니다.

## Validation Triad Requirements (원칙 XIII)

* ERAS 수액전략 KQ에 대해 `validate_single_kq` 수행 결과가 다음과 같은 골든 레벨을 만족해야 합니다.
  * 검증 A (지침 근거 재현): 8 / 9
  * 검증 B (랜드마크 포착): 2 / 2
