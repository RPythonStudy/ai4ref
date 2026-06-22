# Module Specification: common.pubmed (L3)

**File**: [pubmed.py](file:///home/ben/projects/ai4ref/src/common/pubmed.py) | **Spec Version**: 1.0 | **Code Lines Constraint**: <= 300줄 | **Spec Lines Constraint**: <= 200줄

## Summary

PubMed(E-utilities)와 상호작용하는 경량 클라이언트 기능을 정의합니다. 주로 esearch와 efetch_meta API 호출을 래핑하며, 무비용 로컬 운영 정책에 맞추어 API 키 유무에 따른 속도 제한 조율 및 배치 단위 예외 복구를 처리합니다.

## Configuration (Environment Variables)

* `ENTREZ_EMAIL`: 필수/권장 연락용 이메일 (NCBI 정책용). 기본값 `"r.python.ai@gmail.com"`.
* `PUBMED_API_KEY`: 선택사항. 있으면 초당 요청 제한 완화. `"your-"`로 시작하는 플레이스홀더 값은 빈 값으로 처리하여 가드함.

## Rate Limiting Policy

* API 키가 없을 때: 요청 간 간격 `_SLEEP = 0.34`초 (초당 3회 요청 이하 가이드라인 충족).
* API 키가 있을 때: 요청 간 간격 `_SLEEP = 0.11`초 (초당 9회 요청 수준 충족).

## API Specifications

### `esearch`
```python
def esearch(term: str, retmax: int = 300, datetype: str = None,
            reldate: int = None, mindate: str = None, maxdate: str = None) -> list[str]:
```
* **역할**: PubMed에서 특정 검색어를 만족하는 논문의 PMID 목록을 획득합니다.
* **매개변수**:
  * `term`: 검색 쿼리 문자열.
  * `retmax`: 최대 반환 결과 수 (기본 300).
  * `datetype`: 시간 필터링 기준 (예: `"edat"` - 색인일, `"pdat"` - 출판일).
  * `reldate`: 최근 N일 이내 필터링 창 크기 (정수).
  * `mindate`, `maxdate`: 특정 날짜 범위 필터링 (`"YYYY/MM/DD"` 형식).
* **출력**: 매칭된 PMID 문자열의 리스트.
* **설계**: `retmode=json` 방식을 사용하며 결과가 없을 경우 빈 리스트를 반환합니다.

### `efetch_meta`
```python
def efetch_meta(pmids: list[str]) -> dict[str, dict]:
```
* **역할**: 지정된 PMID 목록에 대해 PubMed 메타데이터(제목, 초록, 저널, 출판연도)를 일괄 조회합니다.
* **매개변수**:
  * `pmids`: 조회를 수행할 PMID 문자열 리스트.
* **출력**: `pmid -> {"title": ..., "abstract": ..., "journal": ..., "year": ...}` 구조의 딕셔너리.
* **에러 핸들링 & fail-soft**:
  * 입력받은 PMID들을 **최대 100개 단위의 배치**로 분할하여 개별 XML efetch API 요청을 전송합니다.
  * 특정 배치의 XML 파싱 에러나 HTTP 연결 실패 시, **해당 배치만 스킵하고 오류 로그를 남긴 후 다음 배치 조회를 정상 진행**합니다 (FR-004).
