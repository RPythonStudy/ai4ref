# 논문 중복관리 가이드

## 개요

AI4REF 프로젝트에서 논문 수집 시 중복을 효과적으로 관리하기 위한 3단계 중복점검 시스템 가이드입니다.

## 데이터베이스 스키마

### papers 테이블
- 논문의 메타데이터 저장
- `UNIQUE(source, source_id)`: 기본 중복 방지
- 여러 식별자 필드: `doi`, `pmid`, `pmcid`, `arxiv_id`

### paper_collection_log 테이블
- 논문 수집 이력 추적
- `UNIQUE(paper_id, search_term, collection_name)`: 수집 로그 중복 방지

## 3단계 중복점검 알고리즘

### 1단계: 기본 식별자 중복 체크
**목적**: 동일한 논문의 재수집 방지

**방법**: 
```sql
SELECT id FROM papers 
WHERE source = 'pubmed' AND source_id = %s
```

**결과**:
- 존재하면: 기존 `paper_id` 반환, 3단계로 진행
- 없으면: 2단계로 진행

### 2단계: 교차 식별자 중복 체크
**목적**: 같은 논문이 다른 식별자로 수집되는 것 방지

**방법**: efetch로 수집한 메타데이터의 식별자들을 교차 확인
```sql
-- DOI 확인
SELECT id FROM papers WHERE doi = %s AND doi IS NOT NULL

-- PMCID 확인  
SELECT id FROM papers WHERE pmcid = %s AND pmcid IS NOT NULL

-- Title + Year 조합 (선택사항)
SELECT id FROM papers 
WHERE title = %s AND publication_year = %s
```

**결과**:
- 발견되면: 기존 논문의 식별자 정보 업데이트 후 `paper_id` 반환
- 없으면: 새 논문으로 papers 테이블에 삽입

### 3단계: 수집 로그 중복 체크
**목적**: 동일한 검색어/컬렉션으로 이미 수집된 논문의 재로깅 방지

**방법**:
```sql
SELECT COUNT(*) FROM paper_collection_log 
WHERE paper_id = %s AND search_term = %s AND collection_name = %s
```

**결과**:
- 존재하면: 로그 건너뜀 (`already_logged` 카운트)
- 없으면: 새 로그 레코드 생성 (`inserted` 카운트)

## 메타데이터 일관성 관리

### 우선순위 정책
```python
metadata_priority = {
    'title': 'longest_non_empty',      # 가장 긴 제목 선택
    'abstract': 'longest_non_empty',   # 가장 긴 초록 선택
    'doi': 'first_valid',              # 첫 번째 유효한 DOI 유지
    'publication_year': 'most_reliable', # 신뢰할 만한 소스 우선
    'authors': 'most_complete'         # 가장 완전한 저자 정보
}
```

### 업데이트 전략
- **즉시 업데이트**: 새로운 메타데이터가 더 완전하면 즉시 반영
- **필드별 개선**: 빈 필드는 채우고, 더 나은 정보로 교체
- **이력 추적**: `updated_at` 필드로 마지막 업데이트 시점 기록

### 동시성 제어
- 트랜잭션 레벨에서 중복 체크와 삽입을 원자적으로 처리
- 필요시 SELECT FOR UPDATE로 락 사용

## API 호출 최적화

### efetch 호출 최소화
1. 1단계에서 기존 논문 발견 시 efetch 생략
2. 신규 PMID만 일괄 efetch 처리
3. API 요청 제한 준수 (초당 3회, API 키 있으면 10회)

### 배치 처리
```python
# 권장 구조
new_pmids = filter_existing_pmids(pmids)  # 1단계 필터링
if new_pmids:
    paper_details = efetch_papers(new_pmids)  # 일괄 efetch
    for pmid, details in paper_details.items():
        paper_id = save_or_update_paper(details)  # 2단계 처리
        log_collection(paper_id, search_term, collection_name)  # 3단계 처리
```

## 오류 처리

### 부분 실패 허용
- 일부 PMID의 efetch 실패가 전체 프로세스를 중단하지 않도록 설계
- 실패한 PMID는 별도 로그에 기록하여 재시도 가능

### 데이터 검증
- DOI 형식 검증
- 연도 범위 검증 (1900-현재)
- 필수 필드 존재 여부 확인

## 성능 고려사항

### 인덱스 최적화
```sql
-- 중복 체크용 인덱스
CREATE INDEX idx_papers_pmid ON papers(pmid);
CREATE INDEX idx_papers_doi ON papers(doi);
CREATE INDEX idx_papers_title_year ON papers(title, publication_year);

-- 로그 중복 체크용 인덱스
CREATE INDEX idx_pcl_lookup ON paper_collection_log(paper_id, search_term, collection_name);
```

### 대용량 처리
- PMID 수가 1000개 이상일 때는 청크 단위로 분할 처리
- 메모리 사용량 모니터링

## 모니터링 지표

### 수집 통계
- `total_collected`: 검색으로 수집된 전체 PMID 수
- `total_processed`: 새로 처리된 논문 수 (신규 + 로그 기록)
- `total_logged`: 중복으로 건너뛴 로그 수
- `api_calls_saved`: efetch 호출 절약 수

### 품질 지표
- 중복률: `already_logged / total_collected`
- API 효율성: `api_calls_saved / total_collected`
- 메타데이터 완성도: 필수 필드 채워진 비율

## 구현 체크리스트

- [ ] 3단계 중복점검 로직 구현
- [ ] 메타데이터 병합 룰 적용
- [ ] 트랜잭션 기반 동시성 제어
- [ ] API 호출 최적화
- [ ] 상세한 로깅 및 통계
- [ ] 오류 복구 메커니즘
- [ ] 성능 모니터링
