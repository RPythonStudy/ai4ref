# 논문 중복관리 및 PostgreSQL 운영 가이드

## 1. 논문 중복검색 전략

PubMed 논문 수집 시 다음 3단계 순서로 중복 여부를 체크합니다:

1. **DOI 기준**: 존재할 경우, DB 내 DOI 중복 시 "이미 수집됨"으로 간주
2. **PMID 기준**: DOI가 없으면 PMID 기준 중복 체크  
3. **Title + FirstAuthor + Year**: DOI, PMID 모두 없을 경우 정규화된 3항목이 모두 일치하면 중복

중복체크는 PostgreSQL에 저장된 논문 테이블에서 쿼리로 수행하며, 논문 저장 전 반드시 선행합니다.

### 중복 체크용 SQL 예시
```sql
-- 1. DOI 기준
SELECT id FROM papers WHERE doi = :new_doi;

-- 2. PMID 기준  
SELECT id FROM papers WHERE pmid = :new_pmid;

-- 3. Title + FirstAuthor + Year
SELECT id FROM papers 
WHERE LOWER(TRIM(title)) = LOWER(TRIM(:new_title))
  AND LOWER(TRIM(first_author)) = LOWER(TRIM(:new_first_author))
  AND year = :new_year;
```

## 2. PostgreSQL 도커 운영

### docker-compose.yml
```yaml
version: '3'
services:
  postgres:
    image: postgres:16
    container_name: ai4ref_postgres
    restart: always
    environment:
      POSTGRES_USER: ai4ref
      POSTGRES_PASSWORD: ai4ref_pw
      POSTGRES_DB: ai4ref_db
    ports:
      - "5432:5432"
    volumes:
      - ./db_data:/var/lib/postgresql/data
```

### DB 접속정보
- host: localhost
- port: 5432  
- db: ai4ref_db
- user: ai4ref
- pw: ai4ref_pw

## 3. 개발 지침

### 필수 사항
- DOI/PMID/Title+Author+Year 인덱스 설정 필수
- Collector/Preprocessor 모든 코드에서 3단계 중복체크 적용
- 도커 볼륨(./db_data) 정기 백업 권장

### Python 구현 예시
```python
from common.logger import log_info, log_debug
import psycopg2

def check_duplicate_paper(doi, pmid, title, first_author, year):
    log_debug("논문 중복 체크 시작")
    
    # 1. DOI 체크
    if doi:
        # SQL 쿼리 실행
        log_debug(f"DOI 중복 체크: {doi}")
        
    # 2. PMID 체크  
    if pmid:
        log_debug(f"PMID 중복 체크: {pmid}")
        
    # 3. Title+Author+Year 체크
    log_debug(f"제목/저자/년도 중복 체크")
    
    return False  # 중복 없음
```
