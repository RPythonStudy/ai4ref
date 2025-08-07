# 논문 중복관리 및 PostgreSQL 운영 가이드

## 1. 논문 수집 워크플로우 (필수)

### 단계별 프로세스
```
esearch (검색어 → PMID 목록) 
    ↓
중복 체크 (DB 기존 PMID와 비교)
    ↓
efetch (신규 PMID만 상세정보 수집)
    ↓
DB 저장 (논문 상세정보)
```

### 중복 체크 전략
**efetch 전 필수 체크**: 
- `pubmed_pmids` 테이블에서 기존 PMID 확인
- 중복된 PMID 제외하고 신규 PMID만 efetch 진행
- API 호출 횟수 최소화로 효율성 극대화

**efetch 후 상세 체크**:
1. DOI 기준: 존재할 경우, DB 내 DOI 중복 시 "이미 수집됨"으로 간주
2. PMID 기준: DOI가 없으면 PMID 기준 중복 체크  
3. Title + FirstAuthor + Year: DOI, PMID 모두 없을 경우 정규화된 3항목이 모두 일치하면 중복

**핵심 원칙**: 항상 중복 체크 → efetch 순서로 진행하여 불필요한 API 호출 방지



## 2. PostgreSQL 도커 운영

### docker-compose.yml
```yaml

services:
  postgres:
    image: postgres:16.3-alpine
    container_name: postgres
    restart: unless-stopped
    user: "1000:1000"
    env_file:
      - .env
    volumes:
      - /opt/docker/postgres/data:/var/lib/postgresql/data
    ports:
      - "5432:5432"
    networks:
      - common-infra
networks:
  common-infra:
    driver: bridge
    name: common-infra
```

### DB 접속정보
- host: localhost
- port: 5432  
- db: ai4ref
- user: postgres
- pw: postgres

## 3. 개발 지침

### 필수 사항
- DOI/PMID/Title+Author+Year 인덱스 설정 필수
- Collector/Preprocessor 모든 코드에서 3단계 중복체크 적용
- 도커 볼륨(./db_data) 정기 백업 권장

