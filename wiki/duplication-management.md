# 논문 중복관리 및 PostgreSQL 운영 가이드

## 1. 논문 중복검색 전략

PubMed 논문 수집 시 다음 3단계 순서로 중복 여부를 체크합니다:

1. DOI 기준: 존재할 경우, DB 내 DOI 중복 시 "이미 수집됨"으로 간주
2. PMID 기준: DOI가 없으면 PMID 기준 중복 체크  
3. Title + FirstAuthor + Year: DOI, PMID 모두 없을 경우 정규화된 3항목이 모두 일치하면 중복

중복체크는 PostgreSQL에 저장된 논문 테이블에서 쿼리로 수행하며, 논문 저장 전 반드시 선행합니다.



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

