# 업계 표준 기반 데이터베이스 설계 제안

## 1. 핵심 설계 원칙

### A. 학술 문헌 관리 표준 준수
- **Dublin Core**: title, creator, subject, description, date, identifier
- **CrossRef/DOI**: DOI 중심 중복 관리
- **PubMed NLM**: PMID, PMC, MeSH 용어 표준
- **JATS XML**: 구조화된 논문 메타데이터

### B. 워크플로우 기반 테이블 설계
- **단계별 데이터 추적**: collection → search → filter → enrich
- **중복 방지**: 각 단계별 UNIQUE 제약조건
- **상태 관리**: processing_status, retry_count
- **감사 추적**: created_at, updated_at, processed_at

## 2. 제안하는 테이블 구조

### A. Core Tables (핵심 엔터티)

#### collections (검색 컬렉션)
```sql
CREATE TABLE collections (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,              -- 컬렉션 명
    search_term TEXT NOT NULL,              -- PubMed 검색식
    retmax INTEGER DEFAULT 10000,           -- 최대 검색 결과
    enabled BOOLEAN DEFAULT TRUE,           -- 활성화 여부
    
    -- Zotero 동기화
    zotero_key VARCHAR(12) UNIQUE,          -- Zotero collection key
    zotero_version INTEGER,                 -- API 버전 충돌 방지
    
    -- 메타데이터
    description TEXT,
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW(),
    
    -- 검색 실행 추적
    last_search_at TIMESTAMP,               -- 마지막 검색 시각
    last_search_count INTEGER               -- 마지막 검색 결과 수
);
```

#### papers (논문 마스터)
```sql
CREATE TABLE papers (
    id SERIAL PRIMARY KEY,
    
    -- 다중 소스 지원 (PubMed, arXiv, CrossRef, etc.)
    source VARCHAR(50) NOT NULL,           -- 'pubmed', 'arxiv', 'crossref'
    source_id TEXT NOT NULL,               -- 소스별 ID (PMID, arXiv ID 등)
    
    -- 표준 식별자 (Dublin Core identifier)
    doi TEXT,                              -- Digital Object Identifier
    pmid TEXT,                             -- PubMed ID
    pmcid TEXT,                            -- PMC ID
    arxiv_id TEXT,                         -- arXiv ID
    
    -- 핵심 메타데이터 (Dublin Core 기반)
    title TEXT,                            -- Dublin Core: title
    abstract TEXT,                         -- Dublin Core: description
    authors TEXT,                          -- Dublin Core: creator
    publication_year INTEGER,              -- Dublin Core: date
    journal TEXT,                          -- 저널명
    volume TEXT,                           -- 권
    issue TEXT,                            -- 호
    pages TEXT,                            -- 페이지
    
    -- 주제 분류
    mesh_terms TEXT[],                     -- MeSH 용어 배열
    keywords TEXT[],                       -- 저자 키워드
    
    -- 무료 접근 정보
    is_open_access BOOLEAN,                -- 오픈 액세스 여부
    pdf_url TEXT,                          -- 무료 PDF URL
    pdf_source TEXT,                       -- PDF 출처
    
    -- 품질 지표
    citation_count INTEGER,                -- 인용 횟수
    impact_factor DECIMAL(5,3),            -- 저널 IF
    
    -- 시스템 필드
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW(),
    
    -- 중복 방지 제약조건
    UNIQUE (source, source_id),
    UNIQUE (doi) WHERE doi IS NOT NULL AND doi != ''
);
```

### B. Workflow Tables (워크플로우 단계별)

#### collection_pmids (Step 2: pubmed_esearch 결과)
```sql
CREATE TABLE collection_pmids (
    id SERIAL PRIMARY KEY,
    collection_id INTEGER NOT NULL REFERENCES collections(id) ON DELETE CASCADE,
    pmid TEXT NOT NULL,
    
    -- 검색 메타데이터
    search_term TEXT,                      -- 실제 사용된 검색어
    search_date DATE DEFAULT CURRENT_DATE, -- 검색 날짜
    result_position INTEGER,               -- 검색 결과 내 순서
    esearch_count INTEGER,                 -- 전체 검색 결과 수
    
    created_at TIMESTAMP DEFAULT NOW(),
    
    UNIQUE(collection_id, pmid)            -- 컬렉션별 PMID 중복 방지
);
```

#### unique_pmids (Step 3a: 컬렉션 간 중복 제거)
```sql
CREATE TABLE unique_pmids (
    pmid TEXT PRIMARY KEY,
    first_collection_id INTEGER REFERENCES collections(id),
    first_found_at TIMESTAMP DEFAULT NOW(),
    occurrence_count INTEGER DEFAULT 1,    -- 여러 컬렉션에서 발견 횟수
    
    -- 처리 상태
    processed BOOLEAN DEFAULT FALSE,
    processed_at TIMESTAMP
);
```

#### filtered_pmids (Step 3b: papers 테이블 대비 신규만)
```sql
CREATE TABLE filtered_pmids (
    pmid TEXT PRIMARY KEY,
    ready_for_fetch BOOLEAN DEFAULT TRUE,
    priority INTEGER DEFAULT 0,            -- 우선순위 (최신 논문 우선 등)
    
    -- 필터링 메타데이터
    filter_reason TEXT,                     -- 필터링 사유
    estimated_year INTEGER,                 -- 추정 연도 (정렬용)
    
    created_at TIMESTAMP DEFAULT NOW()
);
```

### C. Relationship Tables (관계)

#### paper_collections (논문-컬렉션 다대다 관계)
```sql
CREATE TABLE paper_collections (
    paper_id INTEGER NOT NULL REFERENCES papers(id) ON DELETE CASCADE,
    collection_id INTEGER NOT NULL REFERENCES collections(id) ON DELETE CASCADE,
    
    -- 연결 메타데이터
    relevance_score DECIMAL(3,2),          -- 관련도 점수 (0.00-1.00)
    search_term TEXT,                      -- 해당 논문을 찾은 검색어
    added_at TIMESTAMP DEFAULT NOW(),
    
    PRIMARY KEY (paper_id, collection_id)
);
```

#### search_history (검색 이력 추적)
```sql
CREATE TABLE search_history (
    id SERIAL PRIMARY KEY,
    collection_id INTEGER REFERENCES collections(id),
    search_term TEXT NOT NULL,
    search_date DATE DEFAULT CURRENT_DATE,
    result_count INTEGER,
    api_response_time INTEGER,             -- API 응답시간 (ms)
    status VARCHAR(20) DEFAULT 'success',  -- success, failed, timeout
    error_message TEXT,
    created_at TIMESTAMP DEFAULT NOW()
);
```

## 3. 중복 관리 전략

### A. 3단계 중복 체크 (업계 표준)
```sql
-- 1. DOI 기준 (최우선)
CREATE UNIQUE INDEX idx_papers_doi ON papers(doi) 
WHERE doi IS NOT NULL AND doi != '';

-- 2. PMID 기준
CREATE UNIQUE INDEX idx_papers_pmid ON papers(pmid) 
WHERE pmid IS NOT NULL AND pmid != '';

-- 3. Title+FirstAuthor+Year 기준 (정규화 함수 필요)
CREATE INDEX idx_papers_fuzzy ON papers(
    LOWER(TRIM(title)), 
    SPLIT_PART(authors, ';', 1), 
    publication_year
) WHERE title IS NOT NULL;
```

### B. 프로세스 상태 추적
```sql
-- workflow_status 테이블로 전체 진행상황 모니터링
CREATE TABLE workflow_status (
    id SERIAL PRIMARY KEY,
    batch_id UUID DEFAULT gen_random_uuid(),
    step_name VARCHAR(50) NOT NULL,        -- 'esearch', 'filter', 'efetch'
    collection_id INTEGER REFERENCES collections(id),
    
    -- 통계
    total_items INTEGER,
    processed_items INTEGER DEFAULT 0,
    success_items INTEGER DEFAULT 0,
    failed_items INTEGER DEFAULT 0,
    
    -- 시간 추적
    started_at TIMESTAMP DEFAULT NOW(),
    completed_at TIMESTAMP,
    estimated_duration INTERVAL,
    
    status VARCHAR(20) DEFAULT 'running'   -- running, completed, failed
);
```

## 4. 성능 최적화 인덱스

```sql
-- 워크플로우 조회 최적화
CREATE INDEX idx_collection_pmids_search ON collection_pmids(search_date, collection_id);
CREATE INDEX idx_unique_pmids_processed ON unique_pmids(processed, first_found_at);
CREATE INDEX idx_filtered_pmids_priority ON filtered_pmids(priority DESC, created_at);

-- 논문 검색 최적화  
CREATE INDEX idx_papers_source_year ON papers(source, publication_year);
CREATE INDEX idx_papers_updated ON papers(updated_at) WHERE title IS NULL;
CREATE INDEX idx_papers_open_access ON papers(is_open_access) WHERE is_open_access = true;

-- 전문검색 지원 (PostgreSQL FTS)
CREATE INDEX idx_papers_fts ON papers USING gin(
    to_tsvector('english', coalesce(title,'') || ' ' || coalesce(abstract,''))
);
```

## 5. 데이터 품질 관리

### A. 제약조건 및 체크
```sql
-- 논문 품질 체크
ALTER TABLE papers ADD CONSTRAINT chk_publication_year 
CHECK (publication_year BETWEEN 1800 AND EXTRACT(YEAR FROM NOW()) + 2);

ALTER TABLE papers ADD CONSTRAINT chk_doi_format 
CHECK (doi ~ '^10\.\d{4,}/.*' OR doi IS NULL);

-- 워크플로우 무결성
ALTER TABLE collection_pmids ADD CONSTRAINT chk_search_date 
CHECK (search_date <= CURRENT_DATE);

ALTER TABLE paper_collections ADD CONSTRAINT chk_relevance_score 
CHECK (relevance_score BETWEEN 0.00 AND 1.00);
```

### B. 자동 정리 정책
```sql
-- 90일 이상 된 임시 테이블 정리
CREATE OR REPLACE FUNCTION cleanup_old_workflow_data() RETURNS void AS $$
BEGIN
    DELETE FROM collection_pmids WHERE created_at < NOW() - INTERVAL '90 days';
    DELETE FROM unique_pmids WHERE processed = true AND processed_at < NOW() - INTERVAL '30 days';
    DELETE FROM filtered_pmids WHERE created_at < NOW() - INTERVAL '7 days';
END;
$$ LANGUAGE plpgsql;
```

## 6. 마이그레이션 계획

### 현재 구조에서 제안 구조로 전환
1. **Phase 1**: 기존 esearch_queue → collection_pmids 변환
2. **Phase 2**: unique_pmids, filtered_pmids 테이블 추가
3. **Phase 3**: papers 테이블 필드 확장 (journal, mesh_terms 등)
4. **Phase 4**: 성능 인덱스 및 제약조건 적용

이 설계는 학술 문헌 관리의 업계 표준을 따르면서도 현재 워크플로우의 요구사항을 충족합니다.
