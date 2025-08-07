CREATE TABLE papers (
    -- 기본 식별자
    id SERIAL PRIMARY KEY,
    
    -- 논문 기본 정보
    title TEXT,
    abstract TEXT,
    publication_year INTEGER,
    doi TEXT UNIQUE,
    
    -- 저자 정보
    authors TEXT,                 -- "김철수; 이영희; 박민수" 형태
    first_author TEXT,
    
    -- 저널 정보
    journal_name TEXT,
    
    -- 수집 관련
    source VARCHAR(50) NOT NULL,  -- 'pubmed', 'arxiv' 등
    source_id TEXT NOT NULL,      -- pmid, arxiv_id 등
    
    -- 검색 추적 필드
    search_term TEXT,             -- 수집에 사용된 검색어
    collection_name TEXT,         -- 해당하는 Zotero 컬렉션명
    
    -- PDF 관리 필드
    pdf_url TEXT,                 -- PDF 다운로드 URL (PMC, 출판사 등)
    pdf_local_path TEXT,          -- 로컬에 저장된 PDF 파일 경로
    pdf_status VARCHAR(20),       -- 'pending', 'downloaded', 'failed', 'unavailable'
    last_pdf_check TIMESTAMP,     -- 마지막 PDF 확인/시도 날짜
    
    -- 확장 가능한 필드 (JSON으로 유연하게)
    extra_metadata JSONB,         -- 나중에 추가될 정보들
    
    -- 시스템 필드
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW(),
    
    -- 중복 방지 및 인덱스
    UNIQUE(source, source_id)
);

-- 인덱스 생성
CREATE INDEX idx_papers_search_term ON papers(search_term);
CREATE INDEX idx_papers_collection_name ON papers(collection_name);
CREATE INDEX idx_papers_pdf_status ON papers(pdf_status);
CREATE INDEX idx_papers_created_at ON papers(created_at);
CREATE INDEX idx_papers_source_id ON papers(source_id);
