CREATE TABLE papers (
    id SERIAL PRIMARY KEY,

    -- 논문 수집 소스 및 고유 식별자
    source VARCHAR(50) NOT NULL,     -- 예: 'pubmed', 'arxiv', 'crossref', 'europepmc', 'semanticscholar'
    source_id TEXT NOT NULL,         -- 해당 소스의 논문ID(PMID, arXiv_id, DOI 등)
    doi TEXT,
    pmid TEXT,
    pmcid TEXT,
    arxiv_id TEXT,

    -- 논문 메타데이터
    title TEXT,
    abstract TEXT,
    publication_year INTEGER,
    authors TEXT,
    journal TEXT,                  -- 저널명 (필수 메타데이터)

    -- 무료 PDF 정보
    pdf_url TEXT,                 -- 대표 PDF 파일(무료/오픈액세스)
    pdf_source TEXT,              -- PDF의 출처(예: 'Europe PMC', 'arXiv', 'PMC', 'Unpaywall' 등)
    pdf_metadata JSONB,           -- 복수 PDF 정보(필요시, [{"source":"PMC","url":"..."}, ...])
    pdf_path VARCHAR(512),         -- PDF 파일의 로컬 경로

    -- 생성/수정 정보
    uploaded BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW(),

    -- 논문 유일성 (중복방지)
    UNIQUE (source, source_id)
);

-- 논문 식별자 인덱스
CREATE INDEX idx_papers_doi         ON papers(doi);
CREATE INDEX idx_papers_pmid        ON papers(pmid);
CREATE INDEX idx_papers_pmcid       ON papers(pmcid);
CREATE INDEX idx_papers_arxiv_id    ON papers(arxiv_id);

-- 무료 PDF 관련 인덱스
CREATE INDEX idx_papers_pdf_url     ON papers(pdf_url);
CREATE INDEX idx_papers_pdf_source  ON papers(pdf_source);

-- 복합 검색/중복검사용(선택)
CREATE INDEX idx_papers_title_year
    ON papers(title, publication_year);

-- PDF 메타데이터(JSONB) 인덱스(필요시)
-- CREATE INDEX idx_papers_pdf_metadata_gin ON papers USING GIN (pdf_metadata);
