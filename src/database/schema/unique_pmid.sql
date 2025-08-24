CREATE TABLE unique_pmid (
    pmid TEXT PRIMARY KEY,
    
    -- 컬렉션 추적
    first_collection_id INTEGER REFERENCES collections(id), -- 최초 발견 컬렉션
    occurrence_count INTEGER DEFAULT 1,    -- 여러 컬렉션에서 발견 횟수
    
    -- 처리 상태
    processed BOOLEAN DEFAULT FALSE,       -- pmid_filter 처리 완료 여부
    processed_at TIMESTAMP,               -- 처리 완료 시각
    
    -- 메타데이터
    first_found_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW()
);

-- 성능 최적화 인덱스  
CREATE INDEX idx_unique_pmid_processed ON unique_pmid(processed, first_found_at);
CREATE INDEX idx_unique_pmid_collection ON unique_pmid(first_collection_id);
CREATE INDEX idx_unique_pmid_occurrence ON unique_pmid(occurrence_count) WHERE occurrence_count > 1;
