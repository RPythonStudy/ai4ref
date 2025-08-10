CREATE TABLE collection_pmid (
    id SERIAL PRIMARY KEY,
    collection_id INTEGER NOT NULL REFERENCES collections(id) ON DELETE CASCADE,
    pmid TEXT NOT NULL,
    
    -- 검색 메타데이터
    collection_name TEXT,                  -- 실제 사용된 검색어
    search_date DATE DEFAULT CURRENT_DATE, -- 검색 날짜
    result_position INTEGER,               -- 검색 결과 내 순서
    esearch_count INTEGER,                 -- 전체 검색 결과 수
    batch_id TEXT,                         -- 동일 검색의 배치 식별자
    
    -- 시스템 필드
    created_at TIMESTAMP DEFAULT NOW(),
    
    -- 중복 방지
    UNIQUE(collection_id, pmid)            -- 컬렉션별 PMID 중복 방지
);

-- 성능 최적화 인덱스
CREATE INDEX idx_collection_pmid_search_date ON collection_pmid(search_date, collection_id);
CREATE INDEX idx_collection_pmid_batch ON collection_pmid(batch_id);
CREATE INDEX idx_collection_pmid_pmid ON collection_pmid(pmid);
