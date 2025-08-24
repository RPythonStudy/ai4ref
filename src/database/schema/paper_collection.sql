CREATE TABLE paper_collection (
    paper_id INTEGER NOT NULL REFERENCES papers(id) ON DELETE CASCADE,
    collection_id INTEGER NOT NULL REFERENCES collections(id) ON DELETE CASCADE,
    
    -- 연결 메타데이터
    relevance_score DECIMAL(3,2),          -- 관련도 점수 (0.00-1.00)
    collection_name TEXT,                  -- 
    added_at TIMESTAMP DEFAULT NOW(),
    
    -- Zotero 동기화 정보
    zotero_item_key VARCHAR(12),           -- Zotero 아이템 키
    zotero_version INTEGER,                -- Zotero API 버전
    synced_at TIMESTAMP,                   -- Zotero 동기화 시각
    
    PRIMARY KEY (paper_id, collection_id)
);

-- 성능 최적화 인덱스
CREATE INDEX idx_paper_collection_relevance ON paper_collection(relevance_score DESC);
CREATE INDEX idx_paper_collection_collection_name ON paper_collection(collection_name);
CREATE INDEX idx_paper_collection_zotero ON paper_collection(zotero_item_key) WHERE zotero_item_key IS NOT NULL;
