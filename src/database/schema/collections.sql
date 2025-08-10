CREATE TABLE collections (
    id SERIAL PRIMARY KEY,                -- 내부 DB 식별자
    name TEXT NOT NULL,                   -- 컬렉션 이름(필수)
    parent_id INTEGER REFERENCES collections(id) ON DELETE CASCADE,  -- 상위 컬렉션(트리구조)
    
    -- 검색 설정 정보
    search_term TEXT,                     -- PubMed 검색어 (search_collections.json의 term)
    retmax INTEGER DEFAULT 10000,        -- 최대 검색 결과 수
    enabled BOOLEAN DEFAULT TRUE,        -- 활성화 여부
    
    -- Zotero 동기화 정보
    zotero_key VARCHAR(12) UNIQUE,        -- Zotero의 collection key (동기화, 업로드 후 저장)
    zotero_version INTEGER,               -- Zotero API version (동기화 충돌 방지용, 선택)
    
    -- 메타데이터
    description TEXT,                     -- 컬렉션 설명 (선택)
    created_at TIMESTAMP DEFAULT NOW(),   -- 생성일(선택)
    updated_at TIMESTAMP DEFAULT NOW()    -- 수정일(선택)
);
