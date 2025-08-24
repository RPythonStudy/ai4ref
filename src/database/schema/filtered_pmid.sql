CREATE TABLE filtered_pmid (
    pmid TEXT PRIMARY KEY,
    
    -- 필터링 결과
    ready_for_fetch BOOLEAN DEFAULT TRUE,  -- efetch 준비 완료 여부
    priority INTEGER DEFAULT 0,            -- 우선순위 (최신 논문 우선 등)
    
    -- 필터링 메타데이터
    filter_reason TEXT,                     -- 필터링 사유 (예: "신규", "기존 논문 존재")
    estimated_year INTEGER,                 -- 추정 연도 (정렬 및 우선순위용)
    
    -- 처리 상태
    fetched BOOLEAN DEFAULT FALSE,          -- efetch 완료 여부
    fetched_at TIMESTAMP,                   -- efetch 완료 시각
    fetch_attempts INTEGER DEFAULT 0,      -- efetch 시도 횟수
    last_error TEXT,                        -- 마지막 오류 메시지
    
    -- 시스템 필드
    created_at TIMESTAMP DEFAULT NOW()
);

-- 성능 최적화 인덱스
CREATE INDEX idx_filtered_pmid_ready ON filtered_pmid(ready_for_fetch, priority DESC);
CREATE INDEX idx_filtered_pmid_fetched ON filtered_pmid(fetched, created_at);
CREATE INDEX idx_filtered_pmid_priority ON filtered_pmid(priority DESC, estimated_year DESC);
CREATE INDEX idx_filtered_pmid_attempts ON filtered_pmid(fetch_attempts) WHERE fetch_attempts > 0;
