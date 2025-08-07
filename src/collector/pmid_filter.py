# pmid_filter.py - PMID 처리 전용 모듈
import psycopg2
import os
from common.logger import log_info, log_error

DB_HOST = os.environ.get("DB_HOST", "localhost")
DB_PORT = os.environ.get("DB_PORT", 5432)
DB_NAME = os.environ.get("POSTGRES_DB", "ai4ref")
DB_USER = os.environ.get("POSTGRES_USER", "postgres")
DB_PASS = os.environ.get("POSTGRES_PASSWORD", "postgres")

def get_connection():
    """데이터베이스 연결 반환"""
    return psycopg2.connect(
        host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
    )

def filter_new_pmids(pmids):
    """
    기존 papers 테이블에서 중복되지 않은 새로운 PMID만 필터링
    """
    try:
        conn = get_connection()
        cur = conn.cursor()
        
        new_pmids = []
        duplicate_count = 0
        
        for pmid in pmids:
            # papers 테이블에서 중복 확인
            cur.execute(
                "SELECT 1 FROM papers WHERE source = 'pubmed' AND source_id = %s",
                (pmid,)
            )
            
            if cur.fetchone():
                duplicate_count += 1
            else:
                new_pmids.append(pmid)
        
        cur.close()
        conn.close()
        
        log_info(f"필터링 완료: 신규 {len(new_pmids)}개, 중복 {duplicate_count}개")
        return {"new_pmids": new_pmids, "duplicates": duplicate_count, "total_checked": len(pmids)}
        
    except Exception as e:
        log_error(f"PMID 필터링 중 오류: {e}")
        return None

def save_pmids_to_papers(pmids, search_term=None, collection_name=None):
    """
    새로운 PMID들을 papers 테이블에 저장 (identifier_only 상태로)
    검색어와 컬렉션 정보도 함께 저장
    """
    try:
        conn = get_connection()
        cur = conn.cursor()
        
        inserted_count = 0
        
        for pmid in pmids:
            try:
                cur.execute("""
                    INSERT INTO papers (source, source_id, search_term, collection_name, created_at) 
                    VALUES ('pubmed', %s, %s, %s, NOW()) 
                    ON CONFLICT (source, source_id) DO NOTHING
                """, (pmid, search_term, collection_name))
                
                if cur.rowcount > 0:
                    inserted_count += 1
                    
            except Exception as e:
                log_error(f"PMID {pmid} 저장 오류: {e}")
        
        conn.commit()
        cur.close()
        conn.close()
        
        log_info(f"PMID 저장 완료: {inserted_count}개")
        return {"inserted": inserted_count}
        
    except Exception as e:
        log_error(f"PMID 저장 중 오류: {e}")
        return None
