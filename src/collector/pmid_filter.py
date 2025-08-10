import sys
import psycopg2
from common.logger import log_info, log_debug, log_error
from common.database import get_db_connection

def create_unique_pmid_from_collection_pmid():
    """collection_pmid에서 중복되지 않은 pmid를 unique_pmid 테이블에 저장합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # 기존 unique_pmid 테이블 초기화
        cur.execute("DELETE FROM unique_pmid")
        log_debug("기존 unique_pmid 테이블 초기화 완료")
        
        # collection_pmid에서 고유한 pmid 조회 및 삽입
        cur.execute("""
            INSERT INTO unique_pmid (pmid, first_collection_id, first_found_at)
            SELECT DISTINCT 
                cp.pmid,
                MIN(cp.collection_id) as first_collection_id,
                NOW() as first_found_at
            FROM collection_pmid cp
            GROUP BY cp.pmid
        """)
        
        unique_count = cur.rowcount
        conn.commit()
        
        log_info(f"unique_pmid 테이블 생성 완료: {unique_count}개")
        
        cur.close()
        conn.close()
        
        return unique_count
        
    except psycopg2.Error as e:
        log_error(f"unique_pmid 생성 오류: {e}")
        return 0

def create_filtered_pmid_from_unique():
    """unique_pmid에서 papers 테이블과 중복되지 않는 pmid를 filtered_pmid에 저장합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # 기존 filtered_pmid 테이블 초기화
        cur.execute("DELETE FROM filtered_pmid")
        log_debug("기존 filtered_pmid 테이블 초기화 완료")
        
        # unique_pmid에서 papers와 중복되지 않는 pmid 조회 및 삽입
        cur.execute("""
            INSERT INTO filtered_pmid (pmid, filter_reason, created_at)
            SELECT 
                up.pmid,
                '신규' as filter_reason,
                NOW() as created_at
            FROM unique_pmid up
            LEFT JOIN papers p ON p.pmid = up.pmid
            WHERE p.pmid IS NULL
        """)
        
        filtered_count = cur.rowcount
        conn.commit()
        
        log_info(f"filtered_pmid 테이블 생성 완료: {filtered_count}개 (papers와 중복 제외)")
        
        cur.close()
        conn.close()
        
        return filtered_count
        
    except psycopg2.Error as e:
        log_error(f"filtered_pmid 생성 오류: {e}")
        return 0

def show_pmid_filter_stats():
    """PMID 필터링 과정의 통계를 표시합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # 각 단계별 통계 조회
        cur.execute("SELECT COUNT(*) FROM collection_pmid")
        collection_pmid_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM unique_pmid")
        unique_pmid_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM filtered_pmid")
        filtered_pmid_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM papers WHERE pmid IS NOT NULL")
        existing_papers_count = cur.fetchone()[0]
        
        log_info("=== PMID 필터링 통계 ===")
        log_info(f"collection_pmid: {collection_pmid_count}개 (전체 수집된 PMID)")
        log_info(f"unique_pmid: {unique_pmid_count}개 (중복 제거된 PMID)")
        log_info(f"기존 papers: {existing_papers_count}개 (이미 저장된 논문)")
        log_info(f"filtered_pmid: {filtered_pmid_count}개 (새로 처리할 PMID)")
        
        if collection_pmid_count > 0:
            dedup_rate = (collection_pmid_count - unique_pmid_count) / collection_pmid_count * 100
            log_info(f"컬렉션 내 중복률: {dedup_rate:.1f}%")
        
        if unique_pmid_count > 0:
            existing_rate = (unique_pmid_count - filtered_pmid_count) / unique_pmid_count * 100
            log_info(f"기존 논문 중복률: {existing_rate:.1f}%")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"통계 조회 오류: {e}")

if __name__ == "__main__":
    log_info("=== PMID 필터링 시작 ===")
    
    # 1단계: collection_pmid에서 unique_pmid 생성
    log_info("1단계: collection_pmid → unique_pmid (중복 제거)")
    unique_count = create_unique_pmid_from_collection_pmid()
    
    if unique_count == 0:
        log_error("unique_pmid 생성 실패. collection_pmid 테이블을 확인하세요")
        sys.exit(1)
    
    # 2단계: unique_pmid에서 filtered_pmid 생성 (papers와 중복 제외)
    log_info("2단계: unique_pmid → filtered_pmid (기존 논문 중복 제외)")
    filtered_count = create_filtered_pmid_from_unique()
    
    # 최종 통계 표시
    show_pmid_filter_stats()
    
    log_info(f"=== PMID 필터링 완료: 최종 {filtered_count}개 PMID가 새로 처리 대상입니다 ===")
