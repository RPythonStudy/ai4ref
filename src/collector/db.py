import psycopg2
import os
from common.logger import log_info, log_error

DB_HOST = os.environ.get("DB_HOST", "localhost")
DB_PORT = os.environ.get("DB_PORT", 5432)
DB_NAME = os.environ.get("POSTGRES_DB", "ai4ref")
DB_USER = os.environ.get("POSTGRES_USER", "postgres")
DB_PASS = os.environ.get("POSTGRES_PASSWORD", "postgres")

def save_pmids_to_db(pmids):
    try:
        conn = psycopg2.connect(
            host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
        )
        
        cur = conn.cursor()
        cur.execute("""
            CREATE TABLE IF NOT EXISTS pubmed_pmids (
                pmid TEXT PRIMARY KEY,
                collected_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        """)
        
        inserted_count = 0
        duplicate_count = 0
        
        for pmid in pmids:
            try:
                cur.execute(
                    "INSERT INTO pubmed_pmids (pmid) VALUES (%s) ON CONFLICT (pmid) DO NOTHING;",
                    (pmid,),
                )
                if cur.rowcount > 0:
                    inserted_count += 1
                else:
                    duplicate_count += 1
            except Exception as e:
                log_error(f"PMID {pmid} 삽입 오류: {e}")
        
        conn.commit()
        
        # 전체 레코드 수 확인
        cur.execute("SELECT COUNT(*) FROM pubmed_pmids;")
        total_count = cur.fetchone()[0]
        
        cur.close()
        conn.close()
        
        log_info(f"저장 완료: 신규 {inserted_count}개, 중복 {duplicate_count}개, 총 {total_count}개")
        return {"inserted": inserted_count, "duplicates": duplicate_count, "total": total_count}
        
    except Exception as e:
        log_error(f"데이터베이스 작업 중 오류: {e}")
        return None

def get_pmids_count():
    """데이터베이스에 저장된 PMID 개수 조회"""
    try:
        conn = psycopg2.connect(
            host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
        )
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM pubmed_pmids;")
        count = cur.fetchone()[0]
        cur.close()
        conn.close()
        log_info(f"데이터베이스 총 PMID 수: {count}")
        return count
    except Exception as e:
        log_error(f"PMID 개수 조회 오류: {e}")
        return None

def get_recent_pmids(limit=10):
    """최근 저장된 PMID 조회"""
    try:
        conn = psycopg2.connect(
            host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
        )
        cur = conn.cursor()
        cur.execute("SELECT pmid, collected_at FROM pubmed_pmids ORDER BY collected_at DESC LIMIT %s;", (limit,))
        results = cur.fetchall()
        cur.close()
        conn.close()
        log_info(f"최근 {len(results)}개 PMID 조회 완료")
        return results
    except Exception as e:
        log_error(f"최근 PMID 조회 오류: {e}")
        return None

def get_all_pmids():
    """모든 PMID 조회"""
    try:
        conn = psycopg2.connect(
            host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
        )
        cur = conn.cursor()
        cur.execute("SELECT pmid FROM pubmed_pmids ORDER BY collected_at;")
        results = [row[0] for row in cur.fetchall()]
        cur.close()
        conn.close()
        log_info(f"전체 PMID {len(results)}개 조회 완료")
        return results
    except Exception as e:
        log_error(f"전체 PMID 조회 오류: {e}")
        return None

def get_connection():
    """데이터베이스 연결 반환"""
    return psycopg2.connect(
        host=DB_HOST, port=DB_PORT, dbname=DB_NAME, user=DB_USER, password=DB_PASS
    )
