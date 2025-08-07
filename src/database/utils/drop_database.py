# drop_database.py - 데이터베이스 삭제 도구 (개발용)
import psycopg2
from common.logger import log_info, log_error

def drop_database():
    """
    ai4ref 데이터베이스를 삭제합니다. (개발용)
    """
    try:
        # postgres 기본 데이터베이스에 연결
        conn = psycopg2.connect(
            host="localhost",
            database="postgres",
            user="postgres",
            password="postgres"
        )
        conn.autocommit = True
        
        cur = conn.cursor()
        
        # ai4ref 데이터베이스 존재 확인
        cur.execute("SELECT 1 FROM pg_catalog.pg_database WHERE datname = 'ai4ref';")
        exists = cur.fetchone()
        
        if not exists:
            log_info("ai4ref 데이터베이스가 존재하지 않습니다")
            cur.close()
            conn.close()
            return True
        
        # 활성 연결 종료 (강제)
        log_info("활성 연결을 종료합니다...")
        cur.execute("""
            SELECT pg_terminate_backend(pid)
            FROM pg_stat_activity
            WHERE datname = 'ai4ref' AND pid <> pg_backend_pid();
        """)
        
        # 데이터베이스 삭제
        log_info("ai4ref 데이터베이스를 삭제합니다...")
        cur.execute("DROP DATABASE ai4ref;")
        log_info("ai4ref 데이터베이스가 삭제되었습니다")
        
        cur.close()
        conn.close()
        return True
        
    except psycopg2.Error as e:
        log_error(f"데이터베이스 삭제 중 오류: {e}")
        return False

def confirm_deletion():
    """
    사용자에게 삭제 확인을 요청합니다.
    """
    log_info("⚠️  주의: ai4ref 데이터베이스와 모든 데이터가 삭제됩니다!")
    log_info("이 작업은 되돌릴 수 없습니다.")
    
    try:
        confirm = input("\n정말로 삭제하시겠습니까? (y/N): ").lower().strip()
        
        if confirm in ['y', 'yes']:
            log_info("삭제를 진행합니다...")
            return True
        else:
            log_info("삭제가 취소되었습니다")
            return False
            
    except KeyboardInterrupt:
        log_info("\n삭제가 취소되었습니다")
        return False

if __name__ == "__main__":
    log_info("=== 데이터베이스 삭제 도구 (개발용) ===")
    
    # 확인 절차
    if confirm_deletion():
        success = drop_database()
        
        if success:
            log_info("데이터베이스 삭제 완료")
            log_info("새로 생성하려면: python src/database/utils/initialize_database.py")
        else:
            log_error("데이터베이스 삭제 실패")
    
    log_info("=== 데이터베이스 삭제 도구 종료 ===")
