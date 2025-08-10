# initialize_database.py - 데이터베이스 초기화 도구
import psycopg2
import sys
import argparse
from pathlib import Path
from common.logger import log_info, log_error, log_debug
from database.utils.analyze_database import show_table_structure, check_column_stats, show_table_columns

def create_database_if_not_exists():
    """ai4ref 데이터베이스 생성"""
    try:
        conn = psycopg2.connect(host="localhost", database="postgres", user="postgres", password="postgres")
        conn.autocommit = True
        cur = conn.cursor()
        
        cur.execute("SELECT 1 FROM pg_catalog.pg_database WHERE datname = 'ai4ref';")
        if cur.fetchone():
            log_info("데이터베이스 ai4ref 존재함")
        else:
            cur.execute("CREATE DATABASE ai4ref;")
            log_info("데이터베이스 ai4ref 생성완료")
        
        cur.close()
        conn.close()
        return True
    except psycopg2.Error as e:
        log_error(f"데이터베이스 생성실패: {e}")
        return False

def drop_database_completely():
    """ai4ref 데이터베이스 완전삭제"""
    log_info("=== 데이터베이스 완전삭제 시작 ===")
    try:
        conn = psycopg2.connect(host="localhost", database="postgres", user="postgres", password="postgres")
        conn.autocommit = True
        cur = conn.cursor()
        
        # 활성연결 종료
        cur.execute("SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = 'ai4ref' AND pid <> pg_backend_pid();")
        # 데이터베이스 삭제
        cur.execute("DROP DATABASE IF EXISTS ai4ref;")
        log_info("데이터베이스 ai4ref 삭제완료")
        
        cur.close()
        conn.close()
        return True
    except psycopg2.Error as e:
        log_error(f"데이터베이스 삭제실패: {e}")
        return False

def create_database_tables():
    """데이터베이스 테이블 생성"""
    # 테이블 생성 순서 (외래키 의존성 고려)
    table_order = ['collections', 'papers', 'collection_pmid', 'unique_pmid', 'filtered_pmid', 'paper_collection']
    
    try:
        with psycopg2.connect(host="localhost", database="ai4ref", user="postgres", password="postgres") as conn:
            with conn.cursor() as cur:
                # 기존 테이블 확인
                existing = []
                for table in table_order:
                    cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = %s);", (table,))
                    if cur.fetchone()[0]:
                        existing.append(table)
                
                # 기존 테이블 처리
                if existing:
                    log_info(f"기존테이블: {', '.join(existing)}")
                    try:
                        confirm = input("모든 테이블을 삭제하고 다시 생성하시겠습니까? (y/n): ")
                        if confirm.lower() != 'y':
                            log_info("테이블생성 취소됨")
                            return True
                    except KeyboardInterrupt:
                        log_info("테이블생성 취소됨")
                        return True
                    
                    # 역순 삭제 (외래키 참조 고려)
                    for table in reversed(table_order):
                        cur.execute(f"DROP TABLE IF EXISTS {table} CASCADE;")
                    log_info("기존테이블 삭제완료")
                
                # 테이블 생성
                schema_dir = Path(__file__).parent.parent / "schema"
                for table in table_order:
                    sql_file = schema_dir / f"{table}.sql"
                    if not sql_file.exists():
                        raise FileNotFoundError(f"스키마파일 없음: {sql_file}")
                    with open(sql_file, 'r', encoding='utf-8') as f:
                        sql_content = f.read()
                        log_debug(f"테이블 생성 SQL 파일: {sql_file}")
                        log_debug(f"SQL 일부: {sql_content[:200]}")
                        cur.execute(sql_content)
                    log_info(f"테이블생성: {table}")
                
                conn.commit()
                log_info("모든테이블 생성완료")
        return True
        
    except psycopg2.Error as e:
        log_error(f"테이블생성실패: {e}")
        return False
    except FileNotFoundError as e:
        log_error(f"파일오류: {e}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ai4ref 데이터베이스 초기화")
    parser.add_argument('--reset', action='store_true', help='기존 데이터베이스 삭제 후 재생성')
    parser.add_argument('--verbose', action='store_true', help='상세 테이블 구조 표시')
    args = parser.parse_args()
    
    log_info("=== 데이터베이스 초기화 시작 ===")
    
    # 1. 리셋 모드: 기존 DB 완전삭제
    if args.reset:
        if not drop_database_completely():
            log_error("데이터베이스 삭제실패")
            sys.exit(1)
    
    # 2. 데이터베이스 생성
    if not create_database_if_not_exists():
        log_error("데이터베이스 생성실패")
        sys.exit(1)
    
    # 3. 테이블 생성
    if not create_database_tables():
        log_error("테이블 생성실패")
        sys.exit(1)
    
    # 4. 완료 및 상태확인
    log_info("=== 초기화 완료 ===")
    log_info("생성완료: collections, papers, collection_pmid, unique_pmid, filtered_pmid, paper_collection")
    
    show_table_columns()
    
    if args.verbose:
        show_table_structure()
        check_column_stats("초기화완료")
    
    log_info("=== 초기화 종료 ===")