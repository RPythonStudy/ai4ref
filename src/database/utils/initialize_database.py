# initialize_database.py - 데이터베이스 초기화 도구
import psycopg2
from pathlib import Path
from common.logger import log_info, log_error

def show_table_structure():
    """
    생성된 papers 테이블의 구조를 표 형태로 보여줍니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()
        
        # 테이블 구조 조회
        cur.execute("""
            SELECT 
                column_name,
                data_type,
                character_maximum_length,
                is_nullable,
                column_default
            FROM information_schema.columns 
            WHERE table_name = 'papers'
            ORDER BY ordinal_position;
        """)
        
        columns = cur.fetchall()
        
        if not columns:
            log_info("papers 테이블이 존재하지 않습니다")
            return
        
        # 테이블 헤더
        print("\n" + "="*80)
        print(f"{'📋 papers 테이블 구조':^80}")
        print("="*80)
        print(f"{'컬럼명':<20} {'타입':<15} {'길이':<8} {'NULL':<8} {'기본값':<20}")
        print("-"*80)
        
        # 테이블 데이터
        for col in columns:
            column_name = col[0]
            data_type = col[1]
            max_length = str(col[2]) if col[2] else ""
            is_nullable = "YES" if col[3] == "YES" else "NO"
            default_val = str(col[4])[:18] + "..." if col[4] and len(str(col[4])) > 18 else str(col[4]) if col[4] else ""
            
            print(f"{column_name:<20} {data_type:<15} {max_length:<8} {is_nullable:<8} {default_val:<20}")
        
        # 인덱스 정보
        cur.execute("""
            SELECT indexname, indexdef 
            FROM pg_indexes 
            WHERE tablename = 'papers'
            ORDER BY indexname;
        """)
        
        indexes = cur.fetchall()
        
        if indexes:
            print("\n" + "-"*80)
            print(f"{'📊 인덱스 정보':^80}")
            print("-"*80)
            for idx_name, idx_def in indexes:
                print(f"• {idx_name}")
                if "UNIQUE" in idx_def:
                    print(f"  └ UNIQUE 제약조건")
                print()
        
        print("="*80 + "\n")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"테이블 구조 조회 중 오류: {e}")

def create_database_if_not_exists():
    """
    ai4ref 데이터베이스가 없으면 생성합니다.
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
        
        if exists:
            log_info("ai4ref 데이터베이스가 이미 존재합니다")
        else:
            cur.execute("CREATE DATABASE ai4ref;")
            log_info("ai4ref 데이터베이스를 생성했습니다")
        
        cur.close()
        conn.close()
        return True
        
    except psycopg2.Error as e:
        log_error(f"데이터베이스 처리 중 오류: {e}")
        return False

def create_papers_table():
    """
    papers 테이블을 생성합니다. 이미 존재하면 삭제 여부를 묻습니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres", 
            password="postgres"
        )
        
        cur = conn.cursor()
        
        # 테이블 존재 확인
        cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'papers');")
        exists = cur.fetchone()[0]
        
        if exists:
            log_info("papers 테이블이 이미 존재합니다")
            try:
                confirm = input("삭제하고 다시 생성하시겠습니까? (y/n): ")
                if confirm.lower() != 'y':
                    log_info("테이블 생성이 취소되었습니다")
                    cur.close()
                    conn.close()
                    return True
            except KeyboardInterrupt:
                log_info("테이블 생성이 취소되었습니다")
                cur.close()
                conn.close()
                return True
            
            cur.execute("DROP TABLE papers CASCADE;")
            log_info("기존 papers 테이블을 삭제했습니다")
        
        # 테이블 생성
        sql_file = Path(__file__).parent.parent / "schema" / "papers_table.sql"
        with open(sql_file, 'r', encoding='utf-8') as f:
            sql = f.read()
        
        cur.execute(sql)
        conn.commit()
        log_info("papers 테이블을 생성했습니다")
        
        cur.close()
        conn.close()
        return True
        
    except psycopg2.Error as e:
        log_error(f"테이블 처리 중 오류: {e}")
        return False
    except FileNotFoundError:
        log_error("papers_table.sql 파일을 찾을 수 없습니다")
        return False

if __name__ == "__main__":
    log_info("=== 데이터베이스 초기화 시작 ===")
    
    # 1. 데이터베이스 생성
    success = create_database_if_not_exists()
    
    if not success:
        log_error("데이터베이스 초기화 실패")
        exit(1)
    
    # 2. 테이블 생성
    success = create_papers_table()
    
    if success:
        log_info("초기화 완료: 데이터베이스와 테이블이 준비되었습니다")
        
        # 3. 테이블 구조 표시
        show_table_structure()
    else:
        log_error("테이블 생성 실패")
    
    log_info("=== 데이터베이스 초기화 종료 ===")