#!/usr/bin/env python3
"""
데이터베이스 테이블 조회 및 관리 도구
개발 단계에서 테이블 상태 확인용
"""

import psycopg2
from tabulate import tabulate
import sys

def connect_db():
    """PostgreSQL 데이터베이스 연결"""
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref", 
            user="postgres",
            password="postgres"
        )
        return conn
    except psycopg2.Error as e:
        print(f"❌ 데이터베이스 연결 실패: {e}")
        return None

def list_all_tables(conn):
    """모든 테이블 목록 조회"""
    print("\n=== 전체 테이블 목록 ===")
    cur = conn.cursor()
    
    # public 스키마의 모든 테이블 조회
    cur.execute("""
        SELECT tablename, tableowner, tablespace 
        FROM pg_tables 
        WHERE schemaname = 'public'
        ORDER BY tablename;
    """)
    
    tables = cur.fetchall()
    
    if tables:
        headers = ['테이블명', '소유자', '테이블스페이스']
        print(tabulate(tables, headers=headers, tablefmt='grid'))
        print(f"\n총 {len(tables)}개의 테이블이 있습니다.")
    else:
        print("❌ 테이블이 없습니다.")
    
    cur.close()

def check_table_exists(conn, table_name):
    """특정 테이블 존재 여부 확인"""
    print(f"\n=== '{table_name}' 테이블 존재 여부 ===")
    cur = conn.cursor()
    
    cur.execute("""
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = %s
        );
    """, (table_name,))
    
    exists = cur.fetchone()[0]
    
    if exists:
        print(f"✅ '{table_name}' 테이블이 존재합니다.")
        return True
    else:
        print(f"❌ '{table_name}' 테이블이 존재하지 않습니다.")
        return False
    
    cur.close()

def describe_table(conn, table_name):
    """테이블 구조 상세 정보"""
    if not check_table_exists(conn, table_name):
        return
    
    print(f"\n=== '{table_name}' 테이블 구조 ===")
    cur = conn.cursor()
    
    # 컬럼 정보 조회
    cur.execute("""
        SELECT 
            column_name,
            data_type,
            character_maximum_length,
            is_nullable,
            column_default
        FROM information_schema.columns 
        WHERE table_schema = 'public' 
        AND table_name = %s
        ORDER BY ordinal_position;
    """, (table_name,))
    
    columns = cur.fetchall()
    
    if columns:
        headers = ['컬럼명', '데이터타입', '최대길이', 'NULL허용', '기본값']
        print(tabulate(columns, headers=headers, tablefmt='grid'))
    
    # 인덱스 정보 조회
    cur.execute("""
        SELECT 
            indexname,
            indexdef
        FROM pg_indexes 
        WHERE schemaname = 'public' 
        AND tablename = %s;
    """, (table_name,))
    
    indexes = cur.fetchall()
    
    if indexes:
        print(f"\n=== '{table_name}' 인덱스 정보 ===")
        for idx_name, idx_def in indexes:
            print(f"📌 {idx_name}: {idx_def}")
    
    cur.close()

def count_table_rows(conn, table_name):
    """테이블 행 개수 조회"""
    if not check_table_exists(conn, table_name):
        return
    
    print(f"\n=== '{table_name}' 데이터 개수 ===")
    cur = conn.cursor()
    
    try:
        cur.execute(f"SELECT COUNT(*) FROM {table_name};")
        count = cur.fetchone()[0]
        print(f"📊 총 {count:,}개의 행이 있습니다.")
    except psycopg2.Error as e:
        print(f"❌ 데이터 조회 실패: {e}")
    
    cur.close()

def drop_table_if_exists(conn, table_name):
    """테이블 삭제 (개발용)"""
    print(f"\n=== '{table_name}' 테이블 삭제 ===")
    
    # 확인 메시지
    confirm = input(f"정말로 '{table_name}' 테이블을 삭제하시겠습니까? (yes/no): ")
    
    if confirm.lower() not in ['yes', 'y']:
        print("❌ 삭제가 취소되었습니다.")
        return
    
    cur = conn.cursor()
    
    try:
        cur.execute(f"DROP TABLE IF EXISTS {table_name} CASCADE;")
        conn.commit()
        print(f"✅ '{table_name}' 테이블이 삭제되었습니다.")
    except psycopg2.Error as e:
        print(f"❌ 테이블 삭제 실패: {e}")
        conn.rollback()
    
    cur.close()

def main():
    """메인 함수"""
    conn = connect_db()
    if not conn:
        return
    
    if len(sys.argv) < 2:
        print("사용법:")
        print("  python db_inspector.py list                    # 모든 테이블 목록")
        print("  python db_inspector.py check <table_name>      # 테이블 존재 확인")
        print("  python db_inspector.py describe <table_name>   # 테이블 구조 조회")
        print("  python db_inspector.py count <table_name>      # 데이터 개수 조회")
        print("  python db_inspector.py drop <table_name>       # 테이블 삭제 (주의!)")
        print("\n예시:")
        print("  python db_inspector.py list")
        print("  python db_inspector.py describe papers")
        conn.close()
        return
    
    command = sys.argv[1].lower()
    
    try:
        if command == "list":
            list_all_tables(conn)
        
        elif command == "check" and len(sys.argv) > 2:
            table_name = sys.argv[2]
            check_table_exists(conn, table_name)
        
        elif command == "describe" and len(sys.argv) > 2:
            table_name = sys.argv[2]
            describe_table(conn, table_name)
        
        elif command == "count" and len(sys.argv) > 2:
            table_name = sys.argv[2]
            count_table_rows(conn, table_name)
        
        elif command == "drop" and len(sys.argv) > 2:
            table_name = sys.argv[2]
            drop_table_if_exists(conn, table_name)
        
        else:
            print("❌ 잘못된 명령어입니다.")
    
    finally:
        conn.close()
        print("\n데이터베이스 연결 종료")

if __name__ == "__main__":
    main()
