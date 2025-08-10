import json
import psycopg2
from pathlib import Path
from common.logger import log_info, log_error

def load_search_collections(json_file_path=None):
    """search_collections.json을 collections 테이블에 덮어쓰기로 로드"""
    log_info("=== collections 테이블 overwrite 시작 ===")
    
    # JSON 파일 경로 설정
    if json_file_path is None:
        project_root = Path(__file__).parent.parent.parent
        json_file_path = project_root / "config" / "search_collections.json"
    
    try:
        # JSON 파일 읽기
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        search_collections = data.get('search_collections', [])
        if not search_collections:
            log_info("JSON 파일에 컬렉션 데이터가 없습니다")
            return True
        
        log_info(f"JSON에서 {len(search_collections)}개 컬렉션 읽음")
        
        # 데이터베이스 연결
        conn = psycopg2.connect(
            host="localhost", database="ai4ref", 
            user="postgres", password="postgres"
        )
        cur = conn.cursor()
        
        # 기존 데이터 확인
        cur.execute("SELECT COUNT(*) FROM collections")
        existing_count = cur.fetchone()[0]
        
        # JSON 컬렉션 이름들 추출
        json_collection_names = {c.get('collection') for c in search_collections if c.get('collection')}
        
        # 기존 컬렉션 이름들 조회
        cur.execute("SELECT name FROM collections")
        existing_names = {row[0] for row in cur.fetchall()}
        
        # 일치/불일치 분석
        matching_count = len(json_collection_names & existing_names)
        new_count = len(json_collection_names - existing_names)
        
        log_info(f"기존 컬렉션: {existing_count}개")
        log_info(f"JSON과 일치하는 것: {matching_count}개, 새로운 것: {new_count}개")
        
        # 전체 삭제 후 재삽입
        cur.execute("DELETE FROM collections")
        log_info(f"기존 데이터 삭제: {existing_count}개")
        
        # 새 데이터 삽입
        inserted_count = 0
        for collection in search_collections:
            name = collection.get('collection')
            term = collection.get('term')
            retmax = collection.get('retmax', 10000)
            enabled = collection.get('enabled', True)
            description = collection.get('description', '')
            
            if name and term:
                cur.execute("""
                    INSERT INTO collections (name, search_term, retmax, enabled, description, created_at)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                """, (name, term, retmax, enabled, description))
                inserted_count += 1
        
        conn.commit()
        log_info(f"overwrite {inserted_count}개 완료")
        
        cur.close()
        conn.close()
        return True
        
    except Exception as e:
        log_error(f"overwrite 실패: {e}")
        if 'conn' in locals():
            conn.rollback()
        return False

def show_collections_table():
    """collections 테이블 현재 상태 조회"""
    try:
        conn = psycopg2.connect(
            host="localhost", database="ai4ref", 
            user="postgres", password="postgres"
        )
        cur = conn.cursor()
        
        cur.execute("SELECT id, name, search_term, retmax, enabled FROM collections ORDER BY id")
        rows = cur.fetchall()
        
        if not rows:
            log_info("collections 테이블이 비어있습니다")
        else:
            log_info(f"=== collections 테이블 현황 ({len(rows)}개) ===")
            for row in rows:
                id, name, search_term, retmax, enabled = row
                status = "활성" if enabled else "비활성"
                log_info(f"[{id}] {name} ({status}) - retmax: {retmax}")
                log_info(f"    검색어: {search_term[:100]}")
        
        cur.close()
        conn.close()
        return True
        
    except Exception as e:
        log_error(f"테이블 조회 오류: {e}")
        return False

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="collections 테이블 로더")
    parser.add_argument('--file', type=str, help='JSON 파일 경로')
    parser.add_argument('--show', action='store_true', help='테이블 상태 조회')
    
    args = parser.parse_args()
    
    if args.show:
        show_collections_table()
    else:
        success = load_search_collections(args.file)
        if success:
            log_info("[완료] collections 테이블 overwrite 성공")
        else:
            log_error("[실패] collections 테이블 overwrite 실패")
