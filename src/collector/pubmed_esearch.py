import os
import sys
import psycopg2
from dotenv import load_dotenv
from Bio import Entrez
from common.logger import log_info, log_debug, log_error
from common.database import get_db_connection


load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL", "r.python.ai@gmail.com")
Entrez.api_key = os.getenv("PUBMED_API_KEY")

def pubmed_esearch(term: str, retmax: int = 10000):
    """PubMed 검색을 수행하여 PMID 리스트를 반환합니다."""
    log_info(f"검색어: {term}")
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        result = Entrez.read(handle)
        handle.close()
        pmids = [str(pmid) for pmid in result.get('IdList', [])]
        log_info(f"총 PMID 개수: {len(pmids)}개")
        return pmids
    except Exception as e:
        log_error(f"검색 오류: {e}")
        return []

def load_collections_from_db():
    """collections 테이블에서 활성화된 컬렉션들을 조회합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        cur.execute("""
            SELECT id, name, search_term, retmax 
            FROM collections 
            WHERE enabled = true
            ORDER BY id
        """)
        
        collections = []
        for row in cur.fetchall():
            id, name, search_term, retmax = row
            collections.append({
                'id': id,
                'name': name,
                'search_term': search_term,
                'retmax': retmax or 10000
            })
        
        cur.close()
        conn.close()
        
        log_debug(f"데이터베이스에서 로드된 컬렉션: {len(collections)}개")
        return collections
        
    except psycopg2.Error as e:
        log_error(f"데이터베이스 조회 오류: {e}")
        return []

def save_pmids_to_collection_pmid(pmids: list, collection_id: int, collection_name: str):
    """수집된 PMID들을 collection_pmid 테이블에 저장합니다."""
    if not pmids:
        return 0
    
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # 기존 중복 PMID 확인
        pmid_list = "(" + ",".join([f"'{pmid}'" for pmid in pmids]) + ")"
        cur.execute(f"""
            SELECT pmid FROM collection_pmid 
            WHERE collection_id = %s AND pmid IN {pmid_list}
        """, (collection_id,))
        
        existing_pmids = set([row[0] for row in cur.fetchall()])
        new_pmids = [pmid for pmid in pmids if pmid not in existing_pmids]
        
        if not new_pmids:
            log_info(f"모든 PMID가 이미 collection_pmid에 존재합니다 (중복 {len(pmids)}개)")
            cur.close()
            conn.close()
            return 0
        
        # 새로운 PMID들을 collection_pmid에 삽입
        insert_data = []
        for i, pmid in enumerate(new_pmids):
            insert_data.append((collection_id, pmid, collection_name, i + 1, len(pmids)))
        
        cur.executemany("""
            INSERT INTO collection_pmid (collection_id, pmid, collection_name, result_position, esearch_count)
            VALUES (%s, %s, %s, %s, %s)
        """, insert_data)
        
        conn.commit()
        cur.close()
        conn.close()
        
        log_info(f"collection_pmid에 저장: 신규 {len(new_pmids)}개, 중복 제외 {len(existing_pmids)}개")
        return len(new_pmids)
        
    except psycopg2.Error as e:
        log_error(f"collection_pmid 저장 오류: {e}")
        return 0

def show_collection_pmid_stats():
    """collection_pmid 테이블의 현재 상태를 조회하여 표시합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # 전체 통계
        cur.execute("SELECT COUNT(*) FROM collection_pmid;")
        total_count = cur.fetchone()[0]
        
        # 컬렉션별 통계
        cur.execute("""
            SELECT c.name, COUNT(cp.pmid) as count
            FROM collections c
            LEFT JOIN collection_pmid cp ON c.id = cp.collection_id
            GROUP BY c.id, c.name
            ORDER BY count DESC;
        """)
        
        collection_stats = cur.fetchall()
        
        log_info(f"=== collection_pmid 테이블 현황 ===")
        log_info(f"총 PMID 수: {total_count}개")
        log_info("컬렉션별 PMID 분포:")
        for name, count in collection_stats:
            log_info(f"  {name}: {count}개")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"통계 조회 오류: {e}")

if __name__ == "__main__":
    log_info("=== PubMed esearch 시작 ===")
    
    # 명령행 인수 처리
    if len(sys.argv) > 1:
        # 수동 검색 모드 (개발/테스트용)
        search_term = sys.argv[1]
        retmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
        
        log_info(f"수동 검색 모드: '{search_term}' (retmax: {retmax})")
        pmids = pubmed_esearch(search_term, retmax)
        
        if pmids:
            log_info(f"수집된 PMID: {len(pmids)}개")
            log_info("수동 검색 모드에서는 collection_pmid에 저장하지 않습니다")
            log_info(f"첫 10개 PMID: {pmids[:10]}")
        else:
            log_info("검색 결과가 없습니다")
        
        sys.exit(0)
    
    # 일반 모드: collections 테이블에서 컬렉션 로드
    collections = load_collections_from_db()
    
    if not collections:
        log_error("활성화된 컬렉션이 없습니다. collections 테이블을 확인하세요")
        sys.exit(1)
    
    log_info(f"수집 대상 컬렉션: {len(collections)}개")
    
    # 수집 통계 초기화
    total_searched = 0
    total_queued = 0
    failed_collections = 0
    
    # 컬렉션별 esearch 실행
    for i, collection in enumerate(collections, 1):
        collection_id = collection['id']
        collection_name = collection['name']
        search_term = collection['search_term']
        retmax = collection['retmax']
        
        log_info(f"검색 {i}/{len(collections)}: '{collection_name}' (retmax: {retmax})")
        
        try:
            # PubMed 검색
            pmids = pubmed_esearch(search_term, retmax)
            total_searched += len(pmids)
            
            if not pmids:
                log_info(f"검색 결과 없음: {collection_name}")
                continue
            
            # collection_pmid에 저장
            queued_count = save_pmids_to_collection_pmid(pmids, collection_id, collection_name)
            total_queued += queued_count
            
        except Exception as e:
            log_error(f"컬렉션 {collection_name} 처리 실패: {e}")
            failed_collections += 1
            continue
    
    # 전체 수집 통계 출력
    log_info("=== 전체 수집 통계 ===")
    log_info(f"총 검색된 PMID: {total_searched}개")
    log_info(f"collection_pmid에 추가된 PMID: {total_queued}개")
    log_info(f"처리 성공한 컬렉션: {len(collections) - failed_collections}/{len(collections)}개")
    
    if failed_collections > 0:
        log_error(f"처리 실패한 컬렉션: {failed_collections}개")
    
    # collection_pmid 테이블 현황 표시
    show_collection_pmid_stats()

    log_info("=== PubMed esearch 완료 ===")
