# create_paper_collection.py - 논문-컬렉션 관계 테이블 생성
from common.database import get_db_connection
from common.logger import log_info, log_error, log_debug

def create_paper_collection_relations():
    """collection_pmid와 papers를 조인하여 paper_collection 관계 생성"""
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_info("=== Paper-Collection 관계 생성 시작 ===")
        
        # 기존 paper_collection 데이터 확인
        cur.execute("SELECT COUNT(*) FROM paper_collection")
        existing_count = cur.fetchone()[0]
        log_info(f"기존 paper_collection 관계: {existing_count}개")
        
        # collection_pmid와 papers 조인으로 관계 생성
        cur.execute("""
            INSERT INTO paper_collection (paper_id, collection_id, collection_name, added_at)
            SELECT DISTINCT 
                p.id as paper_id,
                cp.collection_id,
                cp.collection_name,
                NOW() as added_at
            FROM collection_pmid cp
            JOIN papers p ON p.pmid = cp.pmid
            JOIN collections c ON c.id = cp.collection_id
            WHERE NOT EXISTS (
                SELECT 1 FROM paper_collection pc 
                WHERE pc.paper_id = p.id 
                AND pc.collection_id = cp.collection_id
            )
        """)
        
        new_relations = cur.rowcount
        conn.commit()
        
        # 최종 현황 확인
        cur.execute("SELECT COUNT(*) FROM paper_collection")
        total_relations = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT paper_id) FROM paper_collection")
        unique_papers = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT collection_id) FROM paper_collection")
        unique_collections = cur.fetchone()[0]
        
        log_info(f"신규 관계 생성: {new_relations}개")
        log_info(f"전체 관계: {total_relations}개")
        log_info(f"고유 논문: {unique_papers}개")
        log_info(f"고유 컬렉션: {unique_collections}개")
        
        # 다중 컬렉션 논문 통계
        cur.execute("""
            SELECT COUNT(*) as paper_count, collection_count
            FROM (
                SELECT paper_id, COUNT(*) as collection_count
                FROM paper_collection
                GROUP BY paper_id
            ) stats
            GROUP BY collection_count
            ORDER BY collection_count
        """)
        
        multi_stats = cur.fetchall()
        log_info("논문별 컬렉션 분포:")
        for paper_count, collection_count in multi_stats:
            log_info(f"  {collection_count}개 컬렉션: {paper_count}개 논문")
        
        log_info("=== Paper-Collection 관계 생성 완료 ===")
        
        return new_relations > 0
        
    except Exception as e:
        log_error(f"Paper-Collection 관계 생성 오류: {e}")
        conn.rollback()
        return False
    finally:
        cur.close()
        conn.close()

def show_paper_collection_status():
    """Paper-Collection 관계 현황 조회"""
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_info("=== Paper-Collection 관계 현황 ===")
        
        # 전체 통계
        cur.execute("SELECT COUNT(*) FROM paper_collection")
        total_relations = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT paper_id) FROM paper_collection")
        unique_papers = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT collection_id) FROM paper_collection")
        unique_collections = cur.fetchone()[0]
        
        log_info(f"전체 관계: {total_relations}개")
        log_info(f"고유 논문: {unique_papers}개")
        log_info(f"고유 컬렉션: {unique_collections}개")
        
        if unique_papers > 0:
            avg_collections_per_paper = total_relations / unique_papers
            log_info(f"논문당 평균 컬렉션 수: {avg_collections_per_paper:.1f}개")
        
        # paper_collection 테이블 샘플 디버깅 출력
        cur.execute("SELECT paper_id, collection_id, collection_name FROM paper_collection ORDER BY paper_id, collection_id LIMIT 20;")
        debug_rows = cur.fetchall()
        log_debug("paper_collection 샘플 (상위 20개):")
        for row in debug_rows:
            log_debug(f"  paper_id={row[0]}, collection_id={row[1]}, collection_name={row[2]}")

        # 컬렉션별 논문 수
        cur.execute("""
            SELECT c.name as collection_name, COUNT(pc.paper_id) as paper_count
            FROM collections c
            LEFT JOIN paper_collection pc ON c.id = pc.collection_id
            GROUP BY c.id, c.name
            ORDER BY paper_count DESC
        """)
        
        collection_stats = cur.fetchall()
        log_info("컬렉션별 논문 수:")
        for collection_name, paper_count in collection_stats:
            log_info(f"  {collection_name}: {paper_count}개")
        
        # PDF 보유 논문 통계
        cur.execute("""
            SELECT 
                COUNT(DISTINCT pc.paper_id) as total_papers,
                COUNT(DISTINCT CASE WHEN p.pdf_url IS NOT NULL THEN pc.paper_id END) as papers_with_pdf
            FROM paper_collection pc
            JOIN papers p ON p.id = pc.paper_id
        """)
        
        total_papers, papers_with_pdf = cur.fetchone()
        if total_papers > 0:
            pdf_rate = papers_with_pdf / total_papers * 100
            log_info(f"PDF 보유율: {papers_with_pdf}/{total_papers}개 ({pdf_rate:.1f}%)")
        
        # Zotero 동기화 상태
        cur.execute("""
            SELECT 
                COUNT(*) as total_relations,
                COUNT(CASE WHEN synced_at IS NOT NULL THEN 1 END) as synced_relations
            FROM paper_collection
        """)
        
        total_rel, synced_rel = cur.fetchone()
        if total_rel > 0:
            sync_rate = synced_rel / total_rel * 100
            log_info(f"Zotero 동기화: {synced_rel}/{total_rel}개 ({sync_rate:.1f}%)")
        
    except Exception as e:
        log_error(f"현황 조회 오류: {e}")
    finally:
        cur.close()
        conn.close()

def analyze_multi_collection_papers(limit=10):
    """다중 컬렉션 할당 논문 분석"""
    
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_info("=== 다중 컬렉션 논문 분석 ===")
        
        # 다중 컬렉션 논문 조회
        cur.execute("""
            SELECT 
                p.pmid,
                p.title,
                COUNT(pc.collection_id) as collection_count,
                CASE WHEN p.pdf_url IS NOT NULL THEN 'Y' ELSE 'N' END as has_pdf
            FROM paper_collection pc
            JOIN papers p ON p.id = pc.paper_id
            GROUP BY p.id, p.pmid, p.title, p.pdf_url
            HAVING COUNT(pc.collection_id) > 1
            ORDER BY COUNT(pc.collection_id) DESC, p.pmid
            LIMIT %s
        """, (limit,))
        
        multi_papers = cur.fetchall()
        
        if multi_papers:
            log_info(f"다중 컬렉션 논문 (상위 {len(multi_papers)}개):")
            for pmid, title, count, has_pdf in multi_papers:
                title_short = title[:60] + "..." if len(title) > 60 else title
                log_info(f"  PMID {pmid} ({count}개 컬렉션, PDF:{has_pdf}): {title_short}")
                
                # 해당 논문의 컬렉션 목록
                cur.execute("""
                    SELECT c.name as collection_name
                    FROM paper_collection pc
                    JOIN collections c ON c.id = pc.collection_id
                    JOIN papers p ON p.id = pc.paper_id
                    WHERE p.pmid = %s
                    ORDER BY c.name
                """, (pmid,))
                
                collections = [row[0] for row in cur.fetchall()]
                log_debug(f"    컬렉션: {', '.join(collections)}")
        else:
            log_info("다중 컬렉션 논문이 없습니다")
            
    except Exception as e:
        log_error(f"다중 컬렉션 분석 오류: {e}")
    finally:
        cur.close()
        conn.close()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Paper-Collection 관계 관리")
    parser.add_argument('--status', action='store_true', help='현황 조회')
    parser.add_argument('--analyze', action='store_true', help='다중 컬렉션 분석')
    parser.add_argument('--limit', type=int, default=10, help='분석 대상 논문 수')
    
    args = parser.parse_args()
    
    if args.status:
        show_paper_collection_status()
    elif args.analyze:
        analyze_multi_collection_papers(args.limit)
    else:
        create_paper_collection_relations()
