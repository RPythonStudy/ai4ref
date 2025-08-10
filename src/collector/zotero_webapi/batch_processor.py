"""
Zotero 배치 처리 관련 기능
- 대량 논문 업로드, 진행 상황 관리
"""
from common.database import get_db_connection
from common.logger import log_info, log_error, log_debug
from .item_operations import ZoteroItemManager
from .collection_operations import ZoteroCollectionManager
from .attachment_manager import ZoteroAttachmentManager
from .duplicate_detector import ZoteroDuplicateDetector
from pathlib import Path

class ZoteroBatchProcessor:
    def __init__(self, user_id=None, api_key=None):
        self.item_manager = ZoteroItemManager(user_id, api_key)
        self.collection_manager = ZoteroCollectionManager(user_id, api_key)
        self.attachment_manager = ZoteroAttachmentManager(user_id, api_key)
        self.duplicate_detector = ZoteroDuplicateDetector(user_id, api_key)
    
    def process_paper_collections(self, paper_id):
        """단일 논문의 다중 컬렉션 처리 (중복 감지 포함)"""
        log_debug(f"=== process_paper_collections 시작: Paper ID {paper_id} ===")
        conn = get_db_connection()
        cur = conn.cursor()
        
        try:
            # 논문 정보 조회
            paper = self._get_paper_data(cur, paper_id)
            if not paper:
                return False
            
            # 컬렉션 관계 조회
            collection_rows = self._get_collection_relations(cur, paper_id)
            if not collection_rows:
                return False
            
            log_info(f"논문 처리 시작: PMID {paper['pmid']} ({len(collection_rows)}개 컬렉션)")
            
            # 중복 감지 및 아이템 키 결정
            item_key, creation_type = self.duplicate_detector.resolve_item_conflict(paper)
            if not item_key:
                log_error(f"아이템 처리 실패: PMID {paper['pmid']}")
                return False
            
            # PDF 첨부 (한 번만)
            pdf_attached = self._attach_pdf_if_exists(item_key, paper)
            
            # 모든 컬렉션에 아이템 링크
            success_count = self._link_to_collections(cur, paper_id, item_key, collection_rows)
            
            conn.commit()
            log_debug(f"DB 커밋 완료: Paper ID {paper_id}")
            
            final_result = success_count > 0
            log_info(f"논문 처리 완료: PMID {paper['pmid']} -> {success_count}/{len(collection_rows)} 컬렉션, PDF: {'O' if pdf_attached else 'X'}")
            return final_result
            
        except Exception as e:
            log_error(f"논문 처리 오류 paper_id {paper_id}: {e}")
            import traceback
            log_error(f"스택 트레이스: {traceback.format_exc()}")
            conn.rollback()
            return False
        finally:
            cur.close()
            conn.close()
    
    def upload_all_papers(self, batch_size=10):
        """모든 미동기화 논문 업로드"""
        conn = get_db_connection()
        cur = conn.cursor()
        
        try:
            # 미동기화 논문 조회
            paper_ids = self._get_unsynced_papers(cur)
            
            if not paper_ids:
                log_info("업로드할 논문이 없습니다")
                return
            
            log_info(f"=== Zotero 업로드 시작: {len(paper_ids)}개 논문 ===")
            
            success_count = 0
            failed_count = 0
            
            for i, paper_id in enumerate(paper_ids, 1):
                log_info(f"진행률: {i}/{len(paper_ids)} - Paper ID {paper_id}")
                
                try:
                    if self.process_paper_collections(paper_id):
                        success_count += 1
                    else:
                        failed_count += 1
                        
                except Exception as e:
                    failed_count += 1
                    log_error(f"논문 처리 예외 Paper ID {paper_id}: {e}")
                
                # 배치 단위로 진행 상황 보고
                if i % batch_size == 0:
                    log_debug(f"배치 완료: {i}개 처리됨 (성공: {success_count}, 실패: {failed_count})")
            
            success_rate = success_count / len(paper_ids) * 100 if paper_ids else 0
            log_info(f"=== Zotero 업로드 완료: {success_count}/{len(paper_ids)}개 성공 ({success_rate:.1f}%) ===")
            
        except Exception as e:
            log_error(f"전체 업로드 오류: {e}")
        finally:
            cur.close()
            conn.close()
    
    def show_sync_status(self):
        """동기화 상태 조회"""
        conn = get_db_connection()
        cur = conn.cursor()
        
        try:
            # 전체 현황
            cur.execute("SELECT COUNT(DISTINCT paper_id) FROM paper_collection")
            total_papers = cur.fetchone()[0]
            
            cur.execute("SELECT COUNT(DISTINCT paper_id) FROM paper_collection WHERE synced_at IS NOT NULL")
            synced_papers = cur.fetchone()[0]
            
            cur.execute("SELECT COUNT(*) FROM paper_collection WHERE synced_at IS NOT NULL")
            synced_relations = cur.fetchone()[0]
            
            cur.execute("SELECT COUNT(*) FROM paper_collection")
            total_relations = cur.fetchone()[0]
            
            log_info(f"=== Zotero 동기화 현황 ===")
            log_info(f"논문: {synced_papers}/{total_papers}개 동기화")
            log_info(f"관계: {synced_relations}/{total_relations}개 동기화")
            
            if total_papers > 0:
                paper_rate = synced_papers / total_papers * 100
                log_info(f"논문 동기화율: {paper_rate:.1f}%")
            
            if total_relations > 0:
                relation_rate = synced_relations / total_relations * 100
                log_info(f"관계 동기화율: {relation_rate:.1f}%")
                
        except Exception as e:
            log_error(f"상태 조회 오류: {e}")
        finally:
            cur.close()
            conn.close()
    
    # Helper methods
    def _get_paper_data(self, cur, paper_id):
        """논문 데이터 조회"""
        cur.execute("""
            SELECT pmid, title, abstract, publication_year, authors, doi, journal, pdf_url
            FROM papers WHERE id = %s
        """, (paper_id,))
        
        paper_row = cur.fetchone()
        if not paper_row:
            log_error(f"논문 정보 없음: paper_id {paper_id}")
            return None
        
        return dict(zip(['pmid', 'title', 'abstract', 'publication_year', 'authors', 'doi', 'journal', 'pdf_url'], paper_row))
    
    def _get_collection_relations(self, cur, paper_id):
        """컬렉션 관계 조회"""
        cur.execute("""
            SELECT pc.collection_id, c.name as collection_name, pc.zotero_item_key, pc.synced_at
            FROM paper_collection pc
            JOIN collections c ON pc.collection_id = c.id
            WHERE pc.paper_id = %s
            ORDER BY pc.collection_id
        """, (paper_id,))
        
        collection_rows = cur.fetchall()
        if not collection_rows:
            log_debug(f"컬렉션 관계 없음: paper_id {paper_id}")
            return None
        
        return collection_rows
    
    def _attach_pdf_if_exists(self, item_key, paper):
        """PDF 첨부 처리"""
        if paper['pdf_url'] and Path(paper['pdf_url']).exists():
            return self.attachment_manager.attach_pdf_to_item(item_key, paper['pdf_url'], paper['pmid'])
        else:
            log_debug(f"PDF 첨부 생략: PMID {paper['pmid']} (파일 없음 또는 경로 문제)")
            return False
    
    def _link_to_collections(self, cur, paper_id, item_key, collection_rows):
        """모든 컬렉션에 아이템 링크"""
        success_count = 0
        
        for collection_id, collection_name, _, _ in collection_rows:
            # 컬렉션 확인/생성
            collection_key = self.collection_manager.get_or_create_collection(collection_name)
            if not collection_key:
                log_error(f"컬렉션 키 확보 실패: '{collection_name}'")
                continue
            
            # 컬렉션에 아이템 추가
            if self.collection_manager.add_item_to_collection(item_key, collection_key):
                # 데이터베이스 상태 업데이트
                cur.execute("""
                    UPDATE paper_collection 
                    SET zotero_item_key = %s, synced_at = NOW()
                    WHERE paper_id = %s AND collection_id = %s
                """, (item_key, paper_id, collection_id))
                success_count += 1
            else:
                log_error(f"컬렉션 링크 실패: '{collection_name}'")
        
        return success_count
    
    def _get_unsynced_papers(self, cur):
        """미동기화 논문 목록 조회"""
        cur.execute("""
            SELECT DISTINCT pc.paper_id
            FROM paper_collection pc
            WHERE pc.synced_at IS NULL
            ORDER BY pc.paper_id
        """)
        
        return [row[0] for row in cur.fetchall()]
