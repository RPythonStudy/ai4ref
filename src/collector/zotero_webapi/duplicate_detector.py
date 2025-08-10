"""
Zotero 중복 감지 관련 기능
- PMID 기반 중복 감지, 해결
"""
from .item_operations import ZoteroItemManager
from common.logger import log_info, log_error, log_debug

class ZoteroDuplicateDetector:
    def __init__(self, user_id=None, api_key=None):
        self.item_manager = ZoteroItemManager(user_id, api_key)
    
    def detect_existing_item(self, pmid):
        """PMID로 기존 아이템 존재 여부 확인"""
        return self.item_manager.get_existing_item_by_pmid(pmid)
    
    def resolve_item_conflict(self, paper_data, existing_key=None, stored_key=None):
        """
        아이템 충돌 해결
        - Zotero에 실제 존재하는 아이템 우선 사용
        - 없으면 새로 생성
        """
        pmid = paper_data['pmid']
        log_debug(f"아이템 충돌 해결 시작: PMID {pmid}")
        
        # 1. Zotero에서 실제 존재 여부 확인
        actual_existing = self.detect_existing_item(pmid)
        
        if actual_existing:
            # Zotero에 실제 존재하는 아이템 사용
            log_debug(f"Zotero 기존 아이템 재사용: {actual_existing}")
            return actual_existing, "existing"
        
        # 2. Zotero에 없으므로 새로 생성
        log_debug(f"Zotero에 없음. 새 아이템 생성: PMID {pmid}")
        new_key = self.item_manager.create_item(paper_data)
        
        if new_key:
            log_debug(f"새 아이템 생성 완료: {new_key}")
            return new_key, "created"
        else:
            log_error(f"아이템 생성 실패: PMID {pmid}")
            return None, "failed"
    
    def validate_item_consistency(self, item_key, expected_pmid):
        """아이템의 PMID 일치성 검증"""
        try:
            item = self.item_manager.get_item(item_key)
            if not item:
                return False
            
            extra = item.get('data', {}).get('extra', '')
            return f"PMID: {expected_pmid}" in extra
        except Exception as e:
            log_debug(f"아이템 일치성 검증 오류: {e}")
            return False
