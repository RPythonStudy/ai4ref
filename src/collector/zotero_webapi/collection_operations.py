"""
Zotero 컬렉션 관련 기능
- 컬렉션 생성, 검색, 아이템 추가
"""
import os
from dotenv import load_dotenv
from pyzotero import zotero
from common.logger import log_info, log_error, log_debug

load_dotenv()

class ZoteroCollectionManager:
    def __init__(self, user_id=None, api_key=None):
        self.user_id = user_id or os.getenv("ZOTERO_USER_ID")
        self.api_key = api_key or os.getenv("ZOTERO_API_KEY")
        self.zot = zotero.Zotero(self.user_id, 'user', self.api_key)
    
    def get_or_create_collection(self, collection_name):
        """컬렉션 조회 또는 생성"""
        try:
            collections = self.zot.collections()
            for collection in collections:
                if collection['data']['name'] == collection_name:
                    log_debug(f"기존 컬렉션: '{collection_name}' -> {collection['key']}")
                    return collection['key']
            
            # 컬렉션 생성
            log_info(f"컬렉션 '{collection_name}' 생성 중...")
            result = self.zot.create_collections([{"name": collection_name}])
            
            if result and 'successful' in result:
                collection_key = list(result['successful'].keys())[0]
                log_info(f"컬렉션 생성 성공: '{collection_name}' -> {collection_key}")
                return collection_key
            else:
                log_error(f"컬렉션 생성 실패: '{collection_name}'")
                return None
                
        except Exception as e:
            log_error(f"컬렉션 처리 오류: {e}")
            return None
    
    def add_item_to_collection(self, item_key, collection_key, item_object=None):
        """기존 아이템을 컬렉션에 추가"""
        try:
            # pyzotero의 addto_collection은 전체 아이템 객체가 필요함
            if not item_object:
                # 아이템 객체가 없으면 조회
                from .item_operations import ZoteroItemManager
                item_manager = ZoteroItemManager(self.user_id, self.api_key)
                item_object = item_manager.get_item(item_key)
                
            if not item_object:
                log_debug(f"아이템 객체 없음: {item_key}")
                return False
                
            result = self.zot.addto_collection(collection_key, item_object)
            if result:
                log_debug(f"컬렉션 추가 성공: {item_key} -> {collection_key}")
                return True
            else:
                log_debug(f"컬렉션 추가 실패: {item_key} -> {collection_key}")
                return False
        except Exception as e:
            log_debug(f"컬렉션 추가 오류: {e}")
            return False
    
    def list_collections(self):
        """모든 컬렉션 목록 조회"""
        try:
            collections = self.zot.collections()
            return [(col['key'], col['data']['name']) for col in collections]
        except Exception as e:
            log_error(f"컬렉션 목록 조회 오류: {e}")
            return []
