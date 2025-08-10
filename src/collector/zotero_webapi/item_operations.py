"""
Zotero 아이템 관련 기능
- 아이템 생성, 검색, 변환
"""
import os
from dotenv import load_dotenv
from pyzotero import zotero
from common.logger import log_info, log_error, log_debug

load_dotenv()

class ZoteroItemManager:
    def __init__(self, user_id=None, api_key=None):
        self.user_id = user_id or os.getenv("ZOTERO_USER_ID")
        self.api_key = api_key or os.getenv("ZOTERO_API_KEY")
        self.zot = zotero.Zotero(self.user_id, 'user', self.api_key)
    
    def build_zotero_item(self, paper):
        """논문 데이터를 Zotero 아이템으로 변환"""
        authors = []
        if paper.get("authors"):
            for author in paper["authors"].split(";"):
                parts = author.strip().split(" ", 1)
                if len(parts) == 2:
                    authors.append({
                        "creatorType": "author", 
                        "firstName": parts[0], 
                        "lastName": parts[1]
                    })
        
        item = {
            "itemType": "journalArticle",
            "title": paper["title"] or "",
            "abstractNote": paper.get("abstract", "") or "",
            "creators": authors,
            "date": str(paper["publication_year"]) if paper.get("publication_year") else "",
            "DOI": paper.get("doi", "") or "",
            "extra": f"PMID: {paper['pmid']}"
        }
        
        # 저널 정보 추가
        if paper.get("journal"):
            item["publicationTitle"] = paper["journal"]
            
        return item
    
    def get_existing_item_by_pmid(self, pmid):
        """PMID로 기존 아이템 검색"""
        try:
            # 전체 라이브러리에서 PMID 검색
            items = self.zot.items(q=f"PMID: {pmid}")
            for item in items:
                extra = item['data'].get('extra', '')
                if f"PMID: {pmid}" in extra:
                    log_debug(f"기존 아이템 발견: PMID {pmid} -> {item['key']}")
                    return item['key']
            return None
        except Exception as e:
            log_debug(f"아이템 검색 오류 PMID {pmid}: {e}")
            return None
    
    def create_item(self, paper):
        """새 아이템 생성"""
        try:
            item = self.build_zotero_item(paper)
            result = self.zot.create_items([item])
            
            if result and 'successful' in result and result['successful']:
                item_key = list(result['successful'].keys())[0]
                log_info(f"아이템 생성 성공: PMID {paper['pmid']} -> {item_key}")
                return item_key
            else:
                log_error(f"아이템 생성 실패: PMID {paper['pmid']}")
                return None
                
        except Exception as e:
            log_error(f"아이템 생성 오류 PMID {paper['pmid']}: {e}")
            return None
    
    def get_item(self, item_key):
        """아이템 키로 전체 아이템 조회"""
        try:
            log_debug(f"아이템 조회 시도: {item_key}")
            item = self.zot.item(item_key)
            log_debug(f"아이템 조회 결과 타입: {type(item)}")
            
            if not item:
                log_debug(f"아이템 조회 실패: {item_key}")
                return None
            
            # item이 리스트인 경우 첫 번째 요소 사용
            if isinstance(item, list):
                if not item:
                    log_debug(f"빈 아이템 리스트: {item_key}")
                    return None
                item = item[0]
                log_debug(f"리스트에서 첫 번째 아이템 추출: {type(item)}")
            
            return item
        except Exception as e:
            log_debug(f"아이템 조회 오류: {e}")
            return None
