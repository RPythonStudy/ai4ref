"""
Zotero 클라이언트 연결 및 기본 설정
- pyzotero 인스턴스 관리, 연결 검증
"""
import os
from pyzotero import zotero
from common.logger import log_info, log_error, log_debug

class ZoteroClient:
    def __init__(self, user_id=None, api_key=None):
        """Zotero 클라이언트 초기화"""
        self.user_id = user_id or os.getenv('ZOTERO_USER_ID')
        self.api_key = api_key or os.getenv('ZOTERO_API_KEY')
        
        if not self.user_id or not self.api_key:
            raise ValueError("Zotero 사용자 ID와 API 키가 필요합니다")
        
        self.zot = zotero.Zotero(self.user_id, 'user', self.api_key)
        log_debug(f"Zotero 클라이언트 초기화: User ID {self.user_id}")
    
    def test_connection(self):
        """Zotero API 연결 및 권한 검증 (테스트/운영 분리)"""
        try:
            user_info = self.zot.key_info()
            username = user_info.get('username', 'Unknown')
            log_info(f"Zotero 연결 성공: 사용자명={username}")

            access = user_info.get('access', {})
            user_access = access.get('user', {})
            read_permission = user_access.get('library', False)
            write_permission = user_access.get('write', False)

            if not read_permission:
                log_error("라이브러리 읽기 권한 없음")
            else:
                log_info("라이브러리 읽기 권한 확인")
            if not write_permission:
                log_error("라이브러리 쓰기 권한 없음")
            else:
                log_info("라이브러리 쓰기 권한 확인")

            # 테스트 모드: 실제 데이터 업로드 없이 연결·권한만 검증
            if read_permission and write_permission:
                log_info("연결 및 권한 검증 완료 (테스트 모드)")
                return True
            else:
                log_error("권한 부족: 업로드 불가")
                return False

        except Exception as e:
            log_error(f"Zotero 연결/권한 검증 실패: {e}")
            return False
    
    def get_library_info(self):
        """라이브러리 기본 정보 조회"""
        try:
            # 전체 아이템 수
            items = self.zot.items(limit=1)
            total_items = self.zot.num_items()  # 메서드 호출로 수정
            
            # 컬렉션 수
            collections = self.zot.collections()
            total_collections = len(collections) if collections else 0
            
            log_info(f"라이브러리 현황: 아이템 {total_items}개, 컬렉션 {total_collections}개")
            return {
                'total_items': total_items,
                'total_collections': total_collections,
                'collections': collections
            }
            
        except Exception as e:
            log_error(f"라이브러리 정보 조회 실패: {e}")
            return None
    
    def get_client(self):
        """pyzotero 인스턴스 반환"""
        return self.zot
