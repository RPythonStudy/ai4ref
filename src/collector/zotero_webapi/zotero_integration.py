from collector.zotero.zotero_client import ZoteroClient
from collector.zotero.batch_processor import ZoteroBatchProcessor
from common.logger import log_info, log_error

def main():
    """
    Zotero 전체 파이프라인 실행 (자체 호출용)
    - 연결 및 권한 검증
    - 라이브러리 정보 및 동기화 상태 조회
    - 논문 업로드
    - 최종 동기화 상태 출력
    """
    client = ZoteroClient()
    batch = ZoteroBatchProcessor()

    log_info("Zotero 연결 및 권한 검증 시작")
    if not client.test_connection():
        log_error("Zotero 연결 또는 권한 오류 - 프로그램 종료")
        return

    client.get_library_info()
    batch.show_sync_status()

    log_info("논문 업로드 시작")
    batch.upload_all_papers()

    log_info("최종 동기화 상태")
    batch.show_sync_status()

if __name__ == "__main__":
    main()
