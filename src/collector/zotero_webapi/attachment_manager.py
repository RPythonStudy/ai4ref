"""
Zotero PDF 첨부 관련 기능
- PDF 파일 첨부, 확인
"""
import os
from pathlib import Path
from dotenv import load_dotenv
from pyzotero import zotero
from common.logger import log_info, log_error, log_debug

load_dotenv()

class ZoteroAttachmentManager:
    def __init__(self, user_id=None, api_key=None):
        self.user_id = user_id or os.getenv("ZOTERO_USER_ID")
        self.api_key = api_key or os.getenv("ZOTERO_API_KEY")
        self.zot = zotero.Zotero(self.user_id, 'user', self.api_key)
    
    def attach_pdf_to_item(self, item_key, pdf_path, pmid):
        """아이템에 PDF 첨부"""
        try:
            # PDF 파일 존재 확인
            if not self.check_pdf_exists(pdf_path):
                return False
            
            pdf_file = Path(pdf_path)
            
            # PDF 첨부 업로드
            result = self.zot.attachment_simple([str(pdf_file)], item_key)
            if result:
                log_info(f"PDF 첨부 성공: PMID {pmid} -> {pdf_file.name}")
                return True
            else:
                log_error(f"PDF 첨부 실패: PMID {pmid}")
                return False
                
        except Exception as e:
            log_error(f"PDF 첨부 오류 PMID {pmid}: {e}")
            return False
    
    def check_pdf_exists(self, pdf_path):
        """PDF 파일 존재 확인"""
        try:
            pdf_file = Path(pdf_path)
            if not pdf_file.exists():
                log_debug(f"PDF 파일 없음: {pdf_path}")
                return False
            
            if not pdf_file.suffix.lower() == '.pdf':
                log_debug(f"PDF 파일이 아님: {pdf_path}")
                return False
                
            return True
        except Exception as e:
            log_debug(f"PDF 확인 오류: {e}")
            return False
    
    def get_item_attachments(self, item_key):
        """아이템의 첨부파일 목록 조회"""
        try:
            attachments = self.zot.children(item_key)
            pdf_attachments = [
                att for att in attachments 
                if att['data'].get('contentType') == 'application/pdf'
            ]
            return pdf_attachments
        except Exception as e:
            log_debug(f"첨부파일 조회 오류: {e}")
            return []
