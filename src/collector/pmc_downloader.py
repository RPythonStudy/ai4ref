import asyncio
import aiohttp
import aiofiles
import os
from pathlib import Path
from dotenv import load_dotenv
from common.logger import log_info, log_error, log_debug
from common.database import get_db_connection

# 프로젝트 루트 기준으로 .env 로드
PROJECT_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(dotenv_path=PROJECT_ROOT / ".env")


class PMCDownloader:
    def __init__(self, download_dir=None, batch_size=20, max_concurrent=10):
        # .env의 AI4REF_PDF_DIR 환경변수를 우선 사용
        default_dir = os.getenv("AI4REF_PDF_DIR", "data/pdf")
        self.download_dir = Path(download_dir or default_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)

        self.batch_size = batch_size
        self.semaphore = asyncio.Semaphore(max_concurrent)
        
        # PMC PDF URL 패턴 (우선순위 순)
        self.url_patterns = [
            "https://europepmc.org/articles/{pmcid}?pdf=render",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
        ]
        
    def get_pmc_papers(self):
        """PMCID 있고 PDF 없는 논문들 조회"""
        conn = get_db_connection()
        cur = conn.cursor()
        
        cur.execute("""
            SELECT pmid, pmcid, title 
            FROM papers 
            WHERE pmcid IS NOT NULL 
            AND (pdf_url IS NULL OR pdf_url = '')
            ORDER BY pmid
        """)
        
        papers = cur.fetchall()
        cur.close()
        conn.close()
        
        log_info(f"PMC 다운로드 대상: {len(papers)}개")
        return papers
        
    async def download_single_pdf(self, session, paper):
        """단일 논문 PDF 다운로드"""
        pmid, pmcid, title = paper

        async with self.semaphore:
            try:
                clean_pmcid = pmcid.replace('PMC', '')

                for pattern in self.url_patterns:
                    pdf_url = pattern.format(pmcid=f"PMC{clean_pmcid}")
                    log_debug(f"PDF 다운로드 시도: {pmid} -> {pdf_url}")

                    try:
                        async with session.get(pdf_url) as response:
                            if response.status == 200:
                                content_type = response.headers.get('content-type', '')
                                if 'pdf' in content_type.lower():
                                    success = await self._download_pdf_content(response, pmid)
                                    if success:
                                        await self._update_database(pmid, pmcid, str(response.url))
                                        log_info(f"PDF 다운로드 성공: {pmid}")
                                        return {'pmid': pmid, 'status': 'success', 'url': str(response.url)}
                                else:
                                    log_debug(f"PDF 아님 {pmid}: {content_type}")
                            else:
                                log_debug(f"HTTP 오류 {pmid}: {response.status}")

                    except Exception as e:
                        log_debug(f"URL 시도 실패 {pmid}: {e}")
                        continue

                log_debug(f"PDF 다운로드 실패: {pmid} (모든 URL 실패)")
                return {'pmid': pmid, 'status': 'failed', 'error': 'all_urls_failed'}

            except Exception as e:
                log_error(f"PDF 다운로드 오류 {pmid}: {e}")
                return {'pmid': pmid, 'status': 'error', 'error': str(e)}
                
    async def _download_pdf_content(self, response, pmid):
        """응답에서 PDF 콘텐츠 다운로드"""
        try:
            filename = f"{pmid}_pmc.pdf"
            filepath = self.download_dir / filename
            
            async with aiofiles.open(filepath, 'wb') as f:
                async for chunk in response.content.iter_chunked(8192):
                    await f.write(chunk)
            
            if filepath.stat().st_size > 1024:
                return True
            else:
                filepath.unlink()
                return False
                
        except Exception as e:
            log_debug(f"파일 다운로드 실패 {pmid}: {e}")
            return False
            
    async def _update_database(self, pmid, pmcid, pdf_url):
        """데이터베이스에 PDF 정보 업데이트"""
        try:
            filename = f"{pmid}_pmc.pdf"
            local_path = str(self.download_dir / filename)

            conn = get_db_connection()
            cur = conn.cursor()

            log_debug(f"DB 업데이트 직전 값: pmid={pmid}, pmcid={pmcid}, pdf_url(원본)={pdf_url}, pdf_path(로컬)={local_path}")
            cur.execute("""
                UPDATE papers 
                SET pdf_path = %s, pdf_source = 'PMC', updated_at = NOW()
                WHERE pmid = %s
            """, (local_path, pmid))

            conn.commit()
            cur.close()
            conn.close()

            conn = get_db_connection()
            cur = conn.cursor()
            cur.execute("SELECT pdf_url, pdf_path, pdf_source FROM papers WHERE pmid = %s OR pmcid = %s", (pmid, pmcid))
            rows = cur.fetchall()
            for row in rows:
                log_debug(f"DB 상태: pmid={pmid}, pdf_url={row[0]}, pdf_path={row[1]}, pdf_source={row[2]}")
            cur.close()
            conn.close()

        except Exception as e:
            log_error(f"DB 업데이트 실패 {pmid}: {e}")
            
    async def download_batch(self, papers):
        """논문 배치 다운로드"""
        if not papers:
            log_info("다운로드할 논문이 없습니다")
            return []
            
        log_info(f"PMC 배치 다운로드 시작: {len(papers)}개")
        
        timeout = aiohttp.ClientTimeout(total=30)
        connector = aiohttp.TCPConnector(limit=50, limit_per_host=10)
        headers = {
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        
        async with aiohttp.ClientSession(timeout=timeout, connector=connector, headers=headers) as session:
            tasks = [self.download_single_pdf(session, paper) for paper in papers]
            results = await asyncio.gather(*tasks, return_exceptions=True)
            
        success_count = sum(1 for r in results if isinstance(r, dict) and r.get('status') == 'success')
        failed_count = len(results) - success_count
        
        log_info(f"PMC 배치 완료: 성공 {success_count}, 실패 {failed_count}")
        
        return results
        
    async def download_all(self):
        """모든 PMC 논문 다운로드"""
        log_info("=== PMC PDF 다운로드 시작 ===")
        
        papers = self.get_pmc_papers()
        if not papers:
            log_info("다운로드할 PMC 논문이 없습니다")
            return
            
        all_results = []
        for i in range(0, len(papers), self.batch_size):
            batch = papers[i:i + self.batch_size]
            batch_num = i // self.batch_size + 1
            total_batches = (len(papers) + self.batch_size - 1) // self.batch_size
            
            log_info(f"배치 {batch_num}/{total_batches} 처리 중... ({len(batch)}개)")
            
            batch_results = await self.download_batch(batch)
            all_results.extend(batch_results)
            
            if i + self.batch_size < len(papers):
                await asyncio.sleep(2)
                
        total_success = sum(1 for r in all_results if isinstance(r, dict) and r.get('status') == 'success')
        success_rate = total_success / len(papers) * 100 if papers else 0
        
        log_info(f"=== PMC 다운로드 완료 ===")
        log_info(f"전체: {len(papers)}개, 성공: {total_success}개 ({success_rate:.1f}%)")
        
        return all_results

def show_download_status():
    """현재 PDF 다운로드 상태 조회"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    cur.execute("SELECT COUNT(*) FROM papers")
    total = cur.fetchone()[0]
    
    cur.execute("SELECT COUNT(*) FROM papers WHERE pmcid IS NOT NULL")
    pmcid_count = cur.fetchone()[0]
    
    cur.execute("SELECT COUNT(*) FROM papers WHERE pdf_url IS NOT NULL")
    pdf_count = cur.fetchone()[0]
    
    cur.execute("SELECT COUNT(*) FROM papers WHERE pdf_source = 'PMC'")
    pmc_pdf_count = cur.fetchone()[0]
    
    log_info(f"=== PDF 다운로드 현황 ===")
    log_info(f"전체 논문: {total}개")
    log_info(f"PMCID 보유: {pmcid_count}개")
    log_info(f"PDF 확보: {pdf_count}개")
    log_info(f"PMC PDF: {pmc_pdf_count}개")
    
    if pmcid_count > 0:
        pmc_rate = pmc_pdf_count / pmcid_count * 100
        log_info(f"PMC 성공률: {pmc_rate:.1f}%")
    
    cur.close()
    conn.close()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PMC PDF 다운로더")
    parser.add_argument('--status', action='store_true', help='다운로드 현황 조회')
    parser.add_argument('--batch-size', type=int, default=20, help='배치 크기')
    parser.add_argument('--max-concurrent', type=int, default=10, help='최대 동시 다운로드')
    parser.add_argument('--dir', type=str, help='PDF 저장 디렉토리 (없으면 .env 설정 사용)')
    
    args = parser.parse_args()
    
    if args.status:
        show_download_status()
    else:
        downloader = PMCDownloader(
            download_dir=args.dir,
            batch_size=args.batch_size,
            max_concurrent=args.max_concurrent
        )
        asyncio.run(downloader.download_all())
