# fetcher.py - PubMed 논문 상세 정보 수집
import os
from dotenv import load_dotenv
import requests
import xml.etree.ElementTree as ET
from common.logger import log_info, log_error
from pmid_filter import get_connection

load_dotenv()
PUBMED_API_KEY = os.environ.get("PUBMED_API_KEY")

def fetch_paper_details(pmids, batch_size=200):
    """PMID 목록으로부터 상세 논문 정보 수집"""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    papers = []
    
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
            "rettype": "abstract"
        }
        if PUBMED_API_KEY:
            params["api_key"] = PUBMED_API_KEY
            
        log_info(f"논문 정보 수집: {i+1}-{min(i+batch_size, len(pmids))}/{len(pmids)}")
        
        try:
            r = requests.get(url, params=params)
            if r.status_code == 200:
                root = ET.fromstring(r.text)
                for article in root.findall(".//PubmedArticle"):
                    paper = extract_paper_info(article)
                    if paper:
                        papers.append(paper)
        except Exception as e:
            log_error(f"배치 {i} 처리 오류: {e}")
    
    log_info(f"논문 정보 수집 완료: {len(papers)}개")
    return papers

def extract_paper_info(article):
    """XML에서 논문 정보 추출"""
    try:
        pmid = article.find(".//PMID").text
        title = article.find(".//ArticleTitle")
        title_text = title.text if title is not None else ""
        
        abstract = article.find(".//AbstractText")
        abstract_text = abstract.text if abstract is not None else ""
        
        year = article.find(".//PubDate/Year")
        year_text = year.text if year is not None else ""
        
        # DOI 추출
        doi = ""
        for article_id in article.findall(".//ArticleId"):
            if article_id.get("IdType") == "doi":
                doi = article_id.text
                break
        
        # 저자 정보
        authors = []
        for author in article.findall(".//Author"):
            lastname = author.find("LastName")
            forename = author.find("ForeName")
            if lastname is not None and forename is not None:
                authors.append(f"{forename.text} {lastname.text}")
        
        return {
            "pmid": pmid,
            "title": title_text,
            "abstract": abstract_text,
            "publication_year": int(year_text) if year_text and year_text.isdigit() else None,
            "authors": "; ".join(authors),
            "doi": doi
        }
    except Exception as e:
        log_error(f"논문 정보 추출 오류: {e}")
        return None

def update_papers_with_details(papers):
    """수집한 상세 정보로 papers 테이블 업데이트"""
    try:
        conn = get_connection()
        cur = conn.cursor()
        
        updated_count = 0
        
        for paper in papers:
            try:
                cur.execute("""
                    UPDATE papers 
                    SET title = %s, abstract = %s, publication_year = %s, 
                        authors = %s, doi = %s, updated_at = NOW()
                    WHERE source = 'pubmed' AND source_id = %s
                """, (
                    paper["title"], paper["abstract"], paper["publication_year"],
                    paper["authors"], paper["doi"], paper["pmid"]
                ))
                
                if cur.rowcount > 0:
                    updated_count += 1
                    
            except Exception as e:
                log_error(f"논문 {paper['pmid']} 업데이트 오류: {e}")
        
        conn.commit()
        cur.close()
        conn.close()
        
        log_info(f"논문 상세정보 업데이트 완료: {updated_count}개")
        return {"updated": updated_count}
        
    except Exception as e:
        log_error(f"논문 업데이트 중 오류: {e}")
        return None

def get_pmids_without_details():
    """상세 정보가 없는 PMID 목록 조회 (title이 없는 것들)"""
    try:
        conn = get_connection()
        cur = conn.cursor()
        
        cur.execute("""
            SELECT source_id FROM papers 
            WHERE source = 'pubmed' AND (title IS NULL OR title = '')
        """)
        
        pmids = [row[0] for row in cur.fetchall()]
        
        cur.close()
        conn.close()
        
        log_info(f"상세정보 없는 PMID: {len(pmids)}개")
        return pmids
        
    except Exception as e:
        log_error(f"PMID 조회 중 오류: {e}")
        return []

if __name__ == "__main__":
    # 상세 정보가 없는 PMID들 가져오기
    pmids = get_pmids_without_details()
    
    if pmids:
        # 논문 상세 정보 수집
        papers = fetch_paper_details(pmids)
        
        if papers:
            # papers 테이블 업데이트
            result = update_papers_with_details(papers)
            
            if result:
                log_info(f"작업 완료: {result['updated']}개 논문 정보 업데이트")
        else:
            log_info("수집된 논문 정보가 없습니다")
    else:
        log_info("상세 정보 수집이 필요한 PMID가 없습니다")