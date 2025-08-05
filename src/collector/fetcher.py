import os
from dotenv import load_dotenv
import requests
import xml.etree.ElementTree as ET
import psycopg2
from common.logger import log_info, log_error
import db

load_dotenv()
API_KEY = os.environ.get("PUBMED_API_KEY")

def fetch_paper_details(pmids, batch_size=200):
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
        if API_KEY:
            params["api_key"] = API_KEY
            
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
        
        authors = []
        for author in article.findall(".//Author"):
            lastname = author.find("LastName")
            forename = author.find("ForeName")
            if lastname is not None and forename is not None:
                authors.append(f"{forename.text} {lastname.text}")
        
        # 첫 번째 저자
        first_author = authors[0] if authors else ""
        
        return {
            "pmid": pmid,
            "title": title_text,
            "abstract": abstract_text,
            "year": year_text,
            "authors": "; ".join(authors),
            "doi": doi,
            "first_author": first_author
        }
    except Exception as e:
        log_error(f"논문 정보 추출 오류: {e}")
        return None

def save_papers_to_db(papers):
    try:
        conn = psycopg2.connect(
            host=db.DB_HOST, port=db.DB_PORT, dbname=db.DB_NAME, 
            user=db.DB_USER, password=db.DB_PASS
        )
        cur = conn.cursor()
        
        cur.execute("""
            CREATE TABLE IF NOT EXISTS pubmed_papers (
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                year TEXT,
                authors TEXT,
                doi TEXT,
                first_author TEXT,
                fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        """)
        
        # 기존 테이블에 새 컬럼 추가 (이미 있으면 무시)
        try:
            cur.execute("ALTER TABLE pubmed_papers ADD COLUMN doi TEXT;")
        except:
            pass
        try:
            cur.execute("ALTER TABLE pubmed_papers ADD COLUMN first_author TEXT;")
        except:
            pass
        
        # DOI 고유 인덱스 생성 (NULL 값 제외)
        try:
            cur.execute("""
                CREATE UNIQUE INDEX idx_papers_doi 
                ON pubmed_papers(doi) WHERE doi IS NOT NULL AND doi != '';
            """)
        except:
            pass
        
        # Title + First Author + Year 복합 인덱스 생성
        try:
            cur.execute("""
                CREATE UNIQUE INDEX idx_papers_title_author_year 
                ON pubmed_papers(title, first_author, year) 
                WHERE title IS NOT NULL AND title != '' 
                AND first_author IS NOT NULL AND first_author != ''
                AND year IS NOT NULL AND year != '';
            """)
        except:
            pass
        
        inserted_count = 0
        duplicate_count = 0
        
        for paper in papers:
            try:
                # PMID 중복 검사 (기본 키)
                cur.execute("SELECT pmid FROM pubmed_papers WHERE pmid = %s;", (paper["pmid"],))
                if cur.fetchone():
                    duplicate_count += 1
                    continue
                
                # DOI 중복 검사
                if paper["doi"]:
                    cur.execute("SELECT pmid FROM pubmed_papers WHERE doi = %s;", (paper["doi"],))
                    if cur.fetchone():
                        duplicate_count += 1
                        continue
                
                # Title + First Author + Year 중복 검사
                if paper["title"] and paper["first_author"] and paper["year"]:
                    cur.execute("""
                        SELECT pmid FROM pubmed_papers 
                        WHERE title = %s AND first_author = %s AND year = %s;
                    """, (paper["title"], paper["first_author"], paper["year"]))
                    if cur.fetchone():
                        duplicate_count += 1
                        continue
                
                # 중복이 없으면 삽입
                cur.execute("""
                    INSERT INTO pubmed_papers (pmid, title, abstract, year, authors, doi, first_author)
                    VALUES (%s, %s, %s, %s, %s, %s, %s);
                """, (paper["pmid"], paper["title"], paper["abstract"], 
                      paper["year"], paper["authors"], paper["doi"], paper["first_author"]))
                
                inserted_count += 1
                
            except Exception as e:
                log_error(f"논문 {paper['pmid']} 저장 오류: {e}")
                duplicate_count += 1
                # 트랜잭션 롤백 후 재시작
                conn.rollback()
        
        conn.commit()
        cur.close()
        conn.close()
        
        log_info(f"논문 저장 완료: 신규 {inserted_count}개, 중복 {duplicate_count}개")
        return {"inserted": inserted_count, "duplicates": duplicate_count}
        
    except Exception as e:
        log_error(f"데이터베이스 작업 오류: {e}")
        return None

if __name__ == "__main__":
    # DB에서 PMID 가져오기
    pmids = db.get_all_pmids()
    if pmids:
        papers = fetch_paper_details(pmids)
        if papers:
            save_papers_to_db(papers)
    else:
        log_info("저장된 PMID가 없습니다")
