import os
import sys
import psycopg2
from dotenv import load_dotenv
from Bio import Entrez
import xml.etree.ElementTree as ET
from common.logger import log_info, log_debug, log_error
from common.database import get_db_connection

load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL", "r.python.ai@gmail.com")
Entrez.api_key = os.getenv("PUBMED_API_KEY")

def fetch_pmids_from_filtered():
    """filtered_pmid 테이블에서 처리되지 않은 PMID 목록을 가져옵니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        cur.execute("""
            SELECT pmid FROM filtered_pmid 
            WHERE fetched = false 
            ORDER BY created_at
        """)
        
        pmids = [row[0] for row in cur.fetchall()]
        
        cur.close()
        conn.close()
        
        log_info(f"처리 대상 PMID: {len(pmids)}개")
        return pmids
        
    except psycopg2.Error as e:
        log_error(f"filtered_pmid 조회 오류: {e}")
        return []

def pubmed_efetch(pmids, batch_size=200):
    """PubMed efetch로 논문 상세 정보를 가져옵니다."""
    papers = []
    
    for i in range(0, len(pmids), batch_size):
        batch_pmids = pmids[i:i+batch_size]
        log_info(f"배치 {i//batch_size + 1}: PMID {len(batch_pmids)}개 처리 중")
        
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(batch_pmids),
                rettype="xml",
                retmode="xml"
            )
            
            xml_data = handle.read()
            handle.close()
            
            # XML 파싱
            root = ET.fromstring(xml_data)
            
            for article in root.findall(".//PubmedArticle"):
                paper_info = extract_paper_info(article)
                if paper_info:
                    papers.append(paper_info)
                    
        except Exception as e:
            log_error(f"배치 {i//batch_size + 1} efetch 오류: {e}")
            continue
    
    log_info(f"efetch 완료: {len(papers)}개 논문 정보 수집")
    return papers

def extract_paper_info(article):
    """PubmedArticle XML에서 논문 정보를 추출합니다."""
    try:
        # PMID 추출
        pmid_elem = article.find(".//PMID")
        if pmid_elem is None:
            return None
        pmid = pmid_elem.text
        
        # 제목 추출
        title_elem = article.find(".//ArticleTitle")
        title = title_elem.text if title_elem is not None else ""
        
        # 초록 추출
        abstract_parts = []
        for abstract_text in article.findall(".//AbstractText"):
            if abstract_text.text:
                abstract_parts.append(abstract_text.text)
        abstract = " ".join(abstract_parts) if abstract_parts else ""
        
        # 출간연도 추출
        year_elem = article.find(".//PubDate/Year")
        if year_elem is None:
            year_elem = article.find(".//PubDate/MedlineDate")
            if year_elem is not None and year_elem.text:
                # MedlineDate에서 연도만 추출 (예: "2023 Spring" -> "2023")
                year_text = year_elem.text.split()[0]
            else:
                year_text = None
        else:
            year_text = year_elem.text
        
        publication_year = int(year_text) if year_text and year_text.isdigit() else None
        
        # DOI, PMC ID, arXiv ID 추출
        doi = ""
        pmcid = ""
        arxiv_id = ""
        
        for article_id in article.findall(".//ArticleId"):
            id_type = article_id.get("IdType")
            if id_type == "doi":
                doi = article_id.text
            elif id_type == "pmc":
                pmcid = article_id.text
            elif id_type == "arxiv":
                arxiv_id = article_id.text
        
        # 저자 추출
        authors = []
        for author in article.findall(".//Author"):
            lastname = author.find("LastName")
            forename = author.find("ForeName")
            if lastname is not None and forename is not None:
                authors.append(f"{forename.text} {lastname.text}")
        
        authors_str = "; ".join(authors)
        
        return {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "authors": authors_str,
            "publication_year": publication_year,
            "doi": doi,
            "pmcid": pmcid,
            "arxiv_id": arxiv_id
        }
        
    except Exception as e:
        log_error(f"논문 정보 추출 오류: {e}")
        return None

def check_duplicate_papers(papers):
    """DOI, Title+Author+Year 중복 체크를 수행합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        unique_papers = []
        duplicate_count = 0
        
        for paper in papers:
            pmid = paper["pmid"]
            doi = paper["doi"]
            title = paper["title"]
            authors = paper["authors"]
            year = paper["publication_year"]
            
            is_duplicate = False
            duplicate_reason = ""
            
            # DOI 중복 체크 (DOI가 있는 경우)
            if doi:
                cur.execute("SELECT pmid FROM papers WHERE doi = %s", (doi,))
                existing = cur.fetchone()
                if existing:
                    is_duplicate = True
                    duplicate_reason = f"DOI 중복 (기존 PMID: {existing[0]})"
            
            # Title+Author+Year 중복 체크 (DOI 중복이 없는 경우)
            if not is_duplicate and title and authors and year:
                cur.execute("""
                    SELECT pmid FROM papers 
                    WHERE title = %s AND authors = %s AND publication_year = %s
                """, (title, authors, year))
                existing = cur.fetchone()
                if existing:
                    is_duplicate = True
                    duplicate_reason = f"Title+Author+Year 중복 (기존 PMID: {existing[0]})"
            
            if is_duplicate:
                log_debug(f"중복 논문 제외: PMID {pmid} - {duplicate_reason}")
                duplicate_count += 1
            else:
                unique_papers.append(paper)
        
        cur.close()
        conn.close()
        
        log_info(f"중복 체크 완료: 고유 {len(unique_papers)}개, 중복 제외 {duplicate_count}개")
        return unique_papers
        
    except psycopg2.Error as e:
        log_error(f"중복 체크 오류: {e}")
        return papers

def save_papers_to_db(papers):
    """논문 정보를 papers 테이블에 저장합니다."""
    if not papers:
        return 0
    
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        saved_count = 0
        
        for paper in papers:
            try:
                cur.execute("""
                    INSERT INTO papers (
                        source, source_id, pmid, title, abstract, authors, 
                        publication_year, doi, pmcid, arxiv_id, created_at, updated_at
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                """, (
                    'pubmed',
                    paper["pmid"],
                    paper["pmid"],
                    paper["title"],
                    paper["abstract"],
                    paper["authors"],
                    paper["publication_year"],
                    paper["doi"],
                    paper["pmcid"] if paper["pmcid"] else None,
                    paper["arxiv_id"] if paper["arxiv_id"] else None
                ))
                
                saved_count += 1
                
            except psycopg2.IntegrityError as e:
                # 중복 키 오류 등은 이미 중복 체크에서 걸러져야 하지만, 혹시 모를 경우 처리
                log_error(f"논문 저장 오류 (PMID: {paper['pmid']}): {e}")
                conn.rollback()
                continue
            except Exception as e:
                log_error(f"논문 저장 오류 (PMID: {paper['pmid']}): {e}")
                conn.rollback()
                continue
        
        conn.commit()
        cur.close()
        conn.close()
        
        log_info(f"papers 테이블에 저장 완료: {saved_count}개")
        return saved_count
        
    except psycopg2.Error as e:
        log_error(f"papers 저장 오류: {e}")
        return 0

def update_filtered_pmid_status(pmids, success=True):
    """filtered_pmid 테이블의 처리 상태를 업데이트합니다."""
    if not pmids:
        return
    
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        if success:
            # 성공적으로 처리된 경우
            pmid_list = "(" + ",".join([f"'{pmid}'" for pmid in pmids]) + ")"
            cur.execute(f"""
                UPDATE filtered_pmid 
                SET fetched = true, fetched_at = NOW() 
                WHERE pmid IN {pmid_list}
            """)
        else:
            # 실패한 경우
            pmid_list = "(" + ",".join([f"'{pmid}'" for pmid in pmids]) + ")"
            cur.execute(f"""
                UPDATE filtered_pmid 
                SET fetch_attempts = fetch_attempts + 1,
                    last_error = 'efetch 실패'
                WHERE pmid IN {pmid_list}
            """)
        
        conn.commit()
        cur.close()
        conn.close()
        
        status = "성공" if success else "실패"
        log_debug(f"filtered_pmid 상태 업데이트 ({status}): {len(pmids)}개")
        
    except psycopg2.Error as e:
        log_error(f"filtered_pmid 상태 업데이트 오류: {e}")

def show_efetch_stats():
    """efetch 처리 통계를 표시합니다."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        
        # filtered_pmid 처리 현황
        cur.execute("SELECT COUNT(*) FROM filtered_pmid WHERE fetched = true")
        fetched_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM filtered_pmid WHERE fetched = false")
        pending_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM filtered_pmid WHERE fetch_attempts > 0")
        error_count = cur.fetchone()[0]
        
        # papers 테이블 현황
        cur.execute("SELECT COUNT(*) FROM papers")
        total_papers = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM papers WHERE pmid IS NOT NULL")
        pubmed_papers = cur.fetchone()[0]
        
        log_info("=== efetch 처리 통계 ===")
        log_info(f"filtered_pmid 처리 완료: {fetched_count}개")
        log_info(f"filtered_pmid 처리 대기: {pending_count}개")
        log_info(f"filtered_pmid 오류: {error_count}개")
        log_info(f"papers 테이블 총 논문: {total_papers}개")
        log_info(f"papers 테이블 PubMed 논문: {pubmed_papers}개")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"통계 조회 오류: {e}")

if __name__ == "__main__":
    log_info("=== PubMed efetch 시작 ===")
    
    # filtered_pmid에서 처리할 PMID 목록 가져오기
    pmids = fetch_pmids_from_filtered()
    
    if not pmids:
        log_info("처리할 PMID가 없습니다")
        show_efetch_stats()
        sys.exit(0)
    
    # PubMed efetch로 논문 상세 정보 수집
    log_info(f"efetch 시작: {len(pmids)}개 PMID 처리")
    papers = pubmed_efetch(pmids)
    
    if not papers:
        log_error("efetch에서 논문 정보를 가져오지 못했습니다")
        update_filtered_pmid_status(pmids, success=False)
        sys.exit(1)
    
    # 중복 논문 체크 (DOI, Title+Author+Year)
    log_info("중복 논문 체크 시작")
    unique_papers = check_duplicate_papers(papers)
    
    # papers 테이블에 저장
    if unique_papers:
        log_info(f"papers 테이블에 저장 시작: {len(unique_papers)}개")
        saved_count = save_papers_to_db(unique_papers)
        
        # 성공적으로 처리된 PMID들의 상태 업데이트
        processed_pmids = [paper["pmid"] for paper in unique_papers]
        update_filtered_pmid_status(processed_pmids, success=True)
        
        log_info(f"efetch 완료: {saved_count}개 논문 저장")
    else:
        log_info("저장할 고유 논문이 없습니다 (모두 중복)")
    
    # 전체 통계 표시
    show_efetch_stats()
    
    log_info("=== PubMed efetch 완료 ===")