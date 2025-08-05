import os
import json
import requests
from dotenv import load_dotenv
from db import get_connection
from common.logger import log_info, log_error

load_dotenv()
API_KEY = os.getenv("ZOTERO_API_KEY")
USER_ID = os.getenv("ZOTERO_USER_ID")
ZOTERO_API_URL = f"https://api.zotero.org/users/{USER_ID}/items"

HEADERS = {
    "Authorization": f"Bearer {API_KEY}",
    "Content-Type": "application/json"
}

def build_zotero_item(paper):
    authors = []
    for author in paper["authors"].split(";"):
        parts = author.strip().split(" ", 1)
        if len(parts) == 2:
            authors.append({"creatorType": "author", "firstName": parts[0], "lastName": parts[1]})
    return {
        "itemType": "journalArticle",
        "title": paper["title"],
        "abstractNote": paper["abstract"],
        "creators": authors,
        "date": paper["year"],
        "DOI": paper["doi"],
        "extra": f"PMID: {paper['pmid']}"
    }

def upload_to_zotero_web(papers):
    success = 0
    for paper in papers:
        item = build_zotero_item(paper)
        try:
            r = requests.post(ZOTERO_API_URL, headers=HEADERS, json=[item])
            if r.status_code in [200, 201]:
                log_info(f"Zotero 저장 성공: PMID {paper['pmid']}")
                success += 1
            else:
                log_error(f"Zotero 저장 실패: {paper['pmid']} - {r.status_code}: {r.text}")
        except Exception as e:
            log_error(f"예외 발생: {e}")
    log_info(f"Zotero 저장 완료: 총 {success}개 성공")

def fetch_papers_from_db(limit=None):
    conn = get_connection()
    cur = conn.cursor()
    query = "SELECT pmid, title, abstract, year, authors, doi, first_author FROM pubmed_papers"
    if limit:
        query += " LIMIT %s"
        cur.execute(query, (limit,))
    else:
        cur.execute(query)
    rows = cur.fetchall()
    columns = ["pmid", "title", "abstract", "year", "authors", "doi", "first_author"]
    cur.close()
    conn.close()
    return [dict(zip(columns, row)) for row in rows]

if __name__ == "__main__":
    papers = fetch_papers_from_db()
    if papers:
        upload_to_zotero_web(papers)
    else:
        log_info("Zotero에 업로드할 논문이 없습니다.")
