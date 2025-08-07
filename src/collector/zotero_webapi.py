import os
from dotenv import load_dotenv
from pyzotero import zotero
from pmid_filter import get_connection
from common.logger import log_info, log_error

load_dotenv()
ZOTERO_API_KEY = os.getenv("ZOTERO_API_KEY")
ZOTERO_USER_ID = os.getenv("ZOTERO_USER_ID")
ZOTERO_COLLECTION_NAME = "auto writting"

zot = zotero.Zotero(ZOTERO_USER_ID, 'user', ZOTERO_API_KEY)

def build_zotero_item(paper, collection_key):
    """논문 데이터를 Zotero 아이템으로 변환"""
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
        "date": str(paper["publication_year"]) if paper["publication_year"] else "",
        "DOI": paper["doi"] or "",
        "extra": f"PMID: {paper['source_id']}",
        "collections": [collection_key]
    }

def get_existing_pmids(collection_key):
    """컬렉션 내 기존 PMID 조회"""
    try:
        items = zot.collection_items(collection_key)
        pmids = set()
        for item in items:
            extra = item['data'].get('extra', '')
            if extra.startswith('PMID: '):
                pmids.add(extra.replace('PMID: ', '').strip())
        return pmids
    except Exception as e:
        log_error(f"PMID 조회 오류: {e}")
        return set()

def filter_duplicates(papers, collection_key):
    """중복 논문 제외"""
    existing_pmids = get_existing_pmids(collection_key)
    filtered = [p for p in papers if p['source_id'] not in existing_pmids]
    
    duplicate_count = len(papers) - len(filtered)
    if duplicate_count > 0:
        log_info(f"중복 제외: 전체 {len(papers)}개 → 신규 {len(filtered)}개, 중복 {duplicate_count}개")
    
    return filtered

def get_or_create_collection(name):
    """컬렉션 조회 또는 생성"""
    try:
        collections = zot.collections()
        for collection in collections:
            if collection['data']['name'] == name:
                log_info(f"기존 컬렉션 발견: '{name}' (key: {collection['key']})")
                return collection['key']
        
        # 컬렉션 생성
        log_info(f"컬렉션 '{name}' 생성 중...")
        new_collection = zot.create_collections([{"name": name}])
        
        # 생성 후 다시 조회하여 키 확인
        collections = zot.collections()
        for collection in collections:
            if collection['data']['name'] == name:
                collection_key = collection['key']
                log_info(f"컬렉션 생성 성공: '{name}' (key: {collection_key})")
                return collection_key
        
        log_error(f"컬렉션 생성 후 조회 실패: '{name}'")
        return None
        
    except Exception as e:
        log_error(f"컬렉션 처리 오류: {e}")
        return None

def upload_papers(papers, collection_key):
    """논문을 Zotero에 배치 업로드"""
    papers = filter_duplicates(papers, collection_key)
    
    if not papers:
        log_info("업로드할 신규 논문이 없습니다")
        return
    
    log_info(f"{len(papers)}개 논문 업로드 시작")
    
    try:
        # 50개씩 배치 처리
        batch_size = 50
        total_success = 0
        
        for i in range(0, len(papers), batch_size):
            batch = papers[i:i + batch_size]
            items = [build_zotero_item(paper, collection_key) for paper in batch]
            
            log_info(f"배치 {i//batch_size + 1}: {len(batch)}개 논문 업로드 중")
            created = zot.create_items(items)
            
            # 업로드 성공 확인 (컬렉션은 이미 아이템에 포함됨)
            if created and isinstance(created, dict) and 'successful' in created:
                log_info(f"배치 완료: {len(batch)}개 논문이 컬렉션에 추가됨")
                total_success += len(batch)
            else:
                log_error(f"배치 업로드 실패")
        
        log_info(f"전체 업로드 완료: {total_success}개 성공")
            
    except Exception as e:
        log_error(f"업로드 오류: {e}")

def fetch_papers_from_db(limit=None):
    """DB에서 논문 데이터 조회"""
    conn = get_connection()
    cur = conn.cursor()
    
    query = """
        SELECT source_id, title, abstract, publication_year, authors, doi 
        FROM papers 
        WHERE source = 'pubmed' 
        AND title IS NOT NULL 
        AND title != ''
        ORDER BY created_at DESC
    """
    
    if limit:
        query += " LIMIT %s"
        cur.execute(query, (limit,))
    else:
        cur.execute(query)
        
    rows = cur.fetchall()
    columns = ["source_id", "title", "abstract", "publication_year", "authors", "doi"]
    cur.close()
    conn.close()
    return [dict(zip(columns, row)) for row in rows]

if __name__ == "__main__":
    collection_key = get_or_create_collection(ZOTERO_COLLECTION_NAME)
    
    if not collection_key:
        log_error(f"컬렉션 '{ZOTERO_COLLECTION_NAME}' 처리 실패")
        exit(1)
    
    papers = fetch_papers_from_db()
    if papers:
        upload_papers(papers, collection_key)
    else:
        log_info("Zotero에 업로드할 논문이 없습니다")
