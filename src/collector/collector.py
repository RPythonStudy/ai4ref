import os
from dotenv import load_dotenv
import requests
from common.logger import log_info, log_debug
import db

load_dotenv()
API_KEY = os.environ.get("PUBMED_API_KEY")

def collect_all_pmids(term, retmax=1000):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    pmids, retstart = [], 0
    while True:
        params = {"db":"pubmed", "term":term, "retmax":retmax, "retstart":retstart, "retmode":"json"}
        if API_KEY: params["api_key"] = API_KEY
        log_info(f"API 호출: retstart={retstart}, retmax={retmax}")
        r = requests.get(url, params=params)
        j = r.json() if r.status_code==200 else None
        ids = j["esearchresult"]["idlist"] if j else []
        pmids += ids
        log_debug(f"누적 PMID: {len(pmids)}")
        if not ids or len(pmids) >= int(j["esearchresult"]["count"]): break
        retstart += retmax
    log_info(f"총 PMID 개수: {len(pmids)}")
    return pmids

if __name__ == "__main__":
    import sys
    term = (
        sys.argv[1]
        if len(sys.argv) > 1
        else '("generative AI" OR ChatGPT) AND ("scientific writing" OR "academic writing") AND (productivity OR efficiency)'
    )
    
    # PMID 수집
    pmids = collect_all_pmids(term)
    
    # 데이터베이스에 저장
    result = db.save_pmids_to_db(pmids)
    
    # 결과 확인 (선택사항)
    if result:
        log_info(f"작업 완료 - 수집: {len(pmids)}개, 저장: {result['inserted']}개")
    else:
        log_info("데이터베이스 저장 실패")
