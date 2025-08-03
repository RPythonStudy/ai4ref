import sys
import requests
from common.logger import log_info, log_debug

def collect_papers(term: str):
    retmax = 10
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": retmax,
        "retmode": "json"
    }
    log_debug(f"API 호출: {esearch_url}")
    response = requests.get(esearch_url, params=params)
    
    if response.status_code == 200:
        json_data = response.json()
        log_info(f"검색 완료, 상태코드: {response.status_code}")
        log_debug(f"응답 데이터: {json_data}")
        return json_data
    else:
        log_info(f"API 호출 실패, 상태코드: {response.status_code}")
        return None

if __name__ == "__main__":
    term = sys.argv[1] if len(sys.argv) > 1 else '("generative AI" OR ChatGPT) AND ("scientific writing" OR "academic writing") AND (productivity OR efficiency)'
    collect_papers(term)
