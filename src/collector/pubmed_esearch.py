import os
import json
from dotenv import load_dotenv
from Bio import Entrez
from common.logger import log_info, log_debug
import pmid_filter

load_dotenv()
PUBMED_API_KEY = os.environ.get("PUBMED_API_KEY")

# Entrez 설정
Entrez.email = "ai4ref@example.com"  # NCBI 요구사항
if PUBMED_API_KEY:
    Entrez.api_key = PUBMED_API_KEY

def esearch_biopython(term, retmax=10000):
    """
    BioPython을 사용한 PubMed 검색
    PMID 목록만 반환
    """
    log_info(f"검색어: {term}")
    log_info(f"최대 결과: {retmax}개")
    
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        result = Entrez.read(handle)
        handle.close()
        
        pmids = result['IdList']
        log_info(f"총 PMID 개수: {len(pmids)}")
        return pmids
        
    except Exception as e:
        log_info(f"검색 오류: {e}")
        return []

def load_search_collections():
    """JSON 설정 파일에서 검색 컬렉션 로드"""
    config_path = os.path.join(os.path.dirname(__file__), "../../config/search_collections.json")
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
            return config['search_collections']
    except Exception as e:
        log_info(f"설정 파일 로드 오류: {e}")
        return []

if __name__ == "__main__":
    import sys
    
    # JSON 설정에서 검색 컬렉션 로드
    search_collections = load_search_collections()
    
    if not search_collections:
        log_info("설정 파일에서 검색 컬렉션을 로드할 수 없습니다.")
        sys.exit(1)
    
    # 활성화된 컬렉션만 필터링
    active_collections = [col for col in search_collections if col.get('enabled', True)]
    
    # 명령행 인자로 검색어가 주어지면 해당 검색어만 사용
    if len(sys.argv) > 1:
        # 단일 검색어 모드 (기존 호환성)
        terms = [sys.argv[1]]
        collections = [{"term": sys.argv[1], "collection": "Manual Search"}]
    else:
        # JSON 설정 모드
        collections = active_collections
    
    total_collected = 0
    total_new = 0
    total_saved = 0
    
    for i, collection_config in enumerate(collections, 1):
        term = collection_config['term']
        collection_name = collection_config.get('collection', 'Default Collection')
        retmax = collection_config.get('retmax', 10000)
        
        log_info(f"검색 {i}/{len(collections)} 시작: {collection_name}")
        
        # PMID 수집
        pmids = esearch_biopython(term, retmax)
        
        # 중복 확인 및 필터링
        result = pmid_filter.filter_new_pmids(pmids)
        
        if result and result['new_pmids']:
            # 새로운 PMID들을 데이터베이스에 저장 (검색어와 컬렉션 정보 포함)
            save_result = pmid_filter.save_pmids_to_papers(
                result['new_pmids'], 
                search_term=term, 
                collection_name=collection_name
            )
            saved_count = save_result.get('inserted', 0) if save_result else 0
            
            total_collected += len(pmids)
            total_new += len(result['new_pmids'])
            total_saved += saved_count
            
            log_info(f"검색 {i} 완료 [{collection_name}] - 수집: {len(pmids)}개, 신규: {len(result['new_pmids'])}개, 저장: {saved_count}개")
        else:
            total_collected += len(pmids)
            log_info(f"검색 {i} 완료 [{collection_name}] - 수집: {len(pmids)}개, 신규: 0개 (모두 중복)")
    
    log_info(f"전체 작업 완료 - 총 수집: {total_collected}개, 총 신규: {total_new}개, 총 저장: {total_saved}개")