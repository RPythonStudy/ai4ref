# AI4REF - 의학논문작성 자동화 프로젝트

## 목적
생성형 인공지능, R/Python 및 PostgreSQL 등 오픈소스 도구를 활용하여, 의학연구 논문작성 전 과정의 생산성 향상 및 자동화

## 이 프로젝트 README.md의 특이사항
일반적인 README.md는 프로젝트의 사용자의 이해를 돕고 설치를 안내하는 목적이라면, 이 프로젝트에서는 이 README.md 파일이 Copilot에 대한 지침을 겸하도록 작성하였습니다 (ai4ref/.github/copilot-instructions.md에 symolic link로 동기화). 따라서 일부 문구는 사용자분들의 이해보다는 생성형 인공지능이 이해하기에 적합하도록 작성되었음을 양해해주시길 바랍니다.
그리고 이 프로젝트의 사용자들이라 표현하였지만 이 프로젝트의 공동개발자를 우선 염두에 두고 작성되었음도 양해해 주시길 바랍니다.


## 폴더 구조 설명
   
ai4ref/   
├── src/                               # 핵심 파이썬/R 코드   
│   ├── collector/                     # 논문 수집 모듈 (PubMed, PMC, arXiv 등)   
│   │   ├── load_collections.py        # 컬렉션 검색어 로드   
│   │   ├── pubmed_esearch.py          # PubMed 검색   
│   │   ├── pmid_filter.py             # PMID 중복 필터링   
│   │   ├── pubmed_efetch.py           # 논문 메타데이터 수집   
│   │   ├── pmc_downloader.py          # PDF 다운로드   
│   │   ├── zotero_webapi.py           # Zotero 업로드   
│   │   └── zotero_webapi/             # Zotero 연동 세부 모듈   
│   ├── database/                      # DB 스키마 및 관리   
│   │   ├── schema/                    # 테이블별 SQL 스키마   
│   │   │   ├── collections.sql        # 컬렉션 테이블 스키마   
│   │   │   ├── collection_pmid.sql    # 컬렉션별 PMID 테이블   
│   │   │   ├── unique_pmid.sql        # 중복제거된 PMID 테이블   
│   │   │   ├── filtered_pmid.sql      # 최종 필터링된 PMID 테이블   
│   │   │   ├── papers.sql             # 논문 메타데이터 테이블   
│   │   │   └── paper_collection.sql   # 논문-컬렉션 관계 테이블   
│   │   └── utils/                     # DB 초기화/분석 스크립트   
│   │       ├── initialize_database.py # DB 초기화   
│   │       └── analyze_database.py    # DB 상태/통계 분석   
│   ├── preprocessor/                  # 데이터 전처리 및 필터링   
│   ├── summarizer/                    # AI 기반 논문 요약 및 텍스트 생성   
│   ├── analyzer/                      # 통계 분석, 메타분석, 데이터 해석   
│   ├── reviewer/                      # 논문 검수, 자동 수정, 품질관리   
│   ├── exporter/                      # 결과물 내보내기 (표, 그래프, 보고서 등)   
├── data/   
│   └── pdf/                           # 수집된 논문 PDF 저장 폴더   
├── config/                
│   └── search_collections.json        # 논문 수집용 컬렉션별 검색어/필터 설정   
├── logs/                              # 실행 및 디버깅 로그 (ai4medpaper.log 등)   
├── scripts/               
│   ├── run_full_pipeline.sh           # 전체 워크플로우 자동 실행   
│   ├── daily_collection.sh            # 일일 논문 수집 자동화   
│   ├── python/                        # 파이썬 테스트/유틸 스크립트   
│   │   ├── test_import_logger.py      # 로거 임포트 테스트   
│   │   └── test_syspath.py            # 파이썬 경로 테스트   
│   └── setup/                            
│       └── wsl2_install_postgres.sh   # WSL2용 Postgres 설치   
└── wiki/                              # 상세 개발 문서   
   
## Wiki 문서링크
- 이 프로젝트에서는 README.md 보다 자세한 설명이 필요한 사항들은 Wiki 문서로 작성하였습니다. 링크를 이용하시면 됩니다.

- [프로젝트 개요](wiki/project-overview.md)
- [개발환경구축](wiki/development_enviroment_setup.md)
- [프로젝트 설치](wiki/installation.md)
- [로깅 시스템 가이드](wiki/logging-guide.md)


## 핵심 개발 지침

### 환경 관리
- Python: .venv 가상환경 + direnv로 PYTHONPATH 자동 설정
- R: renv 패키지 환경 격리
- 로그: 통합 로깅 시스템 (JSON 형식)
- DB: PostgreSQL 도커 컨테이너 운영

### 코딩 규칙
- Python import: `from common.logger import log_info` (PYTHONPATH=src 설정됨)
- 모든 프로젝트 스크립트: 로깅 시스템 적용 필수
- 모듈화 설계: 개별 실행 + 연계 파이프라인 지원
- 논문 수집: DOI→PMID→Title+Author+Year 3단계 중복체크




### 터미널 명령어
가상환경 사용을 고려한 명령어 형식:
```bash
# Python 패키지 설치
/home/ben/projects/ai4ref/.venv/bin/pip install requests

# Python 스크립트 실행  
/home/ben/projects/ai4ref/.venv/bin/python src/collector/pubmed_esearch.py
```

### 전체 파이프라인 자동화

#### Bash 스크립트 (권장)
```bash
# 전체 워크플로우 실행
./scripts/run_full_pipeline.sh
```

**워크플로우 단계**:
1. 데이터베이스 초기화 (`initialize_database.py`)
2. 컬렉션 로드 (`load_collections.py`)
3. PubMed 검색 (`pubmed_esearch.py`)
4. PMID 필터링 (`pmid_filter.py`)
5. 논문 메타데이터 수집 (`pubmed_efetch.py`)
6. Paper-Collection 관계 생성 (`create_paper_collection.py`)
7. PDF 다운로드 (`pmc_downloader.py`)
8. Zotero 업로드 (`zotero_webapi_v2.py`)



## 개발/운영 지침
- 각 단계는 src/{collector, preprocessor, summarizer, analyzer, reviewer, exporter} 등 모듈화된 폴더로 분리
- 각 모듈은 독립 실행(CLI, REST API 등) 및 연계 파이프라인 실행 모두 지원
- 입력/출력은 파일 또는 DB를 표준 인터페이스로 사용
- 단계별 파라미터는 config 파일(yaml/json 등) 기반으로 관리
- Docker 등 컨테이너 배포도 지원, 필요시 개별 서비스화 가능
- hybrid(자동+수작업) 검수/수정 구조 적극 지원
- 로그/예외처리 일관성 확보
- 결과물의 reproducibility(재현성)와 추적성(traceability) 확보
- 확장성(인터페이스/API/웹UI 등) 고려


## Copilot 지침
- 이모지/굵은글씨 사용 금지
- Python 가상환경 고려한 스크립트 제안
- 극단적으로 간결하면서도 직관적이고 디버깅 용이한 구조로 스크립트를 제안할 것
- 개발과정에 오류가 발생했던 것은 가급적 디버깅로그로 남기도록 하고, 데이터의 갯수가 변화하는 경우는 인포로그로 출력하도록
- PYTHONPATH=src 설정 인지하여 import 경로 제안
- `python -c "..."` 명령 사용 금지 (멀티라인 코드 파싱 오류, 따옴표 충돌, DB 연결 시간초과 등)
- 대신 다음 방식 사용:
  - 기존 스크립트 활용: `src/database/utils/analyze_database.py`
  - 임시 파일 생성 후 실행: `tests/test_임시작업명.py` 형식으로 생성
  - 간단한 SQL은 psql 직접 사용

