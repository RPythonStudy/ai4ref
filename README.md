# AI4REF - 의학논문작성 자동화 프로젝트

생성형 인공지능, R/Python을 활용한 의료데이터 분석 기반 의학연구 논문작성 전 과정의 생산성 향상 및 자동화 프로젝트입니다.

## 프로젝트 구조
```
src/
├── collector/      # PubMed 논문 수집
├── preprocessor/   # 데이터 전처리/필터링
├── summarizer/     # AI 기반 논문 요약
├── analyzer/       # 통계/메타분석
├── reviewer/       # 검수/수정
├── exporter/       # 결과물 내보내기
├── common/         # 공통 모듈 (로거 등)
└── R/              # R 공용 함수
```

## 핵심 개발 지침

### 환경 관리
- Python: .venv 가상환경 + direnv로 PYTHONPATH 자동 설정
- R: renv 패키지 환경 격리
- 로그: 통합 로깅 시스템 (JSON 형식)
- DB: PostgreSQL 도커 컨테이너 운영

### 코딩 규칙
- Python import: `from common.logger import log_info` (PYTHONPATH=src 설정됨)
- 모든 프로젝트 스크립트: 로깅 시스템 적용 필수
- 교육용 QMD 파일: 로깅 없이 간결하게 작성
- 모듈화 설계: 개별 실행 + 연계 파이프라인 지원
- 논문 수집: DOI→PMID→Title+Author+Year 3단계 중복체크 필수

### 논문 수집 워크플로우
1. load_collections
   - `search_collections.json` 파일에서 collection별 검색어(term) 로드
   - `collections` 테이블에 저장
   - 기본 저장 정책: overwrite (기존 데이터 덮어쓰기 + 변경 추적)
2. pubmed_esearch
   - `collections` term으로 회신받은 pmid를 `collection_pmid`에 저장
3. pmid_filter
   - `collection_pmid`에서 중복되지 않은 pmid를 `unique_pmid`에 저장
   - `unique_pmid`와 `papers` 테이블과 중복되지 않는 pmid 추출하여 `filtered_pmid`로 저장
4. pubmed_efetch
   - `filtered_pmid`로 전체서지정보 회신받음
   - 전체서지정보와 papers 테이블 비교하여 중복을 제외하고 저장
   - 이때, pmid는 이미 중복을 검토했으므로 DOI, Title+Author+Year 중복 체크
5. pdf_downloader
   - papers 테이블에 저장된 정보를 이용하여 오픈 PDF 다운로드
   - 계층적 배치 처리: PMC → Unpaywall → arXiv → Europe PMC 순서
   - 비동기 처리로 대량 논문 고속 다운로드 지원
6. zotero_webapi
   - papers와 collection_pmid 조합으로 paper_collection 테이블 생성
   - 단일 아이템 + 다중 컬렉션 링크 방식으로 Zotero 업로드
   - PDF 첨부파일 자동 업로드 지원

#### Collection 저장 정책
- **Default**: overwrite - 기존 데이터 덮어쓰기 + `updated_at` 갱신
- **Skip**: 중복 시 건너뛰기 (개발/테스트용)
- **Atomic**: 전체 백업 → 삭제 → 재생성 (운영환경 권장)

### PDF 다운로드 전략

#### 대규모 처리 최적화 (수천개 논문 기준)
- **계층적 배치 처리**: 소스별 일괄 처리로 API 효율성 극대화
- **비동기 처리**: asyncio + aiohttp 기반 동시 다운로드
- **적응형 레이트 리미팅**: 실시간 성공률 기반 배치 크기 조정
- **회로 차단기 패턴**: 실패한 소스 자동 우회

#### 다운로드 소스 우선순위
1. **PMC Open Access** (최우선)
   - PMCID 활용한 직접 다운로드
   - 성공률: 60-80%, 처리 대상: 전체 논문의 75%
   - URL: `https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/`

2. **Unpaywall API** (2차)
   - DOI 기반 오픈액세스 논문 탐지
   - 성공률: 50-70%, PMC 실패 논문 대상
   - API: `https://api.unpaywall.org/v2/{doi}`

3. **arXiv 논문** (3차)
   - arXiv DOI 패턴 논문 처리 (`10.48550/arXiv.*`)
   - PDF URL: `https://arxiv.org/pdf/{arxiv_id}.pdf`

4. **Europe PMC/기타** (보완)
   - 추가 소스를 통한 나머지 논문 처리

#### 구현 로드맵
- **Phase 1**: PMC 비동기 배치 구현 (전체 논문의 50-70% 해결)
- **Phase 2**: Unpaywall 통합 (누적 70-85% 해결)
- **Phase 3**: arXiv 및 기타 소스 추가 (90%+ 달성)

#### 성능 목표 (5,000개 논문 기준)
- 처리 시간: 1-2시간 (배치 처리)
- PDF 확보율: 80-90%
- 메모리 사용량: 300-600MB
- API 호출 최적화: 순차 대비 70% 단축

### Zotero 통합 전략

#### 다중 컬렉션 지원 방식
**업계 표준**: 단일 아이템 + 다중 컬렉션 링크 (시나리오 1)

**특징**:
- 하나의 논문 = 하나의 Zotero 아이템
- 여러 컬렉션에 동일 아이템 링크
- `addToCollection()` API 활용
- PDF 첨부파일 단일 저장

**구현 단계**:
1. **Phase 1**: 첫 컬렉션에 아이템 생성 (`create_items()`)
2. **Phase 2**: PDF 첨부파일 업로드 (`attachment_simple()`)
3. **Phase 3**: 추가 컬렉션에 링크 (`addto_collection()`)
4. **Phase 4**: `paper_collection` 테이블 상태 업데이트

**데이터베이스 연동**:
- `paper_collection.zotero_item_key`: 단일 아이템 추적
- `paper_collection.synced_at`: 동기화 완료 시점
- 다중 컬렉션 관계: 동일 `zotero_item_key`로 여러 레코드

**장점**:
- 데이터 정규화 준수 (메타데이터 중복 방지)
- 스토리지 효율성 (PDF 단일 저장)
- 일관성 보장 (수정 시 모든 컬렉션 자동 반영)
- Zotero 네이티브 방식

#### pyzotero 라이브러리 주의사항 (중요!)

**메서드명 규칙**: pyzotero는 camelCase가 아닌 snake_case 사용
- ❌ `addToCollection()` → ✅ `addto_collection()`
- ❌ `createItems()` → ✅ `create_items()`
- ❌ `attachmentSimple()` → ✅ `attachment_simple()`

**addto_collection 사용법**:
```python
# 잘못된 방법
result = zot.addto_collection(collection_key, item_key)  # 오류!

# 올바른 방법
item = zot.item(item_key)  # 전체 아이템 객체 조회
if isinstance(item, list):  # item() 메서드는 리스트 반환
    item = item[0]
result = zot.addto_collection(collection_key, item)  # 전체 객체 전달
```

**주요 이슈들**:
1. **메서드명 오류**: `'Zotero' object has no attribute 'addToCollection'`
   - 해결: `addToCollection` → `addto_collection`

2. **payload 타입 오류**: `string indices must be integers, not 'str'`
   - 해결: 아이템 키 문자열 대신 전체 아이템 객체 전달

3. **리스트 처리 오류**: `list indices must be integers or slices, not str`
   - 해결: `zot.item()` 결과가 리스트인 경우 첫 번째 요소 사용

**검증된 구현 패턴**:
```python
def add_to_collection(self, item_key, collection_key):
    try:
        # 1. 전체 아이템 객체 조회
        item = self.zot.item(item_key)
        
        # 2. 리스트인 경우 첫 번째 요소 추출
        if isinstance(item, list):
            if not item:
                return False
            item = item[0]
        
        # 3. 전체 객체를 payload로 전달
        result = self.zot.addto_collection(collection_key, item)
        return bool(result)
    except Exception as e:
        log_debug(f"컬렉션 추가 오류: {e}")
        return False
```



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

# 특정 단계 건너뛰기
./scripts/run_full_pipeline.sh --skip-init --skip-collection

# 상태만 확인
./scripts/run_full_pipeline.sh --status-only

# 도움말
./scripts/run_full_pipeline.sh --help
```

#### Python 스크립트
```bash
# 전체 워크플로우 실행
.venv/bin/python scripts/run_full_pipeline.py

# 특정 범위만 실행
.venv/bin/python scripts/run_full_pipeline.py --from-step 3 --to-step 6

# 상태만 확인
.venv/bin/python scripts/run_full_pipeline.py --status-only
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

## 상세 문서
- [프로젝트 개요](wiki/project-overview.md)
- [환경 설정 가이드](wiki/environment-setup.md)
- [로깅 시스템 가이드](wiki/logging-guide.md)
- [논문 중복관리 가이드](wiki/duplication-management.md)
- [PDF 다운로드 가이드](wiki/pdf-download-guide.md)

## Copilot 지침
- 이모지/굵은글씨 사용 금지
- Python 가상환경 고려한 스크립트 제안
- 극단적으로 간결하면서도 직관적이고 디버깅 용이한 구조로 스크립트를 제안할 것
- 개발과정에 오류가 발생했던 것은 가급적 디버깅로그로 남기도록 하고, 데이터의 갯수가 변화하는 경우는 인포로그로 출력하도록
- PYTHONPATH=src 설정 인지하여 import 경로 제안

### 터미널 실행 제약사항
- `python -c "..."` 명령 사용 금지 (멀티라인 코드 파싱 오류, 따옴표 충돌, DB 연결 시간초과 등)
- 대신 다음 방식 사용:
  - 기존 스크립트 활용: `src/database/utils/analyze_database.py`
  - 임시 파일 생성 후 실행: `tests/test_임시작업명.py` 형식으로 생성
  - 간단한 SQL은 psql 직접 사용

#### 임시 파일 생성 예시
```bash
# PMC ID 현황 확인용 임시 스크립트
cat > tests/test_pmc_status.py << 'EOF'
from common.database import get_db_connection

conn = get_db_connection()
cur = conn.cursor()
cur.execute("SELECT COUNT(*) FROM papers WHERE pmcid IS NOT NULL")
print(f"PMC ID 있는 논문: {cur.fetchone()[0]}개")
cur.close()
conn.close()
EOF

# 실행
cd /home/ben/projects/ai4ref && .venv/bin/python tests/test_pmc_status.py
```

## 개발 이력 및 해결된 주요 이슈

### Zotero API 통합 (2025-08-09)
**해결된 문제들**:
- ❌ `'Zotero' object has no attribute 'addToCollection'`
- ❌ `string indices must be integers, not 'str'`  
- ❌ `list indices must be integers or slices, not str`

**핵심 수정사항**:
1. **메서드명 정정**: `addToCollection` → `addto_collection`
2. **payload 구조**: 아이템 키 → 전체 아이템 객체
3. **리스트 처리**: `zot.item()` 결과의 첫 번째 요소 추출
4. **디버깅 강화**: 상세 로그로 문제 진단 개선

**영향받은 파일들**:
- `src/collector/zotero_webapi_v2.py`: 메인 Zotero 통합 로직
- `scripts/run_full_pipeline.sh`: Bash 자동화 스크립트
- `scripts/run_full_pipeline.py`: Python 자동화 스크립트
- `README.md`: pyzotero 사용법 가이드 추가

**성과**: 전체 워크플로우 자동화 완성, 산업 표준 Zotero 통합 달성