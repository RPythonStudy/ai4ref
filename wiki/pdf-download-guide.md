# PDF 다운로드 가이드

## 개요

AI4REF 프로젝트의 PDF 다운로드 시스템은 수천개의 의학논문을 효율적으로 수집하기 위해 설계된 계층적 배치 처리 시스템입니다.

## 핵심 전략

### 1. 계층적 배치 처리 (Hierarchical Batch Processing)

소스별로 일괄 처리하여 API 효율성을 극대화합니다:

```
Phase 1: PMC 배치 처리 → 75% 논문 대상
Phase 2: Unpaywall 배치 처리 → 나머지 논문 대상  
Phase 3: 특화 소스 처리 → 최종 보완
```

### 2. 비동기 처리 아키텍처

```python
# 핵심 구현 패턴
import asyncio
import aiohttp

class PDFDownloadPipeline:
    async def process_batch(self, papers, batch_size=50):
        semaphore = asyncio.Semaphore(batch_size)
        async with aiohttp.ClientSession() as session:
            tasks = [self.download_single(session, paper, semaphore) 
                    for paper in papers]
            return await asyncio.gather(*tasks, return_exceptions=True)
```

## 다운로드 소스별 전략

### 1. PMC Open Access (최우선)

**특징:**
- 성공률: 60-80%
- 처리 대상: PMCID 보유 논문 (전체의 75%)
- 안정적인 API, 관대한 레이트 리미트

**URL 패턴:**
```python
pmc_urls = [
    f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/",
    f"https://europepmc.org/articles/{pmcid}?pdf=render"
]
```

**구현 우선순위: 1위** - 즉시 구현하여 대부분의 PDF 확보

### 2. Unpaywall API (2차)

**특징:**
- 성공률: 50-70%
- 처리 대상: DOI 보유 논문 (PMC 실패분)
- 오픈액세스 여부 판단 및 무료 PDF URL 제공

**API 패턴:**
```python
unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email=your@email.com"
# 응답에서 oa_locations[].url_for_pdf 추출
```

### 3. arXiv 논문 (3차)

**특징:**
- arXiv DOI 패턴: `10.48550/arXiv.*`
- 거의 100% 성공률 (해당 논문 대상)
- 직접 PDF URL 생성 가능

**URL 변환:**
```python
# DOI: 10.48550/arXiv.1602.04938
# PDF: https://arxiv.org/pdf/1602.04938.pdf
arxiv_id = doi.replace("10.48550/arXiv.", "")
pdf_url = f"https://arxiv.org/pdf/{arxiv_id}.pdf"
```

### 4. Europe PMC & 기타 (보완)

**특징:**
- 나머지 논문들 대상
- 출판사별 특화 처리
- 보완적 역할

## 성능 최적화

### 배치 크기 최적화

```python
# 소스별 권장 배치 크기
PMC_BATCH_SIZE = 50      # 안정적인 서비스
UNPAYWALL_BATCH_SIZE = 100  # 높은 처리량
ARXIV_BATCH_SIZE = 20    # 보수적 접근
```

### 레이트 리미팅

```python
# 소스별 딜레이 설정
rate_limits = {
    'pmc': 0.33,      # 초당 3회
    'unpaywall': 0.01, # 초당 100회 (관대함)
    'arxiv': 3.0,     # 3초당 1회 (엄격함)
    'europepmc': 1.0  # 초당 1회
}
```

### 회로 차단기 패턴

```python
class CircuitBreaker:
    def __init__(self, failure_threshold=0.1):
        self.failure_threshold = failure_threshold
        self.is_open = False
        
    def check_failure_rate(self, success_count, total_count):
        if total_count > 10 and success_count / total_count < self.failure_threshold:
            self.is_open = True
            log_error(f"회로 차단기 작동: 성공률 {success_count/total_count:.1%}")
```

## 파일 저장 구조

### 디렉토리 구조
```
downloads/
├── pmc/          # PMC 소스 PDF들
├── unpaywall/    # Unpaywall 소스 PDF들
├── arxiv/        # arXiv 논문들
└── others/       # 기타 소스들
```

### 파일명 규칙
```python
# {pmid}_{source}.pdf 형식
filename = f"{paper.pmid}_{source}.pdf"
# 예: 40615394_pmc.pdf
```

### 메타데이터 업데이트
```python
# papers 테이블 업데이트
update_query = """
    UPDATE papers 
    SET pdf_url = %s, pdf_source = %s, updated_at = NOW()
    WHERE pmid = %s
"""
```

## 모니터링 및 로깅

### 진행률 추적
```python
# 단계별 성공률 로깅
log_info(f"PMC 처리: {pmc_success}/{total_pmc} ({pmc_rate:.1%})")
log_info(f"Unpaywall 처리: {unpaywall_success}/{total_unpaywall} ({unpaywall_rate:.1%})")
log_info(f"전체 성공률: {total_success}/{total_papers} ({total_rate:.1%})")
```

### 에러 처리
```python
# 소스별 에러 분류
error_stats = {
    'rate_limit': 0,
    'not_found': 0, 
    'server_error': 0,
    'network_error': 0
}
```

## 구현 로드맵

### Phase 1: PMC 비동기 배치 (1-2주)
- 목표: 전체 논문의 50-70% PDF 확보
- 우선순위: 높음 (즉시 실용성)
- 구현: `src/collector/pdf_downloader.py`

### Phase 2: Unpaywall 통합 (1주)
- 목표: 누적 70-85% PDF 확보
- PMC 실패분 + DOI 전용 논문 처리

### Phase 3: 특화 소스 추가 (선택적)
- arXiv, Europe PMC, 출판사별 처리
- 목표: 90%+ PDF 확보

## 사용법

### 기본 실행
```bash
# PMC 배치 다운로드
cd /home/ben/projects/ai4ref
.venv/bin/python src/collector/pdf_downloader.py --source pmc

# 전체 파이프라인 실행
.venv/bin/python src/collector/pdf_downloader.py --all
```

### 진행상황 확인
```bash
# PDF 다운로드 현황 확인
.venv/bin/python tests/test_pdf_download_status.py
```

## 기대 성능 (5,000개 논문 기준)

- **처리 시간**: 1-2시간 (배치 처리)
- **PDF 확보율**: 80-90%
- **메모리 사용량**: 300-600MB
- **API 호출 최적화**: 순차 처리 대비 70% 단축
- **네트워크 효율성**: Keep-alive 연결 재사용으로 지연시간 최소화

## 문제 해결

### 일반적인 문제들

1. **레이트 리미트 초과**
   ```python
   # 지수 백오프로 재시도
   await asyncio.sleep(2 ** attempt)
   ```

2. **네트워크 타임아웃**
   ```python
   # 타임아웃 설정 조정
   timeout = aiohttp.ClientTimeout(total=30)
   ```

3. **메모리 부족**
   ```python
   # 배치 크기 감소
   batch_size = max(10, current_batch_size // 2)
   ```

4. **디스크 공간 부족**
   ```bash
   # 다운로드 디렉토리 정리
   find downloads/ -name "*.pdf" -size +50M -delete
   ```
