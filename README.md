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

### 논문 수집 워크플로우 (필수 준수)
1. **esearch**: 검색어로 PMID 수집
2. **중복 체크**: DB에서 기존 PMID와 비교하여 중복 제거
3. **efetch**: 중복되지 않는 PMID만 상세 정보 수집
4. **최종 저장**: 수집된 논문 정보를 DB에 저장

**중요**: efetch 전 반드시 DB 중복 체크 수행하여 API 호출 최소화

### 터미널 명령어
가상환경 사용을 고려한 명령어 형식:
```bash
# Python 패키지 설치
/home/ben/projects/ai4ref/.venv/bin/pip install requests

# Python 스크립트 실행  
# Python 스크립트 실행  
/home/ben/projects/ai4ref/.venv/bin/python src/collector/pubmed_esearch.py
```

## 상세 문서
- [프로젝트 개요](wiki/project-overview.md)
- [환경 설정 가이드](wiki/environment-setup.md)
- [로깅 시스템 가이드](wiki/logging-guide.md)
- [논문 중복관리 가이드](wiki/duplication-management.md)

## Copilot 지침
- 이모지/굵은글씨 사용 금지
- Python 가상환경 고려한 스크립트 제안
- 간결하고 디버깅 용이한 구조
- PYTHONPATH=src 설정 인지하여 import 경로 제안