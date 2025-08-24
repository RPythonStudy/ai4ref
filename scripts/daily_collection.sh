#!/bin/bash
# AI4REF 매일 자동 실행 스크립트

# 프로젝트 디렉토리로 이동
cd /home/ben/projects/ai4ref

# 가상환경에서 PubMed 검색 실행
/home/ben/projects/ai4ref/.venv/bin/python src/collector/pubmed_esearch.py

# 로그 확인 (선택적)
echo "=== 최근 로그 ==="
tail -10 logs/ai4ref.log
