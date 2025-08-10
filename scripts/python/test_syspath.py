# scripts/python/path_check.py

import sys
from pathlib import Path
from common.logger import log_info, log_debug

log_info("=== Python 모듈 검색 경로 (sys.path) ===")
for i, path in enumerate(sys.path):
    log_info(f"{i}: {path}")

log_info("=== 현재 작업 디렉토리 ===")
log_info(str(Path.cwd()))
