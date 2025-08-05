#!/usr/bin/env python3

import sys
import os
sys.path.append('/home/ben/projects/ai4ref/src')

from collector.db import save_pmids_to_db, get_pmids_count, get_recent_pmids

# 테스트 데이터
test_pmids = ["12345678", "87654321", "11111111"]

print("=== 데이터베이스 저장 테스트 ===")
result = save_pmids_to_db(test_pmids)
if result:
    print(f"저장 결과: {result}")
else:
    print("저장 실패")

print("\n=== 데이터베이스 조회 테스트 ===")
count = get_pmids_count()
if count is not None:
    print(f"총 PMID 수: {count}")

recent = get_recent_pmids(5)
if recent:
    print("최근 저장된 PMID:")
    for pmid, timestamp in recent:
        print(f"  {pmid} - {timestamp}")
