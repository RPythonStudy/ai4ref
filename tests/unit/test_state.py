# tests/unit/test_state.py
"""SeenStore 단위테스트.

FR-001~FR-007 요구사항 및 fail-soft 예외 복구(SC-001)를 검증합니다.
"""
import json
from pathlib import Path
import pytest

from common.state import SeenStore


def test_seen_store_non_existent_file(tmp_path):
    """지정 경로에 파일이 존재하지 않는 경우 빈 상태로 시작한다 (FR-003)."""
    state_file = tmp_path / "non_existent.jsonl"
    store = SeenStore(str(state_file))
    
    assert len(store) == 0
    assert not store.is_seen("test_key")


def test_seen_store_load_valid_file(tmp_path):
    """정상적인 JSONL 파일로부터 보낸 이력 키 목록을 정확히 로드한다 (FR-002)."""
    state_file = tmp_path / "valid.jsonl"
    records = [
        {"key": "kq1:pmid1", "ts": "2026-06-22T00:00:00Z"},
        {"key": "kq1:pmid2", "ts": "2026-06-22T01:00:00Z"}
    ]
    state_file.write_text("\n".join(json.dumps(r) for r in records) + "\n", encoding="utf-8")
    
    store = SeenStore(str(state_file))
    
    assert len(store) == 2
    assert store.is_seen("kq1:pmid1")
    assert store.is_seen("kq1:pmid2")
    assert not store.is_seen("kq1:pmid3")


def test_seen_store_fail_soft_corrupted_file(tmp_path):
    """손상된 행(구문 에러, 빈 행 등)이 섞여 있어도 예외 없이 무시하고 정상 행만 로드한다 (FR-004, SC-001)."""
    state_file = tmp_path / "corrupted.jsonl"
    
    # 1. 정상행, 2. 빈행, 3. 깨진 JSON, 4. key가 없는 JSON, 5. 정상행
    lines = [
        '{"key": "kq1:pmid1", "ts": "2026-06-22"}',
        '',
        '{"key": "broken_json", ',
        '{"ts": "no_key"}',
        '{"key": "kq1:pmid2", "ts": "2026-06-22"}'
    ]
    state_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    
    store = SeenStore(str(state_file))
    
    # 정상적으로 구문 분석 및 key 추출에 성공한 2건만 메모리에 적재되어야 함
    assert len(store) == 2
    assert store.is_seen("kq1:pmid1")
    assert store.is_seen("kq1:pmid2")
    assert not store.is_seen("broken_json")


def test_seen_store_mark_appends_and_updates_memory(tmp_path):
    """mark 호출 시 디스크 파일 끝에 추가 기록되며 메모리 집합에도 즉시 반영된다 (FR-005)."""
    state_file = tmp_path / "test_mark.jsonl"
    store = SeenStore(str(state_file))
    
    # 처음에는 기록이 없음
    assert not store.is_seen("kq:123")
    
    # mark 호출
    store.mark("kq:123", {"title": "Test Title", "sinks": ["telegram"]})
    
    # 1. 메모리 집합에 반영되었는지 검증
    assert len(store) == 1
    assert store.is_seen("kq:123")
    
    # 2. 디스크 파일에 한 줄의 JSON 객체로 정상 기록되었는지 검증
    content = state_file.read_text(encoding="utf-8").strip()
    data = json.loads(content)
    assert data["key"] == "kq:123"
    assert "ts" in data
    assert data["title"] == "Test Title"
    assert data["sinks"] == ["telegram"]


def test_seen_store_creates_parent_directories(tmp_path):
    """상위 디렉토리가 미존재 상태일 때 mark 수행 시 폴더를 자동 생성한다 (FR-003)."""
    # 존재하지 않는 하위 경로 설계
    nested_dir = tmp_path / "new_folder" / "state"
    state_file = nested_dir / "alerted.jsonl"
    
    assert not nested_dir.exists()
    
    store = SeenStore(str(state_file))
    store.mark("kq:abc")
    
    # 상위 디렉토리가 성공적으로 생성되었는지 검증
    assert nested_dir.exists()
    assert state_file.exists()
