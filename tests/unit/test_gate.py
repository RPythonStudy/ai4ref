# tests/unit/test_gate.py
"""LLM 게이트(judge) 단위테스트 — 배관 결정론 가드.

LLM 비결정성(실 claude_cli 판정 정확도)은 acceptance(gate_baseline.yml) 책임.
여기서는 백엔드를 모의(subprocess.run patch)해 다음을 결정론 검증한다:
  · 출력 스키마 {is_relevant, is_landmark, reason} (FR-006)
  · 2단계 순서 — ① is_relevant=false 면 ② is_landmark 평가 없이 false (FR-001·002·SC-005)
  · fail-soft — 호출/파싱 실패 시 보수적 보류 (FR-009)
  · _extract_json — 래퍼({"result":"..."})·둘러싼 텍스트에서 견고 추출 (FR-010)
  · 백엔드 디스패치 — stub/claude_cli/api(미구현→stub 강등) (FR-007)
  · PROMPT 컨텍스트 — 초록 3000자 절단·question_type→설계라벨 (FR-005·011)
"""
import json
from unittest.mock import patch, MagicMock

from alert.llm_gate import judge, _extract_json


KQ_CTX = {
    "kq": "복부수술 수액전략",
    "comparison": "제한적 vs 자유로운",
    "outcome": "합병증",
    "guideline": "ERAS (2018)",
    "question_type": "intervention",
    "strictness": "loose",
}
PAPER = {"title": "RELIEF trial", "abstract": "restrictive vs liberal fluid in major abdominal surgery"}


def _proc(stdout: str, returncode: int = 0):
    m = MagicMock()
    m.stdout = stdout
    m.returncode = returncode
    return m


# ── 출력 스키마 / 백엔드 디스패치 ────────────────────────────────────────────

def test_stub_schema_keys():
    """stub 백엔드 출력은 {is_relevant, is_landmark, reason} 스키마 (FR-006)."""
    v = judge(PAPER, KQ_CTX, backend="stub")
    assert set(v) == {"is_relevant", "is_landmark", "reason"}
    assert v["is_relevant"] is True and v["is_landmark"] is True


def test_api_backend_falls_back_to_stub():
    """미구현 api 백엔드 선택 시 안전하게 stub 으로 강등 (FR-007)."""
    v = judge(PAPER, KQ_CTX, backend="api")
    assert set(v) == {"is_relevant", "is_landmark", "reason"}
    assert v["is_relevant"] is True  # stub 결과


def test_default_backend_is_stub():
    """기본 백엔드 = stub (키 불필요 개발용)."""
    assert judge(PAPER, KQ_CTX) == judge(PAPER, KQ_CTX, backend="stub")


# ── claude_cli 정상 판정 ────────────────────────────────────────────────────

def test_claude_cli_relevant_landmark():
    """관련·랜드마크 JSON 을 그대로 매핑 (FR-001·003·006)."""
    out = json.dumps({"is_relevant": True, "is_landmark": True, "reason": "practice-changing"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v == {"is_relevant": True, "is_landmark": True, "reason": "practice-changing"}


def test_claude_cli_reject_false_positive():
    """① is_relevant=false 면 ② is_landmark=false (2단계 순서, FR-002·SC-005)."""
    out = json.dumps({"is_relevant": False, "is_landmark": False, "reason": "대상 불일치"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v["is_relevant"] is False and v["is_landmark"] is False


def test_claude_cli_bool_coercion():
    """truthy 값도 bool 로 강제 (스키마 안정)."""
    out = json.dumps({"is_relevant": 1, "is_landmark": 0, "reason": "x"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v["is_relevant"] is True and v["is_landmark"] is False


# ── fail-soft (보수적 보류) ─────────────────────────────────────────────────

def test_claude_cli_failsoft_on_exception():
    """subprocess 예외 → 보수적 보류 (FR-009)."""
    with patch("alert.llm_gate.subprocess.run", side_effect=TimeoutError("boom")):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v["is_relevant"] is False and v["is_landmark"] is False
    assert "판정 실패" in v["reason"]


def test_claude_cli_failsoft_on_is_error_wrapper():
    """claude CLI is_error 래퍼 → 보류 (FR-009)."""
    out = json.dumps({"is_error": True, "result": "rate limit"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v["is_relevant"] is False and v["is_landmark"] is False


def test_claude_cli_failsoft_on_bad_json():
    """비-JSON 출력 → 파싱 실패 → 보류 (FR-009)."""
    with patch("alert.llm_gate.subprocess.run", return_value=_proc("not json at all")):
        v = judge(PAPER, KQ_CTX, backend="claude_cli")
    assert v["is_relevant"] is False and v["is_landmark"] is False


# ── _extract_json 견고 파싱 (FR-010) ────────────────────────────────────────

def test_extract_json_plain():
    assert _extract_json('{"is_landmark": true, "is_relevant": true}')["is_landmark"] is True


def test_extract_json_result_wrapper_str():
    """{"result": "...json 문자열..."} 래퍼에서 내부 판정 추출."""
    inner = '{"is_relevant": true, "is_landmark": false, "reason": "r"}'
    wrapped = json.dumps({"result": inner})
    got = _extract_json(wrapped)
    assert got["is_relevant"] is True and got["is_landmark"] is False


def test_extract_json_surrounding_text():
    """앞뒤 부가 텍스트가 있어도 중괄호 블록 추출."""
    txt = '설명...\n{"is_landmark": true, "is_relevant": true, "reason": "ok"}\n끝'
    assert _extract_json(txt)["is_landmark"] is True


def test_extract_json_garbage_returns_empty():
    assert _extract_json("그냥 텍스트") == {}


# ── PROMPT 컨텍스트 구성 (FR-005·011) ───────────────────────────────────────

def test_prompt_truncates_long_abstract():
    """초록은 3000자로 절단되어 프롬프트에 들어간다 (FR-005)."""
    long_paper = {"title": "T", "abstract": "x" * 5000}
    out = json.dumps({"is_relevant": True, "is_landmark": False, "reason": "ok"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)) as mock_run:
        judge(long_paper, KQ_CTX, backend="claude_cli")
    prompt = mock_run.call_args[0][0][2]   # ["claude","-p",PROMPT,...]
    assert "x" * 3000 in prompt and "x" * 3001 not in prompt


def test_prompt_design_label_fallback_when_qtype_missing():
    """question_type 미지정 → 일반 설계 라벨로 대체 (FR-011)."""
    ctx = dict(KQ_CTX); ctx.pop("question_type")
    out = json.dumps({"is_relevant": True, "is_landmark": False, "reason": "ok"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)) as mock_run:
        judge(PAPER, ctx, backend="claude_cli")
    prompt = mock_run.call_args[0][0][2]
    assert "해당 유형의 적정 연구설계" in prompt


def test_prompt_design_label_for_intervention():
    """intervention → RCT 설계 라벨 주입 (strict 컨텍스트, FR-004)."""
    out = json.dumps({"is_relevant": True, "is_landmark": False, "reason": "ok"})
    with patch("alert.llm_gate.subprocess.run", return_value=_proc(out)) as mock_run:
        judge(PAPER, KQ_CTX, backend="claude_cli")
    prompt = mock_run.call_args[0][0][2]
    assert "무작위 대조시험(RCT)" in prompt
