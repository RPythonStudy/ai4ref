# src/alert/llm_gate.py
"""랜드마크 판정 — 교체형 LLM 백엔드 (sink 처럼 갈아끼움).

backend:
  stub       = 개발용(가짜 판정, 키 불필요)
  claude_cli = 로컬 Claude Code · OAuth 구독 · 추가비용 0   ← 운영 권장
  api        = Anthropic API · 클라우드 · 토큰 과금(예정)

judge() -> {"relevant": bool, "is_landmark": bool, "reason": str}
"""
import json
import subprocess

from common.features import log

PROMPT = """너는 진료지침 감시 보조다. 아래 논문이 현행 지침과 *다른 결론/프로토콜*을 제시하는 '랜드마크'인지 판정하라.
KQ: {kq}
현행 지침: {guideline}
대조(C, 참고): {comparison}
주요 결과(O, 이 결과를 보고하는지 확인): {outcome}
제목: {title}
초록: {abstract}

판정 원칙(중요 · recall 편향): 감시에서 가장 비싼 오류는 '놓침'이다.
- 지침을 바꿀 *가능성*이 보이면 불확실해도 is_landmark=true 로 포함하라(경계선 = 포함).
- 명백히 무관하거나(relevant=false) 지침과 같은 결론을 재확인할 뿐이면 is_landmark=false.

JSON 으로만 답하라:
{{"relevant": true/false, "is_landmark": true/false, "reason": "한 줄 근거(한국어)"}}"""


def judge(paper: dict, kq: dict, backend: str = "stub") -> dict:
    if backend == "claude_cli":
        return _claude_cli(paper, kq)
    if backend == "api":
        log.warning("[llm_gate] backend=api 미구현 → stub 사용")
    return _stub(paper, kq)


def _stub(paper, kq) -> dict:
    # 골격 검증용: 통과시키고 데모 사유 부여. 실판정은 claude_cli 연결 후.
    return {
        "relevant": True,
        "is_landmark": True,
        "reason": "[stub] 실제 판정은 llm_backend=claude_cli 연결 후",
    }


def _claude_cli(paper, kq, model: str = "haiku") -> dict:
    """로컬 Claude Code(OAuth)에 위임. 추가 API 비용 없음.

    분류는 Haiku 로 충분(빠르고 rate 절약). 출력 래퍼:
    {"type":"result", "is_error":bool, "result":"<답변 문자열=우리 JSON>"}.
    """
    prompt = PROMPT.format(
        kq=kq.get("kq", ""),
        guideline=kq.get("guideline", ""),
        comparison=kq.get("comparison", "") or "(개방)",
        outcome=kq.get("outcome", ""),
        title=paper.get("title", ""),
        abstract=(paper.get("abstract", "") or "")[:3000],
    )
    try:
        proc = subprocess.run(
            ["claude", "-p", prompt, "--output-format", "json", "--model", model],
            capture_output=True, text=True, timeout=120,
        )
        raw = (proc.stdout or "").strip()
        wrapper = {}
        try:
            wrapper = json.loads(raw)
        except Exception:
            pass
        if isinstance(wrapper, dict) and wrapper.get("is_error"):
            raise RuntimeError(str(wrapper.get("result") or "claude is_error"))
        data = _extract_json(raw)
        if not data:
            raise RuntimeError(f"판정 JSON 파싱 실패: {raw[:160]}")
        return {
            "relevant": bool(data.get("relevant")),
            "is_landmark": bool(data.get("is_landmark")),
            "reason": str(data.get("reason", "")),
        }
    except Exception as e:
        log.error(f"[llm_gate:claude_cli] 실패 → 보수적 보류: {e}")
        return {"relevant": False, "is_landmark": False, "reason": f"판정 실패: {e}"}


def _extract_json(txt: str) -> dict:
    """claude CLI 출력에서 판정 JSON 추출. 래퍼 {"result": "...json..."} 도 대응."""
    txt = (txt or "").strip()
    try:
        obj = json.loads(txt)
        if isinstance(obj, dict):
            if "is_landmark" in obj:
                return obj
            if isinstance(obj.get("result"), str):
                return _extract_json(obj["result"])
            if isinstance(obj.get("result"), dict):
                return obj["result"]
    except Exception:
        pass
    s, e = txt.find("{"), txt.rfind("}")
    if s != -1 and e > s:
        try:
            return json.loads(txt[s:e + 1])
        except Exception:
            return {}
    return {}
