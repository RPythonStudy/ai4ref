# src/alert/llm_gate.py
"""랜드마크 판정 — 교체형 LLM 백엔드 (sink 처럼 갈아끼움).

backend:
  stub       = 개발용(가짜 판정, 키 불필요)
  claude_cli = 로컬 Claude Code · OAuth 구독 · 추가비용 0   ← 운영 권장
  api        = Anthropic API · 클라우드 · 토큰 과금(예정)

2단계 판정:
  ① 관련성 — KQ 의 대상(P) AND 중재(I) 를 *둘 다* 다루는가 (한쪽만 = false)
  ② 랜드마크 — 지침과 다른 결론·프로토콜의 practice-changing 인가
              (strict = 기대 설계 수준 요구 / loose = recall 편향)

judge() -> {"relevant": bool, "is_landmark": bool, "reason": str}
"""
import json
import subprocess

from common.features import log

# question_type → 기대 연구설계 (게이트가 strict 일 때 요구하는 근거 수준)
DESIGN_LABEL = {
    "intervention": "무작위 대조시험(RCT)",
    "diagnostic":   "진단 정확도 연구",
    "prognostic":   "코호트(예후) 연구",
    "predictive":   "치료반응 예측 연구(RCT 서브그룹/바이오마커)",
}

PROMPT = """너는 진료지침 감시 보조다. 아래 논문을 2단계로 엄격히 판정하라.

[KQ] {kq}
  · 대조(C): {comparison}
  · 주요 결과(O): {outcome}
[현행 지침] {guideline}
[기대 연구설계] {design}
[엄격도] {strictness}

[논문]
  제목: {title}
  초록: {abstract}

판정 단계:
① 관련성 — 이 논문이 KQ 의 **대상(P)과 중재(I)를 둘 다** 직접 다루는가?
   · 한쪽만 맞으면 relevant=false. (대상은 맞으나 다른 중재 / 중재는 맞으나 다른 대상)
   · 예: KQ='복부수술의 수액전략' 인데 논문이 '신장이식의 수액전략'(대상 불일치) 이거나
        '복부수술의 로봇 vs 복강경'(중재 불일치) 이면 relevant=false.
② 랜드마크 — (relevant 일 때만) 현행 지침과 *다른 결론·프로토콜* 을 제시하는 practice-changing 인가?
   · strict: 기대 설계({design}) 수준의 근거여야 is_landmark=true. 그 외 false.
   · loose: 설계 무관, 경계선이면 포함 (recall 편향 — 놓침이 가장 비싼 오류).

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
    """로컬 Claude Code(OAuth)에 위임. 추가 API 비용 없음. 분류는 Haiku 로 충분."""
    qtype = kq.get("question_type", "")
    design = DESIGN_LABEL.get(qtype, "해당 유형의 적정 연구설계")
    prompt = PROMPT.format(
        kq=kq.get("kq", ""),
        comparison=kq.get("comparison", "") or "(개방)",
        outcome=kq.get("outcome", ""),
        guideline=kq.get("guideline", ""),
        design=design,
        strictness=kq.get("strictness", "loose"),
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
