#!/usr/bin/env python3
"""
ai4ref 알림 PoC — guideline-anchored new-treatment alert
=========================================================
RPython 연구회 요청: 특정 임상영역에서 *새 치료법(중재 비교 RCT)이 제안되면*
Telegram 으로 1편 보내기. vault [[2025-08_ai4ref]] §응용 PoC 의 첫 구현.

파이프라인 (검색 → 선별 → 발송):
  1. 검색(recall)  : 진료지침-앵커 query(P·도메인·RCT 게이트, I/C 개방) + rolling 윈도우 → esearch
  2. 선별(top-1)   : esummary 메타 → LLM 게이트("권고 변경 근거?") → 저널티어·최신 랭킹 → top-1
  3. 발송          : Telegram (Bot API; 토큰 없으면 dry-run 출력)

설계 주의 (hub 의 "회고 ≠ 전향" 경고):
  - 저널티어 표는 영향력의 *근사 휴리스틱*. 단일 저널 하드코딩 금지 — 표로 일반화.
  - "새 치료법" 판정의 원칙적 층은 LLM 게이트(현재 heuristic stub; TODO: Claude 연결).
  - P·O 좁게/ I·C 개방 → 기존 중재 목록에 안 묶여 *새* 치료법 포착.

구현 메모:
  - 의존성 0 (stdlib urllib) — DB·biopython 없이 즉시 실행/검증 가능한 PoC.
  - production 전환 시 src/collector.pubmed_esearch(Entrez)+postgres 재사용,
    esummary 대신 efetch 로 abstract 확보해 LLM 게이트 강화.
"""
from __future__ import annotations
import os, sys, json, time, urllib.parse, urllib.request
from datetime import datetime, timedelta, timezone

try:
    import yaml  # config 로드 (가용)
except Exception:
    yaml = None

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = os.getenv("ENTREZ_EMAIL", "r.python.ai@gmail.com")
API_KEY = os.getenv("PUBMED_API_KEY") or None
TOOL = "ai4ref-alert-poc"

# ── 저널 영향력 티어 (단일 하드코딩 X, 표로 일반화 — 휴리스틱 proxy) ───────────
JOURNAL_TIER = {
    100: ["n engl j med", "new england journal", "lancet", "jama:", "jama "],
    60:  ["annals of surgery", "british journal of anaesthesia", "anesthesiology",
          "jama surgery", "intensive care medicine", "critical care medicine",
          "bmj", "british medical journal", "anaesthesia", "circulation"],
}

def _eutils(fcgi: str, params: dict) -> dict:
    q = {**params, "tool": TOOL, "email": EMAIL, "retmode": "json"}
    if API_KEY:
        q["api_key"] = API_KEY
    url = f"{EUTILS}/{fcgi}?{urllib.parse.urlencode(q)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        data = json.load(r)
    time.sleep(0.12 if API_KEY else 0.34)  # E-utilities rate limit (10/s with key, 3/s without)
    return data

def esearch(term: str, retmax: int = 200, mindate: str | None = None,
            maxdate: str | None = None) -> tuple[int, list[str]]:
    """PubMed esearch → (총건수, PMID 리스트). 날짜는 YYYY/MM/DD (pdat)."""
    params = {"db": "pubmed", "term": term, "retmax": str(retmax)}
    if mindate and maxdate:
        params.update({"datetype": "pdat", "mindate": mindate, "maxdate": maxdate})
    d = _eutils("esearch.fcgi", params)["esearchresult"]
    return int(d.get("count", 0)), d.get("idlist", [])

def esummary(pmids: list[str]) -> list[dict]:
    """PMID 리스트 → 선별용 메타(title·journal·year·pubtypes·authors)."""
    if not pmids:
        return []
    d = _eutils("esummary.fcgi", {"db": "pubmed", "id": ",".join(pmids)})["result"]
    out = []
    for pmid in d.get("uids", []):
        it = d[pmid]
        out.append({
            "pmid": pmid,
            "title": it.get("title", "").rstrip("."),
            "journal": it.get("fulljournalname") or it.get("source", ""),
            "pubdate": it.get("sortpubdate") or it.get("pubdate", ""),
            "year": (it.get("sortpubdate") or it.get("pubdate", "") or "")[:4],
            "pubtypes": [p for p in it.get("pubtype", [])],
            "authors": "; ".join(a.get("name", "") for a in it.get("authors", [])[:3]),
            "doi": next((aid.get("value") for aid in it.get("articleids", [])
                         if aid.get("idtype") == "doi"), ""),
        })
    return out

# ── 선별 ─────────────────────────────────────────────────────────────────────
def journal_weight(journal: str) -> int:
    j = (journal or "").lower()
    for w, keys in JOURNAL_TIER.items():
        if any(k in j for k in keys):
            return w
    return 10

def is_rct(paper: dict) -> bool:
    return any("randomized controlled trial" in p.lower() for p in paper["pubtypes"])

def llm_practice_changing_gate(paper: dict) -> tuple[bool, str]:
    """'진료 권고를 바꾸/도전하는 새 치료전략 평가인가?' 판정.
    *** heuristic STUB *** — production: Claude 에 title+abstract 던져 분류.
    현재: RCT + 제목의 비교/시험 신호로 근사."""
    title = paper["title"].lower()
    signals = ["versus", " vs ", " vs.", "compared", "comparison", "randomi", "trial",
             "효과", "비교"]
    rct = is_rct(paper)
    hit = [s for s in signals if s in title]
    ok = rct and bool(hit)
    reason = f"RCT={rct}, 제목신호={hit or '없음'}"
    return ok, reason

def select_top1(papers: list[dict]) -> tuple[dict | None, list[dict]]:
    """게이트 통과분을 (저널티어×10 + 최신연도) 로 랭크 → top-1. 점수표 동반 반환."""
    scored = []
    for p in papers:
        gate_ok, gate_reason = llm_practice_changing_gate(p)
        jw = journal_weight(p["journal"])
        try:
            recency = int(p["year"])
        except ValueError:
            recency = 0
        score = jw * 10 + recency  # 저널티어 지배, 동급내 최신 우선
        scored.append({**p, "gate_ok": gate_ok, "gate_reason": gate_reason,
                       "journal_w": jw, "score": score})
    survivors = [s for s in scored if s["gate_ok"]]
    ranked = sorted(survivors or scored, key=lambda s: s["score"], reverse=True)
    top = ranked[0] if ranked else None
    return top, sorted(scored, key=lambda s: s["score"], reverse=True)

# ── 발송 ─────────────────────────────────────────────────────────────────────
def format_alert(p: dict, anchor: dict) -> str:
    link = f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/"
    return (
        f"🔔 *새 치료 근거 알림* (RPython 연구회)\n"
        f"앵커: {anchor.get('guideline') or anchor.get('id')}\n\n"
        f"*{p['title']}*\n"
        f"{p['authors']} 외 · {p['journal']} ({p['year']})\n"
        f"PMID {p['pmid']} · {link}\n\n"
        f"선택 근거: 저널티어={p['journal_w']} · 게이트({p['gate_reason']}) · score={p['score']}"
    )

def send_telegram(msg: str) -> bool:
    token, chat = os.getenv("TELEGRAM_BOT_TOKEN"), os.getenv("TELEGRAM_CHAT_ID")
    if not (token and chat):
        print("[dry-run] TELEGRAM_BOT_TOKEN/CHAT_ID 미설정 — 발송 생략, 메시지 미리보기:\n")
        print(msg)
        return False
    data = urllib.parse.urlencode({"chat_id": chat, "text": msg,
                                   "parse_mode": "Markdown"}).encode()
    url = f"https://api.telegram.org/bot{token}/sendMessage"
    with urllib.request.urlopen(urllib.request.Request(url, data=data), timeout=20) as r:
        ok = json.load(r).get("ok", False)
    print(f"[telegram] 발송 {'성공' if ok else '실패'}")
    return ok

# ── 오케스트레이션 ───────────────────────────────────────────────────────────
def run_alert(anchor: dict, window_days: int | None = 30,
              mindate: str | None = None, maxdate: str | None = None,
              deliver: bool = False, verbose: bool = True) -> dict | None:
    """앵커 1개에 대해 검색→선별→(발송). 반환 = 선택된 논문(or None)."""
    if window_days and not (mindate and maxdate):
        today = datetime.now(timezone.utc).date()
        mindate = (today - timedelta(days=window_days)).strftime("%Y/%m/%d")
        maxdate = today.strftime("%Y/%m/%d")
    term = " ".join(anchor["query"].split())
    count, pmids = esearch(term, mindate=mindate, maxdate=maxdate)
    if verbose:
        print(f"[검색] 윈도우 {mindate}~{maxdate} | 후보 {count}건 (PMID {len(pmids)} 회수)")
    if not pmids:
        print("[검색] 후보 없음 — 알림 없음")
        return None
    papers = esummary(pmids)
    top, table = select_top1(papers)
    if verbose:
        print(f"\n[선별] 랭킹 (상위 {min(8,len(table))}):")
        print(f"  {'PMID':>9}  {'gate':4}  {'jw':>3}  {'score':>5}  year  journal / title")
        for s in table[:8]:
            mark = "★" if top and s["pmid"] == top["pmid"] else " "
            print(f" {mark}{s['pmid']:>9}  {str(s['gate_ok']):4.4}  {s['journal_w']:>3}  "
                  f"{s['score']:>5}  {s['year']}  {s['journal'][:28]} / {s['title'][:42]}")
    if not top:
        print("\n[선별] 게이트 통과 후보 없음")
        return None
    msg = format_alert(top, anchor)
    print(f"\n[선택] PMID {top['pmid']} — {top['title'][:60]}")
    if deliver:
        send_telegram(msg)
    else:
        print("\n--- 발송 미리보기 (deliver=False) ---\n" + msg)
    return top

def load_anchors(path: str | None = None) -> list[dict]:
    path = path or os.path.join(os.path.dirname(__file__), "..", "..", "config", "alert_anchors.yml")
    if yaml is None:
        raise RuntimeError("pyyaml 필요 (config 로드)")
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f).get("anchors", [])

if __name__ == "__main__":
    # 기본 = RELIEF 검증 모드: ERAS 수액 앵커 + 2018 윈도우 → RELIEF(29742967) top-1 기대.
    # 실서비스 = run_alert(anchor, window_days=30, deliver=True) (rolling).
    anchors = load_anchors()
    anchor = next((a for a in anchors if a["id"] == "eras_fluid_abdominal"), anchors[0])
    print(f"=== ai4ref 알림 PoC — 앵커 '{anchor['id']}' ===\n")

    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if mode == "validate":
        # ── 두-검증 (시간축 분리, T = 지침 제정 시점) ───────────────────────────────
        #    A: guideline_refs (T 이전) 를 재현하나 = 시스템 본질 성능 (recall)
        #    B: post_guideline_landmarks (T 이후) 를 포착하나 = alert 성능
        #    (query) AND PMID[uid] AND <윈도우> 로 포함 판정. 서로 다른 논문.
        q = anchor["query"]
        T = anchor.get("guideline", {}).get("date", "2012-12")
        Ty = int(T[:4])

        refs = anchor.get("guideline_refs", [])
        print(f"=== 검증 A: 지침 근거 재현 (T={T} 이전) — guideline_refs {len(refs)}편 ===")
        a_hit = 0
        for pmid in refs:
            cnt, _ = esearch(f"({q}) AND {pmid}[uid]", retmax=1,
                             mindate="1900/01/01", maxdate=f"{Ty}/12/31")
            ok = cnt >= 1; a_hit += ok
            print(f"  {'✅' if ok else '❌'} {pmid}")
        print(f"  → recall A = {a_hit}/{len(refs)}\n")

        lms = anchor.get("post_guideline_landmarks", [])
        print(f"=== 검증 B: 랜드마크 포착 (T={T} 이후) — landmarks {len(lms)}편 ===")
        b_hit = 0
        for pmid in lms:
            cnt, _ = esearch(f"({q}) AND {pmid}[uid]", retmax=1,
                             mindate=f"{Ty + 1}/01/01", maxdate="2030/12/31")
            ok = cnt >= 1; b_hit += ok
            print(f"  {'✅' if ok else '❌'} {pmid}")
        print(f"  → 포착 B = {b_hit}/{len(lms)}\n")

        print(f"[종합] 검증 A recall {a_hit}/{len(refs)} · 검증 B 포착 {b_hit}/{len(lms)}")
    else:  # rolling (실서비스 시뮬레이션 — top-1 선별·발송)
        run_alert(anchor, window_days=int(sys.argv[2]) if len(sys.argv) > 2 else 30,
                  deliver=os.getenv("ALERT_DELIVER") == "1")
