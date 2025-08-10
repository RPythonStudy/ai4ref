import os, re
import logging
from pyzotero import zotero
from common.database import get_db_connection
from common.logger import log_info, log_error, log_debug

def quiet_http_logs(level=logging.WARNING):
    for name in ("httpx", "httpcore"):
        lg = logging.getLogger(name)
        lg.setLevel(level)      # INFO/DEBUG -> 숨김, WARNING 이상만 출력
        lg.propagate = False    # 루트 로거로 전파 차단

def sync_collection_keys(zot, conn):
    cur = conn.cursor()
    remote = {c["data"]["name"]: c["data"]["key"] for c in zot.collections()}
    cur.execute("SELECT id, name, COALESCE(zotero_key,'') FROM collections ORDER BY id")
    matched = created = 0
    for cid, name, key in cur.fetchall():
        if key: continue
        rkey = remote.get(name)
        if rkey:
            cur.execute("UPDATE collections SET zotero_key=%s, updated_at=NOW() WHERE id=%s", (rkey, cid)); matched += 1
        else:
            res = zot.create_collection([{"name": name, "parentCollection": ""}])
            nkey = res["success"]["0"]
            cur.execute("UPDATE collections SET zotero_key=%s, updated_at=NOW() WHERE id=%s", (nkey, cid)); created += 1
    conn.commit(); cur.close()
    log_info(f"sync_collection_keys: matched={matched}, created={created}")
    return {"matched": matched, "created": created}

# --- 로컬(DB) ↔ 웹(Zotero) 차이 검증: PMID(Extra) 기준 ---
PMID_RX = re.compile(r"(?mi)^\s*PMID:\s*(\d+)\b")

def diff_local_vs_web(check_pdf: bool = True, max_pdf_checks: int = 500):
    """
    로컬(DB)의 pmid→{collection_keys, pdf_path}와
    웹(Zotero)의 pmid→{collection_keys}를 비교.
    반환: only_local / only_web / col_diffs / pdf_needed
    """
    # 1) 로컬: pmid -> {cols, pdf}
    conn = get_db_connection(); cur = conn.cursor()
    cur.execute("""
      SELECT p.pmid,
             COALESCE(p.pdf_path, ''),
             ARRAY_AGG(DISTINCT c.zotero_key ORDER BY c.zotero_key)
      FROM papers p
      JOIN paper_collection pc ON pc.paper_id = p.id
      JOIN collections c       ON c.id = pc.collection_id
      WHERE p.pmid IS NOT NULL
        AND c.zotero_key IS NOT NULL
      GROUP BY p.pmid, p.pdf_path
    """)
    local = {str(pmid): {"cols": set(cols), "pdf": pdf} for pmid, pdf, cols in cur.fetchall()}
    cur.close(); conn.close()

    # 2) 웹: pmid -> {cols}, pmid -> item_key (Extra에서 PMID 추출)
    user_id = os.getenv('ZOTERO_USER_ID'); api_key = os.getenv('ZOTERO_API_KEY')
    zot = zotero.Zotero(user_id, 'user', api_key)

    web, key_by_pmid = {}, {}
    for it in zot.everything(zot.top()):
        m = PMID_RX.search((it["data"].get("extra") or ""))
        if not m: continue
        pmid = m.group(1)
        cols = set(it["data"].get("collections", []))
        web.setdefault(pmid, {"cols": set()})["cols"] |= cols
        key_by_pmid[pmid] = it["key"]

    # 3) 차이 계산
    L, W = set(local), set(web)
    only_local = sorted(L - W)          # 웹에 만들어야 할 아이템
    only_web   = sorted(W - L)          # 참고용

    col_diffs = []
    for pmid in sorted(L & W):
        desired, remote = local[pmid]["cols"], web[pmid]["cols"]
        to_add, to_remove = sorted(desired - remote), sorted(remote - desired)
        if to_add or to_remove: col_diffs.append((pmid, to_add, to_remove))

    # 4) PDF 첨부 필요(선택적으로 children 조회)
    pdf_needed = []
    if check_pdf:
        checks = 0
        for pmid in sorted(L & W):
            if checks >= max_pdf_checks: break
            pdf_path = local[pmid]["pdf"]
            if not pdf_path: continue
            item_key = key_by_pmid.get(pmid)
            if not item_key: continue
            has_pdf = any(
                (ch.get("data", {}).get("itemType") == "attachment" and
                 ch.get("data", {}).get("contentType") == "application/pdf")
                for ch in zot.children(item_key)
            )
            if not has_pdf: pdf_needed.append((pmid, pdf_path, item_key))
            checks += 1

    # 5) 요약
    log_info("=== 로컬↔웹 비교 요약 ===")
    log_info(f"로컬 PMID: {len(local)} / 웹 PMID: {len(web)}")
    log_info(f"웹에 없는 로컬 PMID: {len(only_local)}")
    log_info(f"컬렉션 불일치: {len(col_diffs)}")
    if check_pdf: log_info(f"PDF 첨부 필요(최대 {max_pdf_checks}건 검사): {len(pdf_needed)}")

    return {"only_local": only_local, "only_web": only_web, "col_diffs": col_diffs, "pdf_needed": pdf_needed}

def push_only_local_to_web(limit: int = 100):
    user_id = os.getenv('ZOTERO_USER_ID'); api_key = os.getenv('ZOTERO_API_KEY')
    zot = zotero.Zotero(user_id, 'user', api_key)

    diffs = diff_local_vs_web(check_pdf=False)
    targets = diffs["only_local"][:limit]
    if not targets:
        log_info("웹에 생성할 로컬 PMID가 없습니다")
        return 0

    conn = get_db_connection(); cur = conn.cursor()
    created = 0

    for pmid in targets:
        # 로컬 메타+목표 컬렉션 키 모으기
        cur.execute("""
          SELECT p.id, p.title, COALESCE(p.doi,''), p.publication_year,
                 COALESCE(p.pmcid,''), COALESCE(p.arxiv_id,''), COALESCE(p.pdf_path,'')
          FROM papers p
          WHERE p.pmid = %s
        """, (pmid,))
        row = cur.fetchone()
        if not row: continue
        paper_id, title, doi, year, pmcid, arxiv, pdf_path = row

        cur.execute("""
          SELECT ARRAY_AGG(DISTINCT c.zotero_key ORDER BY c.zotero_key)
          FROM paper_collection pc JOIN collections c ON c.id = pc.collection_id
          WHERE pc.paper_id = %s AND c.zotero_key IS NOT NULL
        """, (paper_id,))
        keys = cur.fetchone()[0] or []
        desired = set(keys)

        # 아이템 생성 (Extra에 PMID와 내부 ID 저장)
        extra_lines = [f"PMID: {pmid}", f"ai4ref_paper_id: {paper_id}"]
        if pmcid: extra_lines.append(f"PMCID: {pmcid}")
        if arxiv: extra_lines.append(f"arXiv: {arxiv}")

        item = {
            "itemType": "journalArticle",
            "title": title or f"(untitled {pmid})",
            "date": str(year or ""),
            "DOI": doi,
            "extra": "\n".join(extra_lines),
            "collections": sorted(list(desired)),
        }
        res = zot.create_items([item])
        item_key = res["success"]["0"]
        created += 1
        log_info(f"[CREATE] PMID {pmid} → {item_key}")

        # PDF 첨부(있으면)
        if pdf_path and os.path.isfile(pdf_path):
            try:
                zot.attachment_simple([pdf_path], parentid=item_key)
                log_info(f"[ATTACH] PMID {pmid} ← {os.path.basename(pdf_path)}")
            except Exception as e:
                log_error(f"[ATTACH-FAIL] PMID {pmid}: {e}")

        # paper_collection 매핑 업데이트
        cur.execute("UPDATE paper_collection SET zotero_item_key=%s WHERE paper_id=%s", (item_key, paper_id))
        conn.commit()

    cur.close(); conn.close()
    log_info(f"웹에 신규 아이템 생성: {created}개")
    return created


if __name__ == "__main__":
    quiet_http_logs()  # <- 먼저 호출
    user_id = os.getenv('ZOTERO_USER_ID'); api_key = os.getenv('ZOTERO_API_KEY')
    if not (user_id and api_key):
        log_error("ZOTERO_USER_ID / ZOTERO_API_KEY 환경변수가 필요합니다"); raise SystemExit(1)
    zot = zotero.Zotero(user_id, 'user', api_key)

    conn = get_db_connection()
    sync_collection_keys(zot, conn)
    conn.close()

    diffs = diff_local_vs_web()
    # 예시 일부만 출력
    if diffs["only_local"][:10]:
        log_info(f"예시 only_local: {', '.join(diffs['only_local'][:10])}")
    for pmid, add, remove in diffs["col_diffs"][:5]:
        log_info(f"PMID {pmid} | to_add={add or '-'} | to_remove={remove or '-'}")
    for pmid, path, _ in diffs["pdf_needed"][:5]:
        log_info(f"PDF 필요: PMID {pmid} -> {os.path.basename(path)}")
        
    push_only_local_to_web(limit=50)
   
