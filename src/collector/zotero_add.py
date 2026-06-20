#!/usr/bin/env python3
"""
ai4ref — PMID/DOI 로 Zotero 에 단건 추가 (경량, postgres 우회)
==============================================================
zotero_webapi.py(local DB ↔ Zotero web 동기)와 별개의 *빠른 단건 적재* 경로.
  .venv/bin/python -m collector.zotero_add 23099039
  .venv/bin/python -m collector.zotero_add --doi 10.1016/j.clnu.2012.08.013 --collection ERAS
의존: pyzotero · python-dotenv · stdlib urllib(efetch). DB·biopython 불요.
.env 필수: ZOTERO_API_KEY, ZOTERO_USER_ID  (선택: ENTREZ_EMAIL, PUBMED_API_KEY)
  → API 키 발급: zotero.org/settings/keys ("Allow write access" 체크). USER_ID 도 같은 페이지.
"""
from __future__ import annotations
import os, sys, json, time, argparse, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def _load_env():
    try:
        from dotenv import load_dotenv
        load_dotenv(os.path.join(os.path.dirname(__file__), "..", "..", ".env"))
    except Exception:
        pass

def _eutils(fcgi: str, params: dict) -> dict:
    email = os.getenv("ENTREZ_EMAIL", "") or "r.python.ai@gmail.com"
    if "your-" in email.lower() or "example.com" in email.lower():  # .env.example placeholder 가드
        email = "r.python.ai@gmail.com"
    api_key = os.getenv("PUBMED_API_KEY") or None
    if api_key and "your-" in api_key.lower():                       # placeholder → 키 없이 호출
        api_key = None
    q = {**params, "tool": "ai4ref-zotero", "email": email, "retmode": "json"}
    if api_key:
        q["api_key"] = api_key
    url = f"{EUTILS}/{fcgi}?{urllib.parse.urlencode(q)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        d = json.load(r)
    time.sleep(0.12 if api_key else 0.34)  # E-utilities rate limit
    return d

def pmid_from_doi(doi: str) -> str | None:
    ids = _eutils("esearch.fcgi", {"db": "pubmed", "term": f"{doi}[doi]"})["esearchresult"]["idlist"]
    return ids[0] if ids else None

def fetch_meta(pmid: str) -> dict:
    """esummary → 서지 메타 (구조 실측 기반: articleids/authors/fulljournalname)."""
    d = _eutils("esummary.fcgi", {"db": "pubmed", "id": pmid})["result"][pmid]
    doi = next((a["value"] for a in d.get("articleids", []) if a.get("idtype") == "doi"), "")
    authors = [a["name"] for a in d.get("authors", []) if a.get("authtype") == "Author"]
    return {
        "pmid": pmid, "doi": doi,
        "title": (d.get("title", "") or "").rstrip("."),
        "authors": authors,
        "journal": d.get("fulljournalname", ""),
        "date": d.get("pubdate", ""),
        "volume": d.get("volume", ""), "issue": d.get("issue", ""), "pages": d.get("pages", ""),
    }

def to_zotero_item(zot, m: dict) -> dict:
    it = zot.item_template("journalArticle")
    it["title"] = m["title"]
    creators = []
    for name in m["authors"]:                       # "Gustafsson UO" → last=Gustafsson, first=UO
        parts = name.split(" ", 1)
        creators.append({"creatorType": "author",
                         "lastName": parts[0],
                         "firstName": parts[1] if len(parts) > 1 else ""})
    it["creators"] = creators
    it["publicationTitle"] = m["journal"]
    it["date"] = m["date"]
    it["volume"] = m["volume"]; it["issue"] = m["issue"]; it["pages"] = m["pages"]
    it["DOI"] = m["doi"]
    it["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/"
    it["extra"] = f"PMID: {m['pmid']}"
    return it

def add(pmid: str | None = None, doi: str | None = None, collection: str | None = None):
    _load_env()
    key, uid = os.getenv("ZOTERO_API_KEY", ""), os.getenv("ZOTERO_USER_ID", "")
    if not key or not uid or "your-" in (key + uid).lower():
        sys.exit("❌ .env 에 ZOTERO_API_KEY / ZOTERO_USER_ID 설정 필요 "
                 "(zotero.org/settings/keys, write access 체크)")
    if doi and not pmid:
        pmid = pmid_from_doi(doi)
        if not pmid:
            sys.exit(f"❌ DOI {doi} → PMID 변환 실패")
    m = fetch_meta(pmid)
    from pyzotero import zotero
    zot = zotero.Zotero(uid, "user", key)
    it = to_zotero_item(zot, m)
    if collection:
        cols = {c["data"]["name"]: c["key"] for c in zot.collections()}
        if collection in cols:
            it["collections"] = [cols[collection]]
        else:
            cres = zot.create_collection([{"name": collection}])
            ckey = cres["successful"]["0"]["key"]
            it["collections"] = [ckey]
            print(f"  + 컬렉션 '{collection}' 신규 생성 (key={ckey})")
    res = zot.create_items([it])
    if res.get("successful"):
        k = res["successful"]["0"]["key"]
        print(f"✅ Zotero 추가: {m['title'][:60]}… (key={k}, PMID={pmid}, DOI={m['doi']})")
    else:
        print("❌ 실패:", json.dumps(res.get("failed"), ensure_ascii=False))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="PMID/DOI → Zotero 단건 추가")
    ap.add_argument("pmid", nargs="?", help="PubMed ID")
    ap.add_argument("--doi", help="DOI (PMID 없을 때)")
    ap.add_argument("--collection", help="Zotero 컬렉션명 (없으면 루트)")
    ap.add_argument("--dry", action="store_true", help="Zotero 미전송, 메타만 출력")
    a = ap.parse_args()
    if not a.pmid and not a.doi:
        sys.exit("PMID 또는 --doi 필요")
    if a.dry:
        _load_env()
        pid = a.pmid or pmid_from_doi(a.doi)
        print(json.dumps(fetch_meta(pid), ensure_ascii=False, indent=2))
    else:
        add(pmid=a.pmid, doi=a.doi, collection=a.collection)
