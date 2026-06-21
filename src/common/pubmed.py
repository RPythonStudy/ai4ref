# src/common/pubmed.py
"""경량 PubMed(E-utilities) 클라이언트 — esearch + efetch 메타.

collector/alert 공용. API 키 없이 동작(rate limit 대응 sleep). 키 있으면 더 빠르게.
.env: ENTREZ_EMAIL(권장), PUBMED_API_KEY(선택)
"""
import os
import time
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from common.features import log

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = os.getenv("ENTREZ_EMAIL", "r.python.ai@gmail.com")
APIKEY = os.getenv("PUBMED_API_KEY", "")
if APIKEY.startswith("your-"):          # .env.example placeholder 가드
    APIKEY = ""

_SLEEP = 0.11 if APIKEY else 0.34


def _params(extra: dict) -> dict:
    q = {"tool": "ai4ref", "email": EMAIL}
    if APIKEY:
        q["api_key"] = APIKEY
    q.update(extra)
    return q


def _get(url: str, timeout: int = 30) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def esearch(term: str, retmax: int = 300, datetype: str = None,
            reldate: int = None, mindate: str = None, maxdate: str = None) -> list:
    q = _params({"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax})
    if datetype:
        q["datetype"] = datetype
    if reldate is not None:
        q["reldate"] = reldate
    if mindate and maxdate:
        q["mindate"] = mindate
        q["maxdate"] = maxdate
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(q)}"
    d = json.loads(_get(url))
    time.sleep(_SLEEP)
    return d.get("esearchresult", {}).get("idlist", [])


def efetch_meta(pmids: list) -> dict:
    """pmid → {title, abstract, journal, year}. 100개씩 배치."""
    out = {}
    for i in range(0, len(pmids), 100):
        batch = pmids[i:i + 100]
        q = _params({"db": "pubmed", "id": ",".join(batch),
                     "rettype": "abstract", "retmode": "xml"})
        url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(q)}"
        try:
            root = ET.fromstring(_get(url))
        except Exception as e:
            log.error(f"[pubmed] efetch 파싱 실패: {e}")
            time.sleep(_SLEEP)
            continue
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//MedlineCitation/PMID") or ""
            te = art.find(".//Article/ArticleTitle")
            title = "".join(te.itertext()).strip() if te is not None else ""
            abst = " ".join("".join(a.itertext())
                            for a in art.findall(".//Abstract/AbstractText")).strip()
            journal = (art.findtext(".//Journal/ISOAbbreviation")
                       or art.findtext(".//Journal/Title") or "")
            year = (art.findtext(".//Journal//PubDate/Year")
                    or (art.findtext(".//Journal//PubDate/MedlineDate") or "")[:4])
            out[pmid] = {"title": title, "abstract": abst,
                         "journal": journal, "year": year}
        time.sleep(_SLEEP)
    return out
