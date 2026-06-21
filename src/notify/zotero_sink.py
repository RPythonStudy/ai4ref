# src/notify/zotero_sink.py
"""Zotero sink — 탐지된 랜드마크를 참고문헌함(Zotero)에 자동 보관.

기존 collector.zotero_add 의 메타·아이템 로직 재사용. pyzotero·자격증명은 *지연 로드* →
키/패키지 없으면 available()=False 로 멀티탭이 자동 skip(fail-soft).

.env: ZOTERO_API_KEY, ZOTERO_USER_ID  (zotero.org/settings/keys, write access)
features.yml: notify.zotero.collection (기본 'ai4ref-landmarks')
"""
import os

from notify.base import Sink, LandmarkItem
from common.features import log


class ZoteroSink(Sink):
    name = "zotero"

    def __init__(self, cfg=None):
        cfg = cfg or {}
        self.collection = cfg.get("collection", "ai4ref-landmarks")
        self.key = os.getenv("ZOTERO_API_KEY", "").strip()
        self.uid = os.getenv("ZOTERO_USER_ID", "").strip()
        self._zot = None
        self._col_key = None

    def available(self) -> tuple[bool, str]:
        if not self.key or not self.uid or "your-" in (self.key + self.uid).lower():
            return False, "ZOTERO_API_KEY/ZOTERO_USER_ID 미설정 (.env, write access)"
        try:
            import pyzotero  # noqa: F401
        except Exception:
            return False, "pyzotero 미설치 (pip install pyzotero)"
        return True, ""

    def _client(self):
        if self._zot is None:
            from pyzotero import zotero
            self._zot = zotero.Zotero(self.uid, "user", self.key)
        return self._zot

    def _collection_key(self) -> str:
        if self._col_key is not None:
            return self._col_key
        zot = self._client()
        cols = {c["data"]["name"]: c["key"] for c in zot.collections()}
        if self.collection in cols:
            self._col_key = cols[self.collection]
        else:
            res = zot.create_collection([{"name": self.collection}])
            self._col_key = res["successful"]["0"]["key"]
            log.info(f"[zotero] 컬렉션 '{self.collection}' 신규 생성")
        return self._col_key

    def send(self, item: LandmarkItem):
        from collector.zotero_add import fetch_meta, to_zotero_item
        zot = self._client()
        m = fetch_meta(item.pmid)
        it = to_zotero_item(zot, m)
        it["collections"] = [self._collection_key()]
        if item.reason:                                   # 한 줄 근거를 Extra 메모에 보강
            it["extra"] = (it.get("extra", "") + f"\nai4ref: {item.reason}").strip()
        res = zot.create_items([it])
        if res.get("successful"):
            return res["successful"]["0"]["key"]
        raise RuntimeError(f"zotero create 실패: {res.get('failed')}")
