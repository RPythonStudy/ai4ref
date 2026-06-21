# src/notify/telegram_sink.py
"""텔레그램 Bot API sendMessage sink.

토큰/chat_id 가 없으면 available()=False 로 멀티탭이 자동 skip(fail-soft).
의존성: 표준 라이브러리(urllib)만 — 별도 설치 불필요.
.env: TELEGRAM_BOT_TOKEN, TELEGRAM_CHAT_ID
"""
import os
import json
import urllib.parse
import urllib.request

from notify.base import Sink, LandmarkItem
from notify.stdout_sink import render

API = "https://api.telegram.org"


class TelegramSink(Sink):
    name = "telegram"

    def __init__(self, cfg=None):
        self.token = os.getenv("TELEGRAM_BOT_TOKEN", "").strip()
        self.chat_id = os.getenv("TELEGRAM_CHAT_ID", "").strip()

    def available(self) -> tuple[bool, str]:
        if not self.token or not self.chat_id:
            return False, "TELEGRAM_BOT_TOKEN/TELEGRAM_CHAT_ID 미설정 (.env)"
        return True, ""

    def send(self, item: LandmarkItem):
        data = urllib.parse.urlencode({
            "chat_id": self.chat_id,
            "text": render(item),
            "disable_web_page_preview": "false",
        }).encode()
        url = f"{API}/bot{self.token}/sendMessage"
        req = urllib.request.Request(url, data=data)
        with urllib.request.urlopen(req, timeout=20) as r:
            res = json.load(r)
        if not res.get("ok"):
            raise RuntimeError(f"telegram API 오류: {res}")
        return str(res["result"]["message_id"])
