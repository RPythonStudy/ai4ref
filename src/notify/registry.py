# src/notify/registry.py
"""멀티탭(분배기): features.yml 에서 켜진 sink 만 모아 랜드마크를 fan-out.

한 sink 가 실패해도 나머지는 계속 — 본체 흐름을 막지 않는다.
새 통로 추가 = SINK_CLASSES 에 한 줄 등록.
"""
from common.features import log
from notify.stdout_sink import StdoutSink
from notify.telegram_sink import TelegramSink
from notify.zotero_sink import ZoteroSink

SINK_CLASSES = {
    "stdout": StdoutSink,
    "telegram": TelegramSink,
    "zotero": ZoteroSink,
}


def build_sinks(features: dict) -> list:
    notify_cfg = (features or {}).get("notify", {}) or {}
    active = []
    for name, cls in SINK_CLASSES.items():
        cfg = notify_cfg.get(name, {}) or {}
        if not cfg.get("enabled"):
            continue
        sink = cls(cfg)
        ok, why = sink.available()
        if not ok:
            log.warning(f"[notify] '{name}' 켜졌으나 prereq 미충족 → skip: {why}")
            continue
        active.append(sink)
        log.info(f"[notify] sink 활성: {name}")
    return active


def fan_out(item, sinks: list) -> dict:
    out = {}
    for s in sinks:
        try:
            out[s.name] = s.send(item)
        except Exception as e:
            out[s.name] = None
            log.error(f"[notify:{s.name}] 전송 실패: {e}")
    return out
