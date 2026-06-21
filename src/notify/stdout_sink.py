# src/notify/stdout_sink.py
"""화면출력 sink — 의존성 0 · 항상 동작. 멀티탭의 기본 구멍(개발·테스트)."""
from notify.base import Sink, LandmarkItem


def render(item: LandmarkItem) -> str:
    head = f"🔬 [{item.kq}] 새 랜드마크"
    meta = " ".join(x for x in [item.journal, item.year] if x)
    title = f"{item.title}" + (f" ({meta})" if meta else "")
    lines = [head, title]
    if item.reason:
        lines.append(f"▸ {item.reason}")
    if item.url:
        lines.append(f"🔗 {item.url}")
    return "\n".join(lines)


class StdoutSink(Sink):
    name = "stdout"

    def __init__(self, cfg=None):
        pass

    def send(self, item: LandmarkItem):
        print("\n" + render(item) + "\n")
        return "stdout"
