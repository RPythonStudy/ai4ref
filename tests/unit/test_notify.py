# tests/unit/test_notify.py
"""Sink 멀티탭(fan-out) 단위테스트 — 배관 결정론 가드.

외부 네트워크(telegram urllib·Zotero pyzotero)는 호출하지 않는다. 검증 대상:
  · render — 제목·메타·근거·링크 포함, 내부 abstract 제외 (FR-006)
  · LandmarkItem — url 자동 파생·key="kq:pmid" (FR-011)
  · build_sinks — enabled 토글 + prereq skip(fail-soft) (FR-002·003)
  · fan_out — 한 통로 전송 실패 격리, 나머지 계속 (FR-005)
"""
import pytest

from notify.base import Sink, LandmarkItem
from notify.stdout_sink import render, StdoutSink
from notify.registry import build_sinks, fan_out


ITEM = LandmarkItem(
    kq="ERAS 수액전략", pmid="29742967", title="RELIEF trial",
    journal="NEJM", year="2018", reason="practice-changing",
    abstract="INTERNAL-ABSTRACT-SHOULD-NOT-RENDER",
)


# ── LandmarkItem ────────────────────────────────────────────────────────────

def test_landmark_url_auto_derived():
    """url 미지정 시 pmid 에서 PubMed 링크 자동 파생 (FR-011)."""
    assert ITEM.url == "https://pubmed.ncbi.nlm.nih.gov/29742967/"


def test_landmark_key():
    """dedup 키 = kq:pmid (§004 연계)."""
    assert ITEM.key == "ERAS 수액전략:29742967"


def test_landmark_explicit_url_kept():
    it = LandmarkItem(kq="k", pmid="1", title="t", url="https://x")
    assert it.url == "https://x"


# ── render (FR-006) ─────────────────────────────────────────────────────────

def test_render_includes_fields():
    out = render(ITEM)
    assert "RELIEF trial" in out and "NEJM" in out and "2018" in out
    assert "practice-changing" in out
    assert ITEM.url in out


def test_render_excludes_internal_abstract():
    """내부 게이트용 abstract 는 렌더에 노출되지 않는다 (FR-006)."""
    assert "INTERNAL-ABSTRACT-SHOULD-NOT-RENDER" not in render(ITEM)


def test_render_without_meta():
    it = LandmarkItem(kq="k", pmid="1", title="제목만")
    out = render(it)
    assert "제목만" in out  # 메타·근거 없어도 동작


# ── build_sinks — 토글 + prereq skip (FR-002·003) ──────────────────────────

def test_build_sinks_stdout_always_active(monkeypatch):
    """stdout 은 의존성 0 → 켜면 항상 활성."""
    feats = {"notify": {"stdout": {"enabled": True}}}
    sinks = build_sinks(feats)
    assert [s.name for s in sinks] == ["stdout"]


def test_build_sinks_telegram_skipped_without_token(monkeypatch):
    """telegram 켜졌으나 토큰 미설정 → prereq 미충족 skip, stdout 은 활성 (FR-003)."""
    monkeypatch.delenv("TELEGRAM_BOT_TOKEN", raising=False)
    monkeypatch.delenv("TELEGRAM_CHAT_ID", raising=False)
    feats = {"notify": {"stdout": {"enabled": True}, "telegram": {"enabled": True}}}
    names = [s.name for s in build_sinks(feats)]
    assert "stdout" in names and "telegram" not in names


def test_build_sinks_disabled_not_instantiated():
    """enabled:false 통로는 활성 목록에 없다 (FR-002)."""
    feats = {"notify": {"stdout": {"enabled": False}, "telegram": {"enabled": False}}}
    assert build_sinks(feats) == []


def test_build_sinks_empty_notify():
    assert build_sinks({}) == []


# ── fan_out — 전송 격리 (FR-005) ────────────────────────────────────────────

class _OkSink(Sink):
    name = "ok"
    def __init__(self): self.called = False
    def send(self, item):
        self.called = True
        return "ref-ok"


class _BoomSink(Sink):
    name = "boom"
    def send(self, item):
        raise RuntimeError("전송 실패")


def test_fan_out_isolates_failure():
    """한 통로 실패해도 다른 통로 계속 — 결과에 실패는 None 으로 기록 (FR-005)."""
    boom, ok = _BoomSink(), _OkSink()
    out = fan_out(ITEM, [boom, ok])
    assert out["boom"] is None        # 실패 격리
    assert out["ok"] == "ref-ok"      # 나머지 정상
    assert ok.called is True          # 실패 뒤 통로도 호출됨


def test_fan_out_all_success():
    a, b = _OkSink(), _OkSink(); b.name = "ok2"
    out = fan_out(ITEM, [a, b])
    assert out == {"ok": "ref-ok", "ok2": "ref-ok"}


def test_fan_out_empty_sinks():
    assert fan_out(ITEM, []) == {}


def test_stdout_sink_send_returns_name(capsys):
    """StdoutSink.send 는 화면 출력 후 'stdout' 반환 (FR-004)."""
    assert StdoutSink().send(ITEM) == "stdout"
    assert "RELIEF trial" in capsys.readouterr().out
