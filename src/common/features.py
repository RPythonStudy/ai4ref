# src/common/features.py
"""기능 플래그 로더 + fail-soft 로거.

config/features.yml 을 읽어 '원하는 기능만' 운영하게 한다.
로거는 프로젝트 logger 가 정상이면 그것을, 아니면 stderr 로 폴백 — 로그 디렉토리
미설정 상태에서도 코어가 죽지 않도록(클론 사용자 fail-soft).
"""
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import yaml
from dotenv import load_dotenv

load_dotenv(override=False)

_ROOT = Path(__file__).resolve().parents[2]          # ai4ref/
DEFAULT_FEATURES = _ROOT / "config" / "features.yml"


def project_root() -> Path:
    return _ROOT


def load_features(path=None) -> dict:
    p = Path(path) if path else DEFAULT_FEATURES
    if not p.exists():
        log.warning(f"[features] {p} 없음 — 모든 기능 비활성으로 간주")
        return {}
    with open(p, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def is_enabled(features: dict, *keys) -> bool:
    """features['a']['b']...['enabled'] 를 안전하게 조회."""
    node = features or {}
    for k in keys:
        if not isinstance(node, dict):
            return False
        node = node.get(k, {})
    if isinstance(node, dict):
        return bool(node.get("enabled", False))
    return bool(node)


# --- fail-soft 로거 -----------------------------------------------------------
def _mk_softlog():
    def make(level):
        def _log(msg):
            try:
                from common.logger import get_logger
                getattr(get_logger(), level)(msg)
            except Exception:
                print(f"[{level.upper()}] {msg}", file=sys.stderr)
        return _log
    return SimpleNamespace(
        debug=make("debug"), info=make("info"),
        warning=make("warning"), error=make("error"),
    )


log = _mk_softlog()
