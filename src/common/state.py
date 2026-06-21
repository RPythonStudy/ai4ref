# src/common/state.py
"""보낸 랜드마크 추적(dedup) — seen-set.

jsonl 백엔드: 한 줄 = 한 전달 기록. 무료·버전관리·클라우드(GitHub Actions 커밋) 친화.
'무엇을·언제·왜 보냈는지'가 영구 로그로 남고, 권위 있는 dedup 기준이 된다.
"""
import json
from pathlib import Path
from datetime import datetime, timezone

from common.features import log, project_root


class SeenStore:
    def __init__(self, path: str):
        p = Path(path)
        self.path = p if p.is_absolute() else (project_root() / p)
        self._seen = self._load()
        log.debug(f"[state] seen-set 로드: {len(self._seen)}건 ({self.path})")

    def _load(self) -> set:
        seen = set()
        if self.path.exists():
            for line in self.path.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    seen.add(json.loads(line)["key"])
                except Exception:
                    continue
        return seen

    def is_seen(self, key: str) -> bool:
        return key in self._seen

    def mark(self, key: str, record: dict | None = None) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        rec = {"key": key, "ts": datetime.now(timezone.utc).isoformat()}
        if record:
            rec.update(record)
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        self._seen.add(key)

    def __len__(self) -> int:
        return len(self._seen)
