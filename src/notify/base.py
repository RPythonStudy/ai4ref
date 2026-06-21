# src/notify/base.py
"""Sink 인터페이스 + 전달 아이템.

모든 출력 통로(stdout·telegram·email·zotero…)는 Sink 를 구현한다 = '같은 플러그 모양'.
새 통로를 추가해도 본체(검색·선별)는 건드리지 않는다.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class LandmarkItem:
    kq: str                 # KQ 라벨/식별자
    pmid: str
    title: str
    journal: str = ""
    year: str = ""
    reason: str = ""        # 한 줄 근거(왜 랜드마크인가)
    url: str = ""
    abstract: str = ""      # 게이트 판정용(내부) — 렌더링엔 안 씀

    def __post_init__(self):
        if not self.url and self.pmid:
            self.url = f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

    @property
    def key(self) -> str:
        return f"{self.kq}:{self.pmid}"


class Sink(ABC):
    name = "base"

    def available(self) -> tuple[bool, str]:
        """prereq 점검 → (사용가능, 사유). 미충족이면 멀티탭이 fail-soft 로 skip."""
        return True, ""

    @abstractmethod
    def send(self, item: LandmarkItem):
        """전달 후 message 참조(또는 None) 반환."""
        ...
