# src/common/query.py
"""PICO → PubMed 검색식 파생 (결정론, stdlib only).

검색식 = `(I 블록) AND (P 블록)`. P·I 검색어를 ` OR ` 로 묶어 각 블록을 만든다.
P·I 만 검색식에 사용(넓은 그물). C·O 는 검색식 제외 — 게이트가 의미 판정(헌법 II).
정련 = P·I OR 리스트에 synonym 항목을 **추가**(리스트 성장 = recall 향상, 헌법 III).
별도 query 필드를 저장하지 않고 매 호출 파생한다(단일 출처).

계약: specs/001-kq-pico-authoring/contracts/build_query.md (C1~C6).
"""
from __future__ import annotations

from typing import Union

ORBlock = Union[list, tuple, str]


def _block(items: ORBlock) -> str:
    """OR 리스트(또는 문자열) → ` OR ` 로 결합한 블록 문자열.

    빈 항목은 제외. 항목 문자열은 strip 외 변형 없이 그대로 사용(C5).
    """
    if isinstance(items, (list, tuple)):
        return " OR ".join(str(t).strip() for t in items if str(t).strip())
    return " ".join((items or "").split())


def build_query(pico: dict) -> str:
    """pico{P, I} → `(I) AND (P)` PubMed 검색식 (C·O 무시, 결정론).

    Args:
        pico: `{P: list[str], I: list[str], C?: str, O?: str}`. P·I 는 OR 리스트.
    Returns:
        `(<I 블록>) AND (<P 블록>)`. 한쪽만 있으면 그 블록만 괄호로 반환.
    """
    P = _block(pico.get("P"))
    I = _block(pico.get("I"))
    if I and P:
        return f"({I}) AND ({P})"
    return f"({I or P})"
