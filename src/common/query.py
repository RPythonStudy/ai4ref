# src/common/query.py
"""PICO → PubMed 검색식 생성.

핵심 로직: P·I 검색어를 OR 로 묶어 블록 → 검색식 = (I) AND (P).
P·I 만 검색식에 사용(넓은 그물, 결정론). C·O 는 검색식 제외, 게이트가 필터(추론).
정련 = P·I 리스트에 synonym 항목을 추가(확장)하는 것 (검증 A 로 6/9→8/9).
"""


def _block(x) -> str:
    """OR 리스트(또는 문자열) → OR 로 join 한 블록."""
    if isinstance(x, (list, tuple)):
        return " OR ".join(str(t).strip() for t in x if str(t).strip())
    return " ".join((x or "").split())


def build_query(pico: dict) -> str:
    """pico{P,I} → "(I) AND (P)". P·I 는 OR 리스트(또는 문자열)."""
    P = _block(pico.get("P"))
    I = _block(pico.get("I"))
    if I and P:
        return f"({I}) AND ({P})"
    return f"({I or P})"
