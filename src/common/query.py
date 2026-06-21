# src/common/query.py
"""PICO → PubMed 검색식 생성 (단일 전략).

P·I 만 검색식에 사용 → 넓은 그물(결정론). C·O 는 검색식에서 제외, 게이트가 필터(추론).
검색식 = (I) AND (P). PICO 가 단일 진실, 검색식은 여기서 파생 (손으로 안 짠다).
"""


def build_query(pico: dict) -> str:
    """pico{P,I} → "(I) AND (P)". 공백·개행 정규화."""
    P = " ".join((pico.get("P") or "").split())
    I = " ".join((pico.get("I") or "").split())
    if I and P:
        return f"({I}) AND ({P})"
    return f"({I or P})"
