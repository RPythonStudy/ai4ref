# Research: KQ·PICO 정의 (Phase 0)

retro-spec 이라 미지수가 거의 없다. 기존 동작 + 헌법에서 결정을 확정한다.

## D1. KQ 레코드 저장 형식

- **Decision**: 단일 YAML 파일 `config/key_questions.yml` 의 `kqs:` 리스트.
- **Rationale**: 사람이 읽고 편집·diff 가능, DB-free(헌법 X), 기존 동작과 동일.
- **Alternatives**: JSON(주석 불가), DB(헌법 X 위배·과함), KQ 당 파일 분리(현 규모 1건엔 과함).

## D2. 검색식 파생 = build_query(pico)

- **Decision**: `(I 항목 OR …) AND (P 항목 OR …)`. 항목은 yaml 의 P·I 리스트 원소를
  그대로 OR 결합, 블록을 AND 결합. C·O 는 제외.
- **Rationale**: 헌법 II·III. acceptance `derived_query` 스냅샷과 일치해야 함.
- **Alternatives**: 정련된 검색식을 별도 저장(III 위배 — 드리프트), C·O 를 검색식 포함
  (II 위배 — 색인 누락으로 recall 붕괴, 검증 A 3/9 실측).

## D3. pt(연구 설계) 처리

- **Decision**: 검색식·기계필터에 `[pt]` 넣지 않음. 게이트(006)의 의미판단으로 위임.
- **Rationale**: 헌법 IV·I. 색인 지연·MeSH 누락이 recall 직접 파괴(검증 A 3/9 실측).
- **Alternatives**: `[pt]` 하드필터 → 거부(실측 붕괴).

## D4. KQ 유효성 검증 범위

- **Decision**: 로드 시 결정론 검증 — P·I 비어있지 않음(검색식 파생 가능), `query` 류
  검색식 필드 부재(III), `design_strictness` 미지정 시 기본 `loose`, `enabled` 기본 처리.
- **Rationale**: FR-002·003·007. 잘못된 KQ 를 조기 표면화(Edge Cases).
- **Alternatives**: 무검증 로드(런타임 오류 지연), 스키마 라이브러리(pydantic 등) 도입
  (현 규모엔 과함 — slim deps·XII).

## D5. 검색식 등가 검증 방법

- **Decision**: 단위테스트로 `build_query(pico)` 출력 == acceptance 의 `derived_query`
  스냅샷 비교(결정론). recall 실측(8/9·B 2/2)은 `validate_kq`(005, PubMed)로.
- **Rationale**: 헌법 XIII·XIV. 재작성 등가를 네트워크 없이 빠르게 확인 + 참값으로 최종.
- **Alternatives**: recall 만으로 검증(네트워크 의존·느림), 검증 생략(등가 보장 불가).

## D6. 수집 트랙(009) 주제의 동일 config 통합 여부

- **Decision**: 본 피처(001) 범위 밖. L2 §009 에서 `type` 필드 통합 vs 분리 확정(헌법 X).
- **Rationale**: 001 은 감시 KQ 스키마에 집중. 조기 통합은 YAGNI.

**결과**: NEEDS CLARIFICATION 없음. Phase 1 진행 가능.
