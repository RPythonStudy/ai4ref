<!--
ai4ref 헌법 초안 (draft) — 2026-06-21
용도: Spec Kit 새출발 세션에서 `specify init` 후 `/speckit.constitution` 입력으로 그대로 투입.
방침(2026-06-21, Dr. Ben 결정): 포팅 예외를 두지 않고 전 코드를 명세-우선 신규 재작성.
  단, 검증 자산(validate_kq golden set, PICO OR-리스트, 게이트 기대판정)은 acceptance 참값으로 보존.
출처: cortex-kit constitution(V~XII·부가섹션) 흡수 + ai4ref HANDOFF 불변 결정 1~10.
-->

# ai4ref Constitution
(진료지침 기반 신규문헌 감시 — alert)

## Part A. 임상 감시 도메인 원칙 (불변 — 재논의 금지)

### I. 결정론 Backbone + 최소 LLM
검색·날짜필터·dedup 은 결정론(재현 가능)으로 처리한다. LLM 은 관련성·랜드마크
판정에만 쓴다. 기계가 할 수 있는 일을 LLM 에 위임하지 않는다.

### II. PICO = OR-블록
`pico.P`/`pico.I` 는 OR 리스트(검색어 항목들)다. 검색식 = `build_query(pico)`
= `(I 항목 OR …) AND (P 항목 OR …)`. `pico.C`/`pico.O` 는 검색식에서 제외하고
게이트 필터 기준으로만 쓴다.

### III. 정련 = OR-리스트 성장
정련은 검증 A 피드백으로 `pico.I`/`pico.P` 에 synonym 항목을 더해 리스트를
키우는 것이다(6/9→8/9 실측). **별도 query 필드를 두지 않는다** — 검색식은 파생물.

### IV. pt(설계)는 게이트 의미판단, 하드필터 금지
검색식·기계필터에 `[pt]` 를 넣지 않는다(넣으면 검증 A 3/9 붕괴·색인 전 신논문
누락, 실측). 설계 판단은 게이트가 *초록에서* 의미로 수행한다.

### V. 게이트 2단계
① 관련성 = 대상(P) AND 중재(I) 동시 매칭(한쪽만 = false → 토픽-인접 FP 차단)
② 랜드마크 = 지침과 다른 결론·practice-changing. `design_strictness`(strict|loose)
는 KQ 별 설정.

### VI. 두-검증 ≠ 게이트
`validate_kq`(검증 A·B)는 *검색 recall* 검증(아는 PMID 를 검색이 찾는가),
게이트는 *선별 precision*(노이즈에서 진짜만). 별개·상보 레이어로 유지한다.

### VII. Recall 편향 (loose 기본)
관련(P·I 일치)이면 경계선도 알린다. 놓침이 가장 비싼 오류다. 기본값 = loose.

### VIII. 모듈 Sink 멀티탭
출력 통로는 `notify/{base,registry,stdout,telegram,zotero}` 의 동일 플러그
모양(`Sink`)을 구현한다. 키 없으면 fail-soft skip. `features.yml` 로 토글한다.
새 통로 추가가 본체(검색·선별)를 건드리지 않는다.

### IX. 추가비용 0
게이트는 로컬 `claude_cli`(OAuth 구독, 모델 Haiku), cron 도 로컬(OpenClaw).
유료 API 호출을 도입하지 않는다.

### X. Config 단순화 (YAML)
운영 config = `key_questions.yml`(감시 KQ·PICO·두-검증) + `features.yml`(토글).
collector 폐기 후 `legacy_collections.yml`·`postgres.yml`·`logging.yml` 는 제거한다.
모든 설정은 YAML.

## Part B. 엔지니어링 규율 (Cortex-Kit 계승 · ai4ref 보정)

### XI. Collocated Spec (3계층 명세 — cortex-kit 차용)
명세는 3계층으로 둔다(cortex-kit 실증 구조):
- **L1 (프로젝트)** = `README.md` 의 기능 카탈로그(번호 기능 + 입력/출력/제약) +
  횡단 설계 `docs/*.md`. rpy-quarto-template fork README 는 **삭제가 아니라 L1 으로 재작성**.
- **L2 (피처)** = `specs/<NNN-feature>/`(Spec Kit 산출물: spec·plan·data-model·
  research·tasks·analysis). `/speckit.*` 명령으로 생성·관리.
- **L3 (컴포넌트)** = 코드 옆 동거 `<모듈>_spec.md`(예: `orchestrator.py` ↔
  `orchestrator_spec.md`). 헤더에 모듈·테스트·L2 참조(FR 번호), 공개 인터페이스 표.
  Spec Kit 이 자동생성하지 않는 **수작업 규율** — 코드와 동일 커밋에 동기화(원칙 XVI).

각 파일은 단일 관심사, 200줄 상한. `.specify/` = 헌법·템플릿(생성 영역). 예외 없음
— 기존 alert 로직도 L3 명세-우선으로 새로 작성한다.

### XII. Procedural Clarity (절차적 명확성)
절차적/함수형 파이프라인을 지향한다. 과도한 Class 상속·Decorator 추상화 금지.
유일한 허용 상속 = `Sink` 플러그 인터페이스(원칙 VIII). 1인 개발자가 흐름을
즉시 읽을 수 있어야 한다.

### XIII. Validation Triad (검증 삼위일체)
구현은 가드레일과 함께 작성한다. 가드레일 = ① `validate_kq`(검색 recall,
검증 A 8/9·B 2/2) ② 게이트 캘리브레이션(데모·백필 FP 0) ③ 스키마/키워드 기반
LLM 출력 검증. **검증 A/B 의 golden PMID, 정련된 PICO OR-리스트, 게이트 기대
판정(RELIEF🎯·FP 거부)은 acceptance 참값으로 보존**하며, 신규 재작성 코드는
이 참값을 재현해야 비로소 등가로 인정한다. LLM 출력은 정확 일치가 아닌 스키마·
키워드로 검증한다.

### XIV. E2E Validation & 이슈 환류
단위가 아닌 파이프라인 데이터 흐름(수집→게이트→sink)을 통합 검증한다(`--demo`).
운영 중 발견 이슈는 KQ 명세·PICO OR-리스트의 입력으로 환류시킨다(정련 = 원칙 III).

### XV. Self-Correction (자가 정정)
코드 수정 전 관련 명세·본 헌법과의 정합성을 사전 점검한다. 구현 세션과
보안/리뷰 세션을 분리한다.

### XVI. No Duplication (중복 금지)
명세에 "기존 구현 참조" 섹션을 둔다. 코드 변경 시 관련 명세를 동일 커밋에서
동기화한다(명세 드리프트 금지). 예: `zotero` 헬퍼는 `common/zotero.py` 단일
출처로 두어 중복 import 를 없앤다.

### XVII. Dependency First (의존성 우선)
단위 구현 전 의존 모듈의 입출력 계약을 확인한다. 미구현 의존성은 인터페이스
계약을 먼저 정의하고 Stub 으로 독립 구현한다.

### XVIII. Language Policy (언어 정책)
명세·주석·커밋·문서는 한글 서술 원칙. 업계 표준 전문용어·코드 식별자는 영어
원어를 그대로 쓰고, 한글 문장 내 영어 혼용을 표준으로 한다.

## Operational Standards (운영 표준)
- **런타임 데이터 격리**: 소스와 실행 데이터(`state/` seen-set jsonl)를 분리하고
  `.gitignore` 로 버전관리에서 제외한다(`.env` 시크릿 포함).
- **설정 불변성**: config 템플릿은 읽기 전용으로 두고 런타임은 값을 주입받는다.
- **환경 자동화**: 경로·의존성은 `pyproject.toml`(uv·src layout·slim deps:
  pyyaml·python-dotenv·requests·pyzotero)로 캡슐화한다. (Makefile 폐기)

## Development Workflow (개발 워크플로)
- **명세 우선**: Spec Kit 흐름으로 헌법(`/speckit.constitution`) → 피처 명세
  (`/speckit.specify`) → 구현. 신규 = KQ 추출 add-on 부터.
- **검증 참값 우선**: 코드는 명세-우선으로 신규 작성하되, 기존 검증 자산
  (validate_kq golden set, PICO OR-리스트, 게이트 기대판정)을 acceptance 참값으로
  먼저 추출·고정한 뒤 구현이 이를 통과하는지로 등가를 판정한다.
- **세션 분리**: 구현과 리뷰/보안을 한 세션에서 동시 수행하지 않는다.
- **커밋 동기화**: 코드 변경 시 관련 `_spec.md`·명세를 동일 커밋에 포함한다.

## Governance (거버넌스)
- 본 헌법은 모든 개발 관행·도구(Spec Kit 등) 사용보다 우선한다.
- Part A(도메인 불변)는 재논의 금지. 변경은 실측 근거 제시 시에만.
- 변경은 시맨틱 버저닝을 따르고, Sync Impact Report 와 템플릿 정합 검토를 거친다.

**Version**: 0.1.0 (draft) | **Ratified**: TODO | **Last Amended**: 2026-06-21
