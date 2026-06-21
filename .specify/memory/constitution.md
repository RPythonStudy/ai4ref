<!--
Sync Impact Report
- Version change: 1.0.0 → 1.1.0 (2026-06-22 개정: 수집 트랙 반영)
- Ratified: 2026-06-21
- Modified principles: X. Config 단순화 — 감시+수집 2트랙(L1 §009~011) 반영: collector
  폐기 → DB-free 부활, 수집 주제 config 추가(통합/분리 L2 확정)
- Added sections:
  - Core Principles I~X (Part A 임상 감시 도메인 — 불변)
  - Core Principles XI~XVIII (Part B 엔지니어링 규율 — cortex-kit 계승)
  - Operational Standards
  - Development Workflow
  - Governance
- 권위 원본: docs/constitution-draft.md (cortex-kit 헌법 V~XII·부가 흡수 + HANDOFF 불변 결정 1~10)
- Templates requiring updates:
  - .specify/templates/plan-template.md ✅ — "Constitution Check" 게이트가 헌법을 동적 참조(범용), 하드코딩 없음
  - .specify/templates/spec-template.md ✅ — 원칙 하드코딩 없음, 충돌 없음
  - .specify/templates/tasks-template.md ✅ — 원칙 하드코딩 없음, 충돌 없음
- Follow-up TODOs: 없음
-->

# ai4ref Constitution
<!-- 진료지침 기반 신규문헌 감시(alert) -->

## Core Principles

### I. 결정론 Backbone + 최소 LLM
검색·날짜필터·dedup 은 결정론(재현 가능) 코드로 처리한다. LLM 은 관련성·랜드마크
판정에만 쓴다. 기계가 결정론으로 수행 가능한 작업을 LLM 에 위임해서는 안 된다.
**근거**: 재현성·디버깅 가능성 확보, LLM 비결정성의 영향 범위 최소화.

### II. PICO = OR-블록
`pico.P`/`pico.I` 는 OR 리스트(검색어 항목들)로 정의한다. 검색식은
`build_query(pico)` = `(I 항목 OR …) AND (P 항목 OR …)` 로 파생한다. `pico.C`/`pico.O`
는 검색식에서 제외하고 게이트 필터 기준으로만 쓴다.
**근거**: recall 우선 검색식 + 의미 필터 분리.

### III. 정련 = OR-리스트 성장
정련은 검증 A 피드백으로 `pico.I`/`pico.P` 에 synonym 항목을 추가해 리스트를
키우는 것이다(6/9→8/9 실측). 별도 `query` 필드를 두어서는 안 된다 — 검색식은 파생물.
**근거**: 단일 출처(PICO 리스트) 유지, 검색식 드리프트 방지.

### IV. pt(설계)는 게이트 의미판단, 하드필터 금지
검색식·기계필터에 `[pt]` 를 넣어서는 안 된다(넣으면 검증 A 3/9 붕괴 + 색인 전
신논문 누락, 실측). 설계 판단은 게이트가 초록에서 의미로 수행한다.
**근거**: 색인 지연·MeSH 누락이 recall 을 직접 파괴하는 실측 결과.

### V. 게이트 2단계
게이트는 ① 관련성 = 대상(P) AND 중재(I) 동시 매칭(한쪽만 매칭 = false → 토픽-인접
FP 차단) ② 랜드마크 = 지침과 다른 결론·practice-changing 의 2단계로 판정한다.
`design_strictness`(strict|loose)는 KQ 별 설정한다.
**근거**: 토픽-인접 false positive 차단(신장이식 P 불일치·로봇 직장절제 I 불일치 실측).

### VI. 두-검증 ≠ 게이트
`validate_kq`(검증 A·B)는 검색 recall 검증(아는 PMID 를 검색이 찾는가), 게이트는
선별 precision(노이즈에서 진짜만)이다. 두 레이어는 별개·상보로 유지하며 한쪽을
다른 쪽으로 대체해서는 안 된다.
**근거**: recall 레이어와 precision 레이어의 책임 분리.

### VII. Recall 편향 (loose 기본)
관련(P·I 일치)이면 경계선 사례도 알린다. 기본값은 loose 다. 놓침(false negative)이
가장 비싼 오류이기 때문이다.
**근거**: 임상 감시에서 누락 비용 > 과알림 비용 (Dr. Ben 선호).

### VIII. 모듈 Sink 멀티탭
출력 통로는 `notify/{base,registry,stdout,telegram,zotero}` 의 동일 플러그 모양
(`Sink`)을 구현한다. 키 부재 시 fail-soft skip 하고 `features.yml` 로 토글한다. 새
통로 추가가 본체(검색·선별)를 건드려서는 안 된다.
**근거**: 출력 확장과 핵심 파이프라인의 결합도 0.

### IX. 추가비용 0
게이트는 로컬 `claude_cli`(OAuth 구독, 모델 Haiku)로, cron 도 로컬(OpenClaw)로
운영한다. 유료 API 호출을 도입해서는 안 된다.
**근거**: 운영 지속가능성 — 종량 과금 의존 제거.

### X. Config 단순화 (YAML)
운영 config 는 감시(`key_questions.yml`: KQ·PICO·두-검증) + 수집(term-based 수집 주제,
L1 §009) + `features.yml`(토글)로 한정한다. 수집 주제를 `key_questions.yml` 에 `type`
필드로 통합할지 별도 config 로 둘지는 L2 에서 확정한다. `postgres.yml`·`logging.yml` 는
제거한다(DB-free·stderr 폴백). collector 는 폐기가 아니라 수집 트랙 backbone 으로
DB-free 재작성한다. 모든 설정은 YAML 로 작성한다.
**근거**: 설정 표면적 최소화·단일 형식. 감시·수집 2트랙(L1) 반영.

### XI. Collocated Spec (3계층 명세 — cortex-kit 차용)
명세는 3계층으로 둔다:
- **L1 (프로젝트)** = `README.md` 의 기능 카탈로그(번호 기능 + 입력/출력/제약) +
  횡단 설계 `docs/*.md`. rpy-quarto-template fork README 는 삭제가 아니라 L1 으로
  재작성한다.
- **L2 (피처)** = `specs/<NNN-feature>/`(Spec Kit 산출물: spec·plan·data-model·
  research·tasks·analysis). `/speckit-*` 명령으로 생성·관리한다.
- **L3 (컴포넌트)** = 코드 옆 동거 `<모듈>_spec.md`(예: `llm_gate.py` ↔
  `llm_gate_spec.md`). 헤더에 모듈·테스트·L2 참조(FR 번호)와 공개 인터페이스 표를
  포함한다. Spec Kit 이 자동생성하지 않는 수작업 규율이며, 코드와 동일 커밋에
  동기화해야 한다(원칙 XVI).

각 파일은 단일 관심사만 담고 200줄을 초과해서는 안 된다. 예외 없음 — 기존 alert
로직도 L3 명세-우선으로 새로 작성한다.
**근거**: 명세-코드 동거로 드리프트 차단, 1인 개발 인지부하 관리.

### XII. Procedural Clarity (절차적 명확성)
절차적/함수형 파이프라인을 지향한다. 과도한 Class 상속·Decorator 추상화를 금지한다.
유일하게 허용하는 상속은 `Sink` 플러그 인터페이스(원칙 VIII)다. 1인 개발자가 흐름을
즉시 읽을 수 있어야 한다.
**근거**: 추상화 비용 > 이득 (1인 유지보수).

### XIII. Validation Triad (검증 삼위일체)
구현은 가드레일과 함께 작성한다. 가드레일 = ① `validate_kq`(검색 recall, 검증 A
8/9·B 2/2) ② 게이트 캘리브레이션(데모·백필 FP 0) ③ 스키마/키워드 기반 LLM 출력
검증. 검증 A/B 의 golden PMID, 정련된 PICO OR-리스트, 게이트 기대판정(RELIEF🎯·FP
거부)은 acceptance 참값으로 보존하며, 신규(재작성) 코드는 이 참값을 재현해야 등가로
인정한다. LLM 출력은 정확 일치가 아닌 스키마·키워드로 검증한다.
**근거**: 재작성 위험(검증·캘리브레이션 소실)을 참값 박제로 회피.

### XIV. E2E Validation & 이슈 환류
단위가 아닌 파이프라인 데이터 흐름(수집→게이트→sink)을 통합 검증한다(`--demo`).
운영 중 발견한 이슈는 KQ 명세·PICO OR-리스트의 입력으로 환류시킨다(정련 = 원칙 III).
**근거**: 결함의 명세 환류로 선순환 유지.

### XV. Self-Correction (자가 정정)
코드 수정 전 관련 명세(L3)·본 헌법과의 정합성을 사전 점검한다. 구현 세션과
보안/리뷰 세션을 분리한다.
**근거**: 명세 이탈·리뷰 편향 최소화.

### XVI. No Duplication (중복 금지)
명세에 "기존 구현 참조" 섹션을 둔다. 코드 변경 시 관련 명세를 동일 커밋에서
동기화한다(명세 드리프트 금지). 예: `zotero` 헬퍼는 `common/zotero.py` 단일
출처로 두어 중복 import 를 없앤다.
**근거**: 단일 출처 원칙, AI 중복 생성 차단.

### XVII. Dependency First (의존성 우선)
단위 구현 전 의존 모듈의 입출력 계약을 확인한다. 미구현 의존성은 인터페이스 계약을
먼저 정의하고 Stub 으로 독립 구현한다.
**근거**: 통합 단계 충돌·재작업 방지.

### XVIII. Language Policy (언어 정책)
명세·주석·커밋·문서는 한글 서술을 원칙으로 한다. 업계 표준 전문용어·코드 식별자는
영어 원어를 그대로 쓰고, 한글 문장 내 영어 혼용을 표준으로 한다.
**근거**: 1인 개발자(Dr. Ben)의 가독성·전문성 동시 확보.

## Operational Standards

- **런타임 데이터 격리**: 소스와 실행 데이터(`state/` seen-set jsonl)를 분리하고
  `.gitignore` 로 버전관리에서 제외한다(`.env` 시크릿 포함).
- **설정 불변성**: config 템플릿은 읽기 전용으로 두고 런타임은 값을 주입받는다.
- **환경 자동화**: 경로·의존성은 `pyproject.toml`(uv·src layout·slim deps:
  pyyaml·python-dotenv·requests·pyzotero)로 캡슐화한다. Makefile 은 폐기한다.

## Development Workflow

- **명세 우선**: Spec Kit 흐름으로 헌법(`/speckit-constitution`) → 피처 명세
  (`/speckit-specify`) → 계획(`/speckit-plan`) → 태스크(`/speckit-tasks`) → 구현
  (`/speckit-implement`) 순으로 진행한다. 신규 작업은 KQ 추출 add-on 부터 시작한다.
- **검증 참값 우선**: 코드는 명세-우선으로 신규 작성하되, 기존 검증 자산
  (validate_kq golden set, PICO OR-리스트, 게이트 기대판정)을 acceptance 참값으로
  먼저 추출·고정한 뒤 구현이 이를 통과하는지로 등가를 판정한다.
- **세션 분리**: 구현과 리뷰/보안을 한 세션에서 동시 수행하지 않는다.
- **커밋 동기화**: 코드 변경 시 관련 `<모듈>_spec.md`·명세를 동일 커밋에 포함한다.

## Governance

- 본 헌법은 모든 개발 관행·도구(Spec Kit 등) 사용보다 우선한다.
- Part A(원칙 I~X, 임상 감시 도메인 불변)는 재논의를 금지한다. 변경은 실측 근거를
  제시할 때에만 허용한다.
- 모든 PR·코드 리뷰·AI 에이전트 작업은 본 원칙 준수 여부를 우선 검증한다.
- 변경은 시맨틱 버저닝을 따른다: MAJOR=원칙 제거/비호환 재정의, MINOR=원칙·섹션
  추가/실질 확장, PATCH=문구·오타·비의미 정정.
- 헌법 수정 절차: 수정 제안 → 영향 분석(Sync Impact Report) → 의존 템플릿 갱신 →
  버전 범프 → 비준.

**Version**: 1.1.0 | **Ratified**: 2026-06-21 | **Last Amended**: 2026-06-22
