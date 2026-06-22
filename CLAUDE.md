# ai4ref — Claude Code 작업 지침

ai4ref = **진료지침 기반 신규문헌 감시(alert)** + **주제 기반 대량수집(collect)**. 2트랙.
현재 집중 = **alert MVP** (감시 파이프라인).

> **세션 시작 시 읽을 것**: `.specify/memory/system-spec.md`(L1 기능 명세) ·
> `.specify/memory/constitution.md`(설계 헌법). 사용자용 소개·실행 = `README.md`.
> 다음 할 일·실측 교훈은 메모리(MEMORY.md)에 자동 로드됨.

## 핵심 불변 결정 (재논의 금지 · 상세 = .specify/memory/constitution.md)

- 결정론 backbone + 최소 LLM
- **PICO = OR-블록**: `pico.P`/`pico.I` = OR 리스트 → 검색식 = `build_query` = `(I) AND (P)`. C·O = 게이트 필터.
- **정련 = OR-리스트 성장** (검증 A 피드백, 6/9→8/9). 별도 query 필드 없음(파생).
- **pt = 게이트 의미판단, 하드필터 ✗** (검색식에 넣으면 검증 A 3/9 붕괴).
- **게이트 2단계**: ① 대상(P) AND 중재(I) 동시 매칭 ② 랜드마크 판정. `design_strictness`(strict|loose, 기본 loose=recall).
- **두-검증(검색 recall) ≠ 게이트(선별 precision)** — 별개 레이어.
- **추가비용 0**: 게이트 = 로컬 `claude` CLI(OAuth·Haiku). 호출 형식 = 메모리 [[claude-gate-cli]].
- config: `key_questions.yml`(감시 KQ) · `features.yml`(토글) · 수집 주제(L2 확정). 모두 YAML.

## 실행 (uv)

```bash
uv run python -m alert.run [--demo | --collect-only]   # 감시 파이프라인
uv run python -m alert.validate_kq                     # 두-검증 8/9 · 2/2
```

## 다음

메모리 [[project-roadmap]] 참조: OpenClaw cron → KQ 추출 add-on → 수집 트랙 → audit.

## Workflow — Issue-First

코드 변경(기능·버그·리팩터링) 커밋 전에 **관련 GitHub 이슈가 존재**해야 하며 커밋에 이슈 번호를 연결한다.

1. **검색** — `gh issue list --state open --search "<키워드>"` 로 기존/유사 이슈 확인.
2. **없으면 생성** — `gh issue create` 로 본문(현상/원인/수정 범위/검증)과 라벨(`enhancement`·`bug`, 필요시 `feature:<NNN-slug>`)을 붙여 먼저 만든다.
3. **연결** — 커밋 메시지 끝에 `Closes #NN`(단독 해결) 또는 `refs #NN`(부분). main 직접 커밋이므로 `Closes` 가 즉시 발동.
4. **예외** — 문서/주석/포매터/설정·메타 작업은 이슈 없이 진행 가능(애매하면 만드는 쪽).

자동화: `scripts/auto-implement.sh <tasks.md>` 가 tasks.md 미완료 묶음을 위 원칙대로 구현→이슈 연결·종료까지 무인 반복(`--dangerously-skip-permissions`). 첫 실행은 지켜볼 것.

<!-- SPECKIT START -->
현재 작업 계획: `specs/007-notification-fanout-sinks/plan.md` (§007 Sink 멀티탭 retro-spec —
랜드마크 fan-out 출력 레이어). 완료: 001·003·004·005·006. 기술 컨텍스트·구조·등가기준은 해당 plan 참조.
<!-- SPECKIT END -->
