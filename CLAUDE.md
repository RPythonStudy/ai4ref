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

<!-- SPECKIT START -->
For additional context about technologies to be used, project structure,
shell commands, and other important information, read the current plan
<!-- SPECKIT END -->
