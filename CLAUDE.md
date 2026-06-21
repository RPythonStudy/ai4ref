# ai4ref — Claude Code 작업 지침

ai4ref = **진료지침 기반 신규문헌 감시(alert)** + 문헌 대량수집(collector, legacy).
현재 집중 = **alert MVP** (감시 파이프라인).

> **세션 시작 시 `HANDOFF.md` 를 먼저 읽을 것** — 현재 상태 · 불변 설계 결정(헌법) · 다음 할 일.
> ⚠️ repo 루트 `README.md` 는 `rpy-quarto-template` 잔재 — 무시.

## 핵심 불변 결정 (재논의 금지 · 상세 = HANDOFF.md)

- 결정론 backbone + 최소 LLM
- **PICO = OR-블록**: `pico.P`/`pico.I` = OR 리스트 → 검색식 = `build_query` = `(I) AND (P)`. C·O = 게이트 필터.
- **정련 = OR-리스트 성장** (검증 A 피드백, 6/9→8/9). 별도 query 필드 없음(파생).
- **pt = 게이트 의미판단, 하드필터 ✗** (검색식에 넣으면 검증 A 3/9 붕괴).
- **게이트 2단계**: ① 대상(P) AND 중재(I) 동시 매칭 ② 랜드마크 판정. `design_strictness`(strict|loose, 기본 loose=recall).
- **두-검증(검색 recall) ≠ 게이트(선별 precision)** — 별개 레이어.
- **추가비용 0**: 게이트 = 로컬 `claude_cli`(OAuth).
- config: `key_questions.yml`(감시 KQ) · `legacy_collections.yml`(수집, 아카이브). 모두 YAML.

## 실행

```bash
.venv/bin/python src/alert/run.py [--demo | --collect-only]   # 감시 파이프라인
.venv/bin/python src/alert/validate_kq.py                     # 두-검증 8/9 · 2/2
```

## 다음

OpenClaw cron(자율 발화) → KQ 추출 add-on(Spec Kit) → audit 정리. (HANDOFF.md §다음)
