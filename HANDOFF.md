# ai4ref alert — 세션 핸드오프 & 설계 헌법 (2026-06-21)

> 이 문서 = ① 지금까지 한 일 + ② 다시 흔들면 안 되는 설계 결정(헌법) + ③ 다음 할 일.
> 새 Claude 세션은 이걸 먼저 읽고 이어가면 됨.
> ⚠️ repo 루트 `README.md` 는 `rpy-quarto-template` 잔재 — 무시할 것 (실제 정체와 불일치).

---

## TL;DR — 현재 상태

**alert(진료지침 기반 신규문헌 감시) MVP 가 end-to-end 작동·검증됨.**
검색 → 게이트(claude_cli) → sink(telegram·zotero). 최근 30일 96건 백필 깨끗(FP 0, 진짜 관련 1건 알림).
**남은 것 = OpenClaw cron(자율 일일 발화) 하나.**

## 빠른 시작

```bash
.venv/bin/python src/alert/run.py --demo          # 배관 검증(가짜 랜드마크)
.venv/bin/python src/alert/run.py --collect-only  # 실 PubMed 후보만(전송 안 함)
.venv/bin/python src/alert/run.py                 # 전체 실행 (reldate 30, 게이트 라이브 → 폰)
.venv/bin/python src/alert/validate_kq.py         # 두-검증 자동채점 (검증 A 8/9 · B 2/2)
```
- config: `config/key_questions.yml`(감시 KQ) · `config/features.yml`(기능 토글) · `.env`(시크릿)
- 게이트 = `claude_cli`(로컬 Claude Code, OAuth 구독 → **추가비용 0**), 모델 Haiku.

## 아키텍처

```
PubMed ──검색·efetch──▶ ai4ref ──판정요청──▶ Claude Code (게이트, OAuth·Haiku)
                          │
              ┌───────────┼───────────┐
        Telegram(알림)  Zotero(보관)  stdout      ← notify sink 멀티탭(features.yml 토글)
                          ▲
                     OpenClaw cron (자율 발화 · 미구현)
```
파이프라인: 수집(esearch edat 증분) → ① 기계필터(날짜·중복·언어) → efetch → ② 게이트(관련·랜드마크) → dedup(seen-set jsonl) → sink fan-out.

---

## 불변 설계 결정 (헌법 — 재논의 금지, 근거 포함)

1. **결정론 backbone + 최소 LLM** — 검색·날짜필터·dedup = 결정론(재현) / 관련·랜드마크 판정만 LLM.
2. **PICO = OR-블록** — `pico.P`/`pico.I` = OR 리스트(검색어 항목들). 검색식 = `build_query(pico)` = `(I 항목 OR …) AND (P 항목 OR …)`. `pico.C`/`pico.O` = 게이트 필터 기준(검색식 제외).
3. **정련 = OR-리스트 성장** — 검증 A 피드백으로 `pico.I`/`pico.P` 에 synonym 항목 추가 → 6/9→8/9. **별도 query 필드 없음**(검색식은 파생). 정련 = 리스트가 자라는 것.
4. **pt(설계)는 게이트 의미판단, 하드필터 ✗** — 검색식에 `[pt]` 넣으면 검증 A 붕괴(근거 9편 중 RCT 3편 → **3/9**, 실측) + 색인 전 신논문 누락. 게이트가 *초록에서* 설계를 의미로 판단.
5. **게이트 2단계** — ① 관련성 = 대상(P) **AND** 중재(I) 둘 다 매칭(한쪽만 = false; **토픽-인접 FP 차단**) ② 랜드마크 = 지침과 다른 결론·practice-changing. `design_strictness`(strict|loose) per-KQ.
6. **두-검증 ≠ 게이트** — `validate_kq`(검증 A 8/9 · B 2/2)는 *검색 recall*(아는 PMID 를 검색이 찾나). 게이트는 *선별 precision*(노이즈서 진짜만). **별개 레이어**, 상보적.
7. **recall 편향 (loose 기본)** — 관련(P·I 일치)이면 경계선도 알림. 놓침이 가장 비싼 오류. (Dr. Ben 선호 = loose)
8. **모듈 sink 멀티탭** — `notify/{base,registry,stdout,telegram,zotero}`. 키 없으면 fail-soft skip. `features.yml` 토글.
9. **추가비용 0** — 게이트 = 로컬 `claude_cli`(OAuth 구독). cron 도 로컬(OpenClaw).
10. **config 2분할** — `key_questions.yml`(임상 감시 KQ, PICO·두-검증) / `legacy_collections.yml`(일반 주제 대량수집, term-based, 아카이브·collector 전용). 모든 config YAML.

---

## 한 일 (커밋)

```
8465c22 두-검증 구조 (검증 A 8/9 · B 2/2)
43f0436 랜드마크 감시 시스템 (sink 멀티탭·claude_cli 게이트·정본 통합)
7c6af09 search_collections.json → YAML
af5a576 kq_anchors 리네임
f274130 PICO operative          ← origin 이 여기까지 push 됨
88f75d7 정련 복원 시도            ─┐
8fc6ba0 P·I OR-블록 모델(최종 PICO) │ 미push
a6dd5c3 key_questions.yml 리네임    │ (origin 보다 앞섬)
8b22959 게이트 튜닝(P·I 매칭+strictness)─┘
```
→ **`git push origin main` 필요** (직접 push 권장: `git -C ~/projects/ai4ref push origin main`).

## 다음 (우선순위)

1. **OpenClaw cron** (자율 일일 발화) — MVP 마지막 조각.
   - 통합 결정 1개: OpenClaw = **Docker 컨테이너** / ai4ref·`claude` = **호스트** → 컨테이너 cron 이 호스트 작업을 어떻게 발화할지 (호스트 명령 위임 vs ai4ref 도 컨테이너 마운트).
   - 운영 시 `reldate=1`(또는 지난 실행 이후) → 하루 3~4건만 → ~1분.
2. **KQ 추출 add-on** — Spec Kit(SDD)으로. 지침 → PICO 추출 → 검증 A 로 OR-블록 정련 → `enabled:false` 토글 레코드. (Spec Kit 은 이 add-on 부터 도입 합의)
3. **정리(audit)** — 템플릿 잔재 청소, `daily_collection.sh` 미완성, 두 Zotero 모듈 응답 파싱 통일(`successful` vs `success`), 하드코딩 DB 접속 등.

## 함정·교훈 (실측 근거)

- **pt 하드필터 = 검증 A 3/9** (실측). 절대 검색식·기계필터에 넣지 말 것.
- **게이트 FP 유형 = 토픽 인접**(P만 또는 I만 매칭): 신장이식 goal-directed(P 불일치), 로봇 vs 복강경 직장절제(I 불일치). → 게이트 ① 의 P·I 동시 매칭으로 차단됨.
- claude 호출: `claude -p "<프롬프트>" --output-format json --model haiku`. 출력 = `{"type":"result","result":"<우리 JSON 문자열>",...}` 래퍼 → `llm_gate._extract_json` 가 처리.
- backfill(96건) ~20분(순차 게이트). 운영 cron(reldate 1)은 ~1분.
- `.env`(텔레그램·Zotero 키)·`state/`(seen-set jsonl) = gitignore. seen-set 현재 3건(RELIEF·OPTIMISE·GDFT메타).
- 발표자료(설계 권위) = vault `knowledge/02_areas/RPythonStudy/2026-06-22_2차_논문수집전략_slidable.md` (슬라이드 6=검색식 도출, 7=선별 필터, 8=두-검증).
