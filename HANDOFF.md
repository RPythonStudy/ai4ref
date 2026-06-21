# ai4ref alert — 세션 핸드오프 & 설계 헌법 (2026-06-21)

> 이 문서 = ① 지금까지 한 일 + ② 다시 흔들면 안 되는 설계 결정(헌법) + ③ 다음 할 일.
> 새 Claude 세션은 이걸 먼저 읽고 이어가면 됨.
> ⚠️ repo 루트 `README.md` 는 `rpy-quarto-template` 잔재 — 무시할 것 (실제 정체와 불일치).

---

## TL;DR — 현재 상태

**alert(진료지침 기반 신규문헌 감시) MVP 가 end-to-end 작동·검증됨.**
검색 → 게이트(claude_cli) → sink(telegram·zotero). 최근 30일 96건 백필 깨끗(FP 0, 진짜 관련 1건 알림).
**남은 것 = OpenClaw cron(자율 일일 발화) 하나.**

## 🔄 새출발 계획 (Spec Kit, 2026-06-21 결정)

**결정 (2026-06-21 갱신)**: ai4ref 의 기존 코드는 대부분 cruft(템플릿 fork·collector·Quarto/R). **cruft 버리고 Spec Kit 으로 새출발**하며, alert MVP 도 *포팅(복사)이 아니라* **전부 명세-우선으로 신규 재작성**한다.
> ⚠️ 방침 전환: 이전엔 "MVP 포팅(복사)·재작성 금지"였으나, 프로젝트별 포팅 예외가 누적 인지부하를 키운다는 판단(Dr. Ben)으로 **단일 규율 = 신규 재작성**으로 통일. 재작성의 본래 위험(검증·캘리브레이션 소실)은 코드 보존이 아니라 **검증 자산을 acceptance 참값으로 박제**해 회피한다 → 헌법 XIII.
> **검증 참값(보존 필수, 재작성 코드가 재현해야 등가)**: `validate_kq` golden PMID(검증 A 8/9·B 2/2), 정련된 PICO OR-리스트(`key_questions.yml`), 게이트 기대판정(RELIEF🎯·토픽-인접 FP 거부).

**백업 완료**: git tag `archive/mvp-pre-speckit`(@039f36c) · `~/projects/ai4ref-archive-2026-06-21.tar.gz` · 원격 `RPythonStudy/ai4ref`.

**헌법 비준 완료**: `.specify/memory/constitution.md` (v1.2.0) — cortex-kit 헌법(V~XII·부가) 흡수 + 본 §불변 설계 결정 1~10. (초안 `docs/constitution-draft.md` 는 비준 후 삭제)

### 명세 3계층 (cortex-kit 차용)
- **L1** = `.specify/memory/system-spec.md`(기능 카탈로그·입력/출력/제약, agent 컨텍스트 정본) + `docs/*.md`. 루트 `README.md` 는 **사용자용**(소개·Quickstart)으로 분리 — agent 컨텍스트엔 README 대신 system-spec·constitution.
- **L2** = `specs/<NNN-feature>/`(Spec Kit `/speckit.*` 산출물).
- **L3** = 코드 옆 `<모듈>_spec.md`(예: `llm_gate.py` ↔ `llm_gate_spec.md`). Spec Kit 자동생성 아님 = 수작업, 코드와 동일 커밋 동기. (헌법 XI)

### REBUILD (명세-우선 신규 재작성 — 기존 코드는 *참조용*, 복사 금지)
아래는 재작성 **범위**(기능 단위)다. 기존 파일을 그대로 옮기지 말고, L3 `<모듈>_spec.md` 먼저 → 구현 → 검증 참값 통과 순으로 새로 짠다.
- config: `key_questions.yml`(PICO OR-블록·두-검증 = 검증 참값) · `features.yml`(토글) — *값은 참값으로 보존, 스키마는 재정의 가능*
- `notify/` (base·registry·stdout·telegram·zotero sink) — Sink 플러그 모양(헌법 VIII·XII 의 유일 허용 상속)
- **감시 트랙** `alert/` (run·llm_gate·validate_kq·classify_qtype) — 결정론 backbone + 게이트 2단계 (L1 §001~008)
- **수집 트랙** `collect/` [신규] — 주제 수집 정의·대량검색·PMC PDF/XML 다운로드·Zotero linked_file·2nd-brain 직접 참조 (L1 §009~011). **기존 `collector/` 의 esearch·efetch·pmc_downloader·zotero_add 를 postgres 제거(DB-free)하여 재활용**
- `common/` (features·pubmed·query·state·**zotero**) — 양 트랙 공유
- `.env`(시크릿)·`.gitignore`·이 문서들(HANDOFF·CLAUDE·README·.specify) = 그대로 이어씀

### 재작성 시 합칠 것 (기존의 부채를 신규에 반영)
- zotero 헬퍼: 기존엔 `zotero_sink` 가 `collector.zotero_add` 의 `fetch_meta`·`to_zotero_item` 를 import. 신규는 **`common/zotero.py` 단일 출처**로 작성(헌법 XVI No Duplication) → 감시(zotero sink)·수집(linked_file 배포) 양 트랙이 공유. `zotero_add.py` 는 이미 DB-free 라 로직 재활용 용이.
- 로거: `pyproject` slim deps 기조 → 표준 `logging` + stderr 폴백으로 단순화(logger.py·logging.yml 재작성 안 함).

### DISCARD (cruft 전량)
> ✅ **실행 완료 (2026-06-22)**: 아래 대부분 삭제됨(template·R·Quarto·빌드산출물 docs/site_libs 1.7M·wiki 서브모듈·Makefile·scripts·DB코드·logger.py·logging.yml·죽은 tests). **`src/collector/*` 는 보류**(수집 트랙 부활 예정). 로깅 = `features.py` stderr 폴백으로 유지. 복구 = `archive/mvp-pre-speckit` 태그·git 이력.
- 템플릿: `_quarto.yml`·`index.qmd`·`posts/`·`wiki/`·`styles/`·`templates/`·`references/`·`utterances.html`·`ai4ref.Rproj`·`renv.lock` (※ `README.md` 는 DISCARD 아님 → L1 명세로 재작성)
- R 트랙: `src/R/`·`src/Rlib/`·`src/preprocessor/`·`src/summerizer/`
- DB·R 인프라: `src/database/*`, `common/database.py`, `config/{postgres.yml, logging.yml}`
- ⚠️ `src/collector/*` 는 **전량 DISCARD 아님** — esearch·efetch·pmc_downloader·zotero_add 는 수집 트랙으로 **DB-free 재작성하여 부활**(위 REBUILD). postgres 결합 코드(get_db_connection·collection_pmid/papers 테이블 워크플로)만 폐기.
- `config/legacy_collections.yml`: term-based 수집 주제 = 수집 트랙(L1 §009)으로 **재정의**. `key_questions.yml` 통합(type 필드) vs 별도 config 는 L2 확정 (헌법 X).
- 오케스트레이션: `Makefile`·`Makefile_project`·`scripts/`
- 패키징: ✅ **완료(2026-06-22, lift-and-shift)** — `requirements.txt` 삭제 → `pyproject.toml`(uv·src 레이아웃·hatchling) + `uv.lock`. 현재 deps = 현행 유지(biopython·aiohttp·psycopg2-binary 포함). 실행 = `uv run python -m alert.run`. **slim 트림(psycopg2 제거·biopython→requests 등)은 DB-free 재작성 때.**

### Spec Kit 실행 순서 (fresh ai4ref 세션에서)
```
0. 검증 참값 추출·고정: validate_kq golden PMID, key_questions.yml 의 PICO OR-리스트,
   게이트 기대판정(RELIEF🎯·FP 거부)을 acceptance 참값 파일로 박제 (재작성 등가 판정 기준)
1. 백업 확인 (위) → cruft DISCARD + pyproject.toml 작성
2. uvx --from git+https://github.com/github/spec-kit.git@v0.11.3 specify init --here --integration claude --script sh --force
   ⚠️ v0.11.3 문법 = `--integration claude` (구 v0.4.1 의 `--ai claude` 아님). `--force` = 비어있지 않은 현재 디렉토리에 merge.
   ⚠️ 버전 고정 필수. 무고정(@없음)은 main HEAD(미릴리스)를 끌어옴. 최신 릴리스 = v0.11.3(2026-06-19).
   (cortex 는 v0.4.1 로 init 됨 → 템플릿 골격 다름. ai4ref 는 v0.11.3 기준으로 통일)
3. /speckit-constitution  — 헌법 비준 (✅ 완료: .specify/memory/constitution.md v1.2.0)
4. /speckit.specify — 현행 alert 감시 파이프라인 retro-spec(기존 동작 = 검증 참값을 명세로)
5. /speckit.specify — KQ 추출 add-on(지침→PICO→검증A 정련→토글) 부터 신규
6. 검증 참값 재현 확인(등가 판정): validate_kq → 8/9·2/2 / 게이트 → RELIEF 🎯·FP 거부
```
> **OpenSpec 은 ai4ref 엔 불필요**(유지할 브라운필드 없음). Dr. Ben 의 *살아있는 다른 프로젝트*(radsafety-pwa 등) 유지보수엔 OpenSpec. (조사 근거: vault `03_resources/AI-agents/2026-06-21_딥리서치_SDD-1인개발자-유지보수전략.md`)

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
[감시 트랙 · alert]
PubMed ──검색·efetch──▶ ai4ref ──판정요청──▶ Claude Code (게이트, OAuth·Haiku)
                          │
              ┌───────────┼───────────┐
        Telegram(알림)  Zotero(보관)  stdout      ← notify sink 멀티탭(features.yml 토글)
                          ▲
                     OpenClaw cron (자율 발화 · 미구현)

[수집 트랙 · collect]
자연어 주제 → config ──term 검색·efetch·PMC다운로드──▶ Gdrive Zotero_attachments (PDF·XML)
                                                          ├─▶ Zotero 서지+컬렉션+linked_file (클라이언트 Zotfile 1단계)
                                                          └─▶ 2nd-brain 직접 참조 (brainify 파싱·LLM, 사본 없음)
```
감시 파이프라인: 수집(esearch edat 증분) → ① 기계필터(날짜·중복·언어) → efetch → ② 게이트(관련·랜드마크) → dedup(seen-set jsonl) → sink fan-out.
수집 파이프라인: 주제(term) esearch → efetch(메타) → PMC PDF/XML 다운로드 → Zotero_attachments 저장(Zotfile 명명) → Zotero linked_file + 2nd-brain 참조.

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
10. **config** (헌법 X·v1.1.0) — `key_questions.yml`(감시 KQ, PICO·두-검증) + `features.yml`(토글) + **수집 주제 config**(term-based, 수집 트랙 L1 §009). 수집 주제를 `key_questions.yml` 에 type 필드로 통합할지 별도 config 로 둘지는 L2 확정. `postgres.yml`·`logging.yml` 제거. collector 는 폐기 아니라 DB-free 부활. 모든 config YAML.

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
3. **수집 트랙** [신규, L1 §009~011] — 자연어 주제 → 대량검색 → PMC PDF/XML → Gdrive `Zotero_attachments` → Zotero linked_file(Zotfile) + 2nd-brain 직접 참조. collector DB-free 재작성. L2 확정 항목: config 통합/분리(헌법 X), Zotfile 링크 반자동화, PMC OA 한계.
4. **정리(audit)** — 템플릿 잔재 청소, `daily_collection.sh` 미완성, 두 Zotero 모듈 응답 파싱 통일(`successful` vs `success`), 하드코딩 DB 접속 등.

## 함정·교훈 (실측 근거)

- **pt 하드필터 = 검증 A 3/9** (실측). 절대 검색식·기계필터에 넣지 말 것.
- **게이트 FP 유형 = 토픽 인접**(P만 또는 I만 매칭): 신장이식 goal-directed(P 불일치), 로봇 vs 복강경 직장절제(I 불일치). → 게이트 ① 의 P·I 동시 매칭으로 차단됨.
- claude 호출: `claude -p "<프롬프트>" --output-format json --model haiku`. 출력 = `{"type":"result","result":"<우리 JSON 문자열>",...}` 래퍼 → `llm_gate._extract_json` 가 처리.
- backfill(96건) ~20분(순차 게이트). 운영 cron(reldate 1)은 ~1분.
- `.env`(텔레그램·Zotero 키)·`state/`(seen-set jsonl) = gitignore. seen-set 현재 3건(RELIEF·OPTIMISE·GDFT메타).
- 발표자료(설계 권위) = vault `knowledge/02_areas/RPythonStudy/2026-06-22_2차_논문수집전략_slidable.md` (슬라이드 6=검색식 도출, 7=선별 필터, 8=두-검증).
