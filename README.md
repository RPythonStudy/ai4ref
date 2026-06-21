# ai4ref (Guideline-anchored New-Literature Surveillance)

진료지침의 Key Question(PICO)을 기준으로 PubMed 신규 문헌을 매일 감시하여,
지침과 다른 결론·practice-changing 랜드마크만 선별해 알림(Telegram)·보관(Zotero)한다.
결정론 backbone 위에 최소 LLM 게이트를 얹은, 추가비용 0의 1인 운영 파이프라인.

> 이 문서는 프로젝트 최상위 명세(L1)다. 피처 명세는 `specs/<NNN-feature>/`(L2),
> 컴포넌트 명세는 코드 옆 `<모듈>_spec.md`(L3)에 둔다. 설계 헌법 =
> `.specify/memory/constitution.md`.

## 핵심 설계 철학 (Core Design Philosophy)

- **결정론 우선 (Deterministic-First):** 검색·날짜필터·dedup 은 재현 가능한 결정론
  코드로 처리하고, LLM 은 관련성·랜드마크 의미판단에만 쓴다. (헌법 I)
- **2층 검증 분리 (Two-Layer Validation):** 검색 recall(두-검증 `validate_kq`)과
  선별 precision(LLM 게이트)을 별개·상보 레이어로 둔다. 한쪽이 다른 쪽을 대체하지
  않는다. (헌법 VI)
- **recall 편향 (Recall-biased):** 관련(P·I 일치)이면 경계선도 알린다. 놓침이 가장
  비싼 오류다. 기본값 loose. (헌법 VII)
- **추가비용 0 (Zero-Marginal-Cost):** 게이트 = 로컬 `claude_cli`(OAuth 구독·Haiku),
  cron = 로컬(OpenClaw). 종량 과금 API 미도입. (헌법 IX)

---

## 주요 기능 (L1 Features)

001. **KQ·PICO 정의 (Key Question & PICO Authoring)** — `key_questions.yml` 에 감시
   대상 지침의 KQ 를 PICO 로 정의한다.
   - **PICO = OR-블록**: `P`(대상)·`I`(중재)는 OR 리스트, `C`·`O`는 게이트 필터 기준.
   - **지침 앵커**: `guideline{name,pmid,date=T}` + Zotero `collection` + `enabled` 토글.
   - **두-검증 앵커**: `guideline_refs`(T 이전 근거 = recall 참값) /
     `post_guideline_landmarks`(T 이후 랜드마크 = precision 참값).
   *입력: 지침·근거문헌 PMID → 출력: KQ 레코드(PICO OR-블록 + 검증 앵커)*
   *제약: 검색식은 P·I 리스트에서 파생, 별도 query 필드 금지. (헌법 II·III)*

002. **KQ 추출 add-on (Guideline → PICO Extraction)** [신규] — 지침 원문에서 PICO 를
   반자동 추출해 001 레코드 초안을 생성한다.
   - 지침 → P·I·C·O 추출 → 검증 A 로 OR-블록 정련(recall 낮으면 synonym 추가) →
     `enabled:false` 레코드로 적재(수동 검토 후 활성화).
   *입력: 진료지침 텍스트 → 출력: enabled:false KQ 초안 + 검증 A recall 점수*
   *제약: Spec Kit(SDD) 신규 피처. 정련 = OR 리스트 성장만(별도 query 없음).*

003. **검색식 빌드 & 증분 수집 (Query Build & Incremental Collection)** — PICO 에서
   검색식을 파생해 PubMed 신규분만 증분 수집한다.
   - `build_query(pico)` = `(I 블록) AND (P 블록)` — 넓은 그물(결정론).
   - `esearch`(edat `reldate` 창 증분) → `efetch`(초록·메타).
   *입력: KQ PICO + reldate(N일) → 출력: 후보 문헌(PMID·제목·초록·저널·연도)*
   *제약: `[pt]` 등 설계 하드필터 검색식 투입 금지(검증 A 붕괴·신논문 누락). (헌법 IV)*

004. **기계 필터 & 중복 차단 (Deterministic Filter & Dedup)** — 후보를 결정론으로
   1차 거른다.
   - 날짜·언어 기계 필터 + `SeenStore`(seen-set jsonl) 중복 차단.
   *입력: 후보 문헌 → 출력: 미본 신규 후보*
   *제약: 결정론·멱등(재실행 안전). seen-set·state 는 gitignore.*

005. **두-검증 (Dual Validation — Recall)** — 검색식이 "아는 PMID"를 실제로 잡는지
   자동 채점한다.
   - 검증 A(`guideline_refs` 재현, 현행 8/9) · 검증 B(`post_guideline_landmarks`
     포착, 2/2).
   *입력: KQ 검증 앵커 → 출력: KQ 별 recall 점수(A/B)*
   *제약: golden PMID·정련 OR-리스트 = acceptance 참값(재작성 등가 기준). (헌법 XIII)*

006. **LLM 게이트 (Relevance & Landmark Gate — Precision)** — 후보 초록을 읽고 2단계로
   판정한다.
   - ① 관련성 = 대상(P) **AND** 중재(I) 동시 매칭(한쪽만 = 거부 → 토픽-인접 FP 차단).
   - ② 랜드마크 = 지침과 다른 결론·practice-changing. `C`·`O`·`question_type`·
     `design_strictness`(strict|loose)를 컨텍스트로 사용.
   - 백엔드 = `claude_cli`(로컬 OAuth·Haiku), 출력 JSON 파싱.
   *입력: 후보 초록 + KQ 게이트 컨텍스트 → 출력: {is_landmark, reason}*
   *사용자: 시스템(자동)*
   *제약: 추가비용 0(로컬). pt = 의미판단(검색식 아님). (헌법 V·IX)*

007. **Sink 멀티탭 (Notification Fan-out)** — 랜드마크를 다중 출력 통로로 fan-out 한다.
   - `notify/{stdout,telegram,zotero}` = 동일 `Sink` 플러그. 키 부재 시 fail-soft skip,
     `features.yml` 토글. Zotero = 랜드마크 영구 보관.
   *입력: 랜드마크 아이템 → 출력: Telegram 알림 · Zotero 보관 · stdout*
   *제약: 새 통로 추가가 본체(검색·선별)를 건드리지 않는다. (헌법 VIII)*

008. **자율 발화 cron (Autonomous Daily Trigger)** [미구현] — 매일 자동으로 파이프라인을
   기동한다.
   - OpenClaw 로컬 cron, `reldate=1`(지난 실행 이후) → 일일 3~4건 증분.
   *입력: 스케줄 → 출력: 일일 감시 실행(003→007)*
   *제약: 로컬 발화(추가비용 0). 컨테이너(OpenClaw)↔호스트(ai4ref·claude) 발화 위임.*

---

## 시스템 운영 및 설계 거버넌스 (Global Policies)

### 1. 결정론 backbone + LLM 경계
검색·날짜필터·dedup = 결정론(재현). 관련성·랜드마크 판정만 LLM. 경계를 넘나들지 않는다.

### 2. recall·precision 2층 분리
두-검증(검색 recall) ≠ 게이트(선별 precision). 별개 레이어로 운영·평가한다.

### 3. 추가비용 0 / 로컬 우선
게이트·cron 모두 로컬(claude_cli OAuth · OpenClaw). 유료 API 미도입.

### 4. config·데이터 거버넌스
운영 config = `key_questions.yml`·`features.yml`(YAML 2종). `state/`(seen-set)·`.env`
(시크릿)은 gitignore. 패키징 = `pyproject.toml`(uv·slim deps).
