# ai4ref — alert (진료지침 기반 신규문헌 감시)

KQ(Key Question)마다 *검색 전략을 설정파일에 미리 구성*해 두고, 지침 제정(T) 이후 새로 색인된
문헌에서 **지침을 바꿀 랜드마크**를 골라 알림·보관한다.

```
[KQ 앵커] kq_anchors.yml (term·guideline·검증셋, 미리 구성)
   └ PubMed 증분 검색(edat reldate, pdat≥T)  ── 결정론
        └ 1차 기계 필터(날짜·중복·언어)          ── 결정론
             └ efetch(제목·초록)
                  └ LLM 게이트(관련·설계·랜드마크)  ── 추론(claude_cli, OAuth)
                       └ dedup(seen-set, jsonl)
                            └ sink fan-out → stdout · Telegram · Zotero
```

## 사용법

```bash
# 배관 검증 (가짜 랜드마크 2건)
python src/alert/run.py --demo

# 실 PubMed 후보만 확인 (판정·전송 없음)
python src/alert/run.py --collect-only --reldate 60 --limit 15

# 전체 실행 (claude_cli 게이트 → 진짜 랜드마크만 알림)
python src/alert/run.py

# KQ 검증 (검증 A 8/9 + 검증 B 2/2 자동 채점)
python src/alert/validate_kq.py
```

## 설정 (단일 정본)

- **`config/kq_anchors.yml`** — KQ 앵커 레코드(진료지침 감시 + 일반 수집 통합).
  `guideline` 있는 레코드 = 감시 KQ(두-검증). 구 `alert_anchors.yml` 은 여기로 통합·폐기(2026-06-21).
  - `kq·question_type·pico·term·guideline{name,pmid,date=T}·guideline_refs(검증A)·post_guideline_landmarks(검증B)·enabled`
- **`config/features.yml`** — 기능/싱크 토글 (alert·notify.{stdout,telegram,zotero}·llm_backend·state).
- **`.env`** — `TELEGRAM_BOT_TOKEN/CHAT_ID`, `ZOTERO_API_KEY/USER_ID`, `ENTREZ_EMAIL`(선택).
  없으면 해당 sink/기능은 fail-soft 로 자동 skip.

## 모듈

| 파일 | 역할 |
|---|---|
| `run.py` | 오케스트레이터 (수집→필터→게이트→dedup→fan-out) |
| `llm_gate.py` | 교체형 LLM 백엔드(stub\|claude_cli\|api) — 관련·랜드마크 판정 |
| `classify_qtype.py` | 0단계 유형 판별 *보조* (권위는 사용자 입력) |
| `validate_kq.py` | 두-검증 자동 채점 (검증A 재현 / 검증B 포착) |
| `../notify/` | sink 멀티탭 (base·registry·stdout·telegram·zotero) |
| `../common/` | features(플래그)·state(seen-set)·pubmed(E-utilities) |

## 두-검증 (시간축 T = 지침 제정)

- **검증 A** (`guideline_refs`, T 이전): 검색식이 지침 근거를 재현하나 = recall. ERAS 2012 → **8/9**.
  (`19602972` Marik=ICU ventilated 간접근거, 검색 한계로 누락 = 정직한 천장)
- **검증 B** (`post_guideline_landmarks`, T 이후): alert 가 지침-이후 랜드마크를 포착하나. RELIEF·OPTIMISE → **2/2**.
- 핵심: 검증 논문 ≠ 시드. 시간축으로 분리 → 순환논증 회피.

## 설계 원칙

- **결정론 backbone + 최소 LLM**: 검색·날짜필터·dedup = 결정론(재현) / 관련·랜드마크 판정 = LLM 의미판단.
- **recall 편향**: 게이트는 "애매하면 포함"(색인 전 신논문·랜드마크 누락 방지). 설계(출판유형)는 하드필터 아닌 의미판단.
- **추가비용 0**: 게이트 = 로컬 Claude Code(OAuth 구독). API 토큰 과금 없음.

## 남은 일

- [ ] OpenClaw cron 등록 (매일/매주 자동 발화; Docker 컨테이너 ↔ 호스트 ai4ref 발화 방식 결정)
- [ ] KQ 추출 add-on (지침 → KQ 자동 추출 → enabled:false 토글 레코드)
