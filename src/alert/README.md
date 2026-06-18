# ai4ref — 신규치료 논문 알림 (guideline-anchored alert) PoC

RPython 연구회 요청: *특정 임상영역에서 새 치료법(중재 비교 RCT)이 제안되면 Telegram 으로 1편* 보내기.
설계 권위 = vault `[[2025-08_ai4ref]]` §응용 PoC (계층형 SR/MA-capable 의 T0 검색 + T1-경량 선별).

## 파이프라인

```
진료지침-앵커 query (P·도메인·RCT 게이트, I/C 개방)
  └ rolling 윈도우 → esearch ──────────────── 검색(recall)
       └ esummary 메타 → LLM 게이트("권고 변경 근거?") → 저널티어·최신 랭킹 → top-1   ── 선별
            └ Telegram (Bot API; 토큰 없으면 dry-run) ── 발송
```

## 사용법

```bash
# 검증 모드 (기본): ERAS 수액 앵커 + 2018 윈도우 → RELIEF(PMID 29742967) top-1 기대
python3 src/alert/guideline_alert.py validate

# 실서비스 시뮬 (rolling, 최근 N일)
python3 src/alert/guideline_alert.py rolling 30
ALERT_DELIVER=1 python3 src/alert/guideline_alert.py rolling 30   # 실제 Telegram 발송
```

검증 출력 예: `[검증] 기대 PMID=29742967 | 선택=29742967 → ✅ PASS`

## 설정

- **앵커**: `config/alert_anchors.yml` — `mode: guideline`(우선) / `topic`(폴백).
  규칙: **P(대상·병태)·O(임상결과)는 좁게, I/C(중재)는 개방** → 기존 중재에 안 묶여 *새* 치료법 포착.
- **환경변수** (`.env`):
  - `ENTREZ_EMAIL`, `PUBMED_API_KEY` (E-utilities; 키 있으면 10 req/s)
  - `TELEGRAM_BOT_TOKEN`, `TELEGRAM_CHAT_ID` (없으면 dry-run 미리보기)

## 설계 주의 (vault hub "회고 ≠ 전향")

- 저널티어 표(`JOURNAL_TIER`)는 영향력의 *근사 휴리스틱* — 단일 저널 하드코딩 X, 표로 일반화.
- "새 치료법" 판정의 원칙적 층 = **LLM 게이트**. 현재 `llm_practice_changing_gate()` 는
  **heuristic STUB**(RCT + 제목 비교/시험 신호). production 은 Claude 에 title+abstract 분류 위임.
- 의존성 0 (stdlib urllib), DB·biopython 불요 — 즉시 실행/검증 가능한 PoC.

## production 전환 TODO

- [ ] LLM 게이트 → Claude 실연결 (efetch 로 abstract 확보 후 분류·PICO 추출)
- [ ] 검색·수집 → `src/collector.pubmed_esearch`(Entrez) + postgres 재사용 (DRY)
- [ ] 발송 → OpenClaw Telegram (cron 자동발화; gateway 경유)
- [ ] 앵커 다건 + dedup(이미 보낸 PMID) + rolling cron
- [ ] precision 튜닝/검증 (known-positive 셋 확장: RELIEF 외)
