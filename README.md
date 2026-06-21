# ai4ref

**진료지침 기반 신규문헌 감시 + 주제 기반 대량수집** — PubMed 문헌을 두 트랙으로 다루는
1인 운영 파이프라인.

- **감시(alert)**: 진료지침의 Key Question(PICO)을 기준으로 신규 문헌을 매일 감시해,
  지침과 다른 결론·practice-changing 랜드마크만 선별해 **Telegram 알림 · Zotero 보관**.
- **수집(collect)**: 자연어 주제를 받아 대량 검색·본문 확보(PDF·XML)하여
  **Zotero · 로컬(Google Drive) · 2nd-brain** 으로 배포.

결정론 backbone 위에 최소 LLM(로컬 OAuth)을 얹어 **추가비용 0** 으로 운영한다.

> **상태**: 감시 MVP 동작(검증 A 8/9·B 2/2, 백필 FP 0) · 수집 트랙 개발 중 ·
> Spec Kit 기반 재작성 진행 중(`uv`/`pyproject` 전환 예정).

## 핵심 설계 철학

- **결정론 우선** — 검색·날짜필터·dedup 은 결정론, LLM 은 관련성·랜드마크 판정에만.
- **2층 검증 분리** — 검색 recall(두-검증)과 선별 precision(LLM 게이트)은 별개 레이어.
- **recall 편향** — 관련이면 경계선도 알린다. 놓침이 가장 비싼 오류(기본 loose).
- **추가비용 0** — 게이트·cron 모두 로컬(claude_cli OAuth · OpenClaw). 유료 API 미도입.

## 기능

상세 명세(입력/출력/제약)는 **[L1 시스템 명세](.specify/memory/system-spec.md)** 참조.

| # | 트랙 | 기능 | 상태 |
|---|---|---|---|
| 001 | 감시 | KQ·PICO 정의 (`key_questions.yml`) | ✅ |
| 002 | 감시 | KQ 추출 add-on (지침→PICO 반자동) | 신규 |
| 003 | 감시 | 검색식 빌드 & 증분 수집 (esearch/efetch) | ✅ |
| 004 | 감시 | 기계 필터 & 중복 차단 (seen-set) | ✅ |
| 005 | 감시 | 두-검증 (검색 recall 자동채점) | ✅ |
| 006 | 감시 | LLM 게이트 (관련성·랜드마크 2단계) | ✅ |
| 007 | 감시 | Sink 멀티탭 (Telegram·Zotero·stdout) | ✅ |
| 008 | 감시 | 자율 발화 cron (일일 자동 기동) | 미구현 |
| 009 | 수집 | 주제 수집 정의 (자연어→config) | 신규 |
| 010 | 수집 | 대량수집 & 본문 확보 (PMC PDF/XML) | 신규 |
| 011 | 수집 | 자산 배포 (Zotero linked + 2nd-brain) | 신규 |

## Quickstart

### 요구사항
- Python 3.11+ 와 [uv](https://docs.astral.sh/uv/) (패키지·가상환경 관리)
- `.env` 시크릿: `ZOTERO_API_KEY`·`ZOTERO_USER_ID`(write 권한), `ENTREZ_EMAIL`,
  `PUBMED_API_KEY`(선택), Telegram 토큰(알림 사용 시). `.env.example` 참고.
- 게이트(라이브 판정)용 로컬 [Claude Code](https://claude.com/claude-code) CLI(`claude`, OAuth).

### 설치
```bash
git clone https://github.com/RPythonStudy/ai4ref.git
cd ai4ref
uv sync                # .venv 생성 + 의존성 설치 (uv.lock 재현)
cp .env.example .env   # 키 채우기
```

### 실행 (감시 트랙)
```bash
uv run python -m alert.run --demo          # 배관 검증(가짜 랜드마크)
uv run python -m alert.run --collect-only  # 실 PubMed 후보만(전송 안 함)
uv run python -m alert.run                 # 전체(reldate 30, 게이트 라이브 → 알림)
uv run python -m alert.validate_kq         # 두-검증 자동채점 (A 8/9 · B 2/2)
```
- 게이트는 `features.yml` 의 `llm_backend: claude_cli` 일 때 라이브 판정(추가비용 0).
- 수집 트랙(009~011)은 개발 중 — 사용법은 구현 후 추가.

## 문서

- [L1 시스템 명세](.specify/memory/system-spec.md) — 기능 카탈로그(입력/출력/제약)
- [설계 헌법](.specify/memory/constitution.md) — 불변 원칙

## 라이선스

미정 (1인 연구 도구).
