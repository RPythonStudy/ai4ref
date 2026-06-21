# Quickstart: KQ·PICO 정의 등가 검증

001 구현이 기존 동작과 등가인지(헌법 XIII) 확인하는 실행 가이드.

## 전제

- `uv sync` 완료, `config/key_questions.yml`(ERAS 2012 KQ 1건) 존재.
- acceptance 참값: `tests/acceptance/recall_baseline.yml`.

## 1. 검색식 파생 등가 (결정론, 네트워크 불요)

```bash
uv run python -c "import yaml; from common.query import build_query; \
r=yaml.safe_load(open('config/key_questions.yml'))['kqs'][0]; print(build_query(r['pico']))"
```
**기대**: 출력이 `recall_baseline.yml` 의 `key_questions[0].derived_query` 와 **문자 단위
일치**(`(I…) AND (P…)`, C·O 미포함).

## 2. KQ 로드·검증 (결정론)

```bash
uv run python -c "from common.kq import load_kqs; ks=load_kqs(); print(len(ks), ks[0]['kq'])"
```
**기대**: 1건 로드, ERAS KQ 제목 출력. 잘못된 KQ(P/I 빈 리스트·query 필드)는 거부됨.

## 3. recall 참값 재현 (PubMed, 005 채점 활용)

```bash
uv run python -m alert.validate_kq
```
**기대**: 검증 A **8/9**(누락 19602972 = 의도된 천장) · 검증 B **2/2** —
`recall_baseline.yml` 과 일치.

## 4. 단위 테스트

```bash
uv run pytest tests/unit/test_query.py -q
```
**기대**: `build_query` 출력 == acceptance `derived_query` (등가 통과).

## 등가 판정

위 1·4(결정론) + 3(recall 참값) 모두 통과 시 001 재작성을 기존 MVP 와 등가로 인정한다.
