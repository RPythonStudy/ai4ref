# Module Specification: alert.llm_gate (L3)

**File**: [llm_gate.py](file:///home/ben/projects/ai4ref/src/alert/llm_gate.py) | **Spec Version**: 1.0 | **Code Lines Constraint**: <= 300줄 | **Spec Lines Constraint**: <= 200줄

## Summary

후보 논문의 제목·초록을 LLM 으로 읽어 **2단계**로 판정하는 정밀(precision) 레이어다.
① 관련성 = KQ 의 대상(P) AND 중재(I) 동시매칭 → ② 랜드마크 = 현행 지침과 다른
practice-changing. 판정 엔진(백엔드)은 sink 처럼 교체형이며 운영 백엔드는 로컬 `claude`
CLI(OAuth·Haiku, 추가비용 0). recall(두-검증 `validate_kq`)과 별개 레이어(헌법 II).

## API Specifications

### `judge`
```python
def judge(paper: dict, kq: dict, backend: str = "stub") -> dict:
```
* **역할**: 후보 1건을 2단계 판정해 게이트 verdict 를 반환.
* **매개변수**:
  * `paper`: `{"title": str, "abstract": str}` — 후보 논문.
  * `kq`: 게이트 컨텍스트 `{kq, comparison(C), outcome(O), guideline, question_type, strictness}`
    (`run.py` 가 `key_questions.yml` 에서 구성).
  * `backend`: `"stub"`(개발용·키불필요) | `"claude_cli"`(운영권장·로컬 OAuth) | `"api"`(미구현→stub 강등).
* **출력** (FR-006, acceptance `gate_baseline.yml` 스키마):
  ```python
  {"is_relevant": bool, "is_landmark": bool, "reason": str}  # reason = 한국어 한 줄
  ```
* **2단계 순서 (FR-001·002·SC-005)**: ② is_landmark 는 ① is_relevant=true 일 때만 의미를 가진다.
  관련 없음이면 is_landmark=false.

### `_claude_cli`
```python
def _claude_cli(paper, kq, model: str = "haiku") -> dict:
```
* 로컬 `claude -p <PROMPT> --output-format json --model haiku` 호출(subprocess, timeout 120s).
* **컨텍스트 주입 (FR-005)**: PROMPT 에 C·O·guideline·기대설계(design)·strictness·title·abstract.
  초록은 **3000자 절단**. `question_type` → `DESIGN_LABEL` 매핑(미지정 시 "해당 유형의 적정 연구설계").
* **strict/loose (FR-004)**: strict = 기대 설계 수준 근거 요구 / loose = 설계 무관·경계선 포함(recall 편향, 기본).
* **fail-soft (FR-009)**: 호출 실패·`is_error` 래퍼·JSON 파싱 실패 시 보수적 보류
  `{is_relevant:false, is_landmark:false, reason:"판정 실패: …"}` + `log.error`, 예외 전파 금지.

### `_stub`
* 골격 검증용 — 항상 `{is_relevant:true, is_landmark:true, reason:"[stub] …"}`. 키 불필요.
* ⚠️ `run.py`: 비-demo 실수집 + stub 백엔드는 모두 랜드마크가 되어 스팸 → stdout sink 로만 제한.

### `_extract_json`
```python
def _extract_json(txt: str) -> dict:
```
* claude CLI 출력에서 판정 JSON 견고 추출 (FR-010):
  * 평문 JSON / `{"result": "...json 문자열..."}` 래퍼(재귀) / `{"result": {...}}` / 앞뒤 텍스트 둘러싼 경우(첫 `{`~마지막 `}`).
  * 추출 실패 시 `{}` 반환(호출부가 보류 처리).

## 상수
* `DESIGN_LABEL`: question_type(intervention/diagnostic/prognostic/predictive) → 기대 연구설계 라벨.
  strict 판정의 컨텍스트로만 사용 — **검색식 하드필터 아님**(헌법 IV·V·IX).

## 보조 모듈 — `classify_qtype.py` (0단계)
* 자연어 주제 → question_type **제안**(heuristic). 권위는 **사용자 `key_questions.yml` 입력**(결정형).
  게이트는 사용자 입력 question_type 을 컨텍스트로 받을 뿐, 자동 분류가 판정을 좌우하지 않는다 (FR-011).

## 등가 기준 (원칙 XIII)
* `tests/acceptance/gate_baseline.yml` 재현:
  * `landmark_true` 2/2 — RELIEF(29742967)·OPTIMISE(24842135): is_relevant=true & is_landmark=true.
  * `reject_false_positive` 2/2 — 신장이식(P 불일치)·로봇수술(I 불일치): is_relevant=false.
* LLM 비결정성은 acceptance 참값 책임. 단위테스트(`tests/unit/test_gate.py`)는 배관(2단계 순서·
  fail-soft·`_extract_json`·백엔드 디스패치·PROMPT 구성)만 백엔드 모의로 결정론 검증.
