# Module Specification: common.state (L3)

**File**: [state.py](file:///home/ben/projects/ai4ref/src/common/state.py) | **Spec Version**: 1.0 | **Code Lines Constraint**: <= 300줄 | **Spec Lines Constraint**: <= 200줄

## Summary

감시 파이프라인에서 이미 발송 완료된 랜드마크 문헌(LandmarkItem)의 이력을 관리하고, 중복 전송을 차단하는 SeenStore를 정의합니다. JSONL 파일 백엔드를 통해 경량화된 영구 로그 기능을 수행합니다.

## SeenStore Class API

```python
class SeenStore:
    def __init__(self, path: str): ...
```
* **역할**: 지정된 경로의 JSONL 파일로부터 보낸 이력 키를 로드하여 메모리 셋(`_seen`)으로 관리합니다.
* **매개변수**:
  * `path`: 기록용 JSONL 파일 경로 (상대 경로 지정 시 프로젝트 루트 기준 절대 경로로 자동 계산).

### `is_seen`
```python
def is_seen(self, key: str) -> bool:
```
* **역할**: 주어진 고유 키(형식: `"{kq}:{pmid}"`)가 이미 발송 이력에 존재하는지 검사합니다.
* **출력**: 중복 여부 (`True` / `False`).

### `mark`
```python
def mark(self, key: str, record: dict | None = None) -> None:
```
* **역할**: 전달 완료된 랜드마크의 키와 메타데이터를 백엔드 JSONL 파일 끝에 추가 기록(append)하고 메모리 셋에 추가합니다.
* **매개변수**:
  * `key`: 고유 식별 키 (예: `"ERAS 수액전략:29742967"`).
  * `record`: 추가적으로 기록할 세부 메타데이터 딕셔너리 (선택사항).
* **설계**: 상위 폴더가 존재하지 않는 경우 자동으로 폴더를 생성(`mkdir(parents=True)`)하며, 기록 시 UTC 타임스탬프(`ts`)를 기본 탑재합니다.

## Fail-Soft Loading Policy

* 파일을 한 줄씩 파싱하여 `json.loads`를 수행하는 도중 예외(손상된 줄, 구문 오류, 빈 행)가 발생하더라도 **에러를 삼키고 로그를 남긴 후 다음 라인을 계속해서 파싱**합니다 (FR-004).
* 지정된 경로에 파일이 존재하지 않는 경우 예외를 내지 않고 빈 상태로 생성합니다.
