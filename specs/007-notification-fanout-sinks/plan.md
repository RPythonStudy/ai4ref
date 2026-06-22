# Implementation Plan: Sink 멀티탭 (Notification Fan-out)

**Branch**: `007-notification-fanout-sinks` | **Date**: 2026-06-22 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/007-notification-fanout-sinks/spec.md`

## Summary

게이트(006)가 선별한 랜드마크를 features.yml 로 켜진 모든 출력 통로(sink)로 fan-out 하는 출력
레이어를 헌법 정합으로 정비한다. 동일 `Sink` 인터페이스(available·send)·토글 활성화·fail-soft
(prereq skip / 전송 격리)·플러그 확장(본체 0줄, 헌법 VIII)·무비용(stdlib + 지연로드, 헌법 III)을
보존하며, 부재했던 단위테스트를 추가하고 L3 명세를 작성한다. 기존 `alert.run` fan-out 동작 보존이
등가 기준.

## Technical Context

**Language/Version**: Python 3.11+ (uv)

**Primary Dependencies**: stdlib (`urllib`, `abc`, `dataclasses`). 선택: `pyzotero`(Zotero 통로
지연 로드, 미설치 시 skip). `common.features`(log·토글), `collector.zotero_add`(Zotero 메타 재사용).

**Storage**: None (전달만 — 중복 차단 seen-set 은 §004, 본 피처 밖)

**Testing**: `pytest` 단위테스트(`tests/unit/test_notify.py`) — 토글 활성화·prereq skip·전송
격리·render(내부필드 제외)·LandmarkItem 링크 파생을 외부 호출 모의로 결정론 검증. 텔레그램/Zotero
실 네트워크는 모의.

**Constraints**:
* 동일 `Sink` 인터페이스 — 본체는 통로 종류 불문 동일 처리 (FR-001).
* fail-soft 2종: 토글 꺼짐/prereq 미충족 → 활성 제외 skip (FR-002·003) / 전송 예외 → 격리·로그·계속 (FR-005).
* 새 통로 추가 = `SINK_CLASSES` 한 줄, 본체 0줄 수정 (헌법 VIII, FR-009·SC-004).
* stdout = 의존성 0·항상 동작, 외부 통로 render 는 stdout `render()` 재사용 (FR-006·007).
* 무비용 — stdlib + pyzotero 지연 로드 (헌법 III, FR-010).

## Constitution Check

| 원칙 | 적용 | 상태 |
|---|---|---|
| III 추가비용 0 / 로컬 | stdout·telegram=stdlib, zotero=pyzotero 지연로드(선택). 유료 없음. | ✅ |
| VIII 플러그 확장 | 새 통로 = SINK_CLASSES 한 줄, 본체(검색·선별) 불변. | ✅ |
| XI 3계층·줄수 | L3 명세(`notify_spec.md`) 작성, 코드 ≤300·`*_spec.md` ≤200줄. | ✅ |
| XII 절차적 명확성 | 명시적 토글·prereq·try 격리. 숨은 마법 없음. | ✅ |
| XVI No Duplication | 외부 통로 render = stdout `render()` 재사용(중복 제거). | ✅ |

## Project Structure

### Documentation (this feature)

```text
specs/007-notification-fanout-sinks/
├── spec.md
├── plan.md                 # 이 파일
├── checklists/requirements.md
└── tasks.md                # [NEXT] /speckit-tasks 산출
```

### Source Code

```text
src/notify/
├── base.py                 # Sink ABC · LandmarkItem (동작 보존)
├── registry.py             # build_sinks(토글·prereq) · fan_out(격리) (동작 보존)
├── stdout_sink.py          # render · 항상동작 (동작 보존)
├── telegram_sink.py        # env 토큰 · urllib (동작 보존)
├── zotero_sink.py          # 지연로드 · 컬렉션 · zotero_add 재사용 (동작 보존)
└── notify_spec.md          # [NEW] notify 패키지 L3 규격서 (Sink 계약·fan-out·fail-soft, ≤200줄)

src/alert/
└── run.py                  # build_sinks/fan_out 호출부 + 안전장치(stub→stdout만) 점검(보존)

tests/unit/
└── test_notify.py          # [NEW] 토글·prereq skip·전송 격리·render·LandmarkItem (모의)
```

## Phase 0 / 1 Notes (retro-spec — 경량)

기존 동작 역기록이라 신규 research·data-model·contracts 산출물은 만들지 않는다(003~006 동일).
- **계약**: `Sink.send(item)->ref|None` · `Sink.available()->(bool,why)` · `build_sinks(features)->list`
  · `fan_out(item, sinks)->{name: ref|None}`. 본문은 `notify_spec.md`(L3)로 일원화(헌법 XI).
- **등가 검증**: 기존 `alert.run` fan-out 동작 보존. 단위테스트는 토글·skip·격리·render 를 결정론 가드.
- **안전장치 경계**: 실수집+stub 백엔드 외부통로 차단은 `run.py` 책임 — 점검·보존만(§007 표면 밖).
